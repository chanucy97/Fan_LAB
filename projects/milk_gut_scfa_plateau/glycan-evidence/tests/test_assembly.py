"""Small assembly adapter tests; mock producers are not biological evidence."""
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from glycan_evidence import assembly


class AssemblyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.r1, self.r2 = self.root / "r1.fq", self.root / "r2.fq"
        self.r1.write_bytes(b"@a/1\nACGT\n+\nIIII\n")
        self.r2.write_bytes(b"@a/2\nACGT\n+\nIIII\n")
        self.ref = self.root / "source.fasta"
        self.ref.write_text(">reference\nACGT\n")
        self.prefix = self.root / "index"
        files = []
        for suffix in ["1", "2", "3", "4", "rev.1", "rev.2"]:
            path = self.root / ("index." + suffix + ".bt2")
            path.write_bytes(("artificial index " + suffix).encode())
            files.append({**assembly.record(path), "suffix": suffix})
        self.manifest = self.root / "index.json"
        assembly.write(self.manifest, {"schema": "glycan-evidence.recruitment-index.v1", "status": "VERIFIED",
                       "reference_fasta": assembly.record(self.ref), "source_genomes": ["toy_reference"],
                       "index_files": files, "index_sequence_verification": {"status": "PASS", "canonical_sequence_sha256": "a" * 64}})
        self.megahit = self.root / "megahit"
        self.megahit.write_text("artificial launcher")
        self.out = self.root / "result"

    def screen(self, state, r1, r2, panel, threads, bowtie, fastp):
        for mate, path in [(1, r1), (2, r2)]:
            with gzip.open(state.out / f"candidate_R{mate}.fastq.gz", "wb") as handle:
                if self.candidates:
                    handle.write(Path(path).read_bytes())
        return {"total_fragments": 1, "candidate_fragments": self.candidates}

    def fake_subprocess(self, argv, stdout=None, **kwargs):
        argv = list(map(str, argv))
        if "--version" in argv:
            stdout.write("MEGAHIT v1.2.9\n")
        else:
            if self.fail:
                return subprocess.CompletedProcess(argv, 7)
            destination = Path(argv[argv.index("-o") + 1])
            destination.mkdir()
            (destination / "final.contigs.fa").write_text(">contig_1\nACGTACGT\n>contig_2\nACGT\n")
        return subprocess.CompletedProcess(argv, 0)

    def run_mock(self, candidates=1, fail=False, **kwargs):
        self.candidates, self.fail = candidates, fail
        def version(state, tool, prefix):
            state.data["versions"][tool] = "synthetic software fixture"
        with patch.object(assembly.frozen, "command_prefix", side_effect=lambda tool, _: [tool]), \
             patch.object(assembly.frozen, "capture_version", side_effect=version), \
             patch.object(assembly.frozen, "screen", side_effect=self.screen), \
             patch.object(assembly.subprocess, "run", side_effect=self.fake_subprocess):
            return assembly.run_assembly(self.r1, self.r2, "TOY", self.prefix, self.ref, self.manifest,
                                         self.out, self.megahit, **kwargs)

    def test_success_has_verified_inputs_and_contigs(self):
        result = self.run_mock(expected_qc_pairs=1)
        self.assertEqual(result["assembly_state"], "assembled")
        self.assertEqual(result["assembly_summary"], {"contigs": 2, "bases": 12, "max_length": 8, "N50": 8})
        self.assertTrue(result["raw_input_validation"]["both_eof"])
        self.assertEqual((self.out / "COMPLETE").read_text().strip(), assembly.record(self.out / "receipt.json")["sha256"])
        self.assertTrue((self.out / "all_contigs.fasta").is_file())

    def test_zero_recruitment_is_completed_unresolved(self):
        result = self.run_mock(candidates=0)
        self.assertEqual(result["assembly_state"], "unresolved_no_recruited_reads")
        self.assertEqual((self.out / "all_contigs.fasta").stat().st_size, 0)
        self.assertEqual(len(result["commands"]), 1)
        self.assertEqual(result["assembly_summary"]["contigs"], 0)

    def test_assembler_failure_does_not_publish(self):
        with self.assertRaises(subprocess.CalledProcessError):
            self.run_mock(fail=True)
        self.assertFalse(self.out.exists())
        failed = list(self.root.glob(".result.running-*/receipt.json"))
        self.assertEqual(len(failed), 1)
        self.assertEqual(json.loads(failed[0].read_text())["status"], "FAILED")
        self.assertFalse((failed[0].parent / "COMPLETE").exists())

    def test_existing_output_is_not_overwritten(self):
        self.out.mkdir()
        (self.out / "user.txt").write_text("preserve")
        with self.assertRaises(FileExistsError):
            self.run_mock()
        self.assertEqual((self.out / "user.txt").read_text(), "preserve")

    def test_expected_complete_qc_count_is_enforced(self):
        with self.assertRaisesRegex(ValueError, "QC count"):
            self.run_mock(expected_qc_pairs=2)
        self.assertFalse(self.out.exists())

    def test_changed_index_content_is_rejected(self):
        (self.root / "index.1.bt2").write_bytes(b"corrupt")
        with self.assertRaisesRegex(ValueError, "index content hash"):
            self.run_mock()

    def test_changed_source_fasta_is_rejected(self):
        self.ref.write_text(">reference\nTGCA\n")
        with self.assertRaisesRegex(ValueError, "FASTA differs"):
            self.run_mock()

    def test_duplicate_manifest_index_member_is_rejected(self):
        data = json.loads(self.manifest.read_text())
        data["index_files"][0]["suffix"] = data["index_files"][1]["suffix"]
        assembly.write(self.manifest, data)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.run_mock()

    def test_extra_index_member_is_rejected(self):
        (self.root / "index.1.bt2l").write_text("unexpected mixed format")
        with self.assertRaisesRegex(ValueError, "exactly six"):
            self.run_mock()

    def test_invalid_fasta_duplicate_or_empty_rejected(self):
        for text in (">a\nACGT\n>a\nACGT\n", ">a\n", "ACGT\n", ">a\nHELLO\n"):
            path = self.root / "bad.fa"
            path.write_text(text)
            with self.assertRaises(ValueError):
                assembly.fasta_metadata(path)

    def test_unsafe_sample_or_invalid_resource_rejected(self):
        with self.assertRaisesRegex(ValueError, "positive"):
            self.run_mock(threads=0)


if __name__ == "__main__":
    unittest.main()
