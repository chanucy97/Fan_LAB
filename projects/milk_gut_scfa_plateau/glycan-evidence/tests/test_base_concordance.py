"""Tiny provenance fixtures and real-BAM regression tests, not a benchmark."""
import copy
import csv
import gzip
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from glycan_evidence.base_concordance import (
    POSITION_FIELDS, SCHEMA_VERSION, _aligned_reference_pairs, count_bam, run_base_concordance,
)
from glycan_evidence.integrity import digest, seal_outputs, verify_outputs
from glycan_evidence.structure import (
    STRUCTURE_FIELDS, file_record, load_config, reverse_complement, sequence_sha, write_tsv,
)

HAVE_PYSAM = importlib.util.find_spec("pysam") is not None


class CigarWalkerTests(unittest.TestCase):
    def test_padding_never_consumes_query_or_reference(self):
        self.assertEqual(list(_aligned_reference_pairs([(0, 1), (6, 1), (0, 3)], 0, 4, 4)),
                         [(0, 0), (1, 1), (2, 2), (3, 3)])

    def test_all_consumption_classes_with_hard_soft_insertion_deletion_skip(self):
        cigar = [(5, 2), (4, 1), (0, 1), (1, 1), (7, 1), (2, 1), (3, 1), (8, 1), (4, 1), (5, 2)]
        self.assertEqual(list(_aligned_reference_pairs(cigar, 1, 6, 6)),
                         [(1, 1), (3, 2), (None, 3), (None, 4), (4, 5)])

    def test_adjacent_match_equal_mismatch_have_same_coordinates(self):
        self.assertEqual(list(_aligned_reference_pairs([(0, 1), (7, 1), (8, 1)], 0, 3, 3)),
                         list(_aligned_reference_pairs([(0, 3)], 0, 3, 3)))

    def test_illegal_operator_length_clipping_and_bounds_are_rejected(self):
        cases = [([(9, 1)], 0, 1, 6), ([(0, 0)], 0, 0, 6), ([(0, -1)], 0, 0, 6),
                 ([(True, 1)], 0, 1, 6), ([(0, True)], 0, 1, 6), ([(0, 1, 2)], 0, 1, 6),
                 ([(0, 1), (5, 1), (0, 1)], 0, 2, 6),
                 ([(0, 1), (4, 1), (0, 1)], 0, 3, 6),
                 ([(0, 4)], 0, 3, 6), ([(0, 4)], 3, 4, 6), ([(0, 1)], -1, 1, 6),
                 ([], 0, 0, 6)]
        for args in cases:
            with self.subTest(args=args), self.assertRaises(ValueError):
                list(_aligned_reference_pairs(*args))


def candidate(reverse=False):
    seq = "ACGTAC"
    return dict(structure_id="structure_000001", contig="ctg", contig_length_bp=len(seq),
                contig_sequence_sha256=sequence_sha(seq), interval_start_0based=0,
                interval_end_0based_exclusive=len(seq), orientation_to_reference="-" if reverse else "+",
                interval_sequence_sha256=sequence_sha(reverse_complement(seq) if reverse else seq),
                left_anchor_boundary_0based=1, right_anchor_boundary_0based=5,
                left_anchor_id="left", right_anchor_id="right", inner_anchor_edges_complete=True,
                left_anchor_inner_edge_complete=True, right_anchor_inner_edge_complete=True)


def alignment(name, seq, flag=0, mapq=60, bq=35, start=0, cigar=None):
    import pysam
    record = pysam.AlignedSegment()
    record.query_name = name
    record.query_sequence = seq
    record.flag = flag
    record.reference_id = 0
    record.reference_start = start
    record.mapping_quality = mapq
    record.cigartuples = cigar or [(0, len(seq))]
    record.query_qualities = [bq] * len(seq)
    return record


class SourceFixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source_dir = self.root / "structure"
        self.source_dir.mkdir()
        self.contigs = self.root / "contigs.fasta"
        self.contigs.write_text(">ctg\nACGTAC\n", encoding="ascii")
        self.table = self.source_dir / "structures.tsv"
        self.receipt_path = self.source_dir / "receipt.json"
        self.out = self.root / "bases"

    def source(self, rows=None, records=None):
        rows = [] if rows is None else rows
        write_tsv(self.table, rows, STRUCTURE_FIELDS)
        value = dict(schema_version="glycan-evidence.structure/1", status="COMPLETE",
                     structure_count=len(rows), contigs=file_record(self.contigs),
                     outputs=[file_record(self.table)],
                     readback=dict(status="NOT_REQUESTED", interpretation="unresolved_without_readback"))
        if records is not None:
            import pysam
            from glycan_evidence.readback import analyze_readback
            path = self.root / "input.bam"
            with pysam.AlignmentFile(str(path), "wb", header={"HD": {"VN": "1.6", "SO": "unsorted"},
                                                             "SQ": [{"SN": "ctg", "LN": 6}]}) as stream:
                for record in records:
                    stream.write(record)
            # Validation evidence is produced by the real existing reader.
            value["alignment_input"] = file_record(path)
            value["readback"] = analyze_readback(path, {"ctg": "ACGTAC"}, rows,
                                                   load_config()["readback"], work_dir=self.root)["receipt"]
        self.save(value)
        (self.source_dir / "COMPLETE").write_text("2026-10-09T00:00:00+00:00\n", encoding="ascii")
        return value

    def save(self, value):
        self.receipt_path.write_text(json.dumps(value), encoding="utf-8")

    def run_bases(self):
        return run_base_concordance(self.receipt_path, self.table, self.out)

    def positions(self):
        with gzip.open(self.out / "structure_000001.positions.tsv.gz", "rt", encoding="utf-8", newline="") as stream:
            return list(csv.DictReader(stream, delimiter="\t"))


class BaseInputTests(SourceFixture):
    def test_zero_candidates_complete_without_importing_alignment_dependency(self):
        self.source()
        result = self.run_bases()
        self.assertEqual(result["schema_version"], SCHEMA_VERSION)
        self.assertEqual(result["status"], "COMPLETE")
        self.assertEqual(result["evidence_status"], "not_applicable_no_candidates")
        self.assertIn("unresolved", result["interpretation"])
        self.assertEqual(result["candidates"], [])
        self.assertEqual(result["files"], {})
        self.assertFalse(result["bam_sha_reverified"])
        self.assertFalse(result["alignment_scan_performed"])
        self.assertEqual(result["alignment_revalidation"], "not_revalidated_no_alignment_requested")
        self.assertEqual(list(self.out.glob("*.positions.tsv.gz")), [])
        self.assertEqual((self.out / "COMPLETE").read_text().strip(), digest(self.out / "receipt.json"))
        seal_outputs(self.out)
        self.assertEqual(verify_outputs(self.out)["status"], "PASS")

    def test_requires_source_complete_marker(self):
        self.source()
        (self.source_dir / "COMPLETE").unlink()
        with self.assertRaisesRegex(ValueError, "COMPLETE marker"):
            self.run_bases()
        self.assertFalse(self.out.exists())

    def test_requires_complete_supported_source_receipt(self):
        source = self.source()
        for field, value in [("status", "RUNNING"), ("schema_version", "other/1")]:
            with self.subTest(field=field):
                broken = copy.deepcopy(source)
                broken[field] = value
                self.save(broken)
                with self.assertRaisesRegex(ValueError, "COMPLETE glycan-evidence.structure/1"):
                    self.run_bases()

    def test_bad_contig_sha_and_structure_content_are_rejected(self):
        source = self.source()
        broken = copy.deepcopy(source)
        broken["contigs"]["sha256"] = "0" * 64
        self.save(broken)
        with self.assertRaisesRegex(ValueError, "SHA/size mismatch: contigs"):
            self.run_bases()
        self.save(source)
        self.table.write_text(self.table.read_text() + "\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "SHA/size mismatch: structures"):
            self.run_bases()

    def test_source_must_bind_table_not_only_candidate_count(self):
        source = self.source()
        source["outputs"] = []
        self.save(source)
        with self.assertRaisesRegex(ValueError, "exactly one structures.tsv"):
            self.run_bases()

    def test_nonzero_candidates_cannot_succeed_without_alignment(self):
        self.source([candidate()])
        with self.assertRaisesRegex(ValueError, "Nonzero candidates require"):
            self.run_bases()
        self.assertFalse(self.out.exists())

    def test_zero_candidates_cannot_hide_missing_requested_readback(self):
        source = self.source()
        for invalid in ({"status": "FAILED"}, None, [], "NOT_REQUESTED"):
            with self.subTest(readback=invalid):
                source["readback"] = invalid
                self.save(source)
                with self.assertRaisesRegex(ValueError, "only zero candidates allow NOT_REQUESTED"):
                    self.run_bases()

    def test_bad_candidate_hash_coordinates_and_duplicate_id_are_rejected(self):
        for field, value, message in [
                ("interval_sequence_sha256", "0" * 64, "interval sequence SHA"),
                ("contig_sequence_sha256", "0" * 64, "contig sequence SHA"),
                ("interval_end_0based_exclusive", 7, "Invalid candidate interval"),
                ("orientation_to_reference", "?", "orientation must"),
                ("structure_id", "../escape", "Unsafe or duplicate")]:
            with self.subTest(field=field):
                row = candidate()
                row[field] = value
                self.source([row])
                with self.assertRaisesRegex(ValueError, message):
                    self.run_bases()
        self.source([candidate(), candidate()])
        with self.assertRaisesRegex(ValueError, "Unsafe or duplicate"):
            self.run_bases()

    def test_zero_candidate_count_must_match_bound_table(self):
        source = self.source()
        source["structure_count"] = 1
        self.save(source)
        with self.assertRaisesRegex(ValueError, "structure_count differs"):
            self.run_bases()

    def test_existing_output_is_preserved(self):
        self.source()
        self.out.mkdir()
        sentinel = self.out / "user.txt"
        sentinel.write_text("keep", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            self.run_bases()
        self.assertEqual(sentinel.read_text(), "keep")


@unittest.skipUnless(HAVE_PYSAM, "Actual BAM fixtures require the structure extra (pysam)")
class BaseBamTests(SourceFixture):
    def test_padding_partial_interval_preserves_nucleotides(self):
        row = candidate()
        row.update(interval_end_0based_exclusive=3, interval_sequence_sha256=sequence_sha("ACG"),
                   right_anchor_boundary_0based=2)
        self.source([row], [alignment("padded", "ACGT", cigar=[(0, 1), (6, 1), (0, 3)])])
        self.run_bases()
        rows = self.positions()
        self.assertEqual([r["assembled_genomic_base"] for r in rows], list("ACG"))
        self.assertEqual([r[b + "_qname_total"] for r, b in zip(rows, "ACG")], ["1"] * 3)
        self.assertEqual([r["assembly_fraction"] for r in rows], ["1.0"] * 3)

    def test_padding_full_interval_and_quality_indices_are_preserved(self):
        row = candidate()
        row.update(interval_end_0based_exclusive=4, interval_sequence_sha256=sequence_sha("ACGT"),
                   right_anchor_boundary_0based=3)
        record = alignment("padded", "ACGT", cigar=[(0, 1), (6, 1), (0, 3)])
        record.query_qualities = [35, 19, 35, 35]
        self.source([row], [record])
        result = self.run_bases()
        rows = self.positions()
        self.assertEqual(result["status"], "COMPLETE")
        self.assertEqual([r[b + "_qname_total"] for r, b in zip(rows, "ACGT")], ["1", "0", "1", "1"])
        self.assertEqual([r["LOW_BQ_read_total"] for r in rows], ["0", "1", "0", "0"])

    def test_mixed_clipping_and_indel_skip_actual_bam(self):
        cigar = [(5, 2), (4, 1), (0, 1), (1, 1), (7, 1), (2, 1), (3, 1), (8, 1), (4, 1), (5, 2)]
        self.source([candidate()], [alignment("mixed", "NATCGN", cigar=cigar)])
        self.run_bases()
        rows = self.positions()
        self.assertEqual((rows[0]["A_qname_total"], rows[1]["C_qname_total"], rows[4]["G_qname_total"]),
                         ("1", "1", "1"))
        self.assertEqual((rows[2]["DEL_OR_SKIP_read_total"], rows[3]["DEL_OR_SKIP_read_total"]), ("1", "1"))
        self.assertEqual(rows[5]["zero_acgt_support"], "True")

    def test_positive_full_interval_and_output_hashes(self):
        self.source([candidate()], [alignment("one", "ACGTAC")])
        result = self.run_bases()
        rows = self.positions()
        self.assertEqual(len(rows), 6)
        self.assertEqual(list(rows[0]), POSITION_FIELDS)
        self.assertEqual([int(row["acgt_qname_total"]) for row in rows], [1] * 6)
        self.assertEqual([float(row["assembly_fraction"]) for row in rows], [1.0] * 6)
        self.assertEqual(result["candidates"][0]["zero_acgt_support_positions"], 0)
        self.assertTrue(result["bam_sha_reverified"])
        item = result["files"]["structure_000001.positions.tsv.gz"]
        path = self.out / "structure_000001.positions.tsv.gz"
        self.assertEqual((item["bytes"], item["sha256"]), (path.stat().st_size, digest(path)))
        self.assertEqual((self.out / "COMPLETE").read_text().strip(), digest(self.out / "receipt.json"))

    def test_agreeing_overlapping_mates_count_one_qname_both_strands(self):
        self.source([candidate()], [alignment("pair", "ACG", 65), alignment("pair", "ACG", 145)])
        self.run_bases()
        first, *_, last = self.positions()
        self.assertEqual((first["A_read_F"], first["A_read_R"], first["A_qname_both"], first["A_qname_total"]),
                         ("1", "1", "1", "1"))
        self.assertEqual(last["zero_acgt_support"], "True")
        self.assertEqual(last["assembly_fraction"], "")

    def test_reverse_conflict_and_n_are_not_acgt_consensus(self):
        self.source([candidate(True)], [alignment("conflict", "A", 65), alignment("conflict", "T", 145),
                                        alignment("n", "N"), alignment("reverse", "A", 16)])
        self.run_bases()
        row = self.positions()[0]
        self.assertEqual((row["genomic_position_0based"], row["full_position_0based"]), ("0", "5"))
        self.assertEqual((row["assembled_genomic_base"], row["assembled_oriented_base"]), ("A", "T"))
        self.assertEqual((row["conflict_qname_both"], row["N_qname_total"], row["acgt_qname_total"]), ("1", "1", "1"))
        self.assertEqual((row["A_qname_R"], row["T_qname_total"]), ("1", "0"))

    def test_mapq_bq_primary_and_duplicate_flag_rules(self):
        self.source([candidate()], [alignment("low", "A", mapq=19), alignment("unknown", "A", mapq=255),
                                   alignment("secondary", "A", flag=256), alignment("supplementary", "A", flag=2048),
                                   alignment("lowbq", "A", bq=19), alignment("edge", "A", mapq=20, bq=20),
                                   alignment("duplicate", "A", flag=1024)])
        result = self.run_bases()
        row = self.positions()[0]
        self.assertEqual((row["A_qname_total"], row["LOW_BQ_read_total"]), ("2", "1"))
        self.assertEqual(result["ledger"]["nonprimary_records_excluded"], 2)
        self.assertEqual(result["ledger"]["mapq_records_excluded"], 2)

    def test_deletion_skip_and_insertion_preserve_reference_positions(self):
        self.source([candidate()], [alignment("del", "AG", cigar=[(0, 1), (2, 1), (0, 1)]),
                                   alignment("skip", "AG", cigar=[(0, 1), (3, 1), (0, 1)]),
                                   alignment("ins", "ATC", cigar=[(0, 1), (1, 1), (0, 1)])])
        self.run_bases()
        rows = self.positions()
        self.assertEqual((rows[1]["DEL_OR_SKIP_read_total"], rows[1]["C_qname_total"]), ("2", "1"))
        self.assertEqual(rows[2]["G_qname_total"], "2")

    def test_alternate_orientation_filter_and_tie(self):
        self.source([candidate(True)], [alignment("one", "CA"), alignment("two", "C", flag=16),
                                       alignment("assembly", "AC")])
        self.run_bases()
        rows = self.positions()
        self.assertEqual(rows[0]["alternate_ge2_qnames_both_orientations"], "C")
        self.assertEqual(rows[0]["alternative_exceeds_assembly"], "True")
        self.assertAlmostEqual(float(rows[0]["assembly_fraction"]), 1 / 3)
        self.assertEqual(rows[1]["alternative_exceeds_assembly"], "False")

    def test_receipt_validation_fields_are_required(self):
        source = self.source([candidate()], [alignment("one", "A")])
        for field in ("qname_validation", "boundary_pair_validation", "sqlite_integrity_check", "header_validation"):
            with self.subTest(field=field):
                broken = copy.deepcopy(source)
                broken["readback"].pop(field)
                self.save(broken)
                with self.assertRaisesRegex(ValueError, "validation missing or failed"):
                    self.run_bases()

    def test_bam_sha_is_rechecked(self):
        source = self.source([candidate()], [alignment("one", "A")])
        source["alignment_input"]["sha256"] = "0" * 64
        self.save(source)
        with self.assertRaisesRegex(ValueError, "SHA/size mismatch: alignment_input"):
            self.run_bases()

    def test_noncontiguous_retained_names_fail_even_with_python_optimized(self):
        self.source([candidate()], [alignment("one", "A"), alignment("two", "A")])
        import pysam
        path = self.root / "repeated.bam"
        with pysam.AlignmentFile(str(path), "wb", header={"SQ": [{"SN": "ctg", "LN": 6}]}) as stream:
            for name in ("one", "two", "one"):
                stream.write(alignment(name, "A"))
        with self.assertRaisesRegex(ValueError, "not contiguous"):
            count_bam(path, {("ctg", 0)})
        env = dict(os.environ)
        env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")
        code = "import sys; from glycan_evidence.base_concordance import count_bam; count_bam(sys.argv[1], {('ctg', 0)})"
        result = subprocess.run([sys.executable, "-O", "-c", code, str(path)], env=env,
                                capture_output=True, text=True, timeout=30)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("ValueError", result.stderr)
        self.assertIn("not contiguous", result.stderr)

    def test_zero_candidates_with_alignment_revalidates_sha_without_position_rows(self):
        self.source([], [alignment("one", "A")])
        result = self.run_bases()
        self.assertTrue(result["bam_sha_reverified"])
        self.assertFalse(result["alignment_scan_performed"])
        self.assertIsNone(result["bam_unchanged_during_read"])
        self.assertEqual(result["files"], {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
