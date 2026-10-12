"""Synthetic orchestration tests; no real reads or installed aligners required."""
import csv
import gzip
import importlib.util
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from glycan_evidence.legacy_v2 import run_sample as pipeline
MODULE_PATH = Path(pipeline.__file__)

FAKE_TOOL = r'''
import gzip, json, os, sys
tool, args = sys.argv[1], sys.argv[2:]
if '--version' in args:
    print(tool + ' synthetic-test-1'); sys.exit(0)
def opt(k): return args[args.index(k)+1]
def read(path):
    fn = gzip.open if path.endswith('.gz') else open
    with fn(path, 'rt') as h: return h.read().splitlines()
def fail(stage):
    if os.environ.get('FAKE_FAIL') == stage:
        sys.stderr.write('Intentional toy process failure\n'); sys.exit(7)
if tool == 'fastp':
    fail('fastp')
    a,b=read(opt('-i')),read(opt('-I'))
    for start in range(0,len(a),4): print('\n'.join(a[start:start+4]+b[start:start+4]))
    open(opt('--json'),'w').write('{}')
    open(opt('--html'),'w').write('<html>synthetic</html>')
    sys.exit(0)
screen = opt('-x').endswith('features')
fail('screen' if screen else 'competitive')
if '--interleaved' in args:
    z=sys.stdin.read().splitlines(); a=[]; b=[]
    for start in range(0,len(z),8): a+=z[start:start+4]; b+=z[start+4:start+8]
else: a,b=read(opt('-1')),read(opt('-2'))
print('@HD\tVN:1.6\tSO:queryname')
comp=str.maketrans('ACGT','TGCA')
for start in range(0,len(a),4):
    name=a[start][1:].split('/')[0]
    for mate,z in ((1,a),(2,b)):
        if screen and os.environ.get('FAKE_MISSING_MATE') and mate==2: continue
        seq,qual=z[start+1],z[start+3]
        unmapped=bool(screen and os.environ.get('FAKE_ALL_UNMAPPED'))
        flag=(77 if mate==1 else 141) if unmapped else (99 if mate==1 else 147)
        if flag & 16: seq=seq.translate(comp)[::-1]; qual=qual[::-1]
        fields=[name,str(flag),'*' if unmapped else 'a','0' if unmapped else str(1 if mate==1 else 101),'0' if unmapped else '42','*' if unmapped else str(len(seq))+'M','=','1','0',seq,qual]
        if not unmapped: fields+=['AS:i:0','YS:i:0']
        print('\t'.join(fields))
'''


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.fake = self.root / "fake.py"
        self.fake.write_text(FAKE_TOOL)
        self.panel = self.root / "panel"
        self.panel.mkdir()
        for prefix in ("features", "genomes"):
            for part in pipeline.INDEX_PARTS:
                (self.panel / f"{prefix}.{part}.bt2").write_bytes(b"synthetic index " + part.encode())
        self.table("panel_contigs.tsv", ["contig", "taxon_group", "accession"], [["a", "A", "acc"]])
        self.table("feature_manifest.tsv", ["feature_id", "channel", "taxon_group", "length", "eligible_identity", "pseudo"], [["id", "identity", "A", 200, 1, 0]])
        self.table("annotations.tsv", ["contig", "feature_id", "subject_start1", "subject_end1", "query_start1", "query_end1", "qseq", "sseq", "pident"], [["a", "id", 1, 200, 1, 200, "A" * 200, "A" * 200, 100]])
        self.r1, self.r2 = self.root / "r1.fq", self.root / "r2.fq"
        self.r1.write_text("@toy/1\n" + "ACGA" * 20 + "\n+\n" + "I" * 80 + "\n")
        self.r2.write_text("@toy/2\n" + "TGCA" * 20 + "\n+\n" + "H" * 80 + "\n")

    def tearDown(self):
        self.temp.cleanup()

    def table(self, name, fields, rows):
        with (self.panel / name).open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(fields)
            writer.writerows(rows)

    def execute(self, skip_qc=True, keep=True, env=None):
        with patch.object(pipeline, "command_prefix", side_effect=lambda tool, _: [sys.executable, str(self.fake), tool]), patch.dict(os.environ, env or {}):
            return pipeline.run(self.r1, self.r2, "toy", self.panel, self.root / "out", 2, skip_qc, None, keep)

    def test_reverse_sequence_and_quality_restored_with_unmapped_partner(self):
        first = "f\t83\ta\t1\t42\t4M\t=\t0\t0\tTCGT\tDCBA\n"
        second = "f\t133\t*\t0\t0\t*\t=\t1\t0\tGCAA\tJKLM\n"
        a, b = self.root / "a.gz", self.root / "b.gz"
        counts = pipeline.select_candidates(io.StringIO(first + second), a, b)
        with gzip.open(a, "rt") as handle:
            self.assertEqual(handle.read(), "@f/1\nACGA\n+\nABCD\n")
        with gzip.open(b, "rt") as handle:
            self.assertEqual(handle.read(), "@f/2\nGCAA\n+\nJKLM\n")
        self.assertEqual(counts["one_mate_mapped"], 1)

    def test_missing_mate_and_missing_sequence_fail(self):
        first = "f\t99\ta\t1\t42\t4M\t=\t1\t0\tACGA\tIIII\n"
        with self.assertRaisesRegex(ValueError, "missing a primary mate"):
            pipeline.select_candidates(io.StringIO(first), self.root / "a.gz", self.root / "b.gz")
        second = "f\t133\t*\t0\t0\t*\t=\t1\t0\t*\t*\n"
        with self.assertRaisesRegex(ValueError, "complete sequence/quality"):
            pipeline.select_candidates(io.StringIO(first + second), self.root / "a.gz", self.root / "b.gz")

    def test_skip_qc_end_to_end_and_candidate_retention(self):
        state = self.execute()
        self.assertEqual(state["status"], "COMPLETE")
        self.assertTrue((self.root / "out/COMPLETE").is_file())
        self.assertTrue((self.root / "out/candidate_R1.fastq.gz").is_file())
        summary = json.loads((self.root / "out/summary.json").read_text())
        self.assertEqual(summary["total_mapped_fragments"], 1)
        self.assertEqual(summary["observed_features_by_margin"]["0"]["identity"], 1)
        screen = next(c for c in state["commands"] if c["stage"] == "feature_screen")["argv"]
        self.assertNotIn("--no-mixed", screen)
        self.assertNotIn("--no-unal", screen)
        self.assertIn("--reorder", screen)
        competitive = next(c for c in state["commands"] if c["stage"] == "competitive_genome_mapping")["argv"]
        self.assertIn("-a", competitive)
        self.assertIn("--no-mixed", competitive)
        self.assertNotIn("--reorder", competitive)
        self.assertIn("--no-discordant", competitive)
        self.assertEqual(competitive[competitive.index("--score-min") + 1], "L,0,-0.4")
        self.assertTrue(all(c["returncode"] == 0 for c in state["commands"]))

    def test_real_qc_stream_uses_same_downstream_and_optional_cleanup(self):
        state = self.execute(skip_qc=False, keep=False)
        self.assertEqual(state["screen"]["candidate_fragments"], 1)
        self.assertTrue((self.root / "out/fastp.json").is_file())
        self.assertFalse((self.root / "out/candidate_R1.fastq.gz").exists())
        self.assertTrue((self.root / "out/competitive.sam.gz").is_file())
        screen = next(c for c in state["commands"] if c["stage"] == "feature_screen")["argv"]
        self.assertIn("--interleaved", screen)
        qc = next(c for c in state["commands"] if c["stage"] == "fastp")["argv"]
        self.assertIn("--detect_adapter_for_pe", qc)
        self.assertEqual(qc[qc.index("--length_required") + 1], "75")

    def test_no_candidates_yields_header_sam_and_all_zero_metrics(self):
        state = self.execute(env={"FAKE_ALL_UNMAPPED": "1"})
        self.assertEqual(state["screen"]["candidate_fragments"], 0)
        self.assertFalse(any(c["stage"] == "competitive_genome_mapping" for c in state["commands"]))
        with (self.root / "out/gene_metrics.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(r["fragments"] == "0" for r in rows))

    def test_process_failure_never_writes_complete(self):
        with self.assertRaisesRegex(RuntimeError, "status 7"):
            self.execute(env={"FAKE_FAIL": "competitive"})
        self.assertFalse((self.root / "out/COMPLETE").exists())
        state = json.loads((self.root / "out/run_manifest.json").read_text())
        self.assertEqual(state["status"], "FAILED")
        self.assertEqual(state["commands"][-1]["returncode"], 7)

    def test_screen_failure_and_missing_mate_are_recorded(self):
        with self.assertRaisesRegex(ValueError, "missing a primary mate"):
            self.execute(env={"FAKE_MISSING_MATE": "1"})
        self.assertFalse((self.root / "out/COMPLETE").exists())
        self.assertEqual(json.loads((self.root / "out/run_manifest.json").read_text())["status"], "FAILED")

    def test_fastp_failure_is_not_hidden_by_successful_empty_screen(self):
        with self.assertRaisesRegex(RuntimeError, "fastp.*status 7"):
            self.execute(skip_qc=False, env={"FAKE_FAIL": "fastp"})
        state = json.loads((self.root / "out/run_manifest.json").read_text())
        self.assertEqual(state["status"], "FAILED")
        self.assertFalse((self.root / "out/COMPLETE").exists())
        self.assertEqual(next(c for c in state["commands"] if c["stage"] == "fastp")["returncode"], 7)

    def test_screen_nonzero_is_not_mistaken_for_zero_candidates(self):
        with self.assertRaisesRegex(RuntimeError, "feature_screen.*status 7"):
            self.execute(env={"FAKE_FAIL": "screen"})
        self.assertFalse((self.root / "out/COMPLETE").exists())

    def test_truncated_query_and_hard_clipping_are_rejected(self):
        fields = "f\t99\ta\t1\t42\t5M\t=\t1\t0\tACGA\tIIII".split("\t")
        with self.assertRaisesRegex(ValueError, "full query-consuming CIGAR"):
            pipeline.restore_read(fields)
        fields[5] = "2H4M"
        with self.assertRaisesRegex(ValueError, "Hard-clipped"):
            pipeline.restore_read(fields)

    def build_cache(self):
        entries = []
        for path in sorted(self.panel.glob("*.bt2")):
            record = pipeline.file_record(path, checksum=True)
            entries.append({key: path.name if key == "name" else record[key] for key in ("name", "bytes", "mtime_ns", "sha256")})
        (self.panel / "index_sha256.json").write_text(json.dumps({"files": entries}))
        return entries

    def test_valid_build_cache_avoids_rehashing_indexes(self):
        entries = self.build_cache()
        original = pipeline.file_record
        def checked(path, checksum=False):
            if str(path).endswith(".bt2") and checksum:
                self.fail("Index contents were rehashed despite a valid build cache")
            return original(path, checksum)
        with patch.object(pipeline, "file_record", side_effect=checked):
            rows, receipt = pipeline.index_records(self.panel)
        self.assertEqual(len(rows), 12)
        self.assertEqual(receipt["checksum_source"], "build_cache_exact_set_stat_verified")
        self.assertEqual({r["sha256"] for r in rows}, {r["sha256"] for r in entries})

    def test_stale_or_incomplete_cache_is_a_failure(self):
        entries = self.build_cache()
        path = self.panel / entries[0]["name"]
        path.write_bytes(path.read_bytes() + b"changed")
        with self.assertRaisesRegex(ValueError, "cache is stale"):
            pipeline.index_records(self.panel)
        entries = self.build_cache()
        (self.panel / "index_sha256.json").write_text(json.dumps({"files": entries[:-1]}))
        with self.assertRaisesRegex(ValueError, "exactly 12"):
            pipeline.index_records(self.panel)

    def test_no_cache_fallback_hashes_full_index_set(self):
        rows, receipt = pipeline.index_records(self.panel)
        self.assertEqual(len(rows), 12)
        self.assertTrue(all(len(row["sha256"]) == 64 for row in rows))
        self.assertEqual(receipt["checksum_source"], "direct_sha256_no_build_cache")
        (self.panel / "genomes.rev.2.bt2").unlink()
        with self.assertRaisesRegex(ValueError, "exactly 12"):
            pipeline.index_records(self.panel)


if __name__ == "__main__":
    unittest.main()
