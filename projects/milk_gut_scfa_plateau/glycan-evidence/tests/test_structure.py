"""Tiny deterministic geometry tests; these are software tests, not validation."""
import copy
import json
import importlib.util
import os
import random
import shutil
import tempfile
import unittest
from contextlib import nullcontext
from pathlib import Path
from unittest.mock import patch

from glycan_evidence.structure import (BLAST_FIELDS, find_structures, load_config, parse_hit,
                                       read_fasta, reverse_complement, run_structure, sequence_sha, summarize_evidence)


def hit(name, start, strand="-", qstart=1, qend=100, qlen=100, contig="c", pident=100):
    length = qend-qstart+1
    return {"qseqid": name, "sseqid": contig, "qlen": qlen, "qstart": qstart, "qend": qend,
            "lo": start, "hi": start+length, "subject_strand": strand, "query_span_fraction": length/qlen,
            "pident": pident, "qseq": "T"*length, "sseq": "T"*length, "hit_id": f"{name}_{start}_{contig}"}


class StructureTests(unittest.TestCase):
    def setUp(self):
        self.config = load_config()
        self.contigs = {"c": "A"*1000, "d": "C"*1000}
        self.a, self.b = "Blon_2174", "Blon_2178"

    def candidates(self, *hits):
        return find_structures(hits, self.contigs, self.config)

    def test_reference_direction_and_physical_boundaries(self):
        rows, seqs, rejected = self.candidates(hit(self.a, 100), hit(self.b, 300))
        row = rows[0]
        self.assertEqual((row["orientation_to_reference"], row["left_anchor_id"], row["right_anchor_id"]), ("+", self.a, self.b))
        self.assertEqual((row["left_anchor_boundary_0based"], row["right_anchor_boundary_0based"]), (200, 300))
        self.assertEqual(row["between_anchor_bp"], 100)
        self.assertTrue(row["inner_anchor_edges_complete"])
        self.assertEqual(row["interval_sequence_sha256"], sequence_sha(seqs[row["structure_id"]]))
        self.assertEqual(row["contig_sequence_sha256"], sequence_sha(self.contigs["c"]))
        self.assertFalse(rejected)

    def test_reverse_contig_swaps_physical_anchor_labels(self):
        rows, sequences, _ = self.candidates(hit(self.a, 300, "+"), hit(self.b, 100, "+"))
        row = rows[0]
        self.assertEqual((row["orientation_to_reference"], row["left_anchor_id"]), ("-", self.b))
        self.assertEqual(row["between_anchor_bp"], 100)
        self.assertEqual(sequences[row["structure_id"]], "T"*300)

    def test_mixed_reference_coding_strands_are_supported(self):
        self.config["anchors"][0]["reference_strand"] = "+"
        rows, _, _ = self.candidates(hit(self.a, 100, "+"), hit(self.b, 300, "-"))
        self.assertEqual(rows[0]["orientation_to_reference"], "+")
        self.assertTrue(rows[0]["inner_anchor_edges_complete"])

    def test_missing_inner_query_start_not_called_exact_boundary(self):
        rows, _, _ = self.candidates(hit(self.a, 100, qstart=2), hit(self.b, 300))
        row = rows[0]
        self.assertFalse(row["inner_anchor_edges_complete"])
        self.assertEqual(row["between_anchor_bp"], "")
        self.assertEqual(row["between_hsp_edges_bp"], 101)
        self.assertIn("unresolved", row["candidate_evidence"])

    def test_outer_truncation_can_preserve_inner_edge(self):
        rows, _, _ = self.candidates(hit(self.a, 100, qend=95), hit(self.b, 300, qstart=6))
        self.assertTrue(rows[0]["inner_anchor_edges_complete"])
        self.assertFalse(rows[0]["both_anchors_full_length"])

    def test_terminal_subject_gap_not_a_complete_inner_edge(self):
        a = hit(self.a, 100)
        a["sseq"] = "-" + a["sseq"][1:]
        rows, _, _ = self.candidates(a, hit(self.b, 300))
        self.assertFalse(rows[0]["left_anchor_inner_edge_complete"])

    def test_partial_anchor_retained_below_screen(self):
        rows, _, _ = self.candidates(hit(self.a, 100, qend=50), hit(self.b, 300))
        self.assertEqual(len(rows), 1)
        self.assertFalse(rows[0]["anchor_screen_pass"])

    def test_wrong_contig_orientation_and_order_are_rejected(self):
        cases = [(hit(self.a, 100), hit(self.b, 300, contig="d"), "anchors_on_different_contigs"),
                 (hit(self.a, 100), hit(self.b, 300, "+"), "inconsistent_anchor_orientation"),
                 (hit(self.a, 300), hit(self.b, 100), "overlap_or_order_conflict")]
        for a, b, reason in cases:
            with self.subTest(reason=reason):
                rows, _, rejected = self.candidates(a, b)
                self.assertFalse(rows)
                self.assertEqual(rejected[0]["reason"], reason)

    def test_config_distance_limit(self):
        self.config["max_inter_anchor_bp"] = 20
        rows, _, rejected = self.candidates(hit(self.a, 100), hit(self.b, 300))
        self.assertFalse(rows)
        self.assertEqual(rejected[0]["reason"], "distance_exceeds_config")

    def test_all_multiple_candidates_are_retained(self):
        rows, _, _ = self.candidates(hit(self.a, 100), hit(self.b, 300), hit(self.b, 500))
        self.assertEqual(len(rows), 2)
        self.assertEqual(summarize_evidence(rows)["observation"], "unresolved_multiple_candidate_structures")

    def test_empty_result_is_unresolved_not_absent(self):
        rows, _, _ = self.candidates(hit(self.a, 100))
        report = summarize_evidence(rows)
        self.assertEqual(report["observation"], "unresolved_no_geometrically_paired_anchors")
        self.assertNotIn("absent", json.dumps(report))

    def test_invalid_config_rejected(self):
        config = copy.deepcopy(self.config)
        config["anchors"][1]["reference_order"] = 0
        with self.assertRaises(ValueError):
            load_config(config)

    def test_sequence_faithful_blast_parse_and_mismatch_rejection(self):
        raw = dict(zip(BLAST_FIELDS, [self.a, "c", 100, 100, 100, 1000, 1, 100, 200, 101, 0, 200, 0, "T"*100, "T"*100]))
        parsed = parse_hit(raw, self.contigs, {self.a: "T"*100})
        self.assertEqual((parsed["lo"], parsed["hi"], parsed["subject_strand"]), (100, 200, "-"))
        raw["sseq"] = "C"*100
        with self.assertRaisesRegex(ValueError, "inconsistent"):
            parse_hit(raw, self.contigs, {self.a: "T"*100})

    def test_duplicate_fasta_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"x.fa"
            path.write_text(">a\nAAA\n>a\nCCC\n")
            with self.assertRaisesRegex(ValueError, "Duplicate"):
                read_fasta(path)

    def test_failure_receipt_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            c, a, out = tmp/"c.fa", tmp/"a.fa", tmp/"out"
            c.write_text(">c\n" + "A"*1000 + "\n")
            a.write_text(f">{self.a}\n" + "T"*100 + f"\n>{self.b}\n" + "T"*100 + "\n")
            with patch("glycan_evidence.structure.executable", side_effect=FileNotFoundError("fixture_missing_tool")):
                with self.assertRaises(FileNotFoundError):
                    run_structure(c, a, out)
            self.assertFalse((out/"COMPLETE").exists())
            self.assertEqual(json.loads((out/"receipt.json").read_text())["status"], "FAILED")
            with self.assertRaises(FileExistsError):
                run_structure(c, a, out)


TEST_BIN_DIR = os.environ.get("GLYCAN_EVIDENCE_TEST_BIN_DIR")
HAS_EXTERNAL_TOOLS = bool(importlib.util.find_spec("pysam")) and all(
    (Path(TEST_BIN_DIR)/name).is_file() if TEST_BIN_DIR else bool(shutil.which(name))
    for name in ("blastn", "bowtie2", "bowtie2-build", "samtools"))


@unittest.skipUnless(HAS_EXTERNAL_TOOLS, "Linux BLAST/Bowtie2/samtools/pysam integration fixture; not installed locally")
class ExternalToolsIntegrationTests(unittest.TestCase):
    def test_full_blast_mapping_forward_reverse_and_shared_decoy(self):
        rng = random.Random(864117)
        sequence = "".join(rng.choice("ACGT") for _ in range(1000))
        artifact_root = os.environ.get("GLYCAN_EVIDENCE_TEST_ARTIFACT_DIR")
        context = (nullcontext(tempfile.mkdtemp(prefix="external_fixture_", dir=artifact_root)) if artifact_root
                   else tempfile.TemporaryDirectory())
        with context as tmp:
            tmp = Path(tmp)
            anchors = tmp/"anchors.fa"
            anchors.write_text(">Blon_2174\n"+reverse_complement(sequence[100:200])+"\n>Blon_2178\n"+reverse_complement(sequence[800:900])+"\n")
            r1, r2 = tmp/"R1.fq", tmp/"R2.fq"
            with r1.open("w") as one, r2.open("w") as two:
                for n, start in enumerate(range(0, 751, 25)):
                    one.write(f"@p{n}/1\n{sequence[start:start+100]}\n+\n"+"I"*100+"\n")
                    two.write(f"@p{n}/2\n{reverse_complement(sequence[start+150:start+250])}\n+\n"+"I"*100+"\n")
            for case in ("forward", "reverse", "identical_decoy", "unique_k5_unknown_mapq"):
                with self.subTest(case=case):
                    fasta = tmp/f"{case}.fa"
                    seq = reverse_complement(sequence) if case == "reverse" else sequence
                    fasta.write_text(">target\n"+seq+"\n"+(">decoy\n"+seq+"\n" if case == "identical_decoy" else ""))
                    out = tmp/case
                    config = load_config()
                    if case == "unique_k5_unknown_mapq":
                        config["readback"]["reporting_mode"] = "k_multi_audit"
                    result = run_structure(fasta, anchors, out, r1=r1, r2=r2, config=config, bin_dir=TEST_BIN_DIR, threads=1)
                    self.assertEqual(result["status"], "COMPLETE")
                    self.assertTrue((out/"COMPLETE").is_file())
                    self.assertEqual(result["structure_count"], 2 if case == "identical_decoy" else 1)
                    candidates = result["summary"]["candidates"]
                    if case in ("identical_decoy", "unique_k5_unknown_mapq"):
                        self.assertTrue(all(c["evidence_level"] == "partial_read_evidence" for c in candidates))
                        self.assertTrue(all(c["breadth_ge_1"] == 0 for c in candidates))
                        if case == "unique_k5_unknown_mapq":
                            self.assertGreater(result["readback"]["counts"].get("mapped_primary_mapq_unknown255", 0), 0)
                    else:
                        self.assertEqual(candidates[0]["evidence_level"], "read_supported_anchor_interval", result["summary"])
                        self.assertEqual(candidates[0]["breadth_ge_1"], 1)
                        self.assertTrue(all(n > 0 for n in candidates[0]["boundary_unique_qnames"]))

    def test_bad_bam_header_fails_before_evidence(self):
        import pysam
        from glycan_evidence.readback import analyze_readback
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp)/"bad.bam"
            with pysam.AlignmentFile(str(bam), "wb", header={"SQ": [{"SN": "c", "LN": 99}]}):
                pass
            with self.assertRaisesRegex(ValueError, "@SQ"):
                analyze_readback(bam, {"c": "A"*100}, [], load_config()["readback"], tmp)


if __name__ == "__main__":
    unittest.main()
