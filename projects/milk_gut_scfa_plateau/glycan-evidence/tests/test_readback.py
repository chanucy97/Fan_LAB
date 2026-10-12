"""Strict SAM/read-evidence fixtures. No fixture establishes field accuracy."""
import copy
import os
import tempfile
import unittest
from pathlib import Path

from glycan_evidence.readback import _blocks, analyze_readback, analyze_records, validate_header, validate_read_pairs
from glycan_evidence.structure import load_config, summarize_evidence


def record(name="pair1", start=100, length=200, flag=99, mapq=40, contig="c", mate_start=100):
    return {"qname": name, "flag": flag, "contig": contig, "start": start, "mapq": mapq,
            "cigar": [(0, length)], "mate_contig": contig, "mate_start": mate_start,
            "template_length": length, "seq": "A"*length, "qualities": bytes([40]*length), "nm": 0}


class ReadbackTests(unittest.TestCase):
    def setUp(self):
        self.contigs = {"c": "A"*500, "decoy": "A"*500}
        self.header = {"HD": {"SO": "unsorted"}, "SQ": [{"SN": k, "LN": len(v)} for k, v in self.contigs.items()]}
        self.config = load_config()["readback"]
        self.structure = {"structure_id": "s", "contig": "c", "interval_start_0based": 100,
                          "interval_end_0based_exclusive": 300, "left_anchor_boundary_0based": 150,
                          "right_anchor_boundary_0based": 250, "left_anchor_inner_edge_complete": True,
                          "right_anchor_inner_edge_complete": True, "inner_anchor_edges_complete": True,
                          "anchor_screen_pass": True,
                          "left_anchor_id": "A", "right_anchor_id": "B", "candidate_evidence": "complete_inner_edges"}

    def analyze(self, records, structures=None):
        # Complete the explicit selected read with a low-MAPQ mate outside the
        # measured interval, unless the fixture already supplies both mates.
        completed = []
        for index, read in enumerate(records):
            read = copy.deepcopy(read)
            other_present = any(r["qname"] == read["qname"] and bool(r["flag"] & 0x40) != bool(read["flag"] & 0x40)
                                and not r["flag"] & (0x100 | 0x800) for r in records)
            if read["flag"] == 99 and not other_present:
                read["mate_start"] = 350
                completed.append(read)
                completed.append(record(read["qname"], start=350, length=100, flag=147, mapq=0, mate_start=read["start"]))
            else:
                completed.append(read)
        return analyze_records(iter(completed), self.header, self.contigs,
                               [self.structure] if structures is None else structures, self.config)

    def quality_core(self, result):
        return next(r for r in result["coverage"] if r["segment"] == "between_anchors" and r["filter"] == "quality_filtered")

    def test_overlapping_mates_are_two_read_depth_one_qname(self):
        result = self.analyze([record(), record(flag=147)])
        self.assertEqual(self.quality_core(result)["min_read_depth"], 2)
        for boundary in result["boundaries"]:
            self.assertEqual(boundary["unique_qnames"], 1)
            self.assertEqual(boundary["forward_primary_mate_alignments"], 1)
            self.assertEqual(boundary["reverse_primary_mate_alignments"], 1)
        self.assertEqual(result["receipt"]["counts"]["alignment_records"], 2)

    def test_identical_decoy_low_mapq_has_no_quality_support(self):
        result = self.analyze([record(mapq=0), record(mapq=0, contig="decoy", flag=99 | 0x100)])
        self.assertEqual(self.quality_core(result)["breadth_ge_1"], 0)
        self.assertTrue(all(r["unique_qnames"] == 0 for r in result["boundaries"]))
        all_primary = next(r for r in result["coverage"] if r["segment"] == "between_anchors" and r["filter"].startswith("all_primary"))
        self.assertEqual(all_primary["min_read_depth"], 1)

    def test_unknown_mapq_255_not_treated_as_high_confidence(self):
        result = self.analyze([record(mapq=255)])
        self.assertEqual(self.quality_core(result)["min_read_depth"], 0)

    def test_each_flank_base_quality_is_checked(self):
        read = record()
        qualities = bytearray(read["qualities"])
        qualities[49] = 0  # Reference 149, within the left boundary window.
        read["qualities"] = qualities
        result = self.analyze([read])
        left, right = result["boundaries"]
        self.assertEqual((left["unique_qnames"], right["unique_qnames"]), (0, 1))
        self.assertEqual(next(r for r in result["coverage"] if r["segment"] == "full_anchored_interval" and r["filter"] == "quality_filtered")["min_read_depth"], 0)

    def test_insertion_breaks_continuous_boundary_support(self):
        read = record(length=201)
        read["cigar"] = [(0, 50), (1, 1), (0, 150)]
        result = self.analyze([read])
        self.assertEqual([r["unique_qnames"] for r in result["boundaries"]], [0, 1])
        self.assertEqual(self.quality_core(result)["min_read_depth"], 1)

    def test_equivalent_aligned_cigar_encodings_preserve_boundary_evidence(self):
        original = self.analyze([record()])
        for cigar in ([(7, 50), (8, 1), (7, 149)], [(7, 50), (7, 150)],
                      [(0, 25), (7, 25), (8, 1), (0, 149)]):
            read = record()
            read["cigar"] = cigar
            with self.subTest(cigar=cigar):
                result = self.analyze([read])
                self.assertEqual(result["coverage"], original["coverage"])
                self.assertEqual(result["boundaries"], original["boundaries"])
                self.assertEqual(result["windows"], original["windows"])

    def test_all_non_aligned_cigar_operators_break_boundary_blocks(self):
        # Coordinate accounting must also preserve query/reference offsets.
        for op, expected in ((1, (150, 200, 51)), (2, (151, 201, 50)),
                             (3, (151, 201, 50)), (4, (150, 200, 51)),
                             (5, (150, 200, 50)), (6, (150, 200, 50))):
            read = record(length=100)
            read["cigar"] = [(7, 50), (op, 1), (8, 50)]
            with self.subTest(op=op):
                self.assertEqual(list(_blocks(read)), [(100, 150, 0), expected])

    def test_deletion_is_neither_depth_nor_boundary_bridge(self):
        read = record(length=199)
        read["cigar"] = [(0, 50), (2, 1), (0, 149)]
        result = self.analyze([read])
        self.assertEqual(result["boundaries"][0]["unique_qnames"], 0)
        self.assertEqual(self.quality_core(result)["min_read_depth"], 0)

    def test_improper_pair_has_depth_but_no_boundary_fragment(self):
        result = self.analyze([record(flag=65)])
        self.assertEqual(self.quality_core(result)["min_read_depth"], 1)
        self.assertEqual(result["boundaries"][0]["unique_qnames"], 0)

    def test_low_depth_and_empty_alignment_remain_unresolved(self):
        for reads in ([], [record(start=100, length=50)]):
            result = self.analyze(reads)
            report = summarize_evidence([self.structure], result)
            self.assertIn("unresolved", report["candidates"][0]["readback_observation"])
            self.assertIn("unresolved", report["structural_classification"])

    def test_no_candidates_still_validates_every_qname(self):
        with self.assertRaisesRegex(ValueError, "non-contiguous"):
            self.analyze([record("x"), record("y"), record("x")], structures=[])

    def test_duplicate_primary_mate_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            self.analyze([record(), record()])

    def test_secondary_and_supplementary_are_not_double_counted(self):
        result = self.analyze([record(), record(flag=99 | 0x100), record(flag=99 | 0x800)])
        self.assertEqual(self.quality_core(result)["min_read_depth"], 1)
        self.assertEqual(result["boundaries"][0]["unique_qnames"], 1)

    def test_boundary_flag_without_actual_mate_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "lacks its primary mate"):
            analyze_records(iter([record()]), self.header, self.contigs, [self.structure], self.config)

    def test_inconsistent_actual_mate_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "inconsistent primary mate"):
            self.analyze([record(), record(flag=147, mate_start=101)])

    def test_missing_or_wrong_decoy_header_is_rejected(self):
        for header in ({"SQ": [{"SN": "c", "LN": 500}]},
                       {"SQ": [{"SN": "c", "LN": 499}, {"SN": "decoy", "LN": 500}]}):
            with self.assertRaisesRegex(ValueError, "@SQ"):
                validate_header(header, self.contigs)

    def test_bad_cigar_sequence_and_missing_qualities_are_rejected(self):
        for key, value in (("seq", "A"), ("qualities", None), ("start", 499)):
            read = record()
            read[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.analyze([read])

    def test_sam_entrypoint_and_bad_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"toy.sam"
            line = "pair1\t99\tc\t101\t40\t200M\t=\t101\t200\t" + "A"*200 + "\t" + "I"*200 + "\tNM:i:0\n"
            line += "pair1\t147\tc\t101\t0\t200M\t=\t101\t-200\t" + "A"*200 + "\t" + "I"*200 + "\tNM:i:0\n"
            path.write_text("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:c\tLN:500\n@SQ\tSN:decoy\tLN:500\n"+line)
            result = analyze_readback(path, self.contigs, [self.structure], self.config, tmp)
            self.assertEqual(self.quality_core(result)["min_read_depth"], 1)
            path.write_text("@SQ\tSN:c\tLN:499\n"+line)
            with self.assertRaisesRegex(ValueError, "@SQ"):
                analyze_readback(path, self.contigs, [self.structure], self.config, tmp)

    def test_fastq_pair_duplicates_truncation_and_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            r1, r2 = Path(tmp)/"r1.fq", Path(tmp)/"r2.fq"
            r1.write_text("@x/1\nACGT\n+\nIIII\n")
            r2.write_text("@x/2\nACGT\n+\nIIII\n")
            self.assertEqual(validate_read_pairs(r1, r2, tmp)["paired_fragments"], 1)
            r1.write_text(r1.read_text()*2)
            r2.write_text(r2.read_text()*2)
            with self.assertRaisesRegex(ValueError, "duplicate FASTQ"):
                validate_read_pairs(r1, r2, tmp)
            r1.write_text("@x/1\nACGT\n+\n")
            with self.assertRaisesRegex(ValueError, "truncated"):
                validate_read_pairs(r1, r2, tmp)

    def test_fastq_same_file_and_hardlink_alias_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            r1, alias = Path(tmp)/"reads.fq", Path(tmp)/"alias.fq"
            r1.write_text("@x\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(ValueError, "distinct files"):
                validate_read_pairs(r1, r1, tmp)
            try:
                os.link(r1, alias)
            except OSError as error:
                self.skipTest(f"Hard links unavailable: {error}")
            with self.assertRaisesRegex(ValueError, "distinct files"):
                validate_read_pairs(r1, alias, tmp)

    def test_fastq_casava_mate_labels_are_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            r1, r2 = Path(tmp)/"r1.fq", Path(tmp)/"r2.fq"
            r1.write_text("@x 1:N:0:A\nACGT\n+\nIIII\n")
            r2.write_text("@x 2:N:0:A\nACGT\n+\nIIII\n")
            self.assertEqual(validate_read_pairs(r1, r2, tmp)["paired_fragments"], 1)
            r2.write_text("@x 1:N:0:A\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(ValueError, "CASAVA mate"):
                validate_read_pairs(r1, r2, tmp)
            r1.write_text("@x/1 2:N:0:A\nACGT\n+\nIIII\n")
            r2.write_text("@x/2 2:N:0:A\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(ValueError, "CASAVA mate"):
                validate_read_pairs(r1, r2, tmp)


if __name__ == "__main__":
    unittest.main()
