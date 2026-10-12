"""Meaningful synthetic evidence tests; all fixture reads are artificial."""
import csv
import importlib.util
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from glycan_evidence.legacy_v2 import map_evidence as engine
MODULE_PATH = Path(engine.__file__)


def sam_record(name, mate, contig, pos, cigar="60M", score=0, secondary=False, reverse=None):
    reverse = mate == 2 if reverse is None else reverse
    flag = 0x1 | 0x2 | (0x40 if mate == 1 else 0x80) | (0x10 if reverse else 0) | (0x100 if secondary else 0)
    return "\t".join(map(str, (name, flag, contig, pos, 42, cigar, "=", 1, 0, "*", "*", f"AS:f:{score}", "YS:i:0")))


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def write_tsv(self, name, fields, rows):
        path = self.root / name
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(fields)
            writer.writerows(rows)
        return path

    def execute(self, records, features, annotations, contigs=(("a", "A", "accA"), ("b", "B", "accB"))):
        cp = self.write_tsv("contigs.tsv", ("contig", "taxon_group", "accession"), contigs)
        fp = self.write_tsv("features.tsv", ("feature_id", "channel", "taxon_group", "length", "eligible_identity", "pseudo"), features)
        ap = self.write_tsv("annotations.tsv", ("contig", "feature_id", "subject_start1", "subject_end1", "query_start1", "query_end1", "qseq", "sseq", "pident"), annotations)
        sam = self.root / "input.sam"
        sam.write_text("@HD\tVN:1.6\tSO:queryname\n" + "\n".join(records) + "\n", encoding="utf-8")
        summary = engine.run(sam, cp, fp, ap, "toy", self.root / "out")
        with (self.root / "out/gene_metrics.tsv").open() as handle:
            rows = {(r["feature_id"], int(r["margin"])): r for r in csv.DictReader(handle, delimiter="\t")}
        return rows, summary

    @staticmethod
    def annotation(feature, start=1, length=300, contig="a"):
        return contig, feature, start, start + length - 1, 1, length, "A" * length, "A" * length, 100

    def test_secondary_and_overlapping_mates_count_once(self):
        records = [sam_record("f", 1, "a", 1), sam_record("f", 2, "a", 31),
                   sam_record("f", 1, "a", 1, secondary=True), sam_record("f", 2, "a", 31, secondary=True)]
        rows, summary = self.execute(records, [("id", "identity", "A", 300, 1, 0), ("unseen", "locus", "", 100, 0, 0)], [self.annotation("id")])
        self.assertEqual(summary["total_mapped_fragments"], 1)
        self.assertEqual(summary["duplicate_position_records_removed"], 2)
        for margin in engine.MARGINS:
            self.assertEqual(rows["id", margin]["fragments"], "1")
            self.assertEqual(rows["id", margin]["covered_bases"], "90")
            self.assertEqual(rows["id", margin]["aligned_bases"], "90")
            self.assertEqual(rows["unseen", margin]["fragments"], "0")

    def test_identity_requires_50_actual_marker_bases(self):
        records = [sam_record("short", 1, "a", 1, "49M"), sam_record("short", 2, "a", 201, "49M"),
                   sam_record("outside", 1, "a", 501), sam_record("outside", 2, "a", 601)]
        rows, _ = self.execute(records, [("id", "identity", "A", 300, 1, 0), ("locus", "locus", "", 300, 0, 0)], [self.annotation("id"), self.annotation("locus")])
        self.assertEqual(rows["id", 0]["fragments"], "0")
        self.assertEqual(rows["locus", 0]["fragments"], "0")
        self.assertEqual(rows["locus", 0]["covered_bases"], "0")
        self.assertEqual(rows["locus", 0]["ambiguous_fragments"], "1")

    def test_fragment_blocks_can_arrive_out_of_input_order_with_all_alternatives(self):
        blocks = {}
        for name, contig in (("a_first", "a"), ("z_last", "a"), ("m_middle", "b")):
            blocks[name] = [sam_record(name, mate, contig, pos, secondary=pos == 501)
                            for pos in (1, 501) for mate in (1, 2)]
        forward = [record for name in ("a_first", "m_middle", "z_last") for record in blocks[name]]
        completed_order = [record for name in ("z_last", "a_first", "m_middle") for record in reversed(blocks[name])]
        grouped = list(engine.sam_groups(io.StringIO("\n".join(completed_order) + "\n")))
        self.assertEqual([name for name, _ in grouped], ["z_last", "a_first", "m_middle"])
        self.assertEqual([len(records) for _, records in grouped], [4, 4, 4])
        features = [("idA", "identity", "A", 1000, 1, 0), ("idB", "identity", "B", 1000, 1, 0)]
        annotations = [self.annotation("idA", length=1000), self.annotation("idB", length=1000, contig="b")]
        expected, summary1 = self.execute(forward, features, annotations)
        actual, summary2 = self.execute(completed_order, features, annotations)
        self.assertEqual(expected, actual)
        self.assertEqual(summary1, summary2)
        self.assertEqual(summary2["total_mapped_fragments"], 3)
        self.assertEqual(actual["idA", 0]["fragments"], "2")
        self.assertEqual(actual["idB", 0]["fragments"], "1")

    def test_cross_group_identity_ambiguous_locus_independent(self):
        records = [sam_record("f", 1, c, 1, secondary=c == "b") for c in ("a", "b")]
        records += [sam_record("f", 2, c, 101, secondary=c == "b") for c in ("a", "b")]
        rows, _ = self.execute(records, [("id", "identity", "A", 300, 1, 0), ("locus", "locus", "A", 300, 0, 0)],
                               [self.annotation("id"), self.annotation("id", contig="b"), self.annotation("locus"), self.annotation("locus", contig="b")])
        self.assertEqual(rows["id", 0]["fragments"], "0")
        self.assertEqual(rows["id", 0]["ambiguous_fragments"], "1")
        self.assertEqual(rows["locus", 0]["fragments"], "1")

    def test_near_best_other_family_and_unknown_are_ambiguous(self):
        records = []
        for name, altpos in (("other", 501), ("unknown", 901)):
            for mate in (1, 2):
                records.extend([sam_record(name, mate, "a", 1, score=0), sam_record(name, mate, "a", altpos, score=-6, secondary=True)])
        rows, _ = self.execute(records, [("A", "locus", "", 200, 0, 0), ("B", "locus", "", 200, 0, 0)], [self.annotation("A", length=200), self.annotation("B", start=501, length=200)])
        self.assertEqual(rows["A", 0]["fragments"], "2")
        self.assertEqual(rows["A", 6]["fragments"], "0")
        self.assertEqual(rows["A", 6]["ambiguous_fragments"], "2")
        self.assertEqual(rows["B", 6]["ambiguous_fragments"], "1")

    def test_two_adjacent_families_and_separate_mates(self):
        records = [sam_record("cross", 1, "a", 51, "100M"), sam_record("cross", 2, "a", 501, "20M"),
                   sam_record("separate", 1, "a", 1, "60M"), sam_record("separate", 2, "a", 121, "60M")]
        rows, _ = self.execute(records, [("A", "locus", "", 100, 0, 0), ("B", "locus", "", 100, 0, 0)], [self.annotation("A", length=100), self.annotation("B", start=101, length=100)])
        self.assertEqual(rows["A", 0]["fragments"], "2")
        self.assertEqual(rows["B", 0]["fragments"], "2")
        self.assertEqual(rows["A", 0]["aligned_bases"], "110")
        self.assertEqual(rows["B", 0]["aligned_bases"], "110")

    def test_reverse_blast_and_gapped_projection(self):
        # Subject 106..101 maps query 1,2,3,4,5,6; q/s gaps each exclude
        # a base. CIGAR positions are reference-forward despite reverse flag.
        annotations = [("a", "rev", 106, 101, 1, 6, "AC-GTAC", "ACG-TAC", 90)]
        features = [("rev", "locus", "", 6, 0, 0)]
        records = [sam_record("rev", 1, "a", 101, "6M", reverse=True), sam_record("rev", 2, "a", 501)]
        rows, _ = self.execute(records, features, annotations)
        self.assertEqual(rows["rev", 0]["covered_bases"], "0")
        self.assertEqual(rows["rev", 0]["ambiguous_fragments"], "1")
        _, _, bins, _ = engine.load_panel(self.root / "contigs.tsv", self.root / "features.tsv", self.root / "annotations.tsv")
        self.assertEqual(engine.project_alignment(engine.Alignment(1, "a", 101, "2M", True, 0), bins)["rev"], {5, 6})

    def test_cigar_does_not_fill_deletions_skips_or_mate_gap(self):
        cigar = "5S30M3I25=4D15X20N30M2H"
        positions = engine.sequenced_reference_positions(1, cigar)
        self.assertEqual(len(positions), 100)
        self.assertFalse(set(range(56, 60)) & set(positions))
        self.assertFalse(set(range(75, 95)) & set(positions))
        rows, _ = self.execute([sam_record("f", 1, "a", 1, cigar), sam_record("f", 2, "a", 201, "50M")],
                               [("gene", "locus", "", 300, 0, 0)], [self.annotation("gene")])
        self.assertEqual(rows["gene", 0]["covered_bases"], "150")

    def test_ties_and_deterministic_projection_do_not_union_locations(self):
        records = [sam_record("f", mate, "a", pos, score=-0.1, secondary=pos == 501) for mate in (1, 2) for pos in (1, 501)]
        annotations = [self.annotation("g", length=100), ("a", "g", 501, 600, 101, 200, "A" * 100, "A" * 100, 100)]
        rows, _ = self.execute(records, [("g", "locus", "", 200, 0, 0)], annotations)
        self.assertEqual(rows["g", 0]["fragments"], "1")
        self.assertEqual(rows["g", 0]["covered_bases"], "60")
        # Changing secondary record order cannot change the selected position.
        rows2, _ = self.execute(list(reversed(records)), [("g", "locus", "", 200, 0, 0)], annotations)
        self.assertEqual(rows, rows2)

    def test_pseudo_and_ineligible_identity_are_never_supporting(self):
        features = [("pseudo", "identity", "A", 300, 1, 1), ("bad", "identity", "A", 300, 0, 0)]
        rows, summary = self.execute([sam_record("f", 1, "a", 1), sam_record("f", 2, "a", 101)], features,
                                     [self.annotation("pseudo"), self.annotation("bad")])
        self.assertEqual(summary["eligible_identity_markers"], 0)
        self.assertTrue(all(r["fragments"] == "0" for r in rows.values()))

    def test_paired_score_uses_ys_and_margin_boundaries(self):
        features = [("A", "locus", "", 200, 0, 0), ("B", "locus", "", 200, 0, 0)]
        annotations = [self.annotation("A", length=200), self.annotation("B", start=501, length=200)]
        records = []
        for mate in (1, 2):
            records.append(sam_record("ys", mate, "a", 1, score=-1).replace("YS:i:0", "YS:i:-11"))
            records.append(sam_record("ys", mate, "a", 501, score=-7, secondary=True))
        rows, _ = self.execute(records, features, annotations)
        self.assertEqual(rows["A", 0]["fragments"], "0")
        self.assertEqual(rows["B", 0]["fragments"], "1")
        self.assertEqual(rows["B", 6]["ambiguous_fragments"], "1")
        for delta in (6, 12, 18):
            records = [sam_record("boundary", mate, "a", pos, score=score, secondary=pos == 501)
                       for mate in (1, 2) for pos, score in ((1, 0), (501, -delta))]
            rows, _ = self.execute(records, features, annotations)
            for margin in engine.MARGINS:
                self.assertEqual(int(rows["A", margin]["fragments"]), int(margin < delta))
                self.assertEqual(int(rows["A", margin]["ambiguous_fragments"]), int(margin >= delta))

    def test_reverse_query_and_conflicting_hsps_do_not_duplicate_bases(self):
        features = [("g", "locus", "", 200, 0, 0)]
        annotations = [("a", "g", 1, 100, 100, 1, "A" * 100, "A" * 100, 100),
                       ("a", "g", 1, 100, 200, 101, "A" * 100, "A" * 100, 99)]
        rows, _ = self.execute([sam_record("f", 1, "a", 1, "60M"), sam_record("f", 2, "a", 501)], features, annotations)
        self.assertEqual(rows["g", 0]["covered_bases"], "60")
        _, _, bins, _ = engine.load_panel(self.root / "contigs.tsv", self.root / "features.tsv", self.root / "annotations.tsv")
        self.assertEqual(engine.project_alignment(engine.Alignment(1, "a", 1, "10M", False, 0), bins)["g"], set(range(91, 101)))

    def test_same_genomic_bases_cannot_resolve_two_locus_homologues(self):
        features = [("A", "locus", "", 300, 0, 0), ("B", "locus", "", 300, 0, 0)]
        rows, summary = self.execute([sam_record("f", 1, "a", 1), sam_record("f", 2, "a", 101)], features,
                                     [self.annotation("A"), self.annotation("B")])
        for margin in engine.MARGINS:
            for feature in ("A", "B"):
                self.assertEqual(rows[feature, margin]["fragments"], "0")
                self.assertEqual(rows[feature, margin]["covered_bases"], "0")
                self.assertEqual(rows[feature, margin]["ambiguous_fragments"], "1")
                self.assertEqual(summary["locus_shared_ambiguity"][str(margin)][feature]["candidate_fragments"], 1)

    def test_only_exclusive_bases_contribute_locus_depth(self):
        features = [("A", "locus", "", 300, 0, 0), ("B", "locus", "", 300, 0, 0)]
        rows, summary = self.execute([sam_record("f", 1, "a", 1, "100M"), sam_record("f", 2, "a", 501)], features,
                                     [self.annotation("A"), self.annotation("B", start=51)])
        self.assertEqual(rows["A", 0]["fragments"], "1")
        self.assertEqual(rows["A", 0]["aligned_bases"], "50")
        self.assertEqual(rows["B", 0]["fragments"], "0")
        self.assertEqual(rows["B", 0]["ambiguous_fragments"], "1")
        self.assertEqual(summary["locus_shared_ambiguity"]["0"]["A"]["candidate_projected_bases"], 50)

    def test_missing_pair_score_fails_explicitly(self):
        record = sam_record("f", 1, "a", 1).replace("\tYS:i:0", "")
        with self.assertRaisesRegex(ValueError, "AS and YS"):
            self.execute([record, sam_record("f", 2, "a", 101)], [("g", "locus", "", 300, 0, 0)], [self.annotation("g")])

    def test_cli_and_unmapped_only_zero_output(self):
        unmapped = "\t".join(map(str, ("none", 77, "*", 0, 0, "*", "*", 0, 0, "*", "*")))
        rows, summary = self.execute([unmapped], [("g", "locus", "", 300, 0, 0)], [self.annotation("g")])
        self.assertEqual(summary["total_fragments"], 1)
        self.assertEqual(summary["total_mapped_fragments"], 0)
        result = subprocess.run([sys.executable, str(MODULE_PATH), "--sam", str(self.root / "input.sam"),
                                 "--contigs", str(self.root / "contigs.tsv"), "--features", str(self.root / "features.tsv"),
                                 "--annotations", str(self.root / "annotations.tsv"), "--sample", "cli", "--out", str(self.root / "cli")], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertFalse(json.loads((self.root / "cli/summary.json").read_text())["biology_calls_made"])
        self.assertEqual(len(rows), 4)


if __name__ == "__main__":
    unittest.main()
