import json
from pathlib import Path
import tempfile
import unittest

from glycan_evidence.structure_inference import (
    UNEXPLAINED, extract_fragment_features, fit_identifiability,
)


def pair(name="q1", ref="A", start=50, mate_start=200, score=0, secondary=False,
         mapq=0, cigar="100M", seq_length=100, nm=0, tlen=None):
    """Starts are zero-based; equal-length paired fixture."""
    span = mate_start+seq_length-start if tlen is None else tlen
    extra = 256 if secondary else 0
    result = []
    for flag, pos, other, length in ((99+extra, start, mate_start, span),
                                    (147+extra, mate_start, start, -span)):
        result.append("\t".join(map(str, [name, flag, ref, pos+1, mapq, cigar, "=", other+1,
            length, "A"*seq_length, "I"*seq_length, "AS:i:"+str(score),
            "YS:i:"+str(score), "NM:i:"+str(nm)])))
    return result


class FeatureTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.meta = {"A": {"length": 500, "junctions": [100, 175]},
                     "B": {"length": 500, "junctions": [100, 175]}}

    def extract(self, lines, config=None):
        path = Path(self.temp.name)/"alignments.sam"
        header = ["@HD\tVN:1.6\tSO:unsorted"] + [f"@SQ\tSN:{c}\tLN:{m['length']}" for c, m in self.meta.items()]
        path.write_text("\n".join(header+lines)+"\n", encoding="utf-8")
        result = extract_fragment_features(path, self.meta, config)
        json.dumps(result, allow_nan=False)
        return result

    def test_one_fragment_mapq_zero_secondary_and_junctions(self):
        result = self.extract(pair()+pair(ref="B", secondary=True))
        self.assertEqual(sum(result["counts"].values()), 1)
        self.assertEqual(result["audit"]["secondary_records"], 2)
        item = result["fragment_ledger"][0]
        self.assertEqual(item["candidates"], ["A", "B"])
        self.assertEqual(item["candidate_evidence"]["A"]["direct_junctions"], [100])
        self.assertEqual(item["candidate_evidence"]["A"]["bridged_junctions"], [175])

    def test_near_best_is_pair_score(self):
        result = self.extract(pair()+pair(ref="B", score=-4, secondary=True))
        self.assertEqual(result["fragment_ledger"][0]["candidates"], ["A"])

    def test_same_candidate_multimap_does_not_double_count_or_invent_junction(self):
        result = self.extract(pair()+pair(start=200, mate_start=350, secondary=True))
        self.assertEqual(sum(result["counts"].values()), 1)
        evidence = result["fragment_ledger"][0]["candidate_evidence"]["A"]
        self.assertEqual(evidence["placements"], 2)
        self.assertEqual(evidence["direct_junctions"], [])
        self.assertEqual(evidence["bridged_junctions"], [])

    def test_same_coordinates_distinct_cigar_variants_are_retained(self):
        result = self.extract(pair()+pair(cigar="50M1I1D49M", nm=2, secondary=True))
        self.assertEqual(sum(result["counts"].values()), 1)
        evidence = result["fragment_ledger"][0]["candidate_evidence"]["A"]
        self.assertEqual(evidence["placements"], 4)
        self.assertEqual(evidence["direct_junctions"], [])
        self.assertEqual(result["audit"]["same_coordinate_variant_groups"], 1)

    def test_technical_summary_has_one_pair_and_two_mates_per_qname(self):
        result = self.extract(pair()+pair(ref="B", secondary=True))
        summary = result["audit"]["technical_summary"]
        self.assertEqual(summary["absolute_template_lengths"]["count"], 1)
        self.assertEqual(summary["absolute_template_lengths"]["mean"], 250)
        self.assertEqual(summary["mate_query_lengths"]["count"], 2)
        self.assertEqual(summary["mate_query_lengths"]["histogram"], {"100": 2})

    def test_paired_ablation_omits_junctions(self):
        result = self.extract(pair(), {"mode": "paired"})
        self.assertEqual(json.loads(next(iter(result["counts"]))), [["A"]])

    def test_missing_mate_fails(self):
        with self.assertRaisesRegex(ValueError, "reciprocal"):
            self.extract(pair()[:1])

    def test_bad_reciprocal_coordinate_fails(self):
        lines = pair()
        fields = lines[1].split("\t")
        fields[7] = "52"
        lines[1] = "\t".join(fields)
        with self.assertRaisesRegex(ValueError, "reciprocal"):
            self.extract(lines)

    def test_bad_tlen_fails(self):
        with self.assertRaisesRegex(ValueError, "TLEN"):
            self.extract(pair(tlen=999))

    def test_noncontiguous_qname_fails(self):
        with self.assertRaisesRegex(ValueError, "noncontiguous"):
            self.extract(pair()+pair(name="q2")+pair())

    def test_duplicate_primary_fails(self):
        with self.assertRaisesRegex(ValueError, "Duplicate primary"):
            self.extract(pair()+pair())

    def test_cigar_insertion_breaks_direct_boundary(self):
        # At boundary 100, an insertion breaks the continuous matched run.
        result = self.extract(pair(cigar="50M1I49M", tlen=249, nm=1))
        self.assertEqual(result["fragment_ledger"][0]["candidate_evidence"]["A"]["direct_junctions"], [])

    def test_rejected_identity_is_local_residual_not_ledger(self):
        result = self.extract(pair(nm=15))
        self.assertEqual(result["counts"], {UNEXPLAINED: 1})
        self.assertEqual(result["fragment_ledger"], [])
        self.assertEqual(result["audit"]["unexplained_fragments"], 1)

    def test_ledger_cap_never_changes_counts(self):
        result = self.extract(pair()+pair(name="q2"), {"max_ledger": 1})
        self.assertEqual(sum(result["counts"].values()), 2)
        self.assertEqual(len(result["fragment_ledger"]), 1)
        self.assertEqual(result["audit"]["ledger_rows_omitted_by_cap"], 1)

    def test_unmapped_retained_only_as_local_residual(self):
        lines = ["\t".join(map(str, ["u", flag, "*", 0, 0, "*", "*", 0, 0, "A"*100, "I"*100])) for flag in (77, 141)]
        result = self.extract(lines)
        self.assertEqual(result["counts"], {UNEXPLAINED: 1})


class FitTest(unittest.TestCase):
    def fit(self, counts, profiles, config=None):
        result = fit_identifiability(counts, profiles, config)
        json.dumps(result, allow_nan=False)
        return result

    def test_two_component_exact_fit(self):
        result = self.fit({"x": 30, "y": 70}, {"A": {"x": 1}, "B": {"y": 1}}, {"tv_slack": 0})
        self.assertAlmostEqual(result["bestfit_tv"], 0)
        self.assertAlmostEqual(result["bestfit_weights"]["A"], .3)
        self.assertAlmostEqual(result["weight_ranges"]["A"][0], .3, places=7)
        self.assertAlmostEqual(result["weight_ranges"]["A"][1], .3, places=7)
        self.assertIsNone(result["biological_call"])

    def test_duplicate_columns_remain_ambiguous_group_identifiable(self):
        result = self.fit({"x": 8, "y": 2}, {"A": {"x": 1}, "A2": {"x": 1}, "B": {"y": 1}},
                          {"tv_slack": 0, "candidate_groups": {"long": ["A", "A2"]}})
        self.assertEqual(result["duplicate_profile_groups"], [["A", "A2"]])
        self.assertAlmostEqual(result["weight_ranges"]["A"][0], 0)
        self.assertAlmostEqual(result["weight_ranges"]["A"][1], .8, places=7)
        self.assertAlmostEqual(result["group_ranges"]["long"][0], .8, places=7)
        self.assertAlmostEqual(result["convex_hull_distance_tv"]["A"], 0)

    def test_convex_hull_not_only_duplicate_columns(self):
        result = self.fit({"x": 5, "y": 5}, {"A": {"x": 1}, "B": {"y": 1}, "C": {"x": 1, "y": 1}})
        self.assertEqual(result["duplicate_profile_groups"], [])
        self.assertAlmostEqual(result["convex_hull_distance_tv"]["C"], 0)
        self.assertAlmostEqual(result["weight_ranges"]["C"][0], 0)
        self.assertAlmostEqual(result["weight_ranges"]["C"][1], 1)

    def test_unmodeled_feature_not_free_background(self):
        result = self.fit({"x": 9, "novel": 1, UNEXPLAINED: 10}, {"A": {"x": 5, UNEXPLAINED: 5}})
        self.assertAlmostEqual(result["bestfit_tv"], .1)
        self.assertEqual(result["accepted_fragment_count"], 10)
        self.assertAlmostEqual(result["unexplained_local_fraction"], .5)
        self.assertAlmostEqual(result["unmodeled_feature_fraction"], .1)

    def test_zero_counts_no_call(self):
        result = self.fit({}, {"A": {"x": 1}})
        self.assertEqual(result["status"], "no_accepted_evidence")
        self.assertIsNone(result["bestfit_weights"])

    def test_invalid_profiles_fail_closed(self):
        for profiles in ({}, {"A": {}}, {"A": {UNEXPLAINED: 2}}, {"A": {"x": float("nan")}}, {"A": {"x": -1}}):
            with self.subTest(profiles=profiles), self.assertRaises(ValueError):
                self.fit({"x": 1}, profiles)

    def test_order_validation(self):
        with self.assertRaises(ValueError):
            self.fit({"x": 1}, {"A": {"x": 1}}, {"candidate_order": ["B"]})
        with self.assertRaises(ValueError):
            self.fit({"x": 1}, {"A": {"x": 1}}, {"candidate_order": []})

    def test_nonintegral_observed_count_rejected(self):
        with self.assertRaisesRegex(ValueError, "integral"):
            self.fit({"x": .5}, {"A": {"x": 1}})

    def test_diagnostics_can_be_explicitly_disabled(self):
        result = self.fit({"x": 1}, {"A": {"x": 1}}, {"diagnose_profiles": False})
        self.assertEqual(result["convex_hull_diagnostics_status"], "not_computed")
        self.assertEqual(result["convex_hull_distance_tv"], {})
        self.assertAlmostEqual(result["bestfit_weights"]["A"], 1)

    def test_ranges_widen_with_budget(self):
        result = self.fit({"x": 3, "y": 7}, {"A": {"x": 1}, "B": {"y": 1}}, {"tv_slack": .1})
        self.assertAlmostEqual(result["weight_ranges"]["A"][0], .2, places=7)
        self.assertAlmostEqual(result["weight_ranges"]["A"][1], .4, places=7)


if __name__ == "__main__":
    unittest.main()
