import copy
import importlib.util
import json
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location("eval_structure", Path(__file__).with_name("evaluate_benchmark.py"))
E = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(E)


def metadata():
    return {
        "bt1": {"length": 8441, "junctions": [2256, 7154], "inner_start": 2256, "inner_end": 7154, "class_id": "ignored_alias"},
        "other_bt1": {"length": 8440, "junctions": [2256, 7153], "inner_start": 2256, "inner_end": 7153},
        "short": {"length": 4367, "junctions": [2256, 3080], "inner_start": 2256, "inner_end": 3080},
    }


def fragment(name, candidates):
    return {"qname": name, "candidates": list(candidates), "candidate_evidence": candidates}


def support(left=True, right=True, last=7154, direct=True):
    points = ([2256] if left else []) + ([last] if right else [])
    return {"direct_junctions": points if direct else [], "bridged_junctions": [] if direct else points}


def extraction(rows):
    meta = metadata()
    return {"config": {"ledger_mode": "accepted"}, "audit": {"accepted_fragments": len(rows), "ledger_rows_omitted_by_cap": 0},
            "fragment_ledger": rows, "candidate_coverage": {
                c: {"length": v["length"], "eligible_base_breadth": .8, "near_best_base_breadth": .5,
                    "eligible_fragments": 20, "near_best_fragments": 10} for c, v in meta.items()}}


def fitted(n, lower=.6, tv=.1):
    return {"status": "fitted_descriptive", "accepted_fragment_count": n, "bestfit_tv": tv,
            "config": {"candidate_groups": {"alt": ["bt1", "other_bt1"], "short": ["short"]}},
            "group_ranges": {"alt": [lower, .9], "short": [0, .2]},
            "bestfit_weights": {"bt1": .8, "other_bt1": 0, "short": .2}}


class EvaluationTests(unittest.TestCase):
    def test_frozen_length_classes_keep_open_set_separate(self):
        self.assertEqual(E.architecture_class(4898), "alternative_4897_4899")
        self.assertEqual(E.architecture_class(4895), "inner_length_4895")
        self.assertEqual(E.architecture_class(4632), "inner_length_4632")
        self.assertEqual(E.architecture_class(4648), "inner_length_4648")
        self.assertNotEqual(E.architecture_class(4632), E.architecture_class(4633))

    def test_distinct_qnames_and_all_candidate_support_required(self):
        rows = [fragment("origin_should_never_be_parsed__1", {"bt1": support(), "other_bt1": support(right=False, last=7153)}),
                fragment("q2", {"bt1": support(), "short": support(last=3080)}),
                fragment("q3", {"bt1": support(direct=False)})]
        result = E.class_evidence(extraction(rows), metadata())["alternative_4897_4899"]
        self.assertEqual(result["left_robust_qnames"], 2)
        self.assertEqual(result["right_robust_qnames"], 1)
        self.assertEqual(result["right_any_qnames"], 3)
        self.assertEqual(result["both_robust_qnames"], 1)
        self.assertEqual(result["right_robust_bridge_qnames"], 1)

    def test_duplicate_or_truncated_ledger_rejected(self):
        row = fragment("q1", {"bt1": support()})
        with self.assertRaises(ValueError):
            E.class_evidence(extraction([row, row]), metadata())
        x = extraction([row]); x["audit"]["ledger_rows_omitted_by_cap"] = 1
        with self.assertRaises(ValueError):
            E.class_evidence(x, metadata())

    def test_baseline_uses_member_max_not_union_and_all_cutoffs(self):
        x = extraction([])
        x["candidate_coverage"]["bt1"]["eligible_base_breadth"] = .49
        x["candidate_coverage"]["bt1"]["near_best_base_breadth"] = .4
        x["candidate_coverage"]["other_bt1"]["eligible_base_breadth"] = .49
        x["candidate_coverage"]["other_bt1"]["near_best_base_breadth"] = .4
        evidence = E.class_evidence(x, metadata())
        baselines = E.reference_decisions(evidence, E.RECOMMENDED_GATE, 0)
        self.assertEqual(len(baselines), 12)
        self.assertNotIn("alternative_4897_4899", baselines["expanded_reference_eligible_0.50"]["supported_classes"])
        self.assertIn("short_824", baselines["expanded_reference_eligible_0.80"]["supported_classes"])
        self.assertNotIn("short_824", baselines["expanded_reference_eligible_0.95"]["supported_classes"])

    def test_matched_support_requires_same_depth_and_both_boundaries(self):
        x = extraction([fragment(str(i), {"bt1": support()}) for i in range(20)])
        evidence = E.class_evidence(x, metadata())
        calls = E.reference_decisions(evidence, E.RECOMMENDED_GATE, 20)
        self.assertEqual(calls["matched_support_eligible_0.80"]["supported_classes"], ["alternative_4897_4899"])
        self.assertEqual(calls["matched_support_near_best_0.80"]["supported_classes"], [])
        self.assertEqual(E.reference_decisions(evidence, E.RECOMMENDED_GATE, 19)["matched_support_eligible_0.80"]["supported_classes"], [])
        weak = copy.deepcopy(evidence)
        weak["alternative_4897_4899"]["right_robust_qnames"] = 2
        self.assertEqual(E.reference_decisions(weak, E.RECOMMENDED_GATE, 20)["matched_support_eligible_0.80"]["supported_classes"], [])

    def test_matched_support_does_not_use_model_tv_or_weight(self):
        x = extraction([fragment(str(i), {"bt1": support()}) for i in range(20)])
        evidence = E.class_evidence(x, metadata())
        model = E.model_decision(fitted(20, lower=.001, tv=.9), evidence, x, E.RECOMMENDED_GATE)
        self.assertFalse(model["supported_classes"])
        matched = E.reference_decisions(evidence, E.RECOMMENDED_GATE, 20)["matched_support_eligible_0.80"]
        self.assertEqual(matched["supported_classes"], ["alternative_4897_4899"])
        self.assertIs(matched["tv_weight_filters"], False)

    def test_gate_boundaries_and_lower_bound_are_fixed(self):
        x = extraction([fragment(str(i), {"bt1": support()}) for i in range(20)])
        evidence = E.class_evidence(x, metadata())
        ok = E.model_decision(fitted(20, lower=.020001, tv=.2), evidence, x, E.RECOMMENDED_GATE)
        self.assertEqual(ok["supported_classes"], ["alternative_4897_4899"])
        exact = E.model_decision(fitted(20, lower=.02), evidence, x, E.RECOMMENDED_GATE)
        self.assertFalse(exact["supported_classes"])
        bad = E.model_decision(fitted(20, tv=.20001), evidence, x, E.RECOMMENDED_GATE)
        self.assertFalse(bad["supported_classes"])
        self.assertIn("poor_or_unavailable_model_fit", bad["global_failures"])

    def test_low_depth_cannot_be_saved_by_large_weight(self):
        x = extraction([fragment(str(i), {"bt1": support()}) for i in range(19)])
        evidence = E.class_evidence(x, metadata())
        decision = E.model_decision(fitted(19), evidence, x, E.RECOMMENDED_GATE)
        self.assertFalse(decision["supported_classes"])
        self.assertIn("alternative_4897_4899", decision["ambiguous_classes"])

    def test_no_accepted_reads_is_no_call_not_absence(self):
        x = extraction([])
        fit = {"status": "no_accepted_evidence", "accepted_fragment_count": 0,
               "bestfit_tv": None, "bestfit_weights": None}
        decision = E.model_decision(fit, E.class_evidence(x, metadata()), x, E.RECOMMENDED_GATE)
        self.assertEqual(decision["supported_classes"], [])
        score = E.score_sets([], [])
        self.assertTrue(score["empty_on_negative_no_absence_claim"])
        self.assertFalse(score["positive_truth_fully_recovered_no_extra"])

    def test_false_additional_and_missed_truth_remain_separate(self):
        score = E.score_sets(["short_824", "inner_length_4648"], ["alternative_4897_4899", "short_824"])
        self.assertEqual(score["recovered_classes"], ["short_824"])
        self.assertEqual(score["additional_false_classes"], ["inner_length_4648"])
        self.assertEqual(score["missed_classes"], ["alternative_4897_4899"])
        self.assertFalse(score["positive_truth_fully_recovered_no_extra"])

    def test_reverse_oriented_inner_coordinates(self):
        sequence = "AAAACCCCGGGGTTT"
        row = {"interval_start_0based": "100", "interval_end_0based_exclusive": "115",
               "left_anchor_boundary_0based": "103", "right_anchor_boundary_0based": "111",
               "orientation_to_reference": "-", "interval_sequence_sha256": E.sha(sequence), "between_anchor_bp": "8"}
        self.assertEqual(E.oriented_inner(row, sequence), sequence[4:12])
        row["interval_sequence_sha256"] = "bad"
        with self.assertRaises(ValueError):
            E.oriented_inner(row, sequence)

    def test_fit_class_must_be_joint_range_not_sum_of_candidate_bounds(self):
        x = extraction([fragment(str(i), {"bt1": support()}) for i in range(20)])
        fit = fitted(20)
        fit["config"]["candidate_groups"]["alt"] = ["bt1"]
        with self.assertRaises(ValueError):
            E.model_decision(fit, E.class_evidence(x, metadata()), x, E.RECOMMENDED_GATE)

    def test_independent_crosscheck_runtime_and_corrupt_baseline(self):
        path = Path(__file__).with_name("analysis_support.py")
        spec = importlib.util.spec_from_file_location("runtime_analysis_support", path)
        runtime = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runtime)
        meta = metadata()
        for value in meta.values():
            value["class_id"] = E.architecture_class(value["inner_end"]-value["inner_start"])
        x = extraction([fragment(str(i), {"bt1": support()}) for i in range(20)])
        fit = fitted(20)
        fit["config"]["candidate_groups"] = {"alternative_4897_4899": ["bt1", "other_bt1"], "short_824": ["short"]}
        fit["group_ranges"] = {"alternative_4897_4899": [.6, .9], "short_824": [0, .2]}
        decision, baseline = runtime.summarize(x, fit, meta, E.RECOMMENDED_GATE)
        result = {"decision": decision, "paired_only_decision": decision, "reference_baseline": baseline}
        evidence = E.class_evidence(x, meta)
        independent = E.model_decision(fit, evidence, x, E.RECOMMENDED_GATE)
        methods = {"paired_junction_model": independent, "junction_profile_ablation": independent,
                   **E.reference_decisions(evidence, E.RECOMMENDED_GATE, 20)}
        self.assertEqual(E.crosscheck_runtime(result, meta, methods)["status"], "PASS")
        result["reference_baseline"]["matched_support"]["eligible_base_breadth"]["0.8"] = []
        with self.assertRaisesRegex(ValueError, "cross-check mismatch"):
            E.crosscheck_runtime(result, meta, methods)


if __name__ == "__main__":
    unittest.main()
