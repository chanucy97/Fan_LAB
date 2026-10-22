import copy
import json
import unittest

import evaluate_assignment as E


def metadata():
    return {"a": {"length": 5000, "inner_start": 100, "inner_end": 4409, "junctions": [100, 4409], "class_id": "canonical_4309"},
            "b": {"length": 2000, "inner_start": 100, "inner_end": 924, "junctions": [100, 924], "class_id": "short_824"}}


def feature(*rows):
    return json.dumps(list(rows), separators=(",", ":"))


UA = feature(["a", [100], []])
UB = feature(["b", [100], []])
SHARED = feature(["a", [4409], []], ["b", [924], []])


def core():
    import glycan_evidence.structure_assignment as module
    return module


class AssignmentEvaluatorTests(unittest.TestCase):
    def test_closed_form_observed_mask_bounds_and_support(self):
        result = E.rebuild_masks({UA: 80, SHARED: 20}, {"a": {UA: 80, SHARED: 20}, "b": {UB: 80, SHARED: 20}}, metadata())
        group = "canonical_4309"
        self.assertEqual(result["accepted_fragments"], 100)
        self.assertEqual(result["observed_support"][group]["informative_qnames"], 80)
        self.assertEqual(result["observed_support"][group]["right"]["possible_qnames"], 20)
        self.assertEqual(result["observed_support"][group]["right"]["exclusive_qnames"], 0)
        nf = result["no_frequency"]["observed_mask_no_profile"]
        self.assertEqual(nf["group_ranges"][group], [.8, 1.])
        self.assertEqual(nf["boundary_bounds"][group]["right"]["lower_fractional_fragments"], 0)
        self.assertEqual(nf["boundary_bounds"][group]["right"]["upper_fractional_fragments"], 20)

    def test_profile_only_edges_change_union_but_create_no_boundary_evidence(self):
        only_b_right = feature(["b", [924], []])
        result = E.rebuild_masks({UA: 80, only_b_right: 20}, {"a": {UA: 80, only_b_right: 20}, "b": {only_b_right: 100}}, metadata())
        self.assertAlmostEqual(result["profile_only_probability_by_candidate"]["a"], .2)
        observed = result["no_frequency"]["observed_mask_no_profile"]
        union = result["no_frequency"]["same_mask_no_profile"]
        self.assertEqual(observed["group_ranges"]["short_824"], [.2, .2])
        self.assertEqual(union["group_ranges"]["short_824"], [0., .2])
        self.assertEqual(union["boundary_bounds"]["canonical_4309"]["right"]["upper_fractional_fragments"], 0)

    def test_hand_solved_shared_boundary_fixture_checks_core(self):
        counts = {UA: 80, SHARED: 20}
        profiles = {"a": {UA: 80, SHARED: 20}, "b": {UB: 80, SHARED: 20}}
        result = E.rebuild_masks(counts, profiles, metadata())
        fit = core().fit_assignment(counts, profiles, metadata(), {"tv_slack": .02, "numerical_tolerance": 0})
        normalized, aliases, checks = E.verify_fit(fit, result)
        self.assertGreater(checks, 30)
        self.assertAlmostEqual(normalized["boundary_bounds"]["canonical_4309"]["right"]["lower_fractional_fragments"], 18.75, places=6)
        calls = E.decision(normalized, result, E.GATE)
        self.assertEqual(calls["supported_classes"], ["canonical_4309"])
        self.assertEqual(E.decision(result["no_frequency"]["same_mask_no_profile"], result, E.GATE, False)["supported_classes"], [])

    def test_high_group_weight_does_not_rescue_rare_shared_boundary(self):
        counts = {UA: 990, SHARED: 10}
        profiles = {"a": {UA: 99, SHARED: 1}, "b": {UB: 99, SHARED: 1}}
        result = E.rebuild_masks(counts, profiles, metadata())
        fit = core().fit_assignment(counts, profiles, metadata(), {"tv_slack": .02})
        normalized, _, _ = E.verify_fit(fit, result)
        self.assertGreater(normalized["group_ranges"]["canonical_4309"][0], .98)
        self.assertAlmostEqual(normalized["boundary_bounds"]["canonical_4309"]["right"]["lower_fractional_fragments"], 0, places=6)
        self.assertEqual(E.decision(normalized, result, E.GATE)["supported_classes"], [])

    def test_possible_support_baseline_does_not_require_exclusive_boundaries(self):
        result = E.rebuild_masks({UA: 80, SHARED: 20}, {"a": {UA: 80, SHARED: 20}, "b": {UB: 80, SHARED: 20}}, metadata())
        evidence = {g: {"eligible_base_breadth_max": .9, "near_best_base_breadth_max": .7} for g in result["groups"]}
        methods = E.possible_baselines(evidence, result, E.GATE)
        self.assertEqual(len(methods), 6)
        self.assertEqual(methods["possible_boundary_eligible_0.80"]["supported_classes"], ["canonical_4309"])
        self.assertEqual(methods["possible_boundary_near_best_0.80"]["supported_classes"], [])
        self.assertNotIn("short_824", methods["possible_boundary_eligible_0.50"]["supported_classes"])

    def test_fractional_lower_bound_not_rounded_to_qname_count(self):
        result = E.rebuild_masks({UA: 80, SHARED: 20}, {"a": {UA: 80, SHARED: 20}, "b": {UB: 80, SHARED: 20}}, metadata())
        bounds = copy.deepcopy(result["no_frequency"]["same_mask_no_profile"])
        bounds["joint_tv"] = .1
        bounds["boundary_bounds"]["canonical_4309"]["right"]["lower_fractional_fragments"] = 2.999999
        self.assertEqual(E.decision(bounds, result, E.GATE)["supported_classes"], [])
        bounds["boundary_bounds"]["canonical_4309"]["right"]["lower_fractional_fragments"] = 3.
        self.assertEqual(E.decision(bounds, result, E.GATE)["supported_classes"], ["canonical_4309"])

    def test_unknown_candidate_or_nonintegral_counts_fail(self):
        profiles = {"a": {UA: 1}, "b": {UB: 1}}
        with self.assertRaises(ValueError):
            E.rebuild_masks({UA: 1.5}, profiles, metadata())
        with self.assertRaises(ValueError):
            E.rebuild_masks({feature(["unknown"]): 1}, profiles, metadata())

    def test_no_evidence_remains_no_call(self):
        profiles = {"a": {UA: 1}, "b": {UB: 1}}
        result = E.rebuild_masks({"__unexplained__": 4}, profiles, metadata())
        fit = core().fit_assignment({"__unexplained__": 4}, profiles, metadata())
        normalized, _, _ = E.verify_fit(fit, result)
        self.assertEqual(E.decision(normalized, result, E.GATE)["supported_classes"], [])
        self.assertEqual(E.decision(None, result, E.GATE, False)["supported_classes"], [])


if __name__ == "__main__":
    unittest.main()
