import json
import unittest
from unittest.mock import patch

from glycan_evidence.structure_assignment import fit_assignment, UNEXPLAINED


def feature(*entries):
    return json.dumps(list(entries), separators=(",", ":"))


UA = feature(["a", [], []])
UB = feature(["b", [], []])
SHARED = feature(["a", [100, 175], []], ["b", [100, 175], []])
META = {"a": dict(length=500, junctions=[100, 175], class_id="A"),
        "b": dict(length=500, junctions=[100, 175], class_id="B")}


class AssignmentTests(unittest.TestCase):
    def fit(self, counts, profiles, config=None, metadata=None):
        result = fit_assignment(counts, profiles, metadata or META, config)
        json.dumps(result, allow_nan=False)
        if result["status"] == "fitted_descriptive":
            self.assertLessEqual(result["marginal_tv_at_joint_bestfit"], result["joint_bestfit_tv"]+1e-8)
            self.assertNotIn("marginal_bestfit_tv", result)
            self.assertLessEqual(result["solver_validation"]["max_equality_residual"], result["solver_validation"]["tolerance"])
            self.assertAlmostEqual(sum(result["bestfit_weights"].values()), 1., places=7)
        self.assertIsNone(result["biological_call"])
        return result

    def test_hand_computed_shared_boundary_positive_lower(self):
        profiles = {"a": {UA: .8, SHARED: .2}, "b": {UB: .8, SHARED: .2}}
        result = self.fit({UA: 80, SHARED: 20}, profiles)
        self.assertAlmostEqual(result["joint_bestfit_tv"], 0., places=7)
        self.assertAlmostEqual(result["group_ranges"]["B"][1], .0125, places=7)
        for side in ("left", "right"):
            bounds = result["boundary_bounds"]["A"][side]
            self.assertAlmostEqual(bounds["lower_mass"], .1875, places=7)
            self.assertAlmostEqual(bounds["lower_fractional_fragments"], 18.75, places=5)
            self.assertEqual(result["observed_support"]["A"][side]["possible_qnames"], 20)
            self.assertEqual(result["observed_support"]["A"][side]["exclusive_qnames"], 0)
            self.assertEqual(result["same_mask_no_profile"]["boundary_bounds"]["A"][side]["lower_mass"], 0)

    def test_fixed_budget_erases_rare_boundary_despite_large_class_lower(self):
        profiles = {"a": {UA: .99, SHARED: .01}, "b": {UB: .99, SHARED: .01}}
        result = self.fit({UA: 990, SHARED: 10}, profiles)
        self.assertAlmostEqual(result["group_ranges"]["A"][0], .99, places=7)
        self.assertAlmostEqual(result["boundary_bounds"]["A"]["right"]["lower_mass"], 0., places=7)
        self.assertEqual(result["observed_support"]["A"]["right"]["possible_qnames"], 10)

    def test_zero_slack_identifies_shared_boundary_in_exact_fixture(self):
        profiles = {"a": {UA: .8, SHARED: .2}, "b": {UB: .8, SHARED: .2}}
        result = self.fit({UA: 80, SHARED: 20}, profiles, {"tv_slack": 0})
        self.assertAlmostEqual(result["boundary_bounds"]["A"]["left"]["lower_mass"], .2, places=7)

    def test_identical_profiles_keep_allocation_ambiguous(self):
        result = self.fit({SHARED: 100}, {"a": {SHARED: 1}, "b": {SHARED: 1}})
        self.assertEqual(result["group_ranges"]["A"], [0., 1.])
        self.assertEqual(result["group_ranges"]["B"], [0., 1.])
        self.assertEqual(result["boundary_bounds"]["A"]["left"]["lower_mass"], 0.)

    def test_mixed_source_hand_fit(self):
        result = self.fit({UA: 24, UB: 56, SHARED: 20},
                          {"a": {UA: .8, SHARED: .2}, "b": {UB: .8, SHARED: .2}}, {"tv_slack": 0})
        self.assertAlmostEqual(result["joint_bestfit_tv"], 0., places=7)
        self.assertAlmostEqual(result["bestfit_weights"]["a"], .3, places=7)
        self.assertAlmostEqual(result["boundary_bounds"]["A"]["left"]["lower_mass"], .06, places=7)

    def test_profile_only_edges_preserved_without_boundary_evidence(self):
        observed_b = feature(["b", [100, 175], []])
        result = self.fit({observed_b: 100}, {"a": {observed_b: 1}, "b": {observed_b: 1}})
        self.assertEqual(result["profile_only_probability_by_candidate"]["a"], 1.)
        self.assertEqual(len(result["profile_only_edges"]), 1)
        self.assertEqual(result["observed_support"]["B"]["left"]["exclusive_qnames"], 100)
        self.assertEqual(result["group_ranges"]["A"], [0., 1.])
        self.assertEqual(result["boundary_bounds"]["A"]["left"]["lower_mass"], 0.)
        self.assertEqual(result["same_mask_no_profile"]["group_ranges"]["B"], [0., 1.])
        self.assertEqual(result["observed_mask_no_profile"]["group_ranges"]["B"], [1., 1.])

    def test_unmodeled_accepted_feature_keeps_mass_and_charges_discrepancy(self):
        new = feature(["a", [100], []])
        result = self.fit({new: 10}, {"a": {UA: 1}, "b": {UB: 1}})
        self.assertAlmostEqual(result["joint_bestfit_tv"], 1.)
        self.assertEqual(result["accepted_fragment_count"], 10)
        self.assertAlmostEqual(result["group_ranges"]["A"][0], 1.)

    def test_literal_unexplained_has_separate_denominator(self):
        result = self.fit({UA: 20, UNEXPLAINED: 80}, {"a": {UA: 1}, "b": {UB: 1}})
        self.assertEqual(result["accepted_fragment_count"], 20)
        self.assertEqual(result["unexplained_local_fragment_count"], 80)
        self.assertEqual(result["unexplained_local_fraction"], .8)

    def test_no_accepted_evidence_is_explicit(self):
        result = self.fit({UNEXPLAINED: 10}, {"a": {UA: 1}, "b": {UB: 1}})
        self.assertEqual(result["status"], "no_accepted_evidence")
        self.assertIsNone(result["joint_bestfit_tv"])
        self.assertIsNone(result["boundary_bounds"])

    def test_unknown_candidates_and_empty_feature_fail_closed(self):
        for bad in (feature(["unknown", [], []]), '[]', '{}', 'not_json'):
            with self.subTest(feature=bad), self.assertRaises(ValueError):
                self.fit({bad: 1}, {"a": {UA: 1}, "b": {UB: 1}})

    def test_invalid_profile_feature_fails_even_without_observed_count(self):
        with self.assertRaises(ValueError):
            self.fit({UA: 1}, {"a": {UA: 1}, "b": {feature(["unknown", [], []]): 1}})

    def test_invalid_junction_not_a_supported_boundary(self):
        with self.assertRaises(ValueError):
            self.fit({feature(["a", [101], []]): 1}, {"a": {UA: 1}, "b": {UB: 1}})

    def test_grouped_variants_do_not_require_unique_candidate(self):
        result = self.fit({SHARED: 30}, {"a": {SHARED: 1}, "b": {SHARED: 1}},
                          {"candidate_groups": {"combined": ["a", "b"]}})
        self.assertEqual(result["group_ranges"]["combined"], [1., 1.])
        self.assertAlmostEqual(result["boundary_bounds"]["combined"]["right"]["lower_fractional_fragments"], 30.)
        self.assertEqual(result["observed_support"]["combined"]["right"]["exclusive_qnames"], 30)

    def test_bridged_evidence_not_relabelled_direct(self):
        bridged = feature(["a", [], [100, 175]])
        result = self.fit({bridged: 30}, {"a": {bridged: 1}, "b": {UB: 1}})
        self.assertAlmostEqual(result["boundary_bounds"]["A"]["left"]["lower_fractional_fragments"], 30.)
        self.assertEqual(result["boundary_bounds"]["A"]["left"]["direct_lower_mass"], 0.)
        self.assertEqual(result["observed_support"]["A"]["left"]["bridged_possible_qnames"], 30)
        self.assertEqual(result["observed_support"]["A"]["left"]["direct_possible_qnames"], 0)

    def test_invalid_counts_profiles_or_order(self):
        for counts in ({UA: .5}, {UA: -1}, {UA: float('nan')}, {UA: True}):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                self.fit(counts, {"a": {UA: 1}, "b": {UB: 1}})
        with self.assertRaises(ValueError):
            self.fit({UA: 1}, {"a": {}, "b": {UB: 1}})
        with self.assertRaises(ValueError):
            self.fit({UA: 1}, {"a": {UA: 1}, "b": {UB: 1}}, {"candidate_order": []})

    def test_invalid_solver_solution_fails_residual_validation(self):
        from scipy.optimize import linprog
        def damaged(*args, **kwargs):
            solution = linprog(*args, **kwargs)
            solution.x[0] += .01
            return solution
        with patch("scipy.optimize.linprog", side_effect=damaged), self.assertRaisesRegex(RuntimeError, "residual"):
            self.fit({UA: 80, SHARED: 20}, {"a": {UA: .8, SHARED: .2}, "b": {UB: .8, SHARED: .2}})


if __name__ == '__main__':
    unittest.main()
