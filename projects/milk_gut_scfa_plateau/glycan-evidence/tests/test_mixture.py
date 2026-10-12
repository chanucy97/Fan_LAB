"""Mathematical and provenance fixtures, not biological benchmark evidence."""
from dataclasses import replace
import hashlib
import json
import math
import unittest

from glycan_evidence.mixture import (QueryMembership, component_from_observations,
                                    fit_mixture)
from glycan_evidence.profile import MarkerPanel, Observation, Source


def observation(counts, name, *, positive=False, panel=None, cluster=None):
    panel = panel or MarkerPanel("target", tuple((f"m{i}", 1000) for i in range(len(counts))), "a" * 64)
    source = Source(name, cluster or name, ("target",) if positive else ())
    return Observation(panel, source, name, tuple(counts), tuple(int(n > 0) for n in counts),
                       tuple(counts), hashlib.sha256(name.encode()).hexdigest())


class MixtureTests(unittest.TestCase):
    def fit(self, counts, b=(8, 2), t=(2, 8), **kwargs):
        background = component_from_observations([observation(b, "B")], role="background")
        target = component_from_observations([observation(t, "T", positive=True)], role="target")
        query = observation(counts, "Q", positive=True)
        membership = QueryMembership((query.source,), (), "d" * 64)
        return fit_mixture(query, background, target, membership, **kwargs)

    def test_binomial_closed_form_and_residual_conservation(self):
        result = self.fit((62, 38))
        self.assertAlmostEqual(result["target_marker_support_weight"], .3, places=10)
        self.assertAlmostEqual(result["h1_js_divergence_nats"], 0, places=12)
        self.assertAlmostEqual(sum(m["count_residual"] for m in result["markers"]), 0, places=10)
        self.assertAlmostEqual(result["target_composite_support_component"], 30, places=9)
        self.assertGreater(result["h1_minus_h0"], 0)
        self.assertIsNone(result["identity_call"])
        self.assertIsNone(result["p_value"])

    def test_nonnegative_boundary_optima(self):
        self.assertEqual(self.fit((95, 5))["target_marker_support_weight"], 0)
        self.assertEqual(self.fit((5, 95))["target_marker_support_weight"], 1)

    def test_disjoint_profiles_have_exact_fraction_not_finite_null_gain(self):
        result = self.fit((3, 7), b=(1, 0), t=(0, 1))
        self.assertAlmostEqual(result["target_marker_support_weight"], .7, places=10)
        self.assertIsNone(result["h0_log_composite_likelihood"])
        self.assertIsNone(result["h1_minus_h0"])
        self.assertIn("h0_impossible", result["fit_improvement_status"])
        json.dumps(result, allow_nan=False)

    def test_global_objective_exceeds_dense_grid(self):
        result = self.fit((13, 7, 19, 2), b=(7, 2, 1, 3), t=(1, 6, 8, 3))
        b, t = [x / 13 for x in (7, 2, 1, 3)], [x / 18 for x in (1, 6, 8, 3)]
        for step in range(1001):
            w = step / 1000
            value = sum(n * math.log((1-w)*p+w*q) for n, p, q in zip((13, 7, 19, 2), b, t))
            self.assertGreaterEqual(result["h1_log_composite_likelihood"] + 1e-10, value)

    def test_scale_invariance_and_count_weighted_objective(self):
        small, large = self.fit((62, 38)), self.fit((6200, 3800))
        self.assertAlmostEqual(small["target_marker_support_weight"], large["target_marker_support_weight"])
        self.assertAlmostEqual(large["h1_minus_h0"], 100 * small["h1_minus_h0"])

    def test_all_zero_retained_without_weight(self):
        result = self.fit((0, 0))
        self.assertEqual(result["status"], "no_support")
        self.assertIsNone(result["target_marker_support_weight"])
        self.assertIsNone(result["h1_js_divergence_nats"])
        self.assertEqual(len(result["markers"]), 2)

    def test_union_support_gap_does_not_fit_supported_subset(self):
        result = self.fit((4, 6, 1), b=(8, 2, 0), t=(2, 8, 0))
        self.assertEqual(result["status"], "reference_support_gap")
        self.assertEqual(result["union_reference_gap_count"], 1)
        self.assertIsNone(result["target_marker_support_weight"])
        self.assertIsNone(result["h1_log_composite_likelihood"])
        self.assertTrue(result["markers"][2]["union_reference_gap"])

    def test_identical_profiles_no_arbitrary_weight(self):
        result = self.fit((3, 7), b=(1, 1), t=(2, 2))
        self.assertEqual(result["status"], "unidentifiable_identical_profiles")
        self.assertIsNone(result["target_marker_support_weight"])
        self.assertEqual(result["h1_minus_h0"], 0)
        self.assertEqual(result["weight_identification_range"], [0., 1.])

    def test_near_collinearity_guard_is_recorded(self):
        result = self.fit((3, 7), b=(500000000, 500000000), t=(500000001, 499999999))
        self.assertEqual(result["status"], "unidentifiable_near_collinear_profiles")
        self.assertEqual(result["near_collinear_tv_guard"], 1e-8)
        self.assertIsNone(result["h1_log_composite_likelihood"])
        for bad in (True, -1, float("nan"), .01):
            with self.assertRaises(ValueError):
                self.fit((3, 7), near_collinear_tv=bad)

    def test_nonidentifiable_observed_subset(self):
        result = self.fit((50, 0, 0), b=(2, 4, 4), t=(2, 7, 1))
        self.assertEqual(result["status"], "unidentifiable_observed_support")
        self.assertIsNone(result["target_marker_support_weight"])
        self.assertEqual(result["h1_minus_h0"], 0)
        self.assertTrue(all(m["fitted_probability"] is None for m in result["markers"]))

    def test_source_weighting_and_component_validation(self):
        a = observation((90, 10), "B1")
        b = observation((1, 9), "B2")
        model = component_from_observations([a, b], role="background")
        self.assertEqual(model.probabilities, (.5, .5))
        for inputs, role in (([], "background"), ([a, a], "background"),
                             ([observation((0, 0), "B")], "background"), ([a], "target")):
            with self.assertRaises(ValueError):
                component_from_observations(inputs, role=role)
        with self.assertRaises(ValueError):
            replace(model, probabilities=(.8, .2))

    def test_target_source_and_cluster_holdout_strict(self):
        b = component_from_observations([observation((8, 2), "B")], role="background")
        t = component_from_observations([observation((2, 8), "T", positive=True)], role="target")
        q = observation((5, 5), "Q", positive=True)
        for source in (t.training_sources[0], Source("renamed_T", "T", ("target",))):
            with self.assertRaises(ValueError):
                fit_mixture(q, b, t, QueryMembership((q.source, source), (), "d" * 64))
        for query in (replace(q, sample="T"), replace(q, input_sha256=t.training_input_sha256[0])):
            with self.assertRaises(ValueError):
                fit_mixture(query, b, t, QueryMembership((q.source,), (), "d" * 64))

    def test_background_seen_and_read_reuse_explicit(self):
        original = observation((8, 2), "B")
        b = component_from_observations([original], role="background")
        t = component_from_observations([observation((2, 8), "T", positive=True)], role="target")
        membership = QueryMembership((original.source,), ("B",), "d" * 64)
        result = fit_mixture(original, b, t, membership)
        self.assertTrue(result["background_source_seen"])
        self.assertTrue(result["background_training_reads_reused"])
        self.assertEqual(result["target_marker_support_weight"], 0)
        with self.assertRaises(ValueError):
            fit_mixture(original, b, t, replace(membership, reused_training_samples=()))
        with self.assertRaises(ValueError):
            fit_mixture(original, b, t, replace(membership, reused_training_samples=("T",)))

    def test_mapping_mismatch_and_missing_membership_rejected(self):
        b = component_from_observations([observation((8, 2), "B")], role="background")
        t = component_from_observations([observation((2, 8), "T", positive=True)], role="target")
        q = observation((5, 5), "Q")
        with self.assertRaises(ValueError):
            fit_mixture(q, b, t, None)
        with self.assertRaises(ValueError):
            fit_mixture(replace(q, panel=replace(q.panel, margin=5)), b, t,
                        QueryMembership((q.source,), (), "d" * 64))

    def test_reuse_cannot_name_a_different_background_training_sample(self):
        b1, b2 = observation((8, 2), "B1"), observation((7, 3), "B2")
        background = component_from_observations([b1, b2], role="background")
        target = component_from_observations([observation((2, 8), "T", positive=True)], role="target")
        for query in (b1, replace(b1, sample="renamed_B1")):
            with self.assertRaises(ValueError):
                fit_mixture(query, background, target,
                            QueryMembership((b1.source,), ("B2",), "d" * 64))
        result = fit_mixture(b1, background, target,
                             QueryMembership((b1.source,), ("B1",), "d" * 64))
        self.assertTrue(result["background_training_reads_reused"])


if __name__ == "__main__":
    unittest.main()
