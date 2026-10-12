"""Scientific and failure-boundary tests for descriptive marker profiles."""
from dataclasses import replace
import hashlib
import math
import unittest

from glycan_evidence.profile import (MarkerPanel, Source, Observation,
                                     observation_from_rows, fit_expected, diagnose)


class ProfileTests(unittest.TestCase):
    def setUp(self):
        self.panel = MarkerPanel("longum", (("a", 100), ("b", 100)), "a" * 64)
        self.train = Source("training", "cluster_a", ("longum",))
        self.test = Source("held_out", "cluster_b", ("longum",))

    def observation(self, counts=(10, 10), source=None, sample="train", panel=None):
        counts = tuple(counts)
        source = source or self.train
        digest = hashlib.sha256(repr((counts, source, sample, panel)).encode()).hexdigest()
        return Observation(panel or self.panel, source or self.train, sample, counts,
                           tuple(int(n > 0) for n in counts), counts, digest)

    def rows(self):
        return [{"feature_id": "b", "length": "100", "fragments": "0",
                 "covered_bases": "0", "aligned_bases": "0"},
                {"feature_id": "a", "length": "100", "fragments": "1",
                 "covered_bases": "50", "aligned_bases": "50"}]

    def parse(self, rows):
        return observation_from_rows(self.panel, self.train, "sample", rows,
                                     input_sha256="c" * 64)

    def test_complete_grid_is_sorted_and_zero_preserved(self):
        observation = self.parse(self.rows())
        self.assertEqual(observation.counts, (1, 0))
        self.assertEqual(observation.covered_bases, (50, 0))

    def test_missing_zero_marker_is_not_observed_zero(self):
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            self.parse(self.rows()[1:])

    def test_duplicate_unknown_or_changed_marker_rejected(self):
        for rows in (self.rows() + self.rows()[:1],
                     [{**self.rows()[0], "feature_id": "other"}, self.rows()[1]],
                     [{**self.rows()[0], "length": 101}, self.rows()[1]]):
            with self.assertRaises(ValueError):
                self.parse(rows)

    def test_strict_numeric_schema(self):
        for value in (-1, 1.2, True, "1.0", "NaN", "", None):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.parse([{**self.rows()[0], "fragments": value}, self.rows()[1]])
        with self.assertRaisesRegex(ValueError, "schema"):
            self.parse([{**self.rows()[0], "extra": 0}, self.rows()[1]])

    def test_support_and_coverage_consistency(self):
        for kwargs in ({"covered_bases": (101, 1)}, {"counts": (0, 1)},
                       {"aligned_bases": (0, 10)}, {"aligned_bases": (1, 1)},
                       {"aligned_bases": (1001, 10)}):
            with self.assertRaises(ValueError):
                replace(self.observation(), **kwargs)

    def test_empty_panel_duplicate_panel_and_wrong_hash(self):
        for markers in ((), (("a", 100), ("a", 100))):
            with self.assertRaises(ValueError):
                MarkerPanel("longum", markers, "a" * 64)
        with self.assertRaises(ValueError):
            replace(self.panel, mapping_contract_sha256="unknown")

    def test_training_needs_exposed_positive_source(self):
        for source in (replace(self.train, positive_targets=()),
                       replace(self.train, exposure="prospective_unseen"),
                       Source("unlabelled", "c", (), "unlabelled_application")):
            with self.assertRaisesRegex(ValueError, "curated target-positive"):
                fit_expected([self.observation(source=source)])

    def test_empty_zero_and_duplicate_training_rejected(self):
        for observations in ([], [self.observation((0, 0))],
                             [self.observation(), self.observation()]):
            with self.assertRaises(ValueError):
                fit_expected(observations)

    def test_each_source_has_equal_weight_despite_technical_depth(self):
        expected = fit_expected([self.observation((90, 10)),
                                 self.observation((90, 10), sample="technical_repeat"),
                                 self.observation((1, 9), self.test, "other_source")])
        self.assertEqual(expected.probabilities, (.5, .5))
        self.assertEqual(len(expected.training_sources), 2)

    def test_alias_cluster_is_also_excluded(self):
        expected = fit_expected([self.observation()])
        for source in (self.train, replace(self.test, cluster_id=self.train.cluster_id)):
            with self.assertRaisesRegex(ValueError, "overlap"):
                diagnose(self.observation(source=source, sample="new_reads_same_source"), expected)

    def test_training_cannot_double_weight_a_nearclone_cluster(self):
        alias = replace(self.test, cluster_id=self.train.cluster_id)
        with self.assertRaisesRegex(ValueError, "share a cluster"):
            fit_expected([self.observation(), self.observation(source=alias, sample="alias")])

    def test_renaming_source_does_not_hide_identical_input(self):
        train = self.observation()
        expected = fit_expected([train])
        renamed = replace(train, sample="renamed", source=self.test)
        with self.assertRaisesRegex(ValueError, "input hash overlap"):
            diagnose(renamed, expected)
        with self.assertRaisesRegex(ValueError, "Duplicate training input hashes"):
            fit_expected([train, renamed])

    def test_expected_mean_cannot_be_silently_replaced(self):
        expected = fit_expected([self.observation()])
        with self.assertRaisesRegex(ValueError, "equal-source mean"):
            replace(expected, probabilities=(.99, .01))

    def test_expected_profile_is_deeply_immutable(self):
        expected = fit_expected([self.observation()])
        values = [.5, .5]
        frozen = replace(expected, probabilities=values, source_probabilities=[values])
        values[0] = .9
        self.assertEqual(frozen.probabilities, (.5, .5))
        self.assertEqual(frozen.source_probabilities, ((.5, .5),))

    def test_mapping_and_panel_changes_are_incompatible(self):
        expected = fit_expected([self.observation()])
        for panel in (replace(self.panel, margin=6),
                      replace(self.panel, mapping_contract_sha256="d" * 64)):
            with self.assertRaisesRegex(ValueError, "mismatch"):
                diagnose(self.observation(source=self.test, sample="held", panel=panel), expected)

    def test_exact_shape_has_zero_divergence_without_p_value(self):
        expected = fit_expected([self.observation()])
        result = diagnose(self.observation((2, 2), self.test, "held"), expected)
        self.assertAlmostEqual(result["js_divergence_nats"], 0)
        self.assertAlmostEqual(result["composite_deviance"], 0)
        self.assertIsNone(result["p_value"])
        self.assertIsNone(result["identity_call"])
        self.assertIn("duplicates_possible", result["count_unit"])
        self.assertEqual(result["training_source_count"], 1)

    def test_low_and_high_count_sparsity_have_different_occupancy_expectations(self):
        expected = fit_expected([self.observation()])
        low = diagnose(self.observation((1, 0), self.test, "low"), expected)
        high = diagnose(self.observation((100, 0), self.test, "high"), expected)
        self.assertEqual(low["js_divergence_nats"], high["js_divergence_nats"])
        self.assertEqual(low["zero_opportunity_mass"], .5)
        self.assertAlmostEqual(low["zero_opportunity_excess_over_iid"], 0)
        self.assertAlmostEqual(high["zero_opportunity_excess_over_iid"], .5)
        self.assertAlmostEqual(high["composite_deviance"], 100 * low["composite_deviance"])

    def test_no_support_yields_no_distribution_test(self):
        result = diagnose(self.observation((0, 0), self.test, "empty"),
                          fit_expected([self.observation()]))
        self.assertEqual(result["status"], "no_support")
        self.assertEqual(result["iid_expected_zero_markers"], 2)
        self.assertEqual(result["zero_opportunity_excess_over_iid"], 0)
        self.assertIsNone(result["js_divergence_nats"])
        self.assertIsNone(result["identity_call"])

    def test_reference_zero_is_not_hidden_with_pseudocount(self):
        expected = fit_expected([self.observation((10, 0))])
        result = diagnose(self.observation((0, 10), self.test, "other"), expected)
        self.assertAlmostEqual(result["js_divergence_nats"], math.log(2))
        self.assertEqual(result["reference_zero_probability_fraction"], 1)
        self.assertEqual(result["deviance_status"], "reference_support_gap")
        self.assertIsNone(result["composite_deviance"])

    def test_single_marker_probability_one_and_zero_total_are_stable(self):
        panel = MarkerPanel("longum", (("a", 100),), "a" * 64)
        expected = fit_expected([self.observation((10,), panel=panel)])
        for count, zeros in ((0, 1), (1, 0)):
            result = diagnose(self.observation((count,), self.test, "test", panel), expected)
            self.assertEqual(result["iid_expected_zero_markers"], zeros)

    def test_serialized_expected_profile_declares_original_units(self):
        value = fit_expected([self.observation()]).to_dict()
        self.assertIn("equal_source", value["estimator"])
        self.assertEqual(value["smoothing"], "none; unsupported reference probability remains zero")
        self.assertEqual(value["panel_sha256"], self.panel.sha256)


if __name__ == "__main__":
    unittest.main()
