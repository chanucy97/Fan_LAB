"""Small synthetic metric fixtures test software, never classifier performance."""

import csv
import hashlib
import itertools
import json
from pathlib import Path
import tempfile
import unittest

from glycan_evidence.legacy_v2 import calibrate_calls as calls


MARGINS = [0, 6, 12, 18]
PLAN = {
    "primary_identity_targets": ["infantis", "longum"],
    "competition_margin_grid": MARGINS,
    "candidate_detection_grid": {"marker_breadth": [0.05, 0.1, 0.2, 0.4],
                                 "minimum_fragments_per_marker": [1, 2, 3],
                                 "supported_identity_markers": [3, 5, 10, 20]},
    "benchmark": {"positive_target_depth_x": [0.05, 0.2, 1, 5], "substitution_rates": [0.005, 0.01]},
}


def make_features():
    result = {}
    for target in PLAN["primary_identity_targets"]:
        for i in range(30):
            key = f"I_{target}_{i:02d}"
            result[key] = {"feature_id": key, "channel": "identity", "taxon_group": target,
                           "length": 1000, "eligible_identity": 1, "pseudo": 0}
    result["excluded"] = {"feature_id": "excluded", "channel": "identity", "taxon_group": "infantis",
                          "length": 1000, "eligible_identity": 0, "pseudo": 0}
    result["pseudo_identity"] = {"feature_id": "pseudo_identity", "channel": "identity", "taxon_group": "infantis",
                                 "length": 1000, "eligible_identity": 1, "pseudo": 1}
    for key, pseudo in [("L_glycosidase", 0), ("L_transposase", 1)]:
        result[key] = {"feature_id": key, "channel": "locus", "taxon_group": "unassigned",
                       "length": 1000, "eligible_identity": 0, "pseudo": pseudo}
    return result


def metric_data(sample, strength):
    """strength(feature, margin) -> (covered proportion, fragments)."""
    features = make_features()
    evidence = {}
    for key, feature in features.items():
        for margin in MARGINS:
            breadth, fragments = strength(feature, margin)
            covered = int(round(breadth * feature["length"]))
            aligned = covered + fragments if covered else 0
            evidence[(key, float(margin))] = {"sample": sample, **feature, "margin": margin,
                                            "fragments": fragments, "covered_bases": covered,
                                            "breadth": covered / feature["length"], "aligned_bases": aligned,
                                            "mean_depth": aligned / feature["length"], "ambiguous_fragments": 0}
    return {"features": features, "evidence": evidence, "schema_sha256": calls.schema_hash(features)}


def feature_annotations():
    return [{**f, "cluster": "H1" if f["channel"] == "locus" else "",
             "product": "transposase pseudogene" if f["feature_id"] == "L_transposase" else "glycosidase" if f["channel"] == "locus" else "identity marker"} for f in make_features().values()]


def target_strength(truth, count=30, breadth=0.4, fragments=3):
    def strength(feature, margin):
        if feature["feature_id"].startswith("I_") and feature["taxon_group"] == truth:
            return (breadth, fragments) if int(feature["feature_id"].rsplit("_", 1)[1]) < count else (0, 0)
        return 0, 0
    return strength


def benchmark_rows(split="calibration", prefix="cal"):
    return [{"case_id": f"{prefix}_{target}_{i}", "split": split, "source_split": split, "truth_taxon": target,
             "accession": f"{prefix}_{target}_genome", "coverage": depth, "error_rate": error,
             "truth_mode": "synthetic_reads_from_real_reference_genome",
             "source_sha256": hashlib.sha256(f"{prefix}_{target}_genome".encode()).hexdigest()}
            for target in PLAN["primary_identity_targets"]
            for i, (depth, error) in enumerate(itertools.product([0.05, 0.2, 1, 5], [0.005, 0.01]))]


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="call_test_")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.metrics = self.root / "metrics"
        self.plan_path = self.root / "plan.json"
        self.plan_path.write_text(json.dumps(PLAN))
        calls.write_table(self.root / "frozen_features.tsv", feature_annotations())
        self.plan, self.targets, self.margins, self.grid = calls.read_plan(self.plan_path)

    def store_metrics(self, sample, data):
        folder = self.metrics / sample
        folder.mkdir(parents=True)
        calls.write_table(folder / "gene_metrics.tsv", list(data["evidence"].values()), calls.METRIC_FIELDS)
        # Explicit synthetic provenance fixtures, not real mapping evidence.
        fixture_record = lambda name: {"path": name, "sha256": hashlib.sha256(name.encode()).hexdigest()}
        run = {"sample": sample, "status": "COMPLETE", "skip_qc": True,
               "commands": [{"stage": "synthetic_fixture", "returncode": 0}], "versions": {"bowtie2": "synthetic_fixture version 0"},
               "panel_metadata": [fixture_record(name) for name in ("panel_contigs.tsv", "feature_manifest.tsv", "annotations.tsv")],
               "panel_index_files": [fixture_record(name) for name in ("features.1.bt2", "genomes.1.bt2")],
               "code": [fixture_record(name) for name in ("run_sample.py", "map_evidence.py")],
               "outputs": [{"path": "gene_metrics.tsv", "sha256": calls.sha256(folder / "gene_metrics.tsv")}]}
        next(r for r in run["panel_metadata"] if r["path"] == "feature_manifest.tsv")["sha256"] = calls.sha256(self.root / "frozen_features.tsv")
        (folder / "run_manifest.json").write_text(json.dumps(run))
        (folder / "COMPLETE").write_text("synthetic software test fixture\n")

    def store_manifest(self, name, rows):
        path = self.root / name
        calls.write_table(path, rows)
        return path

    def calibrate_fixture(self, count=30, breadth=0.4, fragments=3):
        rows = benchmark_rows()
        path = self.store_manifest("cal.tsv", rows)
        for row in rows:
            self.store_metrics(row["case_id"], metric_data(row["case_id"], target_strength(row["truth_taxon"], count, breadth, fragments)))
        out = self.root / "calibration"
        result = calls.main(["--calibration-manifest", str(path), "--metrics-root", str(self.metrics),
                             "--freeze-plan", str(self.plan_path), "--out", str(out)])
        self.assertEqual(result, 0)
        return out / "frozen_rules.json"

    def test_full_grid_and_strict_tie_order(self):
        rows = benchmark_rows()
        loaded = {r["case_id"]: metric_data(r["case_id"], target_strength(r["truth_taxon"])) for r in rows}
        grid, rule, predictions = calls.select_rule(rows, loaded, self.targets, self.grid)
        self.assertEqual(len(grid), 192)
        self.assertEqual(rule, {"margin": 18, "marker_breadth": 0.4,
                                "minimum_fragments_per_marker": 3, "supported_identity_markers": 20})
        self.assertEqual(sum(r["selected"] for r in grid), 1)
        self.assertEqual(len(predictions), 32)
        self.assertTrue(all(r["detected"] == r["truth_positive"] for r in predictions))

    def test_weak_positive_rule_keeps_sensitivity_then_prefers_strict(self):
        rows = benchmark_rows()
        loaded = {r["case_id"]: metric_data(r["case_id"], target_strength(r["truth_taxon"], 8, 0.2, 2)) for r in rows}
        _, rule, _ = calls.select_rule(rows, loaded, self.targets, self.grid)
        self.assertEqual(rule, {"margin": 18, "marker_breadth": 0.2,
                                "minimum_fragments_per_marker": 2, "supported_identity_markers": 5})

    def test_balances_targets_not_case_counts(self):
        rows = benchmark_rows()
        # One longum case, eight infantis cases: macro sensitivity remains 0.5.
        rows = [r for r in rows if r["truth_taxon"] == "infantis"] + [rows[-1]]
        loaded = {r["case_id"]: metric_data(r["case_id"], target_strength("infantis") if r["truth_taxon"] == "infantis" else lambda f, m: (0, 0)) for r in rows}
        grid, _, _ = calls.select_rule(rows, loaded, self.targets, self.grid)
        self.assertEqual(next(r["balanced_target_positive_sensitivity"] for r in grid if r["selected"]), 0.5)

    def test_locus_pseudo_and_excluded_markers_do_not_support_identity(self):
        data = metric_data("s", lambda f, margin: (0.9, 10) if f["feature_id"] in {"excluded", "pseudo_identity", "L_glycosidase", "L_transposase", "I_infantis_00"} else (0, 0))
        rule = {"margin": 0, "marker_breadth": 0.05, "minimum_fragments_per_marker": 1, "supported_identity_markers": 3}
        result = calls.identity_result(data, "infantis", rule)
        self.assertFalse(result["detected"])
        self.assertEqual(result["supported_markers"], 1)
        self.assertEqual(result["eligible_markers"], 30)

    def test_no_zero_false_positive_rule_never_freezes(self):
        rows = benchmark_rows()
        path = self.store_manifest("bad.tsv", rows)
        for row in rows:
            self.store_metrics(row["case_id"], metric_data(row["case_id"], lambda f, margin: (0.8, 10)))
        out = self.root / "failed"
        result = calls.main(["--calibration-manifest", str(path), "--metrics-root", str(self.metrics),
                             "--freeze-plan", str(self.plan_path), "--out", str(out)])
        self.assertEqual(result, 3)
        self.assertFalse((out / "frozen_rules.json").exists())
        self.assertEqual(len(calls.read_table(out / "calibration_grid.tsv")), 192)
        self.assertEqual(json.loads((out / "calibration_status.json").read_text())["status"], "calibration_failed_no_zero_FP_rule")

    def test_test_evaluation_scores_false_positives_without_retuning(self):
        rules = self.calibrate_fixture(count=8, breadth=0.2, fragments=2)
        before = calls.sha256(rules)
        rows = benchmark_rows("test", "heldout")
        path = self.store_manifest("test.tsv", rows)
        for row in rows:
            self.store_metrics(row["case_id"], metric_data(row["case_id"], lambda f, margin: (0.9, 10)))
        out = self.root / "evaluation"
        self.assertEqual(calls.main(["--evaluate-manifest", str(path), "--rules", str(rules), "--metrics-root", str(self.metrics), "--out", str(out)]), 0)
        self.assertEqual(before, calls.sha256(rules))
        confusion = calls.read_table(out / "evaluation_confusion_by_target.tsv")
        self.assertTrue(all(int(row["FP"]) == 8 for row in confusion))
        self.assertEqual(len(calls.read_table(out / "evaluation_confusion_by_target_depth_error.tsv")), 16)
        self.assertFalse((out / "calibration_grid.tsv").exists())

    def test_manifest_split_and_accession_leakage_rejected(self):
        mixed = benchmark_rows()
        mixed[-1]["split"] = "test"
        with self.assertRaisesRegex(ValueError, "refuses source split"):
            calls.load_benchmark_manifest(self.store_manifest("mixed.tsv", mixed), "calibration")
        rules = self.calibrate_fixture()
        rows = benchmark_rows("test", "heldout")
        for row in rows:
            if row["truth_taxon"] == "infantis":
                row["accession"] = "cal_infantis_genome"
        path = self.store_manifest("leak.tsv", rows)
        with self.assertRaisesRegex(ValueError, "also occurs in calibration"):
            calls.main(["--evaluate-manifest", str(path), "--rules", str(rules), "--metrics-root", str(self.metrics), "--out", str(self.root / "never")])
        self.assertFalse((self.root / "never").exists())

    def test_metrics_missing_rows_and_nonfinite_values_are_errors(self):
        data = metric_data("s", target_strength("infantis"))
        data["evidence"].pop(("I_infantis_00", 18.0))
        self.store_metrics("s", data)
        with self.assertRaisesRegex(ValueError, "feature x margin grid"):
            calls.load_metrics(self.metrics / "s/gene_metrics.tsv", "s", MARGINS)
        for value in ["NaN", "inf", "-0.1"]:
            with self.assertRaises(ValueError):
                calls.nonnegative_number(value, "breadth")

    def test_feature_dictionary_drift_rejected(self):
        rules = self.calibrate_fixture()
        rows = benchmark_rows("test", "heldout")
        path = self.store_manifest("drift.tsv", rows)
        for row in rows:
            data = metric_data(row["case_id"], target_strength(row["truth_taxon"]))
            for metric in data["evidence"].values():
                if metric["feature_id"] == "I_infantis_00":
                    metric["eligible_identity"] = 0
            self.store_metrics(row["case_id"], data)
        with self.assertRaisesRegex(ValueError, "Feature dictionary changed"):
            calls.main(["--evaluate-manifest", str(path), "--rules", str(rules), "--metrics-root", str(self.metrics), "--out", str(self.root / "drift")])

    def test_cohort_missing_metrics_retained_and_locus_roles_preserved(self):
        rules = self.calibrate_fixture()
        cohort_path = self.root / "cohort.csv"
        with cohort_path.open("w", newline="") as handle:
            w = csv.DictWriter(handle, fieldnames=["SampleID", "DyadID", "Window"])
            w.writeheader(); w.writerows([{"SampleID": "BF_A", "DyadID": "A", "Window": "D05"}, {"SampleID": "BF_B", "DyadID": "B", "Window": "D14"}])
        data = metric_data("BF_A", lambda f, m: (0.8, 3) if f["channel"] == "locus" or f["taxon_group"] == "infantis" else (0, 0))
        self.store_metrics("BF_A", data)
        annotations = feature_annotations()
        feature_path = self.store_manifest("features.tsv", annotations)
        out = self.root / "cohort_out"
        self.assertEqual(calls.main(["--cohort-manifest", str(cohort_path), "--rules", str(rules), "--metrics-root", str(self.metrics), "--feature-manifest", str(feature_path), "--out", str(out)]), 0)
        identity = calls.read_table(out / "cohort_identity_calls.tsv")
        self.assertEqual(len(identity), 4)
        self.assertEqual({r["identity_state"] for r in identity}, {"detected", "insufficient_evidence"})
        missing = [r for r in identity if r["sample"] == "BF_B"]
        self.assertTrue(all(r["supported_markers"] == "" and r["reason"] == "metrics_missing" for r in missing))
        locus = calls.read_table(out / "cohort_locus_sequence_evidence.tsv")
        self.assertEqual(len(locus), 2 * 2 * 4)
        self.assertTrue(all(r["sequence_evidence_state"] == "well_covered_SEQUENCE_EVIDENCE" for r in locus if r["sample"] == "BF_A"))
        self.assertTrue(all(r["breadth"] == "" and r["sequence_evidence_state"] == "not_assessed_missing_metrics" for r in locus if r["sample"] == "BF_B"))
        self.assertEqual({r["product"] for r in locus if r["feature_id"] == "L_transposase"}, {"transposase pseudogene"})
        self.assertTrue(all(r["function_validated"] == "False" for r in locus))
        summary = calls.read_table(out / "cohort_identity_summary.tsv")
        self.assertTrue(all(int(r["samples_n"]) == 2 and int(r["missing_metrics_n"]) == 1 for r in summary))

    def test_absent_calibration_target_fails(self):
        rows = [r for r in benchmark_rows() if r["truth_taxon"] == "infantis"]
        loaded = {r["case_id"]: metric_data(r["case_id"], target_strength("infantis")) for r in rows}
        with self.assertRaisesRegex(ValueError, "positive and negative"):
            calls.select_rule(rows, loaded, self.targets, self.grid)

    def test_truth_conflict_and_engineered_truth_rejected(self):
        rows = benchmark_rows()
        rows[-1]["truth_taxon"] = "infantis"
        with self.assertRaisesRegex(ValueError, "conflicting truth"):
            calls.load_benchmark_manifest(self.store_manifest("conflict.tsv", rows), "calibration")
        rows = benchmark_rows()
        rows[0]["truth_mode"] = "engineered_fixture"
        with self.assertRaisesRegex(ValueError, "engineered"):
            calls.load_benchmark_manifest(self.store_manifest("engineered.tsv", rows), "calibration")

    def test_calibration_missing_or_duplicate_condition_rejected(self):
        rows = benchmark_rows()[:-1]
        path = self.store_manifest("missing_condition.tsv", rows)
        with self.assertRaisesRegex(ValueError, "complete frozen depth"):
            calls.main(["--calibration-manifest", str(path), "--freeze-plan", str(self.plan_path), "--metrics-root", str(self.metrics), "--out", str(self.root / "missing")])
        rows = benchmark_rows()
        rows.append({**rows[0], "case_id": "duplicated_condition"})
        with self.assertRaisesRegex(ValueError, "Duplicate accession"):
            calls.load_benchmark_manifest(self.store_manifest("duplicate_condition.tsv", rows), "calibration")

    def test_relabelled_source_and_same_source_hash_rejected(self):
        rows = benchmark_rows("test", "heldout")
        rows[0]["source_split"] = "index"
        with self.assertRaisesRegex(ValueError, "original source_split"):
            calls.load_benchmark_manifest(self.store_manifest("relabelled.tsv", rows), "test")
        rules = self.calibrate_fixture()
        rows = benchmark_rows("test", "heldout")
        for row in rows:
            if row["truth_taxon"] == "infantis":
                row["source_sha256"] = hashlib.sha256(b"cal_infantis_genome").hexdigest()
        path = self.store_manifest("source_leak.tsv", rows)
        with self.assertRaisesRegex(ValueError, "source file SHA"):
            calls.main(["--evaluate-manifest", str(path), "--rules", str(rules), "--metrics-root", str(self.metrics), "--out", str(self.root / "leak")])

    def test_pipeline_completion_hash_and_panel_drift_checked(self):
        data = metric_data("sample", target_strength("infantis"))
        self.store_metrics("sample", data)
        folder = self.metrics / "sample"
        with self.assertRaisesRegex(ValueError, "metrics SHA"):
            calls.pipeline_provenance(folder / "gene_metrics.tsv", "sample", "f" * 64)
        (folder / "COMPLETE").unlink()
        with self.assertRaisesRegex(ValueError, "missing COMPLETE"):
            calls.pipeline_provenance(folder / "gene_metrics.tsv", "sample", calls.sha256(folder / "gene_metrics.tsv"))
        rules = self.calibrate_fixture()
        rows = benchmark_rows("test", "heldout")
        for row in rows:
            self.store_metrics(row["case_id"], metric_data(row["case_id"], target_strength(row["truth_taxon"])))
        changed = self.metrics / rows[0]["case_id"] / "run_manifest.json"
        run = json.loads(changed.read_text())
        run["panel_index_files"][0]["sha256"] = "a" * 64
        changed.write_text(json.dumps(run))
        path = self.store_manifest("different_panel.tsv", rows)
        with self.assertRaisesRegex(ValueError, "Mapping panel/code/version changed"):
            calls.main(["--evaluate-manifest", str(path), "--rules", str(rules), "--metrics-root", str(self.metrics), "--out", str(self.root / "different_panel")])


if __name__ == "__main__":
    unittest.main(verbosity=2)
