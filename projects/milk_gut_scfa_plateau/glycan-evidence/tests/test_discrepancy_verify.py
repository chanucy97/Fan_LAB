"""Small, hand-solvable synthetic witnesses; no biological validation claim."""
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from glycan_evidence import discrepancy_verify as verify

try:
    import numpy  # Load once before temporary sys.modules guards below.
except ImportError:
    numpy = None


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


@unittest.skipIf(numpy is None, "NumPy is required for numerical certificate replay")
class DiscrepancyVerifyTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)/"input"
        self.root.mkdir()
        self.out = Path(self.temp.name)/"verified"
        self.sample = "fixture"
        self.folder = self.root/"results"/self.sample
        self.folder.mkdir(parents=True)
        self.markers = ["marker_" + str(i) for i in range(220)]
        self.mapping = {"historical_mapping_contract_sha256": "9"*64,
            "run_panel_metadata_sha256": {"panel.tsv": "8"*64}, "engine_sha256": {}}
        self.raw = {"pairs_n": 192, "both_eof": True, "unique_qnames": True,
            "hash_unit": "complete_decompressed_FASTQ_bytes", "decompressed_sha256": ["a"*64, "b"*64],
            "input_stats": [{"path": "/raw/mate1.gz", "bytes": 100, "mtime_ns": 1},
                            {"path": "/raw/mate2.gz", "bytes": 101, "mtime_ns": 2}]}
        self.screen = {"total_fragments": 160, "candidate_fragments": 40,
            "both_mates_mapped": 30, "one_mate_mapped": 10, "unmapped_fragments": 120}
        self.qc_argv = ["/bin/fastp", "--stdout", "--json", "/old/fastp.json", "--html", "/old/fastp.html"]
        self.screen_argv = ["/bin/bowtie2", "--interleaved", "-"]
        self.old = {"sample": self.sample, "status": "COMPLETE", "skip_qc": False,
            "commands": [{"stage": "fastp", "argv": self.qc_argv},
                {"stage": "feature_screen", "argv": self.screen_argv},
                {"stage": "competitive_genome_mapping", "argv": ["/bin/bowtie2", "--no-unal"]}],
            "screen": self.screen, "versions": {"fastp": "fastp 1.3.7", "bowtie2": "bowtie2 2.5.4"},
            "panel_metadata": [{"path": "/panel/panel.tsv", "sha256": "8"*64}]}
        self.preflight = {"samples": {self.sample: {
            "raw_current_stat": self.raw["input_stats"], "frozen_row": {"expected_qc_pairs": 160},
            "assembly": {"receipt.json": {"sha256": "7"*64, "data": {"raw_input_validation": self.raw}}},
            "old_files": {"run_manifest.json": {"data": self.old}, "screen_counts.json": {"data": self.screen},
                "fastp.json": {"data": {"summary": {"before_filtering": {"total_reads": 384},
                    "after_filtering": {"total_reads": 320}}}}, "summary.json": {"data": {"total_fragments": 32}},
                "candidate_R1.fastq.gz": {"path": "/old/candidate_R1.fastq.gz"},
                "candidate_R2.fastq.gz": {"path": "/old/candidate_R2.fastq.gz"},
                "competitive.sam.gz": {"path": "/old/competitive.sam.gz", "bytes": 1000}}}},
            "bin_versions": {name: {"path": "/bin/" + name, "sha256": str(i)*64, "returncode": 0,
                                    "version": self.old["versions"][name]}
                             for i, name in enumerate(("fastp", "bowtie2"), 1)}}
        write(self.root/"PREFLIGHT.json", self.preflight)
        for name in ("bridge.py", "diagnostic.py", "verify_stage2.py", "frozen/audit/audit_actual.py",
                     "frozen/legacy/run_sample.py", "frozen/legacy/map_evidence.py"):
            path = self.root/name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("raise RuntimeError('INPUT PYTHON MUST NEVER EXECUTE')\n# " + name, encoding="utf-8")
        for name in ("run_sample.py", "map_evidence.py"):
            self.mapping["engine_sha256"][name] = sha(self.root/"frozen/legacy"/name)
        self.refs = []
        for i, (source, role, index) in enumerate((("zero", "background", None), ("bg", "background", 1), ("target", "target", 0))):
            counts = [0]*220
            if index is not None:
                counts[index] = 32
            self.refs.append({"source_id": source, "role": role, "marker_ids": self.markers,
                "input_pairs": 160, "counts": counts, "read_sha256": [str(2*i+1)*64, str(2*i+2)*64],
                "mapping_contract_sha256": "9"*64})
        model = self.root/"frozen/model"
        write(model/"ordered_markers.json", self.markers)
        write(model/"training_profiles.json", self.refs)
        write(model/"mapping_contract.json", self.mapping)
        write(model/"parameters.json", {"support_domain": [19, 5944],
            **{key: sha(model/file) for key, file in (("training_profiles_sha256", "training_profiles.json"),
                ("ordered_markers_sha256", "ordered_markers.json"), ("mapping_contract_sha256", "mapping_contract.json"))}})
        self.z = [i for i in range(220) if i != 1]
        self.bridge_freeze = {"samples": [self.sample], "buckets": 16, "bucket_seed": "fixture-seed",
            "assembly_receipt_sha256": {self.sample: "7"*64},
            "resource_contract": {"maximum_QC_compressed_bytes_per_sample": 10000}}
        self.freeze = {"samples": [self.sample], "bucket_count": 16, "bucket_seed": "fixture-seed",
            "maximum_H_calls_per_sample": 3, "new_q_fits": 0, "simulations": 0,
            "state": "reference_discrepancy_diagnostic", "Z_indices": self.z,
            "numeric_gate": {"feasibility": 1e-12, "convex_gap": 1e-8, "rate_decomposition_mixed_tolerance": 1e-12}}
        self.bridge = {"sample": self.sample, "status": "QC_BRIDGE_PASS_MAPPING_RECOUNT_PENDING",
            "reused_raw_validation": self.raw, "assembly_receipt_sha256": "7"*64,
            "new_H_calls": 0, "new_envelope_fits": 0, "qc_argv": self.qc_argv, "screen_argv": self.screen_argv,
            "qc": {"pairs_n": 160, "both_eof": True, "decompressed_sha256": ["c"*64, "d"*64],
                "ordered_canonical_qname_sha256": "e"*64, "N_buckets": [10]*16, "buckets": 16,
                "bucket_seed_utf8": "fixture-seed", "files": [
                    {"path": "/result/complete_QC_R1.fastq.gz", "bytes": 1000, "sha256": "1"*64},
                    {"path": "/result/complete_QC_R2.fastq.gz", "bytes": 1001, "sha256": "2"*64}],
                "unit": "all_complete_post_QC_paired_fragments_including_no_marker_fragments"},
            "screen": self.screen,
            "candidate_content_checks": [{"mate": m, "original": f"/old/candidate_R{m}.fastq.gz",
                "replayed": f"/new/candidate_R{m}.fastq.gz", "decompressed_sha256": str(m)*64,
                "new_compressed_sha256": str(m+2)*64} for m in (1, 2)],
            "competitive_sam": {"path": "/old/competitive.sam.gz", "bytes": 1000, "sha256": "6"*64}}
        self.make_measurement()
        self.seal()

    def make_measurement(self, zero=False, empty_buckets=False):
        counts = [0]*220
        counts[:3] = [0, 0, 0] if zero else [8, 16, 8]
        denominators = [20]*8+[0]*8 if empty_buckets else [10]*16
        self.bridge["qc"]["N_buckets"] = denominators
        self.query = {"query_id": self.sample, "input_pairs": 160, "counts": counts, "read_sha256": ["c"*64, "d"*64]}
        events = [[0]*220 for _ in range(16)]
        if not zero:
            for i in range(16):
                if empty_buckets:
                    if i < 8:
                        events[i][:3] = [1, 2, 1]
                else:
                    events[i][:3] = [int(i < 8), 1, int(i >= 8)]
        self.recount = {"marker_ids": self.markers, "counts": counts, "N_input_pairs": 160,
            "N_buckets": denominators, "marker_counts_by_bucket": events, "eligible_markers": 220,
            "sum_support": sum(counts), "qc_read_sha256": self.query["read_sha256"],
            "candidate_pairs": 40, "sam_primary_pairs_content_verified": 32, "sam_counters": {"total_fragments": 32},
            "historical_marker_equality": "ALL220_COUNT_COVERAGE_DEPTH_AND_AMBIGUITY_EXACT"}
        fit_specs = [("full", self.refs, list(range(220)), [1., 0., 0.] if zero else [.25, .5, .25]),
                     ("background", self.refs[:2], list(range(220)), [1., 0.] if zero else [.5, .5]),
                     ("local", [self.refs[2], self.refs[0]], self.z, [0., 1.] if zero else [.25, .75])]
        self.fits = {}
        for label, refs, indices, weights in fit_specs:
            support = sum(counts[j] for j in indices); scale = max(support, 1)
            predicted = [sum(r["counts"][j]*w for r, w in zip(refs, weights))/scale for j in indices]
            unsupported = [i for i, j in enumerate(indices) if counts[j] and all(r["counts"][j] == 0 for r in refs)]
            unreachable = sum(counts[indices[i]] for i in unsupported)
            loss = sum((math.sqrt(p)-math.sqrt(counts[j]/scale))**2 for p, j in zip(predicted, indices))
            fit = {"status": "finite_rate_hellinger_projection", "source_input_weights": weights,
                "fitted_scaled_rates": predicted, "loss": loss, "first_order_convex_gap": 0.,
                "unsupported_indices": unsupported, "unsupported_loss_constant": unreachable/scale,
                "raw_solver_fit": {"status": "finite_rate_hellinger_projection", "source_input_weights": weights, "loss": loss}}
            row = {"status": "CERTIFIED", "label": label, "query_id": self.sample, "support": support,
                "N_input_pairs": 160, "reference_source_ids": [r["source_id"] for r in refs], "marker_indices": indices,
                "read_sha256": self.query["read_sha256"], "query_counts_sha256": hashlib.sha256(
                    json.dumps(counts, separators=(",", ":")).encode()).hexdigest(), "new_H_calls": 1, "reused_fit": False,
                "zero_reference_sources": [r["source_id"] for r in refs if not any(r["counts"][j] for j in indices)],
                "unsupported_support": unreachable, "unsupported_fixed_loss": unreachable/scale, "fit": fit}
            self.fits[label] = row
        rows = []
        for j, count in enumerate(counts):
            p = 0. if zero else ([8, 16, 0][j]/160 if j < 3 else 0.)
            rate = count/160
            variation = total = 0.
            same = 0
            for N, events_k in zip(denominators, events):
                if N:
                    rate_k = events_k[j]/N
                    variation += N/160*(rate_k-rate)**2
                    total += N/160*(rate_k-p)**2
                    same += (rate_k-p)*(rate-p) > 0
            row = {"marker_id": self.markers[j], "observed_events": count, "observed_rate": rate,
                "fitted_rate": p, "signed_rate_residual": rate-p, "common_squared_rate_residual": (rate-p)**2,
                "weighted_partition_variation": variation, "weighted_partition_total_residual": total,
                "same_sign_nonempty_buckets": same, "nonempty_buckets": sum(N > 0 for N in denominators),
                "interpretation": "descriptive_technical_partition_no_sampling_CI", "local_marker": j in self.z}
            for label, fitted_key, loss_key in (("full", "fitted_rate", "loss_contribution"),
                    ("background", "background_fitted_rate", "background_loss_contribution"),
                    ("local", "local_fitted_rate", "local_loss_contribution")):
                if label == "local" and j not in self.z:
                    row[fitted_key] = row[loss_key] = None
                else:
                    pos = self.z.index(j) if label == "local" else j
                    scale = max(self.fits[label]["support"], 1)
                    fitted = self.fits[label]["fit"]["fitted_scaled_rates"][pos]*scale/160
                    row[fitted_key] = fitted
                    row[loss_key] = (math.sqrt(fitted*160)-math.sqrt(count))**2/scale
            rows.append(row)
        self.rows = rows
        self.result = {"sample": self.sample, "state": "reference_discrepancy_diagnostic", "N_input_pairs": 160,
            "S_full": sum(counts), "S_local": sum(counts[j] for j in self.z),
            **{"D_"+label: row["fit"]["loss"] for label, row in self.fits.items()},
            "D_background_minus_full": self.fits["background"]["fit"]["loss"]-self.fits["full"]["fit"]["loss"],
            **{key: None for key in verify.NULL_FIELDS}, "new_H_calls": 3, "new_q_fits": 0,
            "full_empirical_support_domain": [19, 5944], "local_empirical_support_domain": [19, 2545],
            "full_outside_domain": zero, "local_outside_domain": True}

    def seal(self):
        """Rebind explicit synthetic trust anchors; used to test semantic rejection."""
        write(self.root/"PREFLIGHT.json", self.preflight)
        model = self.root/"frozen/model"
        write(model/"MANIFEST.json", {"version": "reference-family-compatibility.research.v0.1",
            "identity_call": None, "probability_calibration": False,
            "files_sha256": {p.name: sha(p) for p in model.iterdir() if p.name != "MANIFEST.json"}})
        self.bridge_freeze["files_sha256"] = {name: sha(self.root/name) for name in ("PREFLIGHT.json", "bridge.py")}
        write(self.root/"FREEZE.json", self.bridge_freeze)
        self.freeze["model_manifest_sha256"] = sha(model/"MANIFEST.json")
        names = set(verify.REQUIRED_FROZEN) | {"frozen/legacy/run_sample.py", "frozen/legacy/map_evidence.py"}
        self.freeze["files_sha256"] = {name: sha(self.root/name) for name in sorted(names)}
        write(self.root/"STAGE2_FREEZE.json", self.freeze)
        self.freeze_sha = sha(self.root/"STAGE2_FREEZE.json")
        self.bridge["freeze_sha256"] = sha(self.root/"FREEZE.json")
        write(self.folder/"RECEIPT.json", self.bridge)
        (self.folder/"QC_BRIDGE_PASS").write_text(sha(self.folder/"RECEIPT.json") + "\n", encoding="ascii")
        self.result["input_bridge_receipt_sha256"] = sha(self.folder/"RECEIPT.json")
        for name, value in [("QUERY.json", self.query), ("MARKER_RECOUNT.json", self.recount),
                ("MARKER_RESIDUALS.json", self.rows), ("DESCRIPTIVE_RESULT.json", self.result),
                *[("FIT_"+label+".json", row) for label, row in self.fits.items()]]:
            write(self.folder/name, value)
        write(self.folder/"STAGE2_RESULTS.json", {"sample": self.sample, "stage2_freeze_sha256": self.freeze_sha,
            "status": "AWAITING_INDEPENDENT_WITNESS_REPLAY", "technical_gate": "AWAITING_INDEPENDENT_WITNESS_REPLAY",
            "actual_H_calls": 3, "files_sha256": {name: sha(self.folder/name) for name in sorted(verify.RESULT_FILES)}})
        self.results_sha = sha(self.folder/"STAGE2_RESULTS.json")

    def run_verify(self):
        with patch.dict("sys.modules", {"scipy": None, "scipy.optimize": None}), \
             patch("subprocess.run", side_effect=AssertionError("No subprocess allowed")):
            return verify.run_discrepancy_verify(self.root, self.sample, self.freeze_sha, self.results_sha, self.out)

    def reject(self, pattern=None):
        with self.assertRaisesRegex(ValueError, pattern or ".*"):
            self.run_verify()
        self.assertFalse((self.out/"COMPLETE").exists())

    def test_complete_hand_solved_three_witnesses_without_scipy_or_input_code(self):
        receipt = self.run_verify()
        self.assertEqual(receipt["status"], "COMPLETE")
        self.assertEqual(receipt["optimizer_calls"], 0)
        self.assertEqual(receipt["new_H_calls"], 0)
        self.assertEqual(receipt["historical_H_calls_replayed"], 3)
        for row, expected in zip(receipt["fit_witnesses"], [.25, .5, .5]):
            self.assertAlmostEqual(row["loss_recomputed"], expected, places=14)
        self.assertEqual(receipt["preprocessing"]["post_QC_pairs"], 160)
        self.assertEqual(receipt["preprocessing"]["raw_pairs"], 192)
        self.assertEqual(receipt["preprocessing"]["candidate_groups_omitted_by_no_unal"], 8)
        self.assertEqual((self.out/"COMPLETE").read_text().strip(), sha(self.out/"receipt.json"))
        for name, expected in receipt["files_sha256"].items():
            self.assertEqual(sha(self.out/name), expected)
        self.assertTrue(all(receipt[k] is None for k in verify.NULL_FIELDS))

    def test_zero_support_and_zero_bucket_denominators_are_retained(self):
        self.make_measurement(zero=True, empty_buckets=True); self.seal()
        receipt = self.run_verify()
        self.assertEqual([r["loss_recomputed"] for r in receipt["fit_witnesses"]], [0., 0., 0.])
        self.assertEqual(len(json.loads((self.out/"marker_residuals.json").read_text())), 220)

    def test_nonempty_measurement_with_empty_buckets(self):
        self.make_measurement(empty_buckets=True); self.seal()
        self.assertEqual(self.run_verify()["technical_gate"], "PASS")

    def test_missing_member_rejected(self):
        (self.folder/"QUERY.json").unlink(); self.reject("Missing regular input")

    def test_file_hash_tampering_rejected(self):
        (self.folder/"FIT_full.json").write_text("{}", encoding="utf-8"); self.reject("Hash mismatch")

    def test_wrong_external_anchor_rejected(self):
        self.freeze_sha = "0"*64; self.reject("Hash mismatch")

    def test_required_manifest_member_cannot_be_removed(self):
        stage = json.loads((self.folder/"STAGE2_RESULTS.json").read_text())
        del stage["files_sha256"]["QUERY.json"]
        write(self.folder/"STAGE2_RESULTS.json", stage); self.results_sha = sha(self.folder/"STAGE2_RESULTS.json")
        self.reject("manifest members")

    def test_duplicate_json_key_rejected(self):
        path = self.root/"STAGE2_FREEZE.json"
        path.write_text('{"bucket_count":16,"bucket_count":16}', encoding="utf-8")
        self.freeze_sha = sha(path); self.reject("Duplicate JSON key")

    def test_unsafe_manifest_path_rejected(self):
        self.freeze["files_sha256"]["../outside"] = "1"*64
        write(self.root/"STAGE2_FREEZE.json", self.freeze); self.freeze_sha = sha(self.root/"STAGE2_FREEZE.json")
        self.reject("Unsafe")

    def test_duplicate_case_alias_path_rejected(self):
        self.freeze["files_sha256"]["BRIDGE.PY"] = self.freeze["files_sha256"]["bridge.py"]
        write(self.root/"STAGE2_FREEZE.json", self.freeze); self.freeze_sha = sha(self.root/"STAGE2_FREEZE.json")
        self.reject("Duplicate/case")

    def test_duplicate_physical_manifest_member_rejected(self):
        duplicate = self.root/"duplicate.py"
        os.link(self.root/"bridge.py", duplicate)
        self.freeze["files_sha256"]["duplicate.py"] = sha(duplicate)
        write(self.root/"STAGE2_FREEZE.json", self.freeze); self.freeze_sha = sha(self.root/"STAGE2_FREEZE.json")
        self.reject("Duplicate physical")

    def test_sample_path_traversal_rejected(self):
        self.sample = "../fixture"; self.reject("Unsafe")

    def test_boolean_numeric_value_rejected(self):
        self.query["counts"][0] = True; self.seal(); self.reject("Invalid integer")

    def test_boolean_fit_weight_rejected(self):
        self.fits["full"]["fit"]["source_input_weights"][0] = True; self.seal(); self.reject("Invalid finite number")

    def test_negative_count_rejected(self):
        self.recount["marker_counts_by_bucket"][0][5] = -1; self.seal(); self.reject("Invalid integer")

    def test_nonfinite_json_rejected(self):
        path = self.folder/"QUERY.json"
        path.write_text(json.dumps(self.query)[:-1] + ', "hidden":1e999}', encoding="utf-8")
        stage = json.loads((self.folder/"STAGE2_RESULTS.json").read_text()); stage["files_sha256"]["QUERY.json"] = sha(path)
        write(self.folder/"STAGE2_RESULTS.json", stage); self.results_sha = sha(self.folder/"STAGE2_RESULTS.json")
        self.reject("Non-finite")

    def test_raw_denominator_cannot_replace_post_QC_N(self):
        self.query["input_pairs"] = 192; self.seal(); self.reject("count/N mismatch")

    def test_candidate_denominator_cannot_replace_post_QC_N(self):
        self.query["input_pairs"] = 40; self.recount["N_input_pairs"] = 40; self.seal(); self.reject("bucket N addback")

    def test_reordered_markers_rejected(self):
        self.recount["marker_ids"] = list(reversed(self.markers)); self.seal(); self.reject("marker order")

    def test_marker_bucket_addback_rejected(self):
        self.recount["marker_counts_by_bucket"][0][0] += 1; self.seal(); self.reject("Marker bucket addback")

    def test_marker_event_greater_than_bucket_N_rejected(self):
        self.recount["marker_counts_by_bucket"][0][0] = 11; self.seal(); self.reject("Invalid integer")

    def test_tampered_feasible_but_nonoptimal_certificate_rejected(self):
        self.fits["full"]["fit"]["source_input_weights"] = [.5, .25, .25]
        self.seal(); self.reject("Stored prediction")

    def test_tampered_loss_certificate_rejected(self):
        self.fits["full"]["fit"]["loss"] += .01; self.seal(); self.reject("Arithmetic mismatch")

    def test_consistent_prediction_with_nonoptimal_witness_rejected_by_convex_gap(self):
        fit = self.fits["full"]["fit"]
        fit["source_input_weights"] = [.5, .25, .25]
        fit["fitted_scaled_rates"] = [.25, .25] + [0.]*218
        fit["loss"] = .25 + (math.sqrt(.25)-math.sqrt(.5))**2
        fit["raw_solver_fit"]["source_input_weights"] = fit["source_input_weights"]
        fit["raw_solver_fit"]["loss"] = fit["loss"]
        self.seal(); self.reject("Convex gap")

    def test_negative_feasible_weight_rejected(self):
        self.fits["full"]["fit"]["source_input_weights"][0] = -.1
        self.seal(); self.reject("Invalid finite number")

    def test_changed_reference_order_rejected(self):
        self.fits["full"]["reference_source_ids"] = ["bg", "zero", "target"]
        self.seal(); self.reject("reference order")

    def test_changed_local_coordinate_mask_rejected(self):
        self.freeze["Z_indices"] = self.z[1:]
        self.seal(); self.reject("Local coordinate")

    def test_final_recorded_QC_sizes_respect_cap(self):
        self.bridge_freeze["resource_contract"]["maximum_QC_compressed_bytes_per_sample"] = 2000
        self.seal(); self.reject("QC gzip size")

    def test_duplicate_QC_file_record_rejected(self):
        self.bridge["qc"]["files"][1]["path"] = self.bridge["qc"]["files"][0]["path"]
        self.seal(); self.reject("Duplicate recorded")

    def test_changed_raw_repair_certificate_rejected(self):
        self.fits["full"]["fit"]["raw_solver_fit"]["source_input_weights"] = [.5, .25, .25]
        self.seal(); self.reject("Repaired witness")

    def test_changed_query_hash_in_fit_rejected(self):
        self.fits["full"]["read_sha256"] = ["0"*64, "1"*64]
        self.seal(); self.reject("Fit read/N")

    def test_changed_local_residual_outside_mask_rejected(self):
        self.rows[1]["local_loss_contribution"] = 0.
        self.seal(); self.reject("Nonlocal marker")

    def test_observed_N_buckets_include_fragments_without_markers(self):
        self.bridge["qc"]["N_buckets"] = [2]*16
        self.seal(); self.reject("QC/recount bucket")

    def test_missing_required_field_is_a_value_error(self):
        del self.query["input_pairs"]
        self.seal(); self.reject("Missing or malformed")

    def test_tampered_residual_not_hidden_by_loss_sum(self):
        self.rows[0]["loss_contribution"] += .001
        self.rows[2]["loss_contribution"] -= .001
        self.seal(); self.reject("loss_contribution")

    def test_tampered_partition_variance_rejected(self):
        self.rows[0]["weighted_partition_variation"] += .001; self.seal(); self.reject("weighted_partition_variation")

    def test_q_and_identity_must_remain_null(self):
        self.result["full_q"] = .5; self.seal(); self.reject("must remain null")

    def test_QC_receipt_is_bound_to_descriptor(self):
        self.bridge["qc"]["pairs_n"] = 192
        write(self.folder/"RECEIPT.json", self.bridge)
        (self.folder/"QC_BRIDGE_PASS").write_text(sha(self.folder/"RECEIPT.json"))
        self.reject("completion/descriptor")

    def test_QC_argument_change_rejected_even_with_new_anchor(self):
        self.bridge["qc_argv"] = self.qc_argv + ["--length_required", "5"]
        self.seal(); self.reject("QC argv")

    def test_frozen_tolerance_change_rejected(self):
        self.freeze["numeric_gate"]["convex_gap"] = 1.; self.seal(); self.reject("Changed frozen")

    def test_output_directory_overwrite_rejected(self):
        self.out.mkdir()
        with self.assertRaises(FileExistsError):
            self.run_verify()

    def test_symlink_member_rejected(self):
        link = self.root/"bridge.py"; contents = link.read_bytes()
        target = Path(self.temp.name)/"real_bridge.py"; target.write_bytes(contents); link.unlink()
        try:
            link.symlink_to(target)
        except OSError:
            self.skipTest("Host does not permit creation of symlinks")
        self.reject("Symlink/reparse")

    def test_input_mutation_before_completion_leaves_no_COMPLETE(self):
        original = verify._Inputs.recheck
        calls = 0
        def changed(inputs):
            nonlocal calls
            calls += 1
            if calls == 2:
                (self.folder/"QUERY.json").write_text("{}", encoding="utf-8")
            return original(inputs)
        with patch.object(verify._Inputs, "recheck", changed):
            self.reject("Input changed")


if __name__ == "__main__":
    unittest.main()
