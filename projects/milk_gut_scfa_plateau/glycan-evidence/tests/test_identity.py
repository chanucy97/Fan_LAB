"""Adapter failure and threshold tests; evidence is entirely synthetic."""
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

from glycan_evidence import identity
from glycan_evidence.legacy_v2 import calibrate_calls as legacy


class IdentityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.rules = identity.DEFAULT_RULES
        self.frozen = identity.load_rules(self.rules)
        path = Path(__file__).parent / "legacy_v2/feature_schema.json"
        self.features = json.loads(path.read_text(encoding="utf8"))["features"]
        self.make_metrics()

    def make_metrics(self, sample="TOY", supported_n=20):
        folder = self.root / "metrics" / sample
        folder.mkdir(parents=True, exist_ok=True)
        eligible = [f["feature_id"] for f in self.features if f["channel"] == "identity"
                    and f["taxon_group"] == "infantis" and f["eligible_identity"] == "1" and f["pseudo"] == "0"]
        supported = set(eligible[:supported_n])
        rows = []
        for feature in self.features:
            length = int(feature["length"])
            covered = math.ceil(length * 0.1) if feature["feature_id"] in supported else 0
            for margin in self.frozen["metrics_margins"]:
                rows.append({"sample": sample, **feature, "margin": margin,
                             "fragments": int(bool(covered)), "covered_bases": covered,
                             "breadth": covered / length, "aligned_bases": covered,
                             "mean_depth": covered / length, "ambiguous_fragments": 0})
        legacy.write_table(folder / "gene_metrics.tsv", rows, legacy.METRIC_FIELDS)
        contract = self.frozen["mapping_contract"]
        run = {"sample": sample, "status": "COMPLETE", "skip_qc": True,
               "commands": [{"stage": "synthetic_fixture_not_sequencing_data", "returncode": 0}],
               "versions": {"bowtie2": "bowtie2 version " + contract["bowtie2_version"]},
               "raw_hash_policy": "synthetic fixture; no raw data",
               "outputs": [{"path": "gene_metrics.tsv", "sha256": legacy.sha256(folder / "gene_metrics.tsv")}]}
        for field in ("code", "panel_metadata", "panel_index_files"):
            run[field] = [{"path": name, "sha256": digest} for name, digest in contract[field].items()]
        (folder / "run_manifest.json").write_text(json.dumps(run), encoding="utf8")
        (folder / "COMPLETE").write_text("synthetic fixture\n")
        return folder

    def update_manifest(self, transform):
        path = self.root / "metrics/TOY/run_manifest.json"
        data = json.loads(path.read_text(encoding="utf8"))
        transform(data)
        path.write_text(json.dumps(data), encoding="utf8")

    def apply(self, samples=None, out=None):
        return identity.apply_rules(self.root / "metrics", samples or ["TOY"], self.rules,
                                    out or self.root / "calls")

    def test_frozen_sources_and_rules_are_exact(self):
        self.assertEqual(identity.verify_engine(), identity.ENGINE_SHA256)
        self.assertEqual(legacy.sha256(self.rules), identity.FROZEN_RULES_SHA256)

    def test_threshold_boundary_and_output_units(self):
        receipt = self.apply()
        self.assertEqual(receipt["status"], "COMPLETE")
        result = json.loads((self.root / "calls/identity_calls.json").read_text(encoding="utf8"))
        self.assertEqual([r["identity_state"] for r in result["calls"]], ["detected", "insufficient_evidence"])
        self.assertEqual(result["calls"][0]["supported_markers"], 20)
        self.assertEqual(result["calls"][0]["eligible_markers"], 93)
        self.assertIn("paired_fragment_union", result["coverage_unit"])
        self.assertEqual((self.root / "calls/COMPLETE").read_text().strip(), legacy.sha256(self.root / "calls/receipt.json"))

    def test_below_threshold_is_not_absence(self):
        self.make_metrics("LOW", 19)
        self.apply(["LOW"])
        rows = legacy.read_table(self.root / "calls/identity_calls.tsv")
        self.assertEqual({r["identity_state"] for r in rows}, {"insufficient_evidence"})
        self.assertEqual(rows[0]["supported_markers"], "19")

    def test_missing_complete_rejected_before_output(self):
        (self.root / "metrics/TOY/COMPLETE").unlink()
        with self.assertRaisesRegex(ValueError, "missing COMPLETE"):
            self.apply()
        self.assertFalse((self.root / "calls").exists())

    def test_duplicate_ids_and_unsafe_paths_rejected(self):
        for samples in (["TOY", "TOY"], ["../TOY"], ["TOY\t"], "TOY"):
            with self.assertRaises(ValueError):
                self.apply(samples)

    def test_modified_metrics_hash_rejected(self):
        path = self.root / "metrics/TOY/gene_metrics.tsv"
        with path.open("a", encoding="utf8") as handle:
            handle.write("\n")
        with self.assertRaisesRegex(ValueError, "metrics SHA"):
            self.apply()

    def test_modified_rule_and_rule_schema_rejected(self):
        path = self.root / "rules.json"
        data = json.loads(self.rules.read_text(encoding="utf8"))
        data["selected_rule"]["supported_identity_markers"] = 19
        path.write_text(json.dumps(data), encoding="utf8")
        with self.assertRaisesRegex(ValueError, "modified frozen rules"):
            identity.load_rules(path)
        data["version"] = "2.0"
        path.write_text(json.dumps(data), encoding="utf8")
        with self.assertRaisesRegex(ValueError, "schema/version"):
            identity.load_rules(path)

    def test_duplicate_json_keys_rejected(self):
        path = self.root / "metrics/TOY/run_manifest.json"
        body = path.read_text(encoding="utf8")
        path.write_text('{"status":"FAILED",' + body[1:], encoding="utf8")
        with self.assertRaisesRegex(ValueError, "Duplicate JSON key"):
            self.apply()

    def test_failed_command_or_incomplete_status_rejected(self):
        self.update_manifest(lambda run: run["commands"][0].update(returncode=1))
        with self.assertRaisesRegex(ValueError, "incomplete/failed"):
            self.apply()

    def test_provenance_drift_rejected(self):
        self.update_manifest(lambda run: run["panel_index_files"][0].update(sha256="0" * 64))
        with self.assertRaisesRegex(ValueError, "frozen contract"):
            self.apply()

    def test_duplicate_output_ids_rejected(self):
        self.update_manifest(lambda run: run["outputs"].append(dict(run["outputs"][0])))
        with self.assertRaisesRegex(ValueError, "duplicate or missing output"):
            self.apply()

    def test_duplicate_feature_margin_rejected(self):
        path = self.root / "metrics/TOY/gene_metrics.tsv"
        with path.open("r", encoding="utf8") as handle:
            lines = handle.readlines()
        with path.open("a", encoding="utf8") as handle:
            handle.write(lines[1])
        with self.assertRaisesRegex(ValueError, "duplicate feature x margin"):
            self.apply()

    def test_existing_destination_is_preserved(self):
        self.apply()
        original = legacy.sha256(self.root / "calls/receipt.json")
        with self.assertRaises(FileExistsError):
            self.apply()
        self.assertEqual(original, legacy.sha256(self.root / "calls/receipt.json"))

    def test_write_failure_never_publishes_complete(self):
        with patch.object(identity, "_write_json", side_effect=OSError("disk full")):
            with self.assertRaisesRegex(OSError, "disk full"):
                self.apply()
        self.assertFalse((self.root / "calls").exists())
        self.assertFalse(any(self.root.glob(".calls.running-*/COMPLETE")))

    def mock_engine_output(self, names):
        def run(r1, r2, sample, panel, stage, threads, skip_qc, bin_dir):
            source = self.root / "metrics/TOY"
            for name in ("run_manifest.json", "gene_metrics.tsv", "COMPLETE"):
                shutil.copyfile(source / name, stage / name)
            with gzip.open(stage / "competitive.sam.gz", "wb") as handle:
                for name in names:
                    handle.write(name.encode() + b"\t99\tx\t1\t42\t4M\t=\t10\t0\tACGT\tIIII\n")
            return json.loads((stage / "run_manifest.json").read_text(encoding="utf8"))
        return run

    def wrapper_run(self, names):
        r1, r2 = self.root / "r1.fq", self.root / "r2.fq"
        r1.write_bytes(b"@a/1\nACGT\n+\nIIII\n@b/1\nACGT\n+\nIIII\n")
        r2.write_bytes(b"@a/2\nACGT\n+\nIIII\n@b/2\nACGT\n+\nIIII\n")
        with patch.object(identity, "_verify_panel"), patch.object(identity.engine, "run", side_effect=self.mock_engine_output(names)):
            return identity.run_sample(r1, r2, "TOY", self.root / "fake_panel", self.root / "published", skip_qc=True)

    def test_wrapper_publishes_only_after_full_input_and_sam_checks(self):
        result = self.wrapper_run(["a", "a", "b", "b"])
        self.assertEqual(result["adapter_receipt"]["raw_input_validation"]["pairs_n"], 2)
        self.assertTrue((self.root / "published/COMPLETE").is_file())
        self.assertFalse(list(self.root.glob(".published.running-*")))
        receipt = json.loads((self.root / "published/identity_adapter_receipt.json").read_text(encoding="utf8"))
        self.assertEqual(receipt["published_directory"], str((self.root / "published").resolve()))
        self.assertNotEqual(receipt["execution_directory"], receipt["published_directory"])

    def test_wrapper_rejects_non_contiguous_sam_after_engine_completion(self):
        with self.assertRaisesRegex(ValueError, "Non-contiguous repeated SAM QNAME"):
            self.wrapper_run(["a", "b", "a"])
        self.assertFalse((self.root / "published").exists())
        failed = list(self.root.glob(".published.running-*/ADAPTER_FAILED.json"))
        self.assertEqual(len(failed), 1)
        self.assertFalse((failed[0].parent / "COMPLETE").exists())


class PairedInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def reads(self, names, mate):
        path = self.root / f"R{mate}.fq.gz"
        content = b"".join(b"@" + name.encode() + f"/{mate}\n".encode() + b"ACGT\n+\nIIII\n" for name in names)
        with gzip.open(path, "wb") as handle:
            handle.write(content)
        return path, content

    def test_complete_input_with_exact_decompressed_hashes(self):
        r1, a = self.reads(["a", "b"], 1)
        r2, b = self.reads(["a", "b"], 2)
        result = identity.validate_paired_fastq(r1, r2, self.root)
        self.assertEqual(result["pairs_n"], 2)
        self.assertEqual(result["decompressed_sha256"], [hashlib.sha256(x).hexdigest() for x in (a, b)])
        self.assertFalse(list(self.root.glob("qname-*")))

    def test_duplicate_and_nonconsecutive_names_rejected(self):
        r1, _ = self.reads(["a", "b", "a"], 1)
        r2, _ = self.reads(["a", "b", "a"], 2)
        with self.assertRaisesRegex(ValueError, "Duplicate FASTQ QNAME"):
            identity.validate_paired_fastq(r1, r2, self.root)
        self.assertFalse(list(self.root.glob("qname-*")))

    def test_mate_mismatch_and_different_eof_rejected(self):
        r1, _ = self.reads(["a", "b"], 1)
        r2, _ = self.reads(["a", "c"], 2)
        with self.assertRaisesRegex(ValueError, "QNAME/mate mismatch"):
            identity.validate_paired_fastq(r1, r2, self.root)
        r2, _ = self.reads(["a"], 2)
        with self.assertRaisesRegex(ValueError, "different record counts"):
            identity.validate_paired_fastq(r1, r2, self.root)

    def test_truncated_fastq_is_rejected(self):
        r1 = self.root / "R1.fq"
        r1.write_bytes(b"@a/1\nACGT\n+\n")
        r2, _ = self.reads(["a"], 2)
        with self.assertRaisesRegex(ValueError, "truncated"):
            identity.validate_paired_fastq(r1, r2, self.root)

    def test_illumina_mate_headers_and_quality_are_checked(self):
        r1, r2 = self.root / "R1.fq", self.root / "R2.fq"
        r1.write_bytes(b"@a 2:N:0:INDEX\nACGT\n+\nIIII\n")
        r2.write_bytes(b"@a 1:N:0:INDEX\nACGT\n+\nIIII\n")
        with self.assertRaisesRegex(ValueError, "QNAME/mate mismatch"):
            identity.validate_paired_fastq(r1, r2, self.root)
        r1.write_bytes(b"@a 1:N:0:INDEX\nACGT\n+\nI II\n")
        r2.write_bytes(b"@a 2:N:0:INDEX\nACGT\n+\nIIII\n")
        with self.assertRaisesRegex(ValueError, "characters"):
            identity.validate_paired_fastq(r1, r2, self.root)

    def test_sam_requires_contiguous_groups(self):
        path = self.root / "test.sam"
        records = lambda names: b"".join(name.encode() + b"\t99\tx\t1\t42\t4M\t=\t10\t0\tACGT\tIIII\n" for name in names)
        path.write_bytes(records(["b", "b", "a", "a"]))
        self.assertEqual(identity.validate_sam_qname_groups(path, self.root)["qname_groups_n"], 2)
        path.write_bytes(records(["a", "b", "a"]))
        with self.assertRaisesRegex(ValueError, "Non-contiguous repeated SAM QNAME"):
            identity.validate_sam_qname_groups(path, self.root)


if __name__ == "__main__":
    unittest.main()
