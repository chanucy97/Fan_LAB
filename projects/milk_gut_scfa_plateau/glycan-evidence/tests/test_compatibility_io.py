"""Bridge contracts on explicitly synthetic receipts; no numerical performance claims."""
import hashlib
import json
from pathlib import Path
import subprocess
import unittest
from unittest.mock import patch

from glycan_evidence import compatibility_io as bridge, identity
from glycan_evidence.integrity import digest, seal_outputs
import test_identity


class CompatibilityIOTests(unittest.TestCase):
    def setUp(self):
        fixture = test_identity.IdentityTests(methodName="test_frozen_sources_and_rules_are_exact")
        fixture.setUp()
        self.addCleanup(fixture.doCleanups)
        self.root = fixture.root
        self.run_root = self.root
        self.folder = self.root / "metrics/TOY"
        self.metrics = self.folder / "gene_metrics.tsv"
        self.out = self.root / "prediction"
        self.bundle = self.root / "source_bundle"
        self.bundle.mkdir()
        self.markers = sorted(f["feature_id"] for f in fixture.features if
            f["channel"] == "identity" and f["eligible_identity"] == "1" and f["pseudo"] == "0")
        mapping = dict(engine_sha256=identity.ENGINE_SHA256,
            historical_mapping_contract_sha256=fixture.frozen["mapping_contract_sha256"],
            expected_python_version="fixture Python; synthetic only",
            expected_bowtie2_version="fixture bowtie2 version " + fixture.frozen["mapping_contract"]["bowtie2_version"])
        for name in bridge.PAYLOAD:
            self.save(self.bundle / name, {} if name.endswith(".json") else "synthetic nonexecuted stub")
        self.save(self.bundle / "ordered_markers.json", self.markers)
        self.save(self.bundle / "mapping_contract.json", mapping)
        self.seal_bundle()
        run = self.load(self.folder / "run_manifest.json")
        run["versions"] = {"python": mapping["expected_python_version"], "bowtie2": mapping["expected_bowtie2_version"]}
        stats = [{"path": f"/synthetic/R{i}.fastq", "bytes": 123, "mtime_ns": i} for i in (1, 2)]
        run["raw_inputs"] = stats
        run["screen"] = {"total_fragments": 50, "candidate_fragments": 20}
        self.save(self.folder / "screen_counts.json", run["screen"])
        self.save(self.folder / "summary.json", {"sample": "TOY", "total_fragments": 20})
        self.save(self.folder / "run_manifest.json", run)
        adapter = {"schema": identity.SCHEMA, "status": "COMPLETE", "sample": "TOY",
            "engine_sha256": identity.ENGINE_SHA256, "adapter_sha256": digest(identity.__file__),
            "raw_input_validation": {"pairs_n": 50, "both_eof": True, "unique_qnames": True,
                "hash_unit": "complete_decompressed_FASTQ_bytes", "decompressed_sha256": ["a"*64, "b"*64],
                "input_stats": stats},
            "sam_group_validation": {"qname_groups_n": 20, "non_contiguous_qnames": 0}}
        self.save(self.folder / "identity_adapter_receipt.json", adapter)
        (self.root / "COMPLETE").write_text("synthetic fixture\n")
        self.reseal_run()

    @staticmethod
    def save(path, value):
        path.write_text(json.dumps(value), encoding="utf-8")

    @staticmethod
    def load(path):
        return json.loads(path.read_text(encoding="utf-8"))

    def seal_bundle(self):
        self.save(self.bundle / "MANIFEST.json", {"version": bridge.RELEASE_VERSION,
            "files_sha256": {n: digest(self.bundle / n) for n in sorted(bridge.PAYLOAD)}})
        self.expected = digest(self.bundle / "MANIFEST.json")

    def reseal_run(self):
        run = self.load(self.folder / "run_manifest.json")
        run["outputs"] = [{"path": n, "sha256": digest(self.folder / n)}
                          for n in ("gene_metrics.tsv", "screen_counts.json", "summary.json")]
        self.save(self.folder / "run_manifest.json", run)
        p = self.folder / "identity_adapter_receipt.json"
        if p.exists():
            adapter = self.load(p)
            adapter["run_manifest_sha256"] = digest(self.folder / "run_manifest.json")
            self.save(p, adapter)
        (self.root / "artifact_manifest.json").unlink(missing_ok=True)
        # Keep the synthetic source bundle out of the input run's seal, mirroring separate releases.
        seal_outputs(self.root)

    def change_adapter(self, change):
        p = self.folder / "identity_adapter_receipt.json"
        value = self.load(p)
        change(value)
        self.save(p, value)
        self.reseal_run()

    def fake_api(self, argv, **kwargs):
        self.assertEqual(argv[1:6], ["-E", "-s", "-B", "-X", "utf8"])
        snapshot = Path(argv[argv.index("--bundle") + 1])
        self.assertNotEqual(snapshot, self.bundle)
        self.assertFalse((snapshot / "numpy.py").exists())
        self.assertFalse((snapshot / "__pycache__").exists())
        query = self.load(Path(argv[argv.index("--query") + 1]))
        result = dict(version=bridge.RELEASE_VERSION, package_sha256=self.expected,
            input_binding_sha256=hashlib.sha256(bridge._canonical(query).encode()).hexdigest(),
            input_pairs=query["input_pairs"], composite_support=sum(query["counts"]),
            relation="background_compatible", identity_call=None, p_value=None,
            confidence_interval=None, organism_abundance=None)
        self.save(Path(argv[argv.index("--out") + 1]), result)
        return subprocess.CompletedProcess(argv, 0)

    def run_bridge(self, side_effect=None):
        with patch.object(bridge.subprocess, "run", side_effect=side_effect or self.fake_api) as execute:
            result = bridge.run_compatibility(self.metrics, self.bundle, self.expected, self.out)
        return result, execute

    def rejected(self, pattern):
        with patch.object(bridge.subprocess, "run") as execute:
            with self.assertRaisesRegex(ValueError, pattern):
                bridge.run_compatibility(self.metrics, self.bundle, self.expected, self.out)
            execute.assert_not_called()
        self.assertFalse(self.out.exists())

    def test_receipts_generate_complete_query_without_manual_denominators(self):
        result, _ = self.run_bridge()
        query = self.load(self.out / "query.json")
        self.assertEqual(query["input_pairs"], 50)
        self.assertEqual(len(query["counts"]), 220)
        self.assertEqual(sum(query["counts"]), 20)
        self.assertEqual(query["read_sha256"], ["a"*64, "b"*64])
        self.assertEqual(result["status"], "COMPLETE")
        self.assertIsNone(result["identity_call"])
        self.assertEqual((self.out / "COMPLETE").read_text().strip(), digest(self.out / "receipt.json"))

    def test_wrong_actual_N_cannot_be_substituted(self):
        self.change_adapter(lambda x: x["raw_input_validation"].update(pairs_n=49))
        self.rejected("differs from the no-QC screen denominator")

    def test_bool_N_incomplete_read_audit_and_equal_hashes_rejected(self):
        for key, value in (("pairs_n", True), ("both_eof", False),
                           ("unique_qnames", False), ("decompressed_sha256", ["a"*64, "a"*64])):
            original = self.load(self.folder / "identity_adapter_receipt.json")
            with self.subTest(key=key):
                self.change_adapter(lambda x: x["raw_input_validation"].update({key: value}))
                self.rejected("audited|audit|required")
            self.save(self.folder / "identity_adapter_receipt.json", original)
            self.reseal_run()

    def test_QC_rejected_even_with_unchanged_mapping_digest(self):
        path = self.folder / "run_manifest.json"
        run = self.load(path)
        run["skip_qc"] = False
        self.save(path, run)
        self.reseal_run()
        self.rejected("requires recorded skip_qc=true")

    def test_missing_grid_row_rejected_after_valid_reseal(self):
        lines = self.metrics.read_text().splitlines(True)
        self.metrics.write_text("".join(lines[:-1]), encoding="utf-8")
        self.reseal_run()
        self.rejected("incomplete or unexpected feature")

    def test_missing_adapter_rejected(self):
        (self.folder / "identity_adapter_receipt.json").unlink()
        self.reseal_run()
        self.rejected("Run seal must bind")

    def test_contaminated_bundle_is_rejected_before_any_execution(self):
        (self.bundle / "release_api.py").write_text("raise RuntimeError('must never execute')")
        self.rejected("Release payload changed")

    def test_missing_marker_is_rejected_before_any_execution(self):
        self.save(self.bundle / "ordered_markers.json", self.markers[:-1])
        self.seal_bundle()
        self.rejected("ordered220")

    def test_untrusted_extra_code_and_cache_are_not_copied(self):
        (self.bundle / "numpy.py").write_text("raise RuntimeError('must never execute')")
        (self.bundle / "__pycache__").mkdir()
        (self.bundle / "__pycache__/release_api.cpython-312.pyc").write_bytes(b"untrusted")
        self.reseal_run()
        self.run_bridge()

    def test_failed_subprocess_preserves_failure_without_COMPLETE(self):
        with self.assertRaisesRegex(RuntimeError, "Frozen compatibility API failed"):
            self.run_bridge(lambda argv, **kwargs: subprocess.CompletedProcess(argv, 2))
        self.assertFalse((self.out / "COMPLETE").exists())
        self.assertEqual(self.load(self.out / "receipt.json")["status"], "FAILED")

    def test_existing_output_is_preserved(self):
        self.out.mkdir()
        with self.assertRaises(FileExistsError):
            self.run_bridge()


if __name__ == "__main__":
    unittest.main()
