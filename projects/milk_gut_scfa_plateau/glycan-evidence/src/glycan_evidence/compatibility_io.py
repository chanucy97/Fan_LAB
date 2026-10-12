"""Receipt-bound bridge to the unchanged reference-family research release."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from . import __version__, identity
from .integrity import digest, verify_outputs, write_json

SCHEMA = "glycan-evidence.compatibility-run.v1"
RELEASE_VERSION = "reference-family-compatibility.research.v0.1"
PAYLOAD = frozenset({
    "release_api.py", "frozen_rate_hellinger.py", "frozen_poisson_rate_projection.py",
    "compatibility_core.py", "support_envelope.py", "ordered_markers.json",
    "training_profiles.json", "mapping_contract.json", "parameters.json",
    "calibration78.json", "calibration_binding.json",
})
RELATIONS = frozenset({"target_required_given_reference_family", "background_compatible",
    "reference_family_incompatible", "zero_support", "outside_working_support_domain",
    "numerically_unresolved"})


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _integer(value, minimum=0):
    return type(value) is int and minimum <= value <= 2**53


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _record(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest(path)}


def _verified_bundle(bundle, expected_manifest_sha256):
    """Read data only. No bundle Python is imported before all bytes are verified."""
    bundle = Path(bundle).resolve()
    _require(_sha(expected_manifest_sha256), "An external expected manifest SHA256 is required")
    _require(digest(bundle / "MANIFEST.json") == expected_manifest_sha256, "Release manifest changed")
    manifest = identity._json(bundle / "MANIFEST.json")
    _require(manifest.get("version") == RELEASE_VERSION, "Unsupported compatibility release")
    files = manifest.get("files_sha256")
    _require(type(files) is dict and set(files) == PAYLOAD, "Incomplete or unrecognized release payload")
    for name, expected in files.items():
        path = bundle / name
        _require(_sha(expected) and path.is_file() and not path.is_symlink()
                 and path.resolve().parent == bundle and digest(path) == expected,
                 "Release payload changed: " + name)
    markers = identity._json(bundle / "ordered_markers.json")
    _require(type(markers) is list and len(markers) == 220 and
             all(isinstance(x, str) and x for x in markers) and len(set(markers)) == 220,
             "The complete ordered220 marker universe is required")
    return bundle, manifest, markers, identity._json(bundle / "mapping_contract.json")


def _measured_query(metrics, bundle, manifest, markers, mapping):
    path = Path(metrics).resolve()
    _require(path.name == "gene_metrics.tsv", "Expected a complete-run gene_metrics.tsv")
    folder, sample = path.parent, path.parent.name
    _require(folder.parent.name == "metrics", "Expected sealed run/metrics/sample/gene_metrics.tsv layout")
    run_root = folder.parent.parent
    verify_outputs(run_root)
    seal = identity._json(run_root / "artifact_manifest.json")
    names = {r["path"] for r in seal["files"]}
    needed = ["gene_metrics.tsv", "identity_adapter_receipt.json", "run_manifest.json",
              "screen_counts.json", "summary.json", "COMPLETE"]
    _require({f"metrics/{sample}/{name}" for name in needed} <= names,
             "Run seal must bind metrics, adapter, manifest, screen, summary and completion")
    verified = identity.verify_completed_metrics(folder.parent, sample)
    run = identity._json(folder / "run_manifest.json")
    adapter = identity._json(folder / "identity_adapter_receipt.json")
    _require(adapter.get("adapter_sha256") == digest(identity.__file__), "Unrecognized identity adapter implementation")
    _require(adapter.get("engine_sha256") == mapping.get("engine_sha256"), "Release engine binding mismatch")
    _require(verified["provenance"]["sha256"] == mapping.get("historical_mapping_contract_sha256"),
             "Release historical mapping contract mismatch")
    # QC is NOT included in the historical mapper digest. Verify it separately.
    _require(run.get("skip_qc") is True, "This release bridge requires recorded skip_qc=true; QC-transformed counts are outside its preprocessing contract")
    _require(not any(c.get("stage") == "fastp" for c in run["commands"]), "QC command contradicts skip_qc=true")
    _require(run["versions"].get("python") == mapping.get("expected_python_version") and
             run["versions"].get("bowtie2") == mapping.get("expected_bowtie2_version"), "Release recorded mapper runtime mismatch")
    audit = adapter.get("raw_input_validation", {})
    n = audit.get("pairs_n")
    hashes = audit.get("decompressed_sha256")
    _require(_integer(n, 1), "Actual paired-input N must be a positive audited integer")
    _require(audit.get("both_eof") is True and audit.get("unique_qnames") is True and
             audit.get("hash_unit") == "complete_decompressed_FASTQ_bytes", "Complete paired-input audit required")
    _require(type(hashes) is list and len(hashes) == 2 and all(_sha(x) for x in hashes) and
             len(set(hashes)) == 2, "Two distinct audited decompressed FASTQ hashes required")
    stats = audit.get("input_stats")
    _require(type(stats) is list and len(stats) == 2 and stats == run.get("raw_inputs"), "Raw-input records disagree with mapper input records")
    for item in stats:
        _require(isinstance(item, dict) and isinstance(item.get("path"), str) and item["path"] and
                 _integer(item.get("bytes"), 1) and type(item.get("mtime_ns")) is int,
                 "Malformed audited input file record")
    screen = identity._json(folder / "screen_counts.json")
    summary = identity._json(folder / "summary.json")
    _require(screen == run.get("screen"), "Screen counts disagree with completed mapper manifest")
    _require(type(screen.get("total_fragments")) is int and screen["total_fragments"] == n,
             "Audited input N differs from the no-QC screen denominator")
    candidates = screen.get("candidate_fragments")
    sam = adapter.get("sam_group_validation", {})
    groups = sam.get("qname_groups_n")
    _require(_integer(candidates) and _integer(groups) and groups <= candidates <= n and
             type(sam.get("non_contiguous_qnames")) is int and sam["non_contiguous_qnames"] == 0 and
             summary.get("sample") == sample and type(summary.get("total_fragments")) is int and
             summary["total_fragments"] == groups, "Candidate/SAM audit denominators disagree")
    data = verified["data"]
    selected = {fid for fid, f in data["features"].items() if f["channel"] == "identity"
                and str(f["eligible_identity"]) == "1" and str(f["pseudo"]) == "0"
                and f["taxon_group"] in {"infantis", "longum"}}
    _require(selected == set(markers), "Complete eligible joint220 universe mismatch; never fill missing rows with zero")
    counts = [data["evidence"][(fid, 0)]["fragments"] for fid in markers]
    # The frozen loader represents validated integers as floats; convert exactly.
    _require(all(x == int(x) and 0 <= x <= n for x in counts), "Invalid per-marker fragment counts")
    counts = [int(x) for x in counts]
    _require(sum(counts) <= 2**53, "Composite support exceeds exact integer range")
    query = dict(query_id=sample, marker_ids=markers, counts=counts, input_pairs=n,
                 read_sha256=hashes, mapping_contract_sha256=mapping["historical_mapping_contract_sha256"],
                 release_mapping_sha256=digest(bundle / "mapping_contract.json"))
    records = {name: _record(folder / name) for name in needed}
    records["run_artifact_manifest"] = _record(run_root / "artifact_manifest.json")
    return query, records, {"skip_qc": True, "input_pairs": n, "screen_pairs": screen["total_fragments"],
        "candidate_pairs": candidates, "competitive_qname_groups": groups,
        "raw_input_scope": "verified recorded complete FASTQ audit and sealed mapper provenance; raw FASTQs not reread by this bridge",
        "preprocessing_scope": "no additional QC in the frozen mapping run; original provenance of preprocessed reads is not inferred"}


def run_compatibility(metrics, bundle, expected_manifest_sha256, out):
    """Create a research prediction from sealed metrics; no manual N or read hashes."""
    destination = Path(out).resolve()
    if destination.exists():
        raise FileExistsError(str(destination))
    bundle, manifest, markers, mapping = _verified_bundle(bundle, expected_manifest_sha256)
    query, inputs, preprocessing = _measured_query(metrics, bundle, manifest, markers, mapping)
    destination.mkdir(parents=True, exist_ok=False)
    receipt = dict(schema=SCHEMA, status="RUNNING", package_version=__version__,
        started_utc=datetime.now(timezone.utc).isoformat(), inputs=inputs,
        release_manifest_sha256=expected_manifest_sha256, source_bundle=str(bundle),
        preprocessing=preprocessing, identity_call=None,
        scope="conditional longum reference-family compatibility; not identity, organism abundance or calibrated real-cohort performance",
        code=_record(__file__), commands=[])
    write_json(destination / "receipt.json", receipt)
    try:
        # Execute only a verified copy: ignore source-directory extras and cached bytecode.
        snapshot = destination / "verified_bundle"
        snapshot.mkdir()
        for name in sorted(PAYLOAD | {"MANIFEST.json"}):
            shutil.copyfile(bundle / name, snapshot / name)
        _verified_bundle(snapshot, expected_manifest_sha256)
        write_json(destination / "query.json", query)
        argv = [sys.executable, "-E", "-s", "-B", "-X", "utf8", str(snapshot / "release_api.py"),
                "--bundle", str(snapshot), "--expected-manifest-sha256", expected_manifest_sha256,
                "--query", str(destination / "query.json"), "--out", str(destination / "result.json")]
        with (destination / "api.stdout.log").open("x", encoding="utf-8") as stdout, \
             (destination / "api.stderr.log").open("x", encoding="utf-8") as stderr:
            completed = subprocess.run(argv, cwd=destination, stdout=stdout, stderr=stderr, check=False)
        receipt["commands"].append({"argv": argv, "returncode": completed.returncode})
        if completed.returncode:
            raise RuntimeError(f"Frozen compatibility API failed ({completed.returncode}); inspect api.stderr.log")
        result = identity._json(destination / "result.json")
        _require(result.get("version") == RELEASE_VERSION and result.get("package_sha256") == expected_manifest_sha256
                 and result.get("input_binding_sha256") == hashlib.sha256(_canonical(query).encode()).hexdigest(),
                 "Prediction does not bind the verified query/release")
        _require(result.get("input_pairs") == query["input_pairs"] and result.get("composite_support") == sum(query["counts"])
                 and result.get("relation") in RELATIONS, "Prediction denominator/relationship mismatch")
        _require(all(k in result and result[k] is None for k in
                     ("identity_call", "p_value", "confidence_interval", "organism_abundance")), "Prediction exceeded research scope")
        _verified_bundle(snapshot, expected_manifest_sha256)
        for record in inputs.values():
            _require(_record(record["path"]) == record, "Input report changed during compatibility execution")
        receipt.update(status="COMPLETE", completed_utc=datetime.now(timezone.utc).isoformat(),
                       relation=result["relation"], input_binding_sha256=result["input_binding_sha256"],
                       files={p.relative_to(destination).as_posix(): _record(p) for p in sorted(destination.rglob("*"))
                              if p.is_file() and p.name != "receipt.json"})
        write_json(destination / "receipt.json", receipt)
        (destination / "COMPLETE").write_text(digest(destination / "receipt.json") + "\n", encoding="ascii")
    except BaseException as error:
        receipt.update(status="FAILED", error_type=type(error).__name__, error=str(error))
        (destination / "COMPLETE").unlink(missing_ok=True)
        write_json(destination / "receipt.json", receipt)
        raise
    return receipt
