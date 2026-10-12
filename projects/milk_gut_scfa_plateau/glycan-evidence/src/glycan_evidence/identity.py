"""Strict adapter for the byte-preserved, registered v2 identity engine.

Only the registered pilot rules are accepted. No rule selection, absence call,
functional inference, or reference replacement takes place in this module.
"""
from __future__ import annotations

import csv
from contextlib import closing
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import tempfile
import time
from uuid import uuid4

from .legacy_v2 import calibrate_calls as calls
from .legacy_v2 import run_sample as engine

SCHEMA = "glycan-evidence.identity.v1"
ENGINE_SHA256 = {
    "run_sample.py": "cc92ec525246f2447a734028068fad512af6203a85e528c04a18917f4f033c4c",
    "map_evidence.py": "ca081a202c66fefd0441ae0051c7d6d4fe11a0e8854ec66dddb2a9a696cb7afc",
    "calibrate_calls.py": "92e6b871e173cfe91a5ac29481d50ec71fb5582834745e958242d9fcc255f3be",
}
FROZEN_RULES_SHA256 = "2c2f9650d932b7439873067caf3aad58d8066e6ea86f6228953e7573360dcd25"
DEFAULT_RULES = Path(__file__).with_name("legacy_v2") / "frozen_rules.json"
COVERAGE_UNIT = "reference_feature_bases; depth=paired_fragment_union_per_base"
SCOPE = calls.SCOPE
FIELDS = ["sample", "target", "identity_state", "reason", "supported_markers",
          "eligible_markers", "metrics_available", "flag_scope_panel", *calls.RULE_FIELDS,
          "coverage_unit", "scope"]


def _json(path):
    def unique(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding="utf-8-sig"), object_pairs_hook=unique)


def _write_json(path, data):
    with Path(path).open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(data, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")


def _sample_id(sample):
    if not isinstance(sample, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
        raise ValueError(f"Unsafe or empty sample ID: {sample!r}")
    return sample


def verify_engine():
    """Reject accidental edits to any byte-preserved engine source."""
    folder = Path(engine.__file__).parent
    for name, expected in ENGINE_SHA256.items():
        if calls.sha256(folder / name) != expected:
            raise ValueError(f"Frozen engine hash mismatch: {name}")
    return dict(ENGINE_SHA256)


def load_rules(rules=DEFAULT_RULES):
    """Load the explicitly registered pilot rule artifact, never recalibrate it."""
    verify_engine()
    path = Path(rules)
    data = _json(path)
    if not isinstance(data, dict) or data.get("version") != "1.0.0":
        raise ValueError("Unsupported frozen rule schema/version")
    if calls.sha256(path) != FROZEN_RULES_SHA256:
        raise ValueError("Unregistered or modified frozen rules SHA-256")
    frozen = calls.load_rules(path)
    digest = hashlib.sha256(json.dumps(frozen["mapping_contract"], sort_keys=True,
                                      separators=(",", ":")).encode()).hexdigest()
    if digest != frozen["mapping_contract_sha256"]:
        raise ValueError("Frozen mapping contract hash mismatch")
    return frozen


def verify_completed_metrics(metrics_root, sample, rules=DEFAULT_RULES):
    """Check one complete metrics grid and its frozen provenance; no output writes.

    Historical compact tables retain their original raw-input stat policy. This
    check does not claim to have rehashed or re-read historical FASTQ/SAM files.
    """
    sample = _sample_id(sample)
    frozen = load_rules(rules)
    path = calls.metrics_path(metrics_root, sample)
    run = _json(path.parent / "run_manifest.json")
    if not isinstance(run, dict) or not isinstance(run.get("commands"), list) or not run["commands"]:
        raise ValueError(f"{sample}: missing command provenance")
    if not isinstance(run.get("skip_qc"), bool):
        raise ValueError(f"{sample}: missing boolean QC status")
    output_names = [str(item.get("path", "")).replace("\\", "/").rsplit("/", 1)[-1]
                    for item in run.get("outputs", [])]
    if not output_names or len(output_names) != len(set(output_names)):
        raise ValueError(f"{sample}: duplicate or missing output provenance IDs")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if set(header) != set(calls.METRIC_FIELDS):
            raise ValueError(f"{sample}: unexpected metrics column schema")
        if any(row and len(row) != len(header) for row in reader):
            raise ValueError(f"{sample}: inconsistent metrics row width")
    if len(header) != len(set(header)):
        raise ValueError(f"{sample}: duplicate metrics column names")
    data = calls.load_metrics(path, sample, frozen["metrics_margins"])
    if data["schema_sha256"] != frozen["feature_schema_sha256"]:
        raise ValueError(f"{sample}: feature dictionary differs from frozen schema")
    provenance = calls.pipeline_provenance(path, sample, data["sha256"])
    if provenance["sha256"] != frozen["mapping_contract_sha256"]:
        raise ValueError(f"{sample}: mapping panel/code/version differs from frozen contract")
    for record in run["outputs"]:
        name = str(record.get("path", "")).replace("\\", "/").rsplit("/", 1)[-1]
        candidate = path.parent / name
        if record.get("sha256") and candidate.is_file() and calls.sha256(candidate) != record["sha256"]:
            raise ValueError(f"{sample}: retained output hash mismatch: {name}")
    provenance["raw_input_verification"] = run.get("raw_hash_policy", "not_recorded")
    provenance["historical_qname_uniqueness_rechecked"] = False
    adapter_path = path.parent / "identity_adapter_receipt.json"
    if adapter_path.exists():
        adapter = _json(adapter_path)
        if (adapter.get("schema") != SCHEMA or adapter.get("status") != "COMPLETE"
                or adapter.get("sample") != sample
                or adapter.get("run_manifest_sha256") != calls.sha256(path.parent / "run_manifest.json")
                or adapter.get("engine_sha256") != ENGINE_SHA256):
            raise ValueError(f"{sample}: invalid mapping adapter receipt")
        provenance["adapter_receipt_sha256"] = calls.sha256(adapter_path)
        provenance["recorded_raw_input_validation"] = adapter["raw_input_validation"]
        provenance["recorded_sam_group_validation"] = adapter["sam_group_validation"]
    return {"data": data, "provenance": provenance, "rules": frozen}


def _new_stage(out):
    destination = Path(out).resolve()
    if destination.exists():
        raise FileExistsError(f"Output path must not already exist: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = destination.with_name(f".{destination.name}.running-{uuid4().hex}")
    stage.mkdir()
    return destination, stage


def _failed(stage, error):
    # Preserve all diagnostic files; a failed attempt is never a completed result.
    marker = stage / "COMPLETE"
    if marker.exists():
        marker.unlink()
    _write_json(stage / "ADAPTER_FAILED.json", {"schema": SCHEMA, "status": "FAILED",
                "error_type": type(error).__name__, "error": str(error)})


def apply_rules(metrics_root, samples, rules, out):
    """Apply the registered rules to explicit completed sample IDs.

    Input errors leave no destination. Output I/O errors retain a hidden failed
    staging directory for diagnosis. Existing outputs are never overwritten.
    """
    if isinstance(samples, (str, bytes)):
        raise ValueError("samples must be a sequence of distinct sample IDs")
    samples = [_sample_id(sample) for sample in samples]
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Empty or duplicate sample IDs")
    frozen = load_rules(rules)
    rows, provenance = [], {}
    for sample in samples:
        checked = verify_completed_metrics(metrics_root, sample, rules)
        provenance[sample] = checked["provenance"]
        for target in frozen["targets"]:
            result = calls.identity_result(checked["data"], target, frozen["selected_rule"])
            rows.append({"sample": sample, "target": target,
                         "identity_state": "detected" if result["detected"] else "insufficient_evidence",
                         "reason": "frozen_identity_rule_passed" if result["detected"] else "below_frozen_identity_rule",
                         "supported_markers": result["supported_markers"], "eligible_markers": result["eligible_markers"],
                         "metrics_available": 1, "flag_scope_panel": True, **frozen["selected_rule"],
                         "coverage_unit": COVERAGE_UNIT, "scope": SCOPE})
    destination, stage = _new_stage(out)
    try:
        calls.write_table(stage / "identity_calls.tsv", rows, FIELDS)
        _write_json(stage / "identity_calls.json", {"schema": SCHEMA, "scope": SCOPE,
                    "coverage_unit": COVERAGE_UNIT, "calls": rows})
        receipt = {"schema": SCHEMA, "status": "COMPLETE", "samples_n": len(samples),
                   "calls_n": len(rows), "samples": samples, "engine_sha256": verify_engine(),
                   "adapter_sha256": calls.sha256(__file__), "frozen_rules_sha256": calls.sha256(rules),
                   "selected_rule": frozen["selected_rule"], "scope": SCOPE, "coverage_unit": COVERAGE_UNIT,
                   "nondetection_policy": "insufficient_evidence; no absence call",
                   "input_provenance": provenance,
                   "outputs": [{"name": name, "sha256": calls.sha256(stage / name),
                                "bytes": (stage / name).stat().st_size}
                               for name in ("identity_calls.tsv", "identity_calls.json")]}
        _write_json(stage / "receipt.json", receipt)
        (stage / "COMPLETE").write_text(calls.sha256(stage / "receipt.json") + "\n", encoding="ascii")
        stage.rename(destination)
        return receipt
    except BaseException as error:
        if stage.exists():
            _failed(stage, error)
        raise


def _open_fastq(path):
    return gzip.open(path, "rb") if str(path).endswith(".gz") else Path(path).open("rb")


def _read_fastq(handle, path, pair_no):
    header = handle.readline()
    if not header:
        return None
    sequence, plus, quality = handle.readline(), handle.readline(), handle.readline()
    if not all((sequence, plus, quality)) or not header.startswith(b"@") or not plus.startswith(b"+"):
        raise ValueError(f"Malformed or truncated FASTQ at pair {pair_no}: {path}")
    sequence_b, quality_b = sequence.rstrip(b"\r\n"), quality.rstrip(b"\r\n")
    if not sequence_b or len(sequence_b) != len(quality_b):
        raise ValueError(f"FASTQ sequence/quality length mismatch at pair {pair_no}: {path}")
    if any(value < 33 or value > 126 for value in quality_b) or any(chr(value).isspace() for value in sequence_b):
        raise ValueError(f"Invalid FASTQ sequence/quality characters at pair {pair_no}: {path}")
    header_parts = header[1:].split()
    token = header_parts[0] if header_parts else b""
    if not token:
        raise ValueError(f"Empty FASTQ QNAME at pair {pair_no}: {path}")
    suffix = int(token[-1:]) if token.endswith((b"/1", b"/2")) else None
    name = token[:-2] if suffix else token
    if len(header_parts) > 1 and header_parts[1].startswith((b"1:", b"2:")):
        illumina_mate = int(header_parts[1][:1])
        if suffix is not None and suffix != illumina_mate:
            raise ValueError(f"Conflicting FASTQ mate labels at pair {pair_no}: {path}")
        suffix = illumina_mate
    return name, suffix, header + sequence + plus + quality


def validate_paired_fastq(r1, r2, work_dir):
    """Full-stream pair/EOF/unique-QNAME check with disk-backed exact uniqueness.

    Hashes cover decompressed FASTQ bytes. They are explicitly distinct from
    compressed-file SHA-256 and from legacy file-stat-only provenance.
    """
    started = time.monotonic()
    paths = [Path(r1).resolve(), Path(r2).resolve()]
    if paths[0] == paths[1] or os.path.samefile(*paths):
        raise ValueError("R1 and R2 must be distinct files")
    stats = [engine.file_record(path) for path in paths]
    digests = [hashlib.sha256(), hashlib.sha256()]
    count = 0
    with tempfile.TemporaryDirectory(prefix="qname-", dir=work_dir) as temporary:
        with closing(sqlite3.connect(str(Path(temporary) / "names.sqlite"))) as database:
            database.execute("CREATE TABLE names (name BLOB PRIMARY KEY) WITHOUT ROWID")
            with _open_fastq(paths[0]) as left, _open_fastq(paths[1]) as right:
                while True:
                    a, b = _read_fastq(left, paths[0], count + 1), _read_fastq(right, paths[1], count + 1)
                    if a is None and b is None:
                        break
                    if a is None or b is None:
                        raise ValueError("Paired FASTQs have different record counts")
                    if a[0] != b[0] or a[1] not in (None, 1) or b[1] not in (None, 2):
                        raise ValueError(f"Paired FASTQ QNAME/mate mismatch at pair {count + 1}")
                    try:
                        database.execute("INSERT INTO names VALUES (?)", (a[0],))
                    except sqlite3.IntegrityError as error:
                        raise ValueError(f"Duplicate FASTQ QNAME at pair {count + 1}") from error
                    for digest, record in zip(digests, (a, b)):
                        digest.update(record[2])
                    count += 1
    if not count:
        raise ValueError("Empty paired FASTQ input")
    if stats != [engine.file_record(path) for path in paths]:
        raise ValueError("FASTQ stat changed during full validation")
    return {"pairs_n": count, "both_eof": True, "unique_qnames": True,
            "elapsed_seconds": time.monotonic() - started,
            "scope": "all raw paired FASTQ records before QC; exact disk-backed QNAME uniqueness",
            "input_stats": stats, "hash_unit": "complete_decompressed_FASTQ_bytes",
            "decompressed_sha256": [digest.hexdigest() for digest in digests]}


def validate_sam_qname_groups(path, work_dir):
    """Reject a QNAME reappearing after its group ended, including unmapped SAM."""
    count, current = 0, None
    with tempfile.TemporaryDirectory(prefix="sam-qname-", dir=work_dir) as temporary:
        database = sqlite3.connect(str(Path(temporary) / "names.sqlite"))
        try:
            database.execute("CREATE TABLE names (name BLOB PRIMARY KEY) WITHOUT ROWID")
            with _open_fastq(path) as handle:
                for number, line in enumerate(handle, 1):
                    if line.startswith(b"@") or not line.strip():
                        continue
                    fields = line.rstrip(b"\r\n").split(b"\t")
                    if len(fields) < 11 or not fields[0]:
                        raise ValueError(f"Malformed SAM line {number}")
                    if fields[0] != current:
                        try:
                            database.execute("INSERT INTO names VALUES (?)", (fields[0],))
                        except sqlite3.IntegrityError as error:
                            raise ValueError(f"Non-contiguous repeated SAM QNAME at line {number}") from error
                        current = fields[0]
                        count += 1
        finally:
            database.close()
    return {"qname_groups_n": count, "non_contiguous_qnames": 0}


def _verify_panel(panel, frozen):
    """Validate the fixed reference's metadata and frozen index-cache contract."""
    panel = Path(panel).resolve()
    for name, expected in frozen["mapping_contract"]["panel_metadata"].items():
        if calls.sha256(panel / name) != expected:
            raise ValueError(f"Panel metadata differs from frozen v2: {name}")
    indices, _ = engine.index_records(panel)
    if {Path(row["path"]).name: row["sha256"] for row in indices} != frozen["mapping_contract"]["panel_index_files"]:
        raise ValueError("Panel index hashes differ from frozen v2")


def run_sample(r1, r2, sample, panel, out, threads=1, skip_qc=False, bin_dir=None):
    """Run the unchanged v2 engine, adding validation and atomic publication.

    Production sequencing execution is intended for Linux. Never pass altered
    references to the registered v2 preset; new presets require new validation.
    """
    sample = _sample_id(sample)
    if type(threads) is not int or threads < 1 or not isinstance(skip_qc, bool):
        raise ValueError("threads must be a positive integer and skip_qc a boolean")
    frozen = load_rules()
    panel = Path(panel).resolve()
    _verify_panel(panel, frozen)
    destination, stage = _new_stage(out)
    try:
        inputs = validate_paired_fastq(r1, r2, stage)
        # The engine's completion marker remains hidden inside staging until
        # the wrapper's stricter checks pass and the directory is published.
        run = engine.run(r1, r2, sample, panel, stage, threads, skip_qc, bin_dir)
        groups = validate_sam_qname_groups(stage / "competitive.sam.gz", stage)
        if inputs["input_stats"] != [engine.file_record(path) for path in (r1, r2)]:
            raise ValueError("Raw inputs changed during mapping")
        data = calls.load_metrics(stage / "gene_metrics.tsv", sample, frozen["metrics_margins"])
        provenance = calls.pipeline_provenance(stage / "gene_metrics.tsv", sample, data["sha256"])
        if data["schema_sha256"] != frozen["feature_schema_sha256"] or provenance["sha256"] != frozen["mapping_contract_sha256"]:
            raise ValueError("Completed mapping differs from frozen schema or contract")
        wrapper = {"schema": SCHEMA, "status": "COMPLETE", "sample": sample, "scope": SCOPE,
                   "coverage_unit": COVERAGE_UNIT, "engine_sha256": verify_engine(),
                   "adapter_sha256": calls.sha256(__file__), "raw_input_validation": inputs,
                   "sam_group_validation": groups, "execution_directory": str(stage),
                   "published_directory": str(destination),
                   "run_manifest_sha256": calls.sha256(stage / "run_manifest.json")}
        _write_json(stage / "identity_adapter_receipt.json", wrapper)
        stage.rename(destination)
        return {**run, "adapter_receipt": wrapper}
    except BaseException as error:
        if stage.exists():
            _failed(stage, error)
        raise
