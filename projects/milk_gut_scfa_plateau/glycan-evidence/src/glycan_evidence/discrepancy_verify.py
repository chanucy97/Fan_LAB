"""Data-only replay of sealed complete-QC marker-discrepancy witnesses.

No input Python is imported and no optimizer is called. The floating witness
arithmetic follows frozen audit_actual.audit_fit (FTOL=1e-12, GTOL=1e-8).
This verifies recorded measurements, not the current bytes of large FASTQs/SAMs.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import stat

from . import __version__

FTOL = 1e-12
GTOL = 1e-8
MARKERS = 220
BUCKETS = 16
MAX_INTEGER = 2**53
MAX_FILE_BYTES = 32 * 1024**2
MAX_TOTAL_BYTES = 128 * 1024**2
SCHEMA = "glycan-evidence.discrepancy-verification.v1"
_CODE_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
RESULT_FILES = frozenset({"QUERY.json", "MARKER_RECOUNT.json", "FIT_full.json",
    "FIT_background.json", "FIT_local.json", "MARKER_RESIDUALS.json", "DESCRIPTIVE_RESULT.json"})
REQUIRED_FROZEN = frozenset({"FREEZE.json", "PREFLIGHT.json", "bridge.py", "diagnostic.py",
    "verify_stage2.py", "frozen/audit/audit_actual.py", "frozen/model/MANIFEST.json",
    "frozen/model/training_profiles.json", "frozen/model/ordered_markers.json",
    "frozen/model/mapping_contract.json", "frozen/model/parameters.json"})
NULL_FIELDS = ("identity_call", "p_value", "confidence_interval", "organism_abundance", "full_q", "local_q")


def _need(condition, message):
    if not condition:
        raise ValueError(message)


def _integer(value, name, minimum=0, maximum=MAX_INTEGER):
    _need(type(value) is int and minimum <= value <= maximum, "Invalid integer: " + name)
    return value


def _number(value, name, minimum=0):
    _need(type(value) in (int, float) and math.isfinite(value) and value >= minimum,
          "Invalid finite number: " + name)
    return value


def _close(actual, saved, name):
    _number(saved, name, -math.inf)
    _need(math.isfinite(actual) and abs(actual-saved) <= FTOL*max(1., abs(actual), abs(saved)),
          "Arithmetic mismatch: " + name)


def _sha(value):
    _need(type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value), "Invalid SHA256")
    return value


def _hashes(value, name):
    _need(type(value) is list and len(value) == 2, "Two read hashes required: " + name)
    for item in value:
        _sha(item)
    _need(value[0] != value[1], "Duplicate read hashes: " + name)
    return value


def _ints(value, length, name, maximum=MAX_INTEGER):
    _need(type(value) is list and len(value) == length, "Invalid vector shape: " + name)
    for item in value:
        _integer(item, name, maximum=maximum)
    _need(sum(value) <= MAX_INTEGER, "Integer sum outside exact range: " + name)
    return value


def _numbers(value, length, name, minimum=0):
    _need(type(value) is list and len(value) == length, "Invalid vector shape: " + name)
    for item in value:
        _number(item, name, minimum)
    return value


def _relative(name):
    _need(type(name) is str and name and not any(c in name for c in "\\:\x00\r\n"), "Unsafe path")
    parts = name.split("/")
    _need(not name.startswith("/") and all(p not in ("", ".", "..") and p == p.rstrip(" .")
          and not re.fullmatch(r"(?i)(con|prn|aux|nul|com[1-9]|lpt[1-9])(?:\..*)?", p)
          for p in parts), "Unsafe or noncanonical path: " + name)
    return parts


def _no_links(path):
    for part in [path, *path.parents]:
        if part.exists() or part.is_symlink():
            info = part.lstat()
            _need(not part.is_symlink() and not (getattr(info, "st_file_attributes", 0) & 0x400),
                  "Symlink/reparse path is not allowed: " + str(part))


def _external_path(path):
    path = Path(path).absolute()
    _need(".." not in path.parts, "Parent traversal in path")
    _no_links(path)
    return path.resolve()


def _unique_pairs(pairs):
    result = {}
    for key, value in pairs:
        _need(key not in result, "Duplicate JSON key: " + key)
        result[key] = value
    return result


def _finite_json(value):
    if type(value) is float:
        _need(math.isfinite(value), "Non-finite JSON number")
    elif type(value) is dict:
        for item in value.values():
            _finite_json(item)
    elif type(value) is list:
        for item in value:
            _finite_json(item)


class _Inputs:
    """Read bounded files once; bind parsed bytes and recheck before COMPLETE."""
    def __init__(self, root):
        self.root = root
        self.bytes = {}
        self.records = {}

    def load(self, name, expected=None):
        parts = _relative(name)
        path = self.root.joinpath(*parts)
        _no_links(path)
        _need(path.resolve().is_relative_to(self.root), "Input escaped root")
        _need(path.is_file() and stat.S_ISREG(path.stat().st_mode), "Missing regular input: " + name)
        if name not in self.bytes:
            _need(path.stat().st_size <= MAX_FILE_BYTES, "Input exceeds small-file limit: " + name)
            with path.open("rb") as handle:
                data = handle.read(MAX_FILE_BYTES + 1)
            _need(len(data) <= MAX_FILE_BYTES and sum(map(len, self.bytes.values())) + len(data) <= MAX_TOTAL_BYTES,
                  "Input exceeds small-file budget")
            self.bytes[name] = data
            self.records[name] = {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
        if expected is not None:
            _need(self.records[name]["sha256"] == _sha(expected), "Hash mismatch: " + name)
        return self.bytes[name]

    def json(self, name, expected=None):
        value = json.loads(self.load(name, expected).decode("utf-8"), object_pairs_hook=_unique_pairs)
        _finite_json(value)
        return value

    def manifest(self, files, prefix="", required=frozenset(), exact=False):
        _need(type(files) is dict and 0 < len(files) <= 128, "Invalid file manifest")
        _need(set(files) == required if exact else required <= set(files), "Missing or unrecognized required manifest members")
        names, identities = set(), set()
        for name, expected in files.items():
            _relative(name)
            _need(name.casefold() not in names, "Duplicate/case-aliased manifest path")
            names.add(name.casefold())
            self.load(prefix + name, expected)
            info = (self.root/(prefix + name)).stat()
            identity = (info.st_dev, info.st_ino)
            _need(not info.st_ino or identity not in identities, "Duplicate physical manifest path")
            identities.add(identity)

    def recheck(self):
        for name, saved in self.records.items():
            path = self.root/name
            _no_links(path)
            _need(path.is_file() and path.stat().st_size == saved["bytes"], "Input changed during verification: " + name)
            with path.open("rb") as handle:
                value = hashlib.file_digest(handle, "sha256").hexdigest()
            _need(value == saved["sha256"], "Input changed during verification: " + name)


def _markers(value):
    _need(type(value) is list and len(value) == MARKERS and
          all(type(v) is str and v for v in value) and len(set(value)) == MARKERS,
          "Complete ordered 220-marker universe required")
    return value


def _command(commands, stage):
    matches = [c for c in commands if c.get("stage") == stage]
    _need(len(matches) == 1, "Missing or duplicate command: " + stage)
    return _argv(matches[0]["argv"])


def _argv(value):
    _need(type(value) is list and value and all(type(v) is str and v and "\x00" not in v for v in value),
          "Malformed recorded argv")
    return value


def _qc_command(value):
    result = list(_argv(value))
    for flag in ("--json", "--html"):
        _need(result.count(flag) == 1 and result.index(flag)+1 < len(result), "Missing/duplicate QC report flag")
        result[result.index(flag)+1] = "REPORT_PATH"
    _need("--stdout" in result, "Complete QC stdout contract missing")
    return result


def _file_records(value, count, name):
    _need(type(value) is list and len(value) == count, "Invalid file records: " + name)
    seen = set()
    for row in value:
        path = row["path"]
        _need(type(path) is str and path and not any(x in path for x in ("\x00", "\n", "\r")), "Unsafe recorded path")
        normalized = path.replace("\\", "/")
        parts = normalized.split("/")
        _need((normalized.startswith("/") or re.match(r"^[A-Za-z]:/", normalized)) and
              all(p not in ("", ".", "..") and p == p.rstrip(" .") for p in parts[1:]), "Unsafe recorded path")
        _need(normalized.casefold() not in seen, "Duplicate recorded file path")
        seen.add(normalized.casefold())
        _integer(row["bytes"], name + " bytes", 1)
        if "mtime_ns" in row:
            _integer(row["mtime_ns"], name + " mtime_ns", maximum=2**63-1)
        if "sha256" in row:
            _sha(row["sha256"])
    return value


def _bridge(inputs, prefix, sample, freeze, bridge_freeze, preflight, mapping, query, recount, result):
    bridge = inputs.json(prefix + "RECEIPT.json")
    bridge_sha = inputs.records[prefix + "RECEIPT.json"]["sha256"]
    _need(inputs.load(prefix + "QC_BRIDGE_PASS").decode("ascii").strip() == bridge_sha
          == result["input_bridge_receipt_sha256"], "QC completion/descriptor binding mismatch")
    _need(bridge["sample"] == sample and bridge["status"] == "QC_BRIDGE_PASS_MAPPING_RECOUNT_PENDING",
          "Incomplete/wrong-sample QC bridge")
    _need(bridge["freeze_sha256"] == inputs.records["FREEZE.json"]["sha256"], "QC freeze binding mismatch")
    _need(bridge["new_H_calls"] == _integer(bridge["new_H_calls"], "bridge H") == 0 and
          _integer(bridge["new_envelope_fits"], "bridge envelope") == 0, "QC bridge has unexpected fitting")
    entry = preflight["samples"][sample]
    old = entry["old_files"]["run_manifest.json"]["data"]
    asm = entry["assembly"]["receipt.json"]
    _need(bridge["assembly_receipt_sha256"] == bridge_freeze["assembly_receipt_sha256"][sample]
          == asm["sha256"], "Assembly receipt binding mismatch")
    _sha(asm["sha256"])
    raw = bridge["reused_raw_validation"]
    _need(raw == asm["data"]["raw_input_validation"] and raw["both_eof"] is True and
          raw["unique_qnames"] is True and raw["hash_unit"] == "complete_decompressed_FASTQ_bytes",
          "Complete raw validation record missing")
    raw_n = _integer(raw["pairs_n"], "raw pairs", 1)
    _hashes(raw["decompressed_sha256"], "raw")
    _file_records(raw["input_stats"], 2, "raw")
    _need(raw["input_stats"] == entry["raw_current_stat"], "Preflight/raw-stat binding mismatch")
    qc = bridge["qc"]
    n = _integer(qc["pairs_n"], "post-QC pairs", 1)
    _need(n <= raw_n and qc["both_eof"] is True and
          qc["unit"] == "all_complete_post_QC_paired_fragments_including_no_marker_fragments",
          "Full post-QC input audit missing")
    _need(_integer(qc["buckets"], "QC buckets") == BUCKETS and
          qc["bucket_seed_utf8"] == freeze["bucket_seed"] == bridge_freeze["bucket_seed"], "Bucket contract mismatch")
    _need(_integer(bridge_freeze["buckets"], "frozen QC buckets") == BUCKETS, "QC bucket count mismatch")
    _sha(qc["ordered_canonical_qname_sha256"])
    _hashes(qc["decompressed_sha256"], "post-QC")
    _file_records(qc["files"], 2, "QC")
    for row in qc["files"]:
        _sha(row["sha256"])
    cap = _integer(bridge_freeze["resource_contract"]["maximum_QC_compressed_bytes_per_sample"], "QC bytes cap", 1)
    _need(sum(r["bytes"] for r in qc["files"]) <= cap, "Final recorded QC gzip size exceeds frozen cap")
    _need(query["input_pairs"] == recount["N_input_pairs"] == n ==
          _integer(entry["frozen_row"]["expected_qc_pairs"], "preflight post-QC N", 1),
          "Actual post-QC denominator mismatch")
    _need(query["read_sha256"] == recount["qc_read_sha256"] == qc["decompressed_sha256"], "Post-QC read hash mismatch")
    _need(recount["N_buckets"] == qc["N_buckets"], "QC/recount bucket denominator mismatch")
    _ints(qc["N_buckets"], BUCKETS, "QC buckets", n)
    _need(sum(qc["N_buckets"]) == n, "QC bucket N addback mismatch")
    _need(old["sample"] == sample and old["status"] == "COMPLETE" and old["skip_qc"] is False,
          "Historical complete QC run missing")
    _need(_qc_command(bridge["qc_argv"]) == _qc_command(_command(old["commands"], "fastp")), "QC argv changed")
    _need(_argv(bridge["screen_argv"]) == _command(old["commands"], "feature_screen") and
          bridge["screen_argv"][-2:] == ["--interleaved", "-"], "Screen argv changed")
    for tool, argv in (("fastp", bridge["qc_argv"]), ("bowtie2", bridge["screen_argv"])):
        binary = preflight["bin_versions"][tool]
        _need(binary["path"] == argv[0] and binary["version"] == old["versions"][tool] and
              _integer(binary["returncode"], "tool version returncode") == 0, "Frozen tool/version mismatch")
        _sha(binary["sha256"])
    historical_fastp = entry["old_files"]["fastp.json"]["data"]["summary"]
    _need(_integer(historical_fastp["before_filtering"]["total_reads"], "raw total reads") == 2*raw_n and
          _integer(historical_fastp["after_filtering"]["total_reads"], "QC total reads") == 2*n,
          "Raw/post-QC fastp denominator mismatch")
    screen = bridge["screen"]
    _need(screen == old["screen"] == entry["old_files"]["screen_counts.json"]["data"], "Full historical screen differs")
    for key, value in screen.items():
        _integer(value, "screen " + key)
    candidates = screen["candidate_fragments"]
    _need(screen["total_fragments"] == n and candidates <= n and
          screen["both_mates_mapped"] + screen["one_mate_mapped"] == candidates and
          candidates + screen["unmapped_fragments"] == n, "Screen denominator mismatch")
    checks = bridge["candidate_content_checks"]
    _need(type(checks) is list and len(checks) == 2, "Candidate content checks missing")
    for mate, check in enumerate(checks, 1):
        _need(_integer(check["mate"], "candidate mate") == mate and
              check["original"] == entry["old_files"][f"candidate_R{mate}.fastq.gz"]["path"], "Candidate source mismatch")
        _sha(check["decompressed_sha256"])
        _sha(check["new_compressed_sha256"])
    _need(len({c["replayed"] for c in checks}) == 2 and len({c["decompressed_sha256"] for c in checks}) == 2,
          "Duplicate candidate content records")
    sam = bridge["competitive_sam"]
    _file_records([sam], 1, "competitive SAM")
    historical_sam = entry["old_files"]["competitive.sam.gz"]
    _need(sam["path"] == historical_sam["path"] and sam["bytes"] ==
          _integer(historical_sam["bytes"], "historical SAM bytes", 1), "SAM source binding mismatch")
    _sha(sam["sha256"])
    groups = _integer(recount["sam_primary_pairs_content_verified"], "SAM paired content")
    _need(_integer(recount["candidate_pairs"], "candidate pairs") == candidates and groups <= candidates and
          groups == _integer(recount["sam_counters"]["total_fragments"], "SAM groups") ==
          _integer(entry["old_files"]["summary.json"]["data"]["total_fragments"], "historical SAM groups"),
          "Candidate/SAM groups mismatch")
    _need("--no-unal" in _command(old["commands"], "competitive_genome_mapping"), "Unmapped-candidate omission contract missing")
    for name, expected in mapping["run_panel_metadata_sha256"].items():
        records = [r for r in old["panel_metadata"] if PurePosixPath(r["path"]).name == name]
        _need(len(records) == 1 and records[0]["sha256"] == _sha(expected), "Training/measurement panel mismatch")
    for name in ("run_sample.py", "map_evidence.py"):
        inputs.load("frozen/legacy/" + name, mapping["engine_sha256"][name])
    return {"raw_pairs": raw_n, "post_QC_pairs": n, "candidate_pairs": candidates,
            "mapped_QNAME_groups": groups, "candidate_groups_omitted_by_no_unal": candidates-groups,
            "post_QC_decompressed_sha256": qc["decompressed_sha256"], "receipt_sha256": bridge_sha,
            "large_file_bytes_rehashed": False,
            "scope": "hash-bound complete QC and candidate-content audit records; FASTQ/SAM bytes not reread"}


def _audit_fit(row, query, refs, indices, label, np):
    """Frozen floating certificate replay, with strict types before NumPy coercion."""
    _need(row["status"] == "CERTIFIED" and row["label"] == label and row["query_id"] == query["query_id"], "Fit status/label/query mismatch")
    _need(row["reference_source_ids"] == [r["source_id"] for r in refs], "Fit reference order mismatch")
    _ints(row["marker_indices"], len(indices), "fit marker indices", MARKERS-1)
    _need(row["marker_indices"] == indices, "Fit marker indices mismatch")
    _need(row["read_sha256"] == query["read_sha256"] and
          _integer(row["N_input_pairs"], "fit N", 1) == query["input_pairs"], "Fit read/N binding mismatch")
    canonical = json.dumps(query["counts"], sort_keys=True, separators=(",", ":"), allow_nan=False)
    _need(row["query_counts_sha256"] == hashlib.sha256(canonical.encode()).hexdigest(), "Fit counts hash mismatch")
    _need(row["reused_fit"] is False and _integer(row["new_H_calls"], "fit historical H calls") == 1, "Fit execution kind mismatch")
    n = np.array([query["counts"][i] for i in indices], dtype=float)
    support = sum(query["counts"][i] for i in indices)
    scale = max(support, 1)
    y = n/scale
    matrix = np.array([r["counts"] for r in refs], dtype=float).T[indices, :] / np.array([r["input_pairs"] for r in refs])
    matrix *= query["input_pairs"]/scale
    matrix = np.ascontiguousarray(matrix)
    _need(_integer(row["support"], "fit support") == support, "Fit support mismatch")
    fit = row["fit"]
    _need(fit["status"] == "finite_rate_hellinger_projection", "Uncertified fit status")
    a = np.array(_numbers(fit["source_input_weights"], len(refs), "witness weights"), dtype=float)
    _need(abs(a.sum()-1) <= FTOL, "Witness simplex sum mismatch")
    z = matrix @ a
    _need(np.isfinite(z).all() and (z >= 0).all(), "Invalid fitted prediction")
    saved = _numbers(fit["fitted_scaled_rates"], len(indices), "fitted scaled rates")
    _need(np.allclose(z, saved, rtol=FTOL, atol=FTOL), "Stored prediction mismatch")
    available = matrix.sum(axis=1) > 0
    positive = (y > 0) & available
    _need((z[positive] > 0).all(), "Gradient boundary")
    unsupported = np.flatnonzero((y > 0) & (~available)).tolist()
    unreachable = sum(int(n[i]) for i in unsupported)
    constant = unreachable/scale
    loss = float(np.sum((np.sqrt(z)-np.sqrt(y))**2))
    ratio = np.zeros(len(y)); ratio[positive] = np.sqrt(y[positive]/z[positive])
    gradient = matrix.T @ (1-ratio)
    raw_gap = float(gradient@a-gradient.min()); gap = max(0., raw_gap)
    _need(math.isfinite(raw_gap) and -FTOL <= raw_gap <= GTOL+FTOL, "Convex gap outside frozen tolerance")
    for key, value in (("loss", loss), ("first_order_convex_gap", gap), ("unsupported_loss_constant", constant)):
        _number(fit[key], key)
        _close(value, fit[key], label + " " + key)
    _ints(fit["unsupported_indices"], len(unsupported), "unsupported indices", len(indices)-1)
    _need(fit["unsupported_indices"] == unsupported and
          _integer(row["unsupported_support"], "unsupported support") == unreachable, "Unsupported count/index mismatch")
    _number(row["unsupported_fixed_loss"], "unsupported fixed loss")
    _close(constant, row["unsupported_fixed_loss"], "unsupported fixed loss")
    zeros = [refs[i]["source_id"] for i in np.flatnonzero(matrix.sum(axis=0) == 0)]
    _need(zeros and row["zero_reference_sources"] == zeros, "Genuine zero reference mismatch")
    _need(constant <= loss+FTOL and loss-gap <= (1. if support else 0.)+FTOL, "Analytic minimum bounds failed")
    raw = fit["raw_solver_fit"]
    _need(raw["status"] == "finite_rate_hellinger_projection", "Raw witness status mismatch")
    raw_a = np.array(_numbers(raw["source_input_weights"], len(refs), "raw weights", -FTOL), dtype=float)
    _need(abs(raw_a.sum()-1) <= FTOL, "Raw witness feasibility mismatch")
    feasible = np.maximum(raw_a, 0.); feasible /= feasible.sum()
    _need(np.allclose(a, feasible, rtol=FTOL, atol=FTOL), "Repaired witness mismatch")
    _number(raw["loss"], "raw loss")
    _close(loss, raw["loss"], "raw/repaired loss")
    # Raw weights may be slightly negative within the original repair tolerance.
    # The accepted witness is strictly nonnegative; no new fit is made here.
    return {"label": label, "support": support, "loss_recomputed": loss, "gap_recomputed": gap,
            "unsupported_support": unreachable, "unsupported_fraction": constant,
            "raw_min_weight": float(raw_a.min()), "fitted_scaled_rates": z.tolist(),
            "historical_H_calls": 1, "optimizer_calls_in_verification": 0}


def _residuals(rows, query, recount, fits, z_indices, markers, np):
    _need(type(rows) is list and len(rows) == MARKERS, "Residual table denominator mismatch")
    n = query["input_pairs"]
    counts = np.array(query["counts"], dtype=float)
    bucket_n = np.array(recount["N_buckets"], dtype=float)
    events = np.array(recount["marker_counts_by_bucket"], dtype=float)
    weights = bucket_n/n
    use = bucket_n > 0
    rates = counts/n
    partition = events[use]/bucket_n[use, None]
    # Match original stage2: use saved rates only after independent certification.
    predictions = {label: np.array(row["fit"]["fitted_scaled_rates"])*max(row["support"], 1)/n
                   for label, row in fits.items()}
    predicted = predictions["full"]
    variation = (weights[use, None]*(partition-rates)**2).sum(axis=0)
    common = (rates-predicted)**2
    total = (weights[use, None]*(partition-predicted)**2).sum(axis=0)
    _need(np.allclose(total, common+variation, rtol=FTOL, atol=FTOL), "Weighted decomposition identity failed")
    nonempty = int(use.sum())
    signs = (((partition-predicted)*(rates-predicted)) > 0).sum(axis=0)
    local_positions = {j: i for i, j in enumerate(z_indices)}
    verified = []
    for j, row in enumerate(rows):
        _need(row["marker_id"] == markers[j] and _integer(row["observed_events"], "residual events") == query["counts"][j],
              "Residual marker order/events mismatch")
        _need(type(row["local_marker"]) is bool and row["local_marker"] == (j in local_positions), "Local marker flag mismatch")
        expected = {"observed_rate": float(rates[j]), "fitted_rate": float(predicted[j]),
            "signed_rate_residual": float(rates[j]-predicted[j]), "common_squared_rate_residual": float(common[j]),
            "weighted_partition_variation": float(variation[j]), "weighted_partition_total_residual": float(total[j]),
            "background_fitted_rate": float(predictions["background"][j])}
        for label, key in (("full", "loss_contribution"), ("background", "background_loss_contribution")):
            expected[key] = (math.sqrt(predictions[label][j]*n)-math.sqrt(counts[j]))**2/max(fits[label]["support"], 1)
        if j in local_positions:
            rate = float(predictions["local"][local_positions[j]])
            expected.update(local_fitted_rate=rate,
                local_loss_contribution=(math.sqrt(rate*n)-math.sqrt(counts[j]))**2/max(fits["local"]["support"], 1))
        else:
            _need(row["local_fitted_rate"] is None and row["local_loss_contribution"] is None, "Nonlocal marker has local results")
        for key, value in expected.items():
            _number(row[key], key, -math.inf if key == "signed_rate_residual" else 0)
            _close(value, row[key], markers[j] + " " + key)
        _need(_integer(row["nonempty_buckets"], "nonempty buckets") == nonempty and
              _integer(row["same_sign_nonempty_buckets"], "same-sign buckets") == int(signs[j]), "Partition sign/count mismatch")
        _need(row["interpretation"] == "descriptive_technical_partition_no_sampling_CI", "Partition interpretation mismatch")
        verified.append({"marker_id": markers[j], "observed_events": query["counts"][j], **expected})
    for label, key in (("full", "loss_contribution"), ("background", "background_loss_contribution"), ("local", "local_loss_contribution")):
        _close(sum(row.get(key, 0.) for row in verified), fits[label]["fit"]["loss"], label + " marker-loss addback")
    return verified, float(np.max(abs(total-common-variation)))


def _write(path, value):
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run_discrepancy_verify(root, sample, expected_freeze_sha256, expected_results_sha256, out):
    """Verify one sealed 220-marker/16-bucket result and publish an independent receipt.

    Both SHA256 arguments must come from an external trusted handoff. ``out``
    must not exist. Only bounded files under ``root`` are read, never raw reads
    or executables. Returns the COMPLETE receipt; invalid data raises ValueError
    (missing files/OSError can also propagate). NumPy is needed; SciPy is not.
    """
    try:
        return _run_discrepancy_verify(root, sample, expected_freeze_sha256, expected_results_sha256, out)
    except (KeyError, TypeError, IndexError, OverflowError) as error:
        raise ValueError("Missing or malformed discrepancy field: " + str(error)) from error


def _run_discrepancy_verify(root, sample, expected_freeze_sha256, expected_results_sha256, out):
    _need(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == _CODE_SHA256,
          "Verifier source changed after import; start a fresh process")
    _relative(sample)
    _need("/" not in sample and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", sample), "Unsafe sample name")
    _sha(expected_freeze_sha256); _sha(expected_results_sha256)
    root = _external_path(root)
    destination = _external_path(out)
    _need(root.is_dir(), "Missing source root")
    if destination.exists():
        raise FileExistsError(str(destination))
    _need(not root.is_relative_to(destination) and not destination.is_relative_to(root/"results")
          and not destination.is_relative_to(root/"frozen"), "Output overlaps protected input tree")
    inputs = _Inputs(root)
    freeze = inputs.json("STAGE2_FREEZE.json", expected_freeze_sha256)
    inputs.manifest(freeze["files_sha256"], required=REQUIRED_FROZEN)
    _need(type(freeze["samples"]) is list and sample in freeze["samples"] and
          len(set(freeze["samples"])) == len(freeze["samples"]), "Sample not in unique frozen sample set")
    _need(_integer(freeze["bucket_count"], "frozen bucket count") == BUCKETS and
          type(freeze["bucket_seed"]) is str and freeze["bucket_seed"], "Invalid frozen bucket contract")
    _need(_integer(freeze["maximum_H_calls_per_sample"], "frozen H calls") == 3 and
          _integer(freeze["new_q_fits"], "frozen q fits") == 0 and
          _integer(freeze["simulations"], "frozen simulations") == 0 and
          freeze["state"] == "reference_discrepancy_diagnostic", "Unsupported discrepancy contract")
    for key, value in (("feasibility", FTOL), ("convex_gap", GTOL), ("rate_decomposition_mixed_tolerance", FTOL)):
        _number(freeze["numeric_gate"][key], "frozen tolerance")
        _need(freeze["numeric_gate"][key] == value, "Changed frozen numerical tolerance")
    prefix = "results/" + sample + "/"
    stage = inputs.json(prefix + "STAGE2_RESULTS.json", expected_results_sha256)
    _need(stage["sample"] == sample and stage["stage2_freeze_sha256"] == expected_freeze_sha256 and
          _integer(stage["actual_H_calls"], "historical H calls") == 3 and
          stage["status"] == "AWAITING_INDEPENDENT_WITNESS_REPLAY" and
          stage["technical_gate"] == "AWAITING_INDEPENDENT_WITNESS_REPLAY", "Stage result/freeze binding mismatch")
    inputs.manifest(stage["files_sha256"], prefix, RESULT_FILES, exact=True)
    bridge_freeze = inputs.json("FREEZE.json")
    inputs.manifest(bridge_freeze["files_sha256"], required={"PREFLIGHT.json", "bridge.py"})
    _need(sample in bridge_freeze["samples"], "Sample absent from QC freeze")
    manifest = inputs.json("frozen/model/MANIFEST.json", freeze["model_manifest_sha256"])
    _need(manifest["version"] == "reference-family-compatibility.research.v0.1" and
          manifest["identity_call"] is None and manifest["probability_calibration"] is False,
          "Unsupported frozen reference model")
    inputs.manifest(manifest["files_sha256"], "frozen/model/",
                    {"ordered_markers.json", "training_profiles.json", "mapping_contract.json", "parameters.json"})
    markers = _markers(inputs.json("frozen/model/ordered_markers.json"))
    mapping = inputs.json("frozen/model/mapping_contract.json")
    parameters = inputs.json("frozen/model/parameters.json")
    for key, file in (("training_profiles_sha256", "training_profiles.json"), ("ordered_markers_sha256", "ordered_markers.json"),
                      ("mapping_contract_sha256", "mapping_contract.json")):
        _need(parameters[key] == inputs.records["frozen/model/" + file]["sha256"], "Model parameter binding mismatch")
    query = inputs.json(prefix + "QUERY.json")
    recount = inputs.json(prefix + "MARKER_RECOUNT.json")
    result = inputs.json(prefix + "DESCRIPTIVE_RESULT.json")
    n = _integer(query["input_pairs"], "query N", 1)
    _ints(query["counts"], MARKERS, "query counts", n)
    _hashes(query["read_sha256"], "query")
    _need(query["query_id"] == sample and recount["marker_ids"] == markers, "Query/marker order mismatch")
    _ints(recount["counts"], MARKERS, "recount counts", n)
    _need(query["counts"] == recount["counts"] and _integer(recount["N_input_pairs"], "recount N", 1) == n,
          "Query/recount count/N mismatch")
    _need(_integer(recount["eligible_markers"], "marker count") == MARKERS and
          _integer(recount["sum_support"], "recount support") == sum(query["counts"]) and
          recount["historical_marker_equality"] == "ALL220_COUNT_COVERAGE_DEPTH_AND_AMBIGUITY_EXACT",
          "Complete historical marker recount missing")
    bucket_n = _ints(recount["N_buckets"], BUCKETS, "bucket N", n)
    _need(sum(bucket_n) == n, "Full QC bucket N addback mismatch")
    events = recount["marker_counts_by_bucket"]
    _need(type(events) is list and len(events) == BUCKETS, "Invalid bucket event shape")
    for i, row in enumerate(events):
        _ints(row, MARKERS, "bucket counts", bucket_n[i])
    _need([sum(row[j] for row in events) for j in range(MARKERS)] == query["counts"], "Marker bucket addback mismatch")
    preprocessing = _bridge(inputs, prefix, sample, freeze, bridge_freeze, inputs.json("PREFLIGHT.json"),
                            mapping, query, recount, result)
    refs = inputs.json("frozen/model/training_profiles.json")
    _need(type(refs) is list and refs, "Training profiles missing")
    seen = set()
    for ref in refs:
        source = ref["source_id"]
        _need(type(source) is str and source and source not in seen, "Duplicate/invalid reference source")
        seen.add(source)
        _need(ref["role"] in {"target", "background"} and ref["marker_ids"] == markers, "Reference role/marker order mismatch")
        ref_n = _integer(ref["input_pairs"], "reference N", 1)
        _ints(ref["counts"], MARKERS, "reference counts", ref_n)
        _hashes(ref["read_sha256"], "reference")
        _need(not set(query["read_sha256"]) & set(ref["read_sha256"]), "Query/training read leakage")
        _need(ref["mapping_contract_sha256"] == mapping["historical_mapping_contract_sha256"], "Reference mapping contract mismatch")
    bg = [r for r in refs if r["role"] == "background"]
    targets = [r for r in refs if r["role"] == "target"]
    _need(bg and targets, "Both reference roles required")
    z = freeze["Z_indices"]
    _ints(z, len(z), "local indices", MARKERS-1)
    _need(z and z == [j for j in range(MARKERS) if all(r["counts"][j] == 0 for r in bg)], "Local coordinate set mismatch")
    import numpy as np
    fits, verified = {}, []
    for label, references, indices in (("full", refs, list(range(MARKERS))),
            ("background", bg, list(range(MARKERS))), ("local", targets+[bg[0]], z)):
        fits[label] = inputs.json(prefix + "FIT_" + label + ".json")
        verified.append(_audit_fit(fits[label], query, references, indices, label, np))
    full, back = fits["full"]["fit"], fits["background"]["fit"]
    _need(back["loss"]+back["first_order_convex_gap"]+GTOL >= full["loss"]-full["first_order_convex_gap"], "Nested loss reversal")
    residuals, identity_error = _residuals(inputs.json(prefix + "MARKER_RESIDUALS.json"), query, recount, fits, z, markers, np)
    _need(result["sample"] == sample and result["state"] == "reference_discrepancy_diagnostic" and
          _integer(result["N_input_pairs"], "descriptive N", 1) == n and
          _integer(result["S_full"], "full support") == fits["full"]["support"] and
          _integer(result["S_local"], "local support") == fits["local"]["support"], "Descriptive denominator mismatch")
    _need(all(key in result and result[key] is None for key in NULL_FIELDS), "Classification/q/probability/abundance must remain null")
    _need(_integer(result["new_H_calls"], "historical H budget") == 3 and
          _integer(result["new_q_fits"], "historical q budget") == 0, "Descriptive fitting budget mismatch")
    for label in fits:
        _number(result["D_" + label], "descriptive loss")
        _close(fits[label]["fit"]["loss"], result["D_" + label], "descriptive " + label)
    _close(back["loss"]-full["loss"], result["D_background_minus_full"], "descriptive loss difference")
    for label, domain in (("full", parameters["support_domain"]), ("local", [19, 2545])):
        _ints(domain, 2, "support domain")
        _need(domain[0] <= domain[1] and result[label+"_empirical_support_domain"] == domain and
              type(result[label+"_outside_domain"]) is bool and result[label+"_outside_domain"] ==
              (not domain[0] <= fits[label]["support"] <= domain[1]), "Support domain annotation mismatch")
    inputs.recheck()
    receipt = {"schema": SCHEMA, "status": "COMPLETE", "technical_gate": "PASS", "sample": sample,
        "package_version": __version__, "verified_utc": datetime.now(timezone.utc).isoformat(),
        "stage2_freeze_sha256": expected_freeze_sha256, "stage2_results_sha256": expected_results_sha256,
        "inputs": inputs.records, "preprocessing": preprocessing, "fit_witnesses": verified,
        "marker_count": MARKERS, "bucket_count": BUCKETS, "max_weighted_identity_error": identity_error,
        "tolerances": {"feasibility": FTOL, "convex_gap": GTOL, "rate_arithmetic": FTOL},
        "optimizer_calls": 0, "new_H_calls": 0, "new_q_fits": 0, "historical_H_calls_replayed": 3,
        "numpy_version": np.__version__, "verifier_sha256": _CODE_SHA256,
        "coverage_unit": "complete_post_QC_paired_QNAME_fragments; multi-marker event counts",
        "scope": "saved floating certificate and measured-input record verification; descriptive reference discrepancy; no new scientific validation",
        **{key: None for key in NULL_FIELDS}}
    _no_links(destination)
    destination.mkdir(parents=True, exist_ok=False)
    try:
        _write(destination/"replayed_fits.json", verified)
        _write(destination/"marker_residuals.json", residuals)
        receipt["files_sha256"] = {name: hashlib.sha256((destination/name).read_bytes()).hexdigest()
                                  for name in ("replayed_fits.json", "marker_residuals.json")}
        inputs.recheck()
        _need(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == _CODE_SHA256,
              "Verifier source changed during verification")
        _write(destination/"receipt.json", receipt)
        with (destination/"COMPLETE").open("x", encoding="ascii") as handle:
            handle.write(hashlib.sha256((destination/"receipt.json").read_bytes()).hexdigest() + "\n")
    except BaseException:
        (destination/"COMPLETE").unlink(missing_ok=True)
        raise
    return receipt
