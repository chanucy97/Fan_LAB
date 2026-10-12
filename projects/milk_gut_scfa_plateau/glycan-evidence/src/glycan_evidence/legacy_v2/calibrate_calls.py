#!/usr/bin/env python3
"""Calibrate identity calls, evaluate frozen rules, or describe cohort evidence.

Reads gene_metrics.tsv only; never reads FASTQ/SAM or tunes on evaluation data.
All non-detections in cohort output remain insufficient_evidence.
"""

import argparse
import csv
from collections import defaultdict
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import sys


VERSION = "1.0.0"
META_FIELDS = ("feature_id", "channel", "taxon_group", "length", "eligible_identity", "pseudo")
METRIC_FIELDS = ("sample", *META_FIELDS, "margin", "fragments", "covered_bases", "breadth", "aligned_bases", "mean_depth", "ambiguous_fragments")
RULE_FIELDS = ("margin", "marker_breadth", "minimum_fragments_per_marker", "supported_identity_markers")
SCOPE = "competitive_reference_panel; no sample-specific detectability or host-background validation"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_table(path, required=(), delimiter="\t"):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not set(required) <= set(reader.fieldnames or ()):
            raise ValueError(f"{path}: missing columns {sorted(set(required) - set(reader.fieldnames or ()))}")
        rows = list(reader)
    return rows


def write_table(path, rows, fields=None):
    if fields is None:
        fields = list(rows[0]) if rows else []
    with Path(path).open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False, allow_nan=False)
        handle.write("\n")


def nonnegative_number(value, label, integer=False):
    try:
        parsed = float(value)
    except (ValueError, TypeError) as exc:
        raise ValueError(f"Invalid {label}: {value!r}") from exc
    if not math.isfinite(parsed) or parsed < 0 or (integer and parsed != int(parsed)):
        raise ValueError(f"Invalid {label}: {value!r}")
    return int(parsed) if integer else parsed


def metadata(row):
    result = {k: row[k] for k in META_FIELDS}
    for key in ("length", "eligible_identity", "pseudo"):
        result[key] = nonnegative_number(row[key], key, integer=True)
    if not result["feature_id"] or result["channel"] not in {"identity", "locus"}:
        raise ValueError("Invalid feature ID or channel")
    if result["length"] == 0 or result["eligible_identity"] not in {0, 1} or result["pseudo"] not in {0, 1}:
        raise ValueError("Invalid feature length or binary eligibility/pseudo flag")
    return result


def schema_hash(features):
    content = [features[key] for key in sorted(features)]
    return hashlib.sha256(json.dumps(content, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def pipeline_provenance(path, sample, metrics_sha):
    folder = Path(path).parent
    manifest = folder / "run_manifest.json"
    if not (folder / "COMPLETE").is_file() or not manifest.is_file():
        raise ValueError(f"{sample}: missing COMPLETE/run_manifest.json from the mapping pipeline")
    run = json.loads(manifest.read_text(encoding="utf-8-sig"))
    if run.get("status") != "COMPLETE" or run.get("sample") != sample:
        raise ValueError(f"{sample}: incomplete or mismatched mapping run")
    if any(command.get("returncode") != 0 for command in run.get("commands", [])):
        raise ValueError(f"{sample}: mapping commands have an incomplete/failed status")

    def hash_records(field, required=()):
        records = {}
        for item in run.get(field, []):
            name = str(item.get("path", "")).replace("\\", "/").rsplit("/", 1)[-1]
            digest = item.get("sha256", "")
            if not name or name in records or not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
                raise ValueError(f"{sample}: missing/invalid {field} SHA-256 records")
            records[name] = digest.lower()
        if not records or not set(required) <= set(records):
            raise ValueError(f"{sample}: incomplete {field} provenance")
        return records

    # Large retained SAM output intentionally lacks SHA; gene_metrics must have one.
    output_records = {str(item.get("path", "")).replace("\\", "/").rsplit("/", 1)[-1]: item for item in run.get("outputs", [])}
    if output_records.get("gene_metrics.tsv", {}).get("sha256") != metrics_sha:
        raise ValueError(f"{sample}: metrics SHA does not match completed pipeline output")
    indices = hash_records("panel_index_files")
    if not all(any(name.startswith(prefix + ".") for name in indices) for prefix in ("features", "genomes")):
        raise ValueError(f"{sample}: both feature and genome index provenance are required")
    version = run.get("versions", {}).get("bowtie2")
    if not version:
        raise ValueError(f"{sample}: Bowtie2 version is missing")
    match = re.search(r"\bversion\s+(\S+)", version)
    contract = {"panel_metadata": hash_records("panel_metadata", ("panel_contigs.tsv", "feature_manifest.tsv", "annotations.tsv")),
                "panel_index_files": indices, "code": hash_records("code", ("run_sample.py", "map_evidence.py")),
                "bowtie2_version": match.group(1) if match else version}
    digest = hashlib.sha256(json.dumps(contract, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return {"contract": contract, "sha256": digest, "run_manifest_sha256": sha256(manifest), "skip_qc": run.get("skip_qc")}


def metrics_path(root, sample):
    # Case IDs/sample IDs name one direct child directory, never arbitrary paths.
    if not sample or sample in {".", ".."} or any(x in sample for x in ("/", "\\", ":")):
        raise ValueError(f"Unsafe sample/case ID: {sample!r}")
    root = Path(root).resolve()
    path = (root / sample / "gene_metrics.tsv").resolve()
    if not path.is_relative_to(root):
        raise ValueError("Metrics path escapes metrics root")
    return path


def load_metrics(path, sample, margins):
    rows = read_table(path, METRIC_FIELDS)
    features, evidence = {}, {}
    for row in rows:
        if row["sample"] != sample:
            raise ValueError(f"{path}: sample field disagrees with manifest ID {sample}")
        feature = metadata(row)
        feature_id = feature["feature_id"]
        if feature_id in features and features[feature_id] != feature:
            raise ValueError(f"{path}: inconsistent feature metadata")
        features[feature_id] = feature
        parsed = {**feature, "sample": sample}
        for key in ("margin", "breadth", "mean_depth"):
            parsed[key] = nonnegative_number(row[key], key)
        for key in ("fragments", "covered_bases", "aligned_bases", "ambiguous_fragments"):
            parsed[key] = nonnegative_number(row[key], key, integer=True)
        if parsed["breadth"] > 1 or parsed["covered_bases"] > feature["length"]:
            raise ValueError(f"{path}: breadth exceeds feature bounds")
        if not math.isclose(parsed["breadth"], parsed["covered_bases"] / feature["length"], abs_tol=1e-9):
            raise ValueError(f"{path}: breadth and covered bases disagree")
        if not math.isclose(parsed["mean_depth"], parsed["aligned_bases"] / feature["length"], abs_tol=1e-9):
            raise ValueError(f"{path}: mean depth and aligned bases disagree")
        if parsed["aligned_bases"] < parsed["covered_bases"]:
            raise ValueError(f"{path}: aligned bases smaller than covered bases")
        key = feature_id, parsed["margin"]
        if key in evidence:
            raise ValueError(f"{path}: duplicate feature x margin")
        evidence[key] = parsed
    if not features:
        raise ValueError(f"{path}: no feature rows")
    expected = {(f, float(m)) for f in features for m in margins}
    if set(evidence) != expected:
        raise ValueError(f"{path}: incomplete or unexpected feature x margin grid")
    return {"features": features, "evidence": evidence, "schema_sha256": schema_hash(features), "sha256": sha256(path)}


def load_benchmark_manifest(path, required_split):
    rows = read_table(path, ("case_id", "split", "source_split", "truth_taxon", "accession", "coverage", "error_rate", "truth_mode", "source_sha256"))
    if not rows:
        raise ValueError("Empty benchmark manifest")
    if len({r["case_id"] for r in rows}) != len(rows):
        raise ValueError("Duplicate benchmark case IDs")
    accession_truth, accession_hash, conditions = {}, {}, set()
    for row in rows:
        if row["split"] != required_split:
            raise ValueError(f"{required_split} mode refuses source split {row['split']}; supply a separate manifest")
        if row["source_split"] != row["split"]:
            raise ValueError("Benchmark split disagrees with original source_split")
        if row["truth_mode"] != "synthetic_reads_from_real_reference_genome":
            raise ValueError("Identity benchmark refuses engineered or unrecognized truth mode")
        if not row["truth_taxon"] or not row["accession"]:
            raise ValueError("Missing benchmark truth taxon or source accession")
        previous = accession_truth.setdefault(row["accession"], row["truth_taxon"])
        if previous != row["truth_taxon"]:
            raise ValueError("One source accession has conflicting truth taxon labels")
        if not re.fullmatch(r"[0-9a-fA-F]{64}", row["source_sha256"]):
            raise ValueError("Benchmark source SHA-256 is missing or malformed")
        row["source_sha256"] = row["source_sha256"].lower()
        if accession_hash.setdefault(row["accession"], row["source_sha256"]) != row["source_sha256"]:
            raise ValueError("One source accession has conflicting source hashes")
        row["coverage"] = nonnegative_number(row["coverage"], "coverage")
        row["error_rate"] = nonnegative_number(row["error_rate"], "error rate")
        if row["coverage"] == 0:
            raise ValueError("Coverage must be positive")
        if row["error_rate"] >= 1:
            raise ValueError("Error rate must be below one")
        condition = row["accession"], row["coverage"], row["error_rate"]
        if condition in conditions:
            raise ValueError("Duplicate accession x coverage x error-rate condition")
        conditions.add(condition)
    return rows


def load_case_metrics(rows, root, margins, expected_schema=None, expected_pipeline=None):
    loaded, hashes = {}, {}
    for row in rows:
        sample = row["case_id"]
        path = metrics_path(root, sample)
        data = load_metrics(path, sample, margins)
        if expected_schema is None:
            expected_schema = data["schema_sha256"]
        if data["schema_sha256"] != expected_schema:
            raise ValueError(f"Feature dictionary changed for case {sample}")
        pipeline = pipeline_provenance(path, sample, data["sha256"])
        if expected_pipeline is None:
            expected_pipeline = pipeline["sha256"]
        if pipeline["sha256"] != expected_pipeline:
            raise ValueError(f"Mapping panel/code/version changed for case {sample}")
        data["pipeline"] = pipeline
        loaded[sample], hashes[sample] = data, {"gene_metrics_sha256": data["sha256"], "run_manifest_sha256": pipeline["run_manifest_sha256"]}
    return loaded, expected_schema, hashes


def read_plan(path):
    plan = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    targets = plan["primary_identity_targets"]
    if targets != ["infantis", "longum"]:
        raise ValueError("This implementation requires frozen targets infantis and longum")
    grid = plan["candidate_detection_grid"]
    values = {"margin": plan["competition_margin_grid"], **{k: grid[k] for k in RULE_FIELDS[1:]}}
    for name, items in values.items():
        if not items or len(items) != len(set(items)):
            raise ValueError(f"Empty or duplicate threshold grid: {name}")
        for item in items:
            number = nonnegative_number(item, name, integer=name in {"minimum_fragments_per_marker", "supported_identity_markers"})
            if name != "margin" and number <= 0:
                raise ValueError("Detection thresholds must be positive")
            if name == "marker_breadth" and number > 1:
                raise ValueError("Marker breadth exceeds one")
    rules = [dict(zip(RULE_FIELDS, combination)) for combination in itertools.product(*(values[k] for k in RULE_FIELDS))]
    return plan, targets, values["margin"], rules


def identity_result(data, target, rule):
    eligible, supported = [], []
    for feature_id, f in data["features"].items():
        if f["channel"] != "identity" or f["taxon_group"] != target or not f["eligible_identity"] or f["pseudo"]:
            continue
        eligible.append(feature_id)
        metric = data["evidence"][(feature_id, float(rule["margin"]))]
        if metric["breadth"] >= rule["marker_breadth"] and metric["fragments"] >= rule["minimum_fragments_per_marker"]:
            supported.append(feature_id)
    return {"detected": len(supported) >= rule["supported_identity_markers"],
            "supported_markers": len(supported), "eligible_markers": len(eligible)}


def predict(rows, loaded, targets, rule):
    predictions = []
    for row in rows:
        for target in targets:
            result = identity_result(loaded[row["case_id"]], target, rule)
            predictions.append({k: row[k] for k in ("case_id", "split", "truth_taxon", "accession", "coverage", "error_rate")} | {
                "target": target, "truth_positive": int(row["truth_taxon"] == target), "detected": int(result["detected"]),
                "identity_state": "detected" if result["detected"] else "insufficient_evidence",
                "supported_markers": result["supported_markers"], "eligible_markers": result["eligible_markers"],
                "flag_scope_panel": True})
    return predictions


def confusion(rows):
    tp = sum(r["truth_positive"] and r["detected"] for r in rows)
    fp = sum(not r["truth_positive"] and r["detected"] for r in rows)
    fn = sum(r["truth_positive"] and not r["detected"] for r in rows)
    tn = sum(not r["truth_positive"] and not r["detected"] for r in rows)
    positives, negatives = tp + fn, fp + tn
    return {"cases_n": len(rows), "source_genomes_n": len({r["accession"] for r in rows}),
            "positive_cases_n": positives, "negative_cases_n": negatives,
            "positive_source_genomes_n": len({r["accession"] for r in rows if r["truth_positive"]}),
            "negative_source_genomes_n": len({r["accession"] for r in rows if not r["truth_positive"]}),
            "TP": tp, "FP": fp, "FN": fn, "TN": tn,
            "sensitivity": tp / positives if positives else "",
            "observed_false_positive_fraction": fp / negatives if negatives else "",
            "insufficient_cases_n": fn + tn,
            "call_fraction": (tp + fp) / len(rows) if rows else ""}


def summaries(predictions, keys):
    groups = defaultdict(list)
    for row in predictions:
        groups[tuple(row[k] for k in keys)].append(row)
    return [{**dict(zip(keys, key)), **confusion(group)} for key, group in sorted(groups.items())]


def select_rule(rows, loaded, targets, rules):
    grid_rows, candidates = [], []
    for target in targets:
        if not any(r["truth_taxon"] == target for r in rows) or not any(r["truth_taxon"] != target for r in rows):
            raise ValueError(f"Calibration needs positive and negative source cases for {target}")
    for index, rule in enumerate(rules, 1):
        predictions = predict(rows, loaded, targets, rule)
        result = {"rule_id": f"rule_{index:03d}", **rule}
        sensitivity, false_positives = Fraction(0), 0
        for target in targets:
            counts = confusion([r for r in predictions if r["target"] == target])
            result.update({f"{target}_{k}": v for k, v in counts.items()})
            sensitivity += Fraction(counts["TP"], counts["positive_cases_n"])
            false_positives += counts["FP"]
        sensitivity /= len(targets)
        result.update({"balanced_target_positive_sensitivity": float(sensitivity), "cross_target_false_detections": false_positives,
                       "zero_observed_FP": int(false_positives == 0), "selected": 0})
        grid_rows.append(result)
        if false_positives == 0:
            key = (sensitivity, rule["supported_identity_markers"], rule["marker_breadth"], rule["minimum_fragments_per_marker"], rule["margin"])
            candidates.append((key, index - 1, rule))
    if not candidates:
        return grid_rows, None, []
    _, selected_index, chosen = max(candidates, key=lambda item: item[0])
    grid_rows[selected_index]["selected"] = 1
    return grid_rows, chosen, predict(rows, loaded, targets, chosen)


def new_output(path):
    out = Path(path).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Output directory must be new or empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    return out


def provenance(mode, manifest, metrics_hashes):
    return {"version": VERSION, "mode": mode, "script_sha256": sha256(__file__), "python_version": sys.version,
            "manifest_sha256": sha256(manifest), "metrics_file_sha256": metrics_hashes,
            "scope": SCOPE, "inference_unit": "source reference genomes; depth/error cases are technical conditions, not additional strains",
            "benchmark_confusion_definition": "Binary detected vs not detected. FN/TN include insufficient-evidence outputs; they are not biological absence calls.",
            "nondetection_policy": "insufficient_evidence; no not_detected call until sample-specific detectability is independently validated"}


def calibrate(args):
    plan, targets, margins, candidates = read_plan(args.freeze_plan)
    rows = load_benchmark_manifest(args.calibration_manifest, "calibration")
    expected_conditions = set(itertools.product(plan["benchmark"]["positive_target_depth_x"], plan["benchmark"]["substitution_rates"]))
    for accession in {r["accession"] for r in rows}:
        actual = {(r["coverage"], r["error_rate"]) for r in rows if r["accession"] == accession}
        if actual != expected_conditions:
            raise ValueError(f"{accession}: calibration conditions do not match the complete frozen depth x error grid")
    loaded, fingerprint, hashes = load_case_metrics(rows, args.metrics_root, margins)
    mapping = next(iter(loaded.values()))["pipeline"]
    grid, selected, predictions = select_rule(rows, loaded, targets, candidates)
    run = provenance("calibration", args.calibration_manifest, hashes)
    run.update({"freeze_plan_sha256": sha256(args.freeze_plan), "feature_schema_sha256": fingerprint,
                "mapping_contract_sha256": mapping["sha256"], "mapping_contract": mapping["contract"],
                "source_genomes_n": len({r["accession"] for r in rows}), "technical_cases_n": len(rows),
                "source_accessions": sorted({r["accession"] for r in rows}), "lineage_independence_certified_by_this_program": False,
                "selection": "Zero observed calibration cross-target false detections; maximize mean target positive sensitivity across all calibration depth/error cases; ties prefer marker count, breadth, fragments and margin, all descending.",
                "grid_rules_n": len(grid), "population_specificity_guarantee": False})
    out = new_output(args.out)
    write_table(out / "calibration_grid.tsv", grid)
    if selected is None:
        run["status"] = "calibration_failed_no_zero_FP_rule"
        write_json(out / "calibration_status.json", run)
        write_json(out / "run_metadata.json", run)
        return 3
    counts = summaries(predictions, ["target"])
    selected_row = next(r for r in grid if r["selected"])
    frozen = {**run, "status": "pilot_calibrated", "targets": targets, "metrics_margins": margins, "selected_rule": selected,
              "selected_rule_id": selected_row["rule_id"], "balanced_target_positive_sensitivity": selected_row["balanced_target_positive_sensitivity"],
              "observed_calibration_false_detections": selected_row["cross_target_false_detections"],
              "calibration_counts": counts, "calibration_source_accessions": run["source_accessions"],
              "calibration_source_sha256": sorted({r["source_sha256"] for r in rows}),
              "freeze_plan": plan, "claim": "Pilot calibration on the listed reference genomes; not a guarantee of 99% specificity or population-wide sensitivity."}
    write_table(out / "calibration_selected_predictions.tsv", predictions)
    write_table(out / "calibration_counts.tsv", counts)
    write_json(out / "frozen_rules.json", frozen)
    run.update({"status": "pilot_calibrated", "frozen_rules_sha256": sha256(out / "frozen_rules.json")})
    write_json(out / "run_metadata.json", run)
    return 0


def load_rules(path):
    rules = json.loads(Path(path).read_text(encoding="utf-8-sig"))
    if rules.get("status") != "pilot_calibrated" or rules.get("observed_calibration_false_detections") != 0:
        raise ValueError("Rules are not successfully frozen pilot calibration rules")
    if rules.get("targets") != ["infantis", "longum"]:
        raise ValueError("Unexpected rule targets")
    if set(rules.get("selected_rule", {})) != set(RULE_FIELDS):
        raise ValueError("Incomplete frozen detection rule")
    for name, value in rules["selected_rule"].items():
        numeric = nonnegative_number(value, name, integer=name in {"minimum_fragments_per_marker", "supported_identity_markers"})
        if (name != "margin" and numeric <= 0) or (name == "marker_breadth" and numeric > 1):
            raise ValueError("Invalid frozen detection threshold")
    if rules["selected_rule"]["margin"] not in rules["metrics_margins"]:
        raise ValueError("Frozen margin missing from metrics grid")
    return rules


def evaluate(args):
    frozen = load_rules(args.rules)
    rows = load_benchmark_manifest(args.evaluate_manifest, "test")
    if {r["accession"] for r in rows} & set(frozen["calibration_source_accessions"]):
        raise ValueError("Evaluation genome accession also occurs in calibration")
    if {r["source_sha256"] for r in rows} & set(frozen["calibration_source_sha256"]):
        raise ValueError("Evaluation source file SHA also occurs in calibration")
    loaded, _, hashes = load_case_metrics(rows, args.metrics_root, frozen["metrics_margins"], frozen["feature_schema_sha256"], frozen["mapping_contract_sha256"])
    predictions = predict(rows, loaded, frozen["targets"], frozen["selected_rule"])
    out = new_output(args.out)
    write_table(out / "evaluation_predictions.tsv", predictions)
    write_table(out / "evaluation_confusion_by_target.tsv", summaries(predictions, ["target"]))
    write_table(out / "evaluation_confusion_by_target_depth_error.tsv", summaries(predictions, ["target", "coverage", "error_rate"]))
    run = provenance("evaluation_only_no_threshold_selection", args.evaluate_manifest, hashes)
    run.update({"frozen_rules_sha256": sha256(args.rules), "selected_rule": frozen["selected_rule"],
                "technical_cases_n": len(rows), "source_genomes_n": len({r["accession"] for r in rows}),
                "lineage_independence_certified_by_this_program": False, "population_specificity_guarantee": False})
    write_json(out / "run_metadata.json", run)
    return 0


def cohort(args):
    frozen = load_rules(args.rules)
    rows = read_table(args.cohort_manifest, delimiter=",")
    if not rows:
        raise ValueError("Empty cohort manifest")
    sample_column = args.sample_column or next((k for k in ("SampleID", "sample", "case_id") if k in rows[0]), None)
    if sample_column is None or sample_column not in rows[0]:
        raise ValueError("Cohort needs SampleID/sample/case_id, or --sample-column")
    if len({r[sample_column] for r in rows}) != len(rows):
        raise ValueError("Duplicate cohort sample IDs")
    feature_rows = read_table(args.feature_manifest, META_FIELDS)
    if sha256(args.feature_manifest) != frozen["mapping_contract"]["panel_metadata"]["feature_manifest.tsv"]:
        raise ValueError("Cohort feature annotation file differs from the calibrated mapping manifest")
    annotations = {r["feature_id"]: r for r in feature_rows}
    if len(annotations) != len(feature_rows):
        raise ValueError("Duplicate feature-manifest IDs")
    features = {key: metadata(row) for key, row in annotations.items()}
    if schema_hash(features) != frozen["feature_schema_sha256"]:
        raise ValueError("Cohort feature manifest differs from calibrated feature dictionary")
    identity, locus, hashes, missing = [], [], {}, []
    margins, rule = frozen["metrics_margins"], frozen["selected_rule"]
    for row in rows:
        sample = row[sample_column]
        path = metrics_path(args.metrics_root, sample)
        data = load_metrics(path, sample, margins) if path.exists() else None
        if data is not None and data["schema_sha256"] != frozen["feature_schema_sha256"]:
            raise ValueError(f"Cohort metrics feature dictionary changed: {sample}")
        if data is None:
            missing.append(sample)
        else:
            pipeline = pipeline_provenance(path, sample, data["sha256"])
            if pipeline["sha256"] != frozen["mapping_contract_sha256"]:
                raise ValueError(f"Cohort mapping panel/code/version changed: {sample}")
            hashes[sample] = {"gene_metrics_sha256": data["sha256"], "run_manifest_sha256": pipeline["run_manifest_sha256"], "skip_qc": pipeline["skip_qc"]}
        context = {key: row[key] for key in ("MotherID", "DyadID", "Window", "NominalDay", "eligible_frozen_BF76", "identity_note") if key in row}
        for target in frozen["targets"]:
            result = identity_result(data, target, rule) if data else None
            detected = result["detected"] if result else False
            identity.append({"sample": sample, **context, "target": target,
                             "identity_state": "detected" if detected else "insufficient_evidence",
                             "reason": "frozen_identity_rule_passed" if detected else "below_frozen_identity_rule" if result else "metrics_missing",
                             "supported_markers": result["supported_markers"] if result else "",
                             "eligible_markers": result["eligible_markers"] if result else sum(f["channel"] == "identity" and f["taxon_group"] == target and f["eligible_identity"] and not f["pseudo"] for f in features.values()),
                             "metrics_available": int(data is not None), "flag_scope_panel": True, **rule})
        for feature_id, f in features.items():
            if f["channel"] != "locus":
                continue
            for margin in margins:
                metric = data["evidence"][(feature_id, float(margin))] if data else None
                well_covered = metric is not None and metric["breadth"] >= 0.8 and metric["fragments"] >= 3
                locus.append({"sample": sample, **context, **annotations[feature_id], "margin": margin,
                              **{key: metric[key] if metric else "" for key in ("fragments", "covered_bases", "breadth", "aligned_bases", "mean_depth", "ambiguous_fragments")},
                              "sequence_evidence_state": "well_covered_SEQUENCE_EVIDENCE" if well_covered else "below_descriptive_coverage_rule" if metric else "not_assessed_missing_metrics",
                              "at_identity_calibrated_margin": int(margin == rule["margin"]),
                              "metrics_available": int(data is not None), "flag_scope_panel": True,
                              "function_validated": False, "locus_subspecies_linkage_established": False})
    out = new_output(args.out)
    write_table(out / "cohort_identity_calls.tsv", identity)
    write_table(out / "cohort_locus_sequence_evidence.tsv", locus)
    write_table(out / "cohort_identity_summary.tsv", [
        {"target": target, "samples_n": len(rows), "metrics_available_n": len(rows) - len(missing),
         "detected_n": sum(r["target"] == target and r["identity_state"] == "detected" for r in identity),
         "insufficient_evidence_n": sum(r["target"] == target and r["identity_state"] == "insufficient_evidence" for r in identity),
         "missing_metrics_n": len(missing)} for target in frozen["targets"]])
    run = provenance("cohort_frozen_identity_and_descriptive_locus_evidence", args.cohort_manifest, hashes)
    run.update({"frozen_rules_sha256": sha256(args.rules), "feature_manifest_sha256": sha256(args.feature_manifest),
                "samples_n": len(rows), "metrics_missing_samples": missing,
                "feature_annotation_policy": "Retain every locus, including pseudo/transposase annotations; no module-level function inferred.",
                "locus_descriptive_rule": {"breadth_min": 0.8, "fragments_min": 3, "applies_to": "each locus at each reported margin", "validation_claim": "none"},
                "identity_margin_flag": "Reference annotation only: identity threshold calibration is not locus detection-limit or function calibration.",
                "cohort_selection": "All rows in the supplied cohort manifest retained; no new biological eligibility filter applied."})
    write_json(out / "run_metadata.json", run)
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--calibration-manifest", type=Path)
    mode.add_argument("--evaluate-manifest", type=Path)
    mode.add_argument("--cohort-manifest", type=Path)
    parser.add_argument("--metrics-root", required=True, type=Path)
    parser.add_argument("--freeze-plan", type=Path)
    parser.add_argument("--rules", type=Path)
    parser.add_argument("--feature-manifest", type=Path)
    parser.add_argument("--sample-column")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args(argv)
    if args.calibration_manifest:
        if args.freeze_plan is None or args.rules is not None:
            parser.error("Calibration requires --freeze-plan and does not accept --rules")
        return calibrate(args)
    if args.rules is None or args.freeze_plan is not None:
        parser.error("Evaluation/cohort require --rules and do not accept --freeze-plan")
    if args.cohort_manifest:
        if args.feature_manifest is None:
            parser.error("Cohort mode requires --feature-manifest to preserve locus annotations")
        return cohort(args)
    return evaluate(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError) as error:
        print(f"calibrate_calls: {error}", file=sys.stderr)
        sys.exit(2)
