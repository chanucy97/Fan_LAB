#!/usr/bin/env python3
"""Independent, fixed-rule evaluation of existing local-structure fixtures.

Truth files are read only by this evaluator. All input cases remain in the
denominator. Reference-profile calls and exact assembly sequence recovery are
different endpoints and are reported separately.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
import math
from pathlib import Path

SCHEMA = "local-structure-independent-benchmark.v1"
RECOMMENDED_GATE = {
    "min_accepted_fragments": 20,
    "max_bestfit_tv": 0.20,
    "min_group_lower_weight_exclusive": 0.02,
    "min_class_exclusive_boundary_qnames": 3,
    "baseline_breadth_cutoffs": [0.5, 0.8, 0.95],
}


def sha(text):
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def file_record(path):
    path = Path(path)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def architecture_class(inner_length):
    """Prespecified length classes, not species, functional or allele labels."""
    if isinstance(inner_length, bool) or int(inner_length) != float(inner_length):
        raise ValueError("Inner length must be an integer")
    value = int(inner_length)
    if value < 0:
        raise ValueError("Negative inner length")
    if value == 4309:
        return "canonical_4309"
    if 4897 <= value <= 4899:
        return "alternative_4897_4899"
    if value == 824:
        return "short_824"
    return "inner_length_" + str(value)


def bool_value(value):
    if value in (True, "True", "true"):
        return True
    if value in (False, "False", "false", ""):
        return False
    raise ValueError("Invalid historical Boolean: " + str(value))


def finite(value, label, low=0.0, high=1.0):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("Non-numeric " + label)
    if not math.isfinite(value) or not low <= value <= high:
        raise ValueError("Out-of-range " + label)
    return float(value)


def normalize_metadata(metadata):
    if not isinstance(metadata, dict) or not metadata:
        raise ValueError("Nonempty candidate metadata required")
    result = {}
    for candidate, item in metadata.items():
        if not isinstance(candidate, str) or not candidate:
            raise ValueError("Invalid candidate ID")
        left = item.get("inner_start", item.get("inner_start0"))
        right = item.get("inner_end", item.get("inner_end0"))
        length = item["length"]
        if any(type(x) is not int for x in (left, right, length)):
            raise ValueError("Integer interval coordinates required")
        if not 0 < left < right < length:
            raise ValueError("Invalid inner boundaries")
        if sorted(item["junctions"]) != [left, right]:
            raise ValueError("Junctions differ from inner boundaries")
        result[candidate] = dict(item, inner_start=left, inner_end=right,
                                 evaluation_class=architecture_class(right-left))
    return result


def normalize_gate(config):
    """Accept a direct gate or an explicitly named nested evaluation_gate."""
    raw = config.get("evaluation_gate", config)
    gate = dict(raw)
    required = set(RECOMMENDED_GATE)
    if not required <= set(gate):
        raise ValueError("Frozen gate missing keys: " + str(sorted(required-set(gate))))
    for key in ("min_accepted_fragments", "min_class_exclusive_boundary_qnames"):
        if type(gate[key]) is not int or gate[key] < 1:
            raise ValueError("Invalid gate " + key)
    finite(gate["max_bestfit_tv"], "max_bestfit_tv")
    finite(gate["min_group_lower_weight_exclusive"], "min_group_lower_weight_exclusive")
    cutoffs = gate["baseline_breadth_cutoffs"]
    if not isinstance(cutoffs, list) or len(cutoffs) != len(set(cutoffs)) or not cutoffs:
        raise ValueError("Invalid baseline cutoff grid")
    for value in cutoffs:
        finite(value, "breadth cutoff")
    return gate


def class_evidence(extraction, metadata):
    """Count QNAMEs once; require every compatible placement to support a boundary.

    Each extraction candidate_evidence already intersects a candidate's own
    near-best placements. We intersect again across every near-best candidate
    for class-specific, robust support. No origin information is parsed from
    QNAME strings.
    """
    meta = normalize_metadata(metadata)
    audit = extraction["audit"]
    if audit.get("ledger_rows_omitted_by_cap", 0):
        raise ValueError("Truncated fragment ledger cannot establish boundary counts")
    if extraction.get("config", {}).get("ledger_mode") == "none":
        raise ValueError("Fragment ledger required")
    groups = defaultdict(list)
    for candidate, item in meta.items():
        groups[item["evaluation_class"]].append(candidate)
    sets = {group: {name: set() for name in (
        "near_best_qnames", "class_exclusive_qnames", "left_any_qnames", "right_any_qnames",
        "left_robust_qnames", "right_robust_qnames", "left_robust_direct_qnames",
        "right_robust_direct_qnames", "left_robust_bridge_qnames", "right_robust_bridge_qnames",
        "both_robust_qnames")}
        for group in groups}
    seen, accepted = set(), 0
    for row in extraction["fragment_ledger"]:
        name = row["qname"]
        if not isinstance(name, str) or not name or name in seen:
            raise ValueError("Duplicate or invalid ledger QNAME")
        seen.add(name)
        evidence = row["candidate_evidence"]
        if set(row["candidates"]) != set(evidence) or len(row["candidates"]) != len(evidence):
            raise ValueError("Ledger candidate list disagrees with evidence")
        if not evidence:
            continue
        if not set(evidence) <= set(meta):
            raise ValueError("Unknown ledger candidate")
        accepted += 1
        classes = {meta[c]["evaluation_class"] for c in evidence}
        for group in classes:
            sets[group]["near_best_qnames"].add(name)
        exclusive = len(classes) == 1
        group = next(iter(classes)) if exclusive else None
        if exclusive:
            sets[group]["class_exclusive_qnames"].add(name)
        robust_both = []
        for side, field in (("left", "inner_start"), ("right", "inner_end")):
            flags = []
            for candidate, item in evidence.items():
                junction = meta[candidate][field]
                direct, bridge = set(item["direct_junctions"]), set(item["bridged_junctions"])
                if not direct <= set(meta[candidate]["junctions"]) or not bridge <= set(meta[candidate]["junctions"]):
                    raise ValueError("Unexpected junction in ledger")
                if junction in direct | bridge:
                    sets[meta[candidate]["evaluation_class"]][side+"_any_qnames"].add(name)
                flags.append((junction in direct, junction in bridge))
            robust = exclusive and all(a or b for a, b in flags)
            robust_both.append(robust)
            if robust:
                sets[group][side+"_robust_qnames"].add(name)
                if all(a for a, _ in flags):
                    sets[group][side+"_robust_direct_qnames"].add(name)
                if all(b for _, b in flags):
                    sets[group][side+"_robust_bridge_qnames"].add(name)
        if all(robust_both):
            sets[group]["both_robust_qnames"].add(name)
    if accepted != int(audit["accepted_fragments"]):
        raise ValueError("Accepted fragment count and complete ledger disagree")
    coverage = extraction["candidate_coverage"]
    if set(coverage) != set(meta):
        raise ValueError("Candidate coverage dictionary differs from panel")
    result = {}
    for group, members in groups.items():
        eligible, near = [], []
        for candidate in members:
            row = coverage[candidate]
            if row["length"] != meta[candidate]["length"]:
                raise ValueError("Candidate coverage length mismatch")
            e = finite(row["eligible_base_breadth"], "eligible breadth")
            n = finite(row["near_best_base_breadth"], "near-best breadth")
            if n > e+1e-12:
                raise ValueError("Near-best breadth exceeds eligible breadth")
            if row["near_best_fragments"] > row["eligible_fragments"]:
                raise ValueError("Near-best fragments exceed eligible fragments")
            eligible.append(e)
            near.append(n)
        result[group] = dict(members=sorted(members), **{key: len(value) for key, value in sets[group].items()},
                             eligible_base_breadth_max=max(eligible),
                             near_best_base_breadth_max=max(near))
    return result


def _ranges_by_class(fit, evidence):
    if fit["status"] == "no_accepted_evidence":
        return {group: None for group in evidence}
    if fit["status"] != "fitted_descriptive":
        raise ValueError("Unexpected fit status")
    named = fit["config"]["candidate_groups"]
    ranges = fit["group_ranges"]
    output = {}
    for group, row in evidence.items():
        keys = [key for key, members in named.items() if set(members) == set(row["members"])]
        if len(keys) != 1 or keys[0] not in ranges:
            raise ValueError("Fit must contain exactly one range for each full architecture class: " + group)
        pair = ranges[keys[0]]
        if not isinstance(pair, list) or len(pair) != 2:
            raise ValueError("Invalid group range")
        lo, hi = [finite(x, "group range") for x in pair]
        if lo > hi+1e-8:
            raise ValueError("Reversed group range")
        output[group] = [lo, hi]
    return output


def model_decision(fit, evidence, extraction, gate):
    gate = normalize_gate(gate)
    ranges = _ranges_by_class(fit, evidence)
    n = int(extraction["audit"]["accepted_fragments"])
    if float(fit["accepted_fragment_count"]) != n:
        raise ValueError("Fit denominator and extraction denominator disagree")
    tv = fit["bestfit_tv"]
    if n:
        finite(tv, "TV")
    elif tv is not None:
        raise ValueError("No-evidence fit has numeric TV")
    global_failures = []
    if n < gate["min_accepted_fragments"]:
        global_failures.append("insufficient_accepted_fragments")
    if tv is None or tv > gate["max_bestfit_tv"]:
        global_failures.append("poor_or_unavailable_model_fit")
    supported, ambiguous, details = [], [], {}
    for group, row in evidence.items():
        bounds = ranges[group]
        failures = list(global_failures)
        lower = bounds[0] if bounds else None
        upper = bounds[1] if bounds else None
        if lower is None or not lower > gate["min_group_lower_weight_exclusive"]:
            failures.append("group_lower_bound_not_above_threshold")
        for side in ("left", "right"):
            if row[side+"_robust_qnames"] < gate["min_class_exclusive_boundary_qnames"]:
                failures.append("insufficient_"+side+"_boundary_fragments")
        if not failures:
            supported.append(group)
        elif upper is not None and upper > gate["min_group_lower_weight_exclusive"]:
            ambiguous.append(group)
        weights = fit["bestfit_weights"]
        point = sum(finite(weights[c], "point weight") for c in row["members"]) if weights else None
        details[group] = dict(row, weight_lower=lower, weight_upper=upper,
                             point_weight=point, supported=not failures, failure_reasons=failures)
    return {"supported_classes": sorted(supported), "ambiguous_classes": sorted(ambiguous),
            "global_failures": global_failures, "accepted_fragments": n,
            "bestfit_tv": tv, "unmodeled_feature_fraction": fit.get("unmodeled_feature_fraction"),
            "unexplained_local_fraction": fit.get("unexplained_local_fraction"),
            "class_evidence": details}


def reference_decisions(evidence, gate, accepted_fragments):
    gate = normalize_gate(gate)
    if type(accepted_fragments) is not int or accepted_fragments < 0:
        raise ValueError("Nonnegative accepted-fragment count required")
    output = {}
    for field in ("eligible_base_breadth_max", "near_best_base_breadth_max"):
        for cutoff in gate["baseline_breadth_cutoffs"]:
            suffix = field.removesuffix("_base_breadth_max") + "_" + format(cutoff, ".2f")
            name = "expanded_reference_" + suffix
            output[name] = {"supported_classes": sorted(k for k, v in evidence.items() if v[field] >= cutoff),
                            "cutoff": cutoff, "field": field, "ambiguous_classes": [],
                            "mapping_rule": "quality-eligible reciprocal paired placements; MAPQ not filtered; class member maximum, not union breadth"}
            failures, supported, ambiguous = {}, [], []
            for group, row in evidence.items():
                reasons = []
                if accepted_fragments < gate["min_accepted_fragments"]:
                    reasons.append("insufficient_accepted_fragments")
                if row[field] < cutoff:
                    reasons.append("insufficient_reference_breadth")
                for side in ("left", "right"):
                    if row[side+"_robust_qnames"] < gate["min_class_exclusive_boundary_qnames"]:
                        reasons.append("insufficient_"+side+"_boundary_fragments")
                failures[group] = reasons
                if not reasons:
                    supported.append(group)
                elif row[field] >= cutoff:
                    ambiguous.append(group)
            output["matched_support_"+suffix] = {
                "supported_classes": sorted(supported), "ambiguous_classes": sorted(ambiguous),
                "accepted_fragments": accepted_fragments, "cutoff": cutoff, "field": field,
                "class_failure_reasons": failures,
                "mapping_rule": "same depth and robust class-exclusive two-boundary gates as the model plus reference breadth; no TV or weight gate",
                "tv_weight_filters": False,
            }
    return output


def score_sets(predicted, truth):
    predicted, truth = set(predicted), set(truth)
    return {"truth_classes": sorted(truth), "predicted_classes": sorted(predicted),
            "recovered_classes": sorted(predicted & truth), "missed_classes": sorted(truth-predicted),
            "additional_false_classes": sorted(predicted-truth),
            "positive_truth_fully_recovered_no_extra": bool(truth) and predicted == truth,
            "no_class_call": not predicted,
            "negative_case_any_false_class": not truth and bool(predicted),
            "empty_on_negative_no_absence_claim": not truth and not predicted}


def crosscheck_runtime(result, metadata, methods):
    """Compare runtime support summaries with independently recomputed values."""
    meta = normalize_metadata(metadata)
    aliases = {}
    for item in meta.values():
        alias = item.get("class_id", item["evaluation_class"])
        target = item["evaluation_class"]
        if alias in aliases and aliases[alias] != target:
            raise ValueError("Runtime class spans distinct frozen evaluation classes")
        aliases[alias] = target
    checks = 0

    def names(values):
        if not isinstance(values, list) or len(values) != len(set(values)) or not set(values) <= set(aliases):
            raise ValueError("Invalid runtime class list")
        return sorted(aliases[x] for x in values)

    def equal(actual, expected, label):
        nonlocal checks
        if isinstance(expected, float):
            if actual is None or not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError("Runtime cross-check mismatch: " + label)
        elif actual != expected:
            raise ValueError("Runtime cross-check mismatch: " + label)
        checks += 1

    for key, method in (("decision", "paired_junction_model"), ("paired_only_decision", "junction_profile_ablation")):
        raw, expected = result[key], methods[method]
        equal(names(raw["supported_classes"]), expected["supported_classes"], key+" supported classes")
        equal(names(raw["ambiguous_classes"]), expected["ambiguous_classes"], key+" ambiguous classes")
        equal(raw["accepted_fragments"], expected["accepted_fragments"], key+" accepted count")
        equal(raw["bestfit_tv"], expected["bestfit_tv"], key+" TV")
        equal(raw["global_gate_pass"], not expected["global_failures"], key+" global gate")
        if set(raw["classes"]) != set(aliases):
            raise ValueError("Runtime class summary differs from panel")
        for alias, values in raw["classes"].items():
            row = expected["class_evidence"][aliases[alias]]
            equal(sorted(values["members"]), sorted(row["members"]), key+" members")
            for source, target in (("left", "left_robust_qnames"), ("right", "right_robust_qnames"),
                                   ("left_direct", "left_robust_direct_qnames"), ("right_direct", "right_robust_direct_qnames"),
                                   ("lower_weight", "weight_lower"), ("upper_weight", "weight_upper"),
                                   ("supported", "supported"), ("eligible_base_breadth", "eligible_base_breadth_max"),
                                   ("near_best_base_breadth", "near_best_base_breadth_max")):
                equal(values[source], row[target], key+" "+alias+" "+source)
    baselines = result["reference_baseline"]
    for method, values in methods.items():
        if not method.startswith(("expanded_reference_", "matched_support_")):
            continue
        field = values["field"].removesuffix("_max")
        source = baselines["matched_support"] if method.startswith("matched_support_") else baselines
        equal(names(source[field][str(values["cutoff"])]), values["supported_classes"], method)
    return {"status": "PASS", "independently_equal_values": checks,
            "scope": "runtime decisions, boundary counts, group ranges, max breadth and all twelve baseline calls"}


def fasta(text):
    rows, name = {}, None
    for line in text.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in rows:
                raise ValueError("Duplicate FASTA name")
            rows[name] = ""
        elif line.strip():
            if name is None:
                raise ValueError("FASTA sequence before header")
            rows[name] += line.strip().upper()
    if any(set(seq)-set("ACGTRYSWKMBDHVN") for seq in rows.values()):
        raise ValueError("Invalid nucleotide in historical FASTA")
    return rows


def rc(sequence):
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def canonical(sequence):
    return min(sequence, rc(sequence))


def read_table_text(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def _history_text(history_case, key, required=True):
    item = history_case["files"].get(key)
    if not item or not item["exists"]:
        if required:
            raise ValueError("Missing required historical file: " + key)
        return ""
    data = item["text"]
    if hashlib.sha256(data.encode("utf-8")).hexdigest() != item["sha256"]:
        raise ValueError("Historical text hash mismatch")
    return data


def oriented_inner(row, sequence):
    start, end = int(row["interval_start_0based"]), int(row["interval_end_0based_exclusive"])
    left, right = int(row["left_anchor_boundary_0based"]), int(row["right_anchor_boundary_0based"])
    if len(sequence) != end-start or sha(sequence) != row["interval_sequence_sha256"]:
        raise ValueError("Historical interval sequence/coordinate/hash mismatch")
    if not start <= left <= right <= end:
        raise ValueError("Invalid historical boundary geometry")
    if row["orientation_to_reference"] == "+":
        lo, hi = left-start, right-start
    elif row["orientation_to_reference"] == "-":
        lo, hi = end-right, end-left
    else:
        raise ValueError("Invalid historical orientation")
    inner = sequence[lo:hi]
    if row.get("between_anchor_bp") and len(inner) != int(row["between_anchor_bp"]):
        raise ValueError("Historical inner length mismatch")
    return inner


def truth_and_assembly(case, history_case):
    """Use archived sequence truth; never infer truth from new method output."""
    truth_full = fasta(_history_text(history_case, "truth_full"))
    if case["fixture_set"] == "represented_stress24":
        truth_inner = fasta(_history_text(history_case, "truth_inner"))
        if set(truth_full) != set(truth_inner):
            raise ValueError("Historical full/inner truth IDs differ")
        expected = sum(bool(c["positive"]) for c in case["components"])
        if len(truth_full) != expected:
            raise ValueError("Stress positive truth denominator mismatch")
    else:
        rows = read_table_text(_history_text(history_case, "truth_structure_table"))
        if len(rows) != 1 or len(truth_full) != 1:
            raise ValueError("Later-source fixture must retain one source truth interval")
        truth_inner = {r["structure_id"]: oriented_inner(r, truth_full[r["structure_id"]]) for r in rows}
    structures_text = _history_text(history_case, "structure_table", required=False)
    rows = read_table_text(structures_text) if structures_text else []
    assembled = fasta(_history_text(history_case, "anchored_fasta", required=bool(rows)))
    if set(assembled) != {r["structure_id"] for r in rows}:
        raise ValueError("Historical assembly structures/FASTA IDs differ")
    if not structures_text:
        receipt = json.loads(_history_text(history_case, "task_receipt"))
        if receipt.get("candidate_count") != 0 or not receipt.get("geometry_state", "").startswith("NOT_APPLICABLE"):
            raise ValueError("Missing historical geometry not justified by retained receipt")
    full_geometry, inner_geometry, candidate_rows = set(), set(), []
    matched_full, matched_inner = set(), set()
    for row in rows:
        sequence = assembled[row["structure_id"]]
        inner = oriented_inner(row, sequence)
        complete_inner = bool_value(row["inner_anchor_edges_complete"])
        complete_outer = bool_value(row["both_anchors_full_length"])
        group = architecture_class(len(inner))
        if complete_inner:
            inner_geometry.add(group)
        if complete_inner and complete_outer:
            full_geometry.add(group)
        full_matches = [key for key, val in truth_full.items() if canonical(val) == canonical(sequence)]
        inner_matches = [key for key, val in truth_inner.items() if canonical(val) == canonical(inner)]
        if complete_outer and complete_inner:
            matched_full.update(full_matches)
        if complete_inner:
            matched_inner.update(inner_matches)
        candidate_rows.append({"structure_id": row["structure_id"], "class": group,
                               "full_length": len(sequence), "inner_length": len(inner),
                               "complete_outer": complete_outer, "complete_inner": complete_inner,
                               "exact_full_truth_ids": full_matches, "exact_inner_truth_ids": inner_matches})
    true_classes = sorted({architecture_class(len(seq)) for seq in truth_inner.values()})
    return {"truth_classes": true_classes, "truth_alleles": len(truth_full),
            "truth_inner_lengths": sorted(len(seq) for seq in truth_inner.values()),
            "assembly_full_geometry_classes": sorted(full_geometry),
            "assembly_inner_geometry_classes": sorted(inner_geometry),
            "assembly_full_exact_alleles": len(matched_full),
            "assembly_inner_exact_alleles": len(matched_inner),
            "assembly_full_exact_ids": sorted(matched_full), "assembly_inner_exact_ids": sorted(matched_inner),
            "assembly_candidates": candidate_rows}


def evaluate_case(case, result, metadata, gate, history_case):
    if result["case_id"] != case["case_id"] or history_case["case_id"] != case["case_id"]:
        raise ValueError("Case ID mismatch")
    evidence = class_evidence(result["extraction"], metadata)
    main = model_decision(result["fit"], evidence, result["extraction"], gate)
    paired = model_decision(result["paired_only_fit"], evidence, result["extraction"], gate)
    paired["ablation_scope"] = "junction categories removed from fitted profiles only; same boundary evidence retained in decision gate"
    methods = {"paired_junction_model": main, "junction_profile_ablation": paired}
    methods.update(reference_decisions(evidence, gate, int(result["extraction"]["audit"]["accepted_fragments"])))
    truth = truth_and_assembly(case, history_case)
    methods["ordinary_megahit_inner_geometry"] = {"supported_classes": truth["assembly_inner_geometry_classes"], "ambiguous_classes": []}
    methods["ordinary_megahit_full_geometry"] = {"supported_classes": truth["assembly_full_geometry_classes"], "ambiguous_classes": []}
    crosscheck = crosscheck_runtime(result, metadata, methods)
    for item in methods.values():
        item.update(score_sets(item["supported_classes"], truth["truth_classes"]))
    represented = set(truth["truth_classes"]) <= set(evidence)
    state = "negative_whole_interval_deletion" if not truth["truth_classes"] else (
        "represented_length_class" if represented else "unrepresented_length_class")
    # Additional refusal is never added to the recovery count.
    comparisons = {}
    for name, item in methods.items():
        if name == "paired_junction_model":
            continue
        main_hits, baseline_hits = set(main["recovered_classes"]), set(item["recovered_classes"])
        comparisons[name] = {"new_true_classes": sorted(main_hits-baseline_hits),
                             "lost_true_classes": sorted(baseline_hits-main_hits),
                             "new_false_classes": sorted(set(main["additional_false_classes"])-set(item["additional_false_classes"])),
                             "removed_false_classes": sorted(set(item["additional_false_classes"])-set(main["additional_false_classes"])),
                             "abstention_without_new_recovery": main["no_class_call"] and not item["no_class_call"] and not main_hits}
    return {"case_id": case["case_id"], "fixture_set": case["fixture_set"], "truth_stratum": state,
            "source_coverages": [c["coverage"] for c in case["components"]],
            "error_rate": case["error_rate"], "truth_and_assembly": truth,
            "runtime_crosscheck": crosscheck, "methods": methods, "new_method_vs": comparisons}


def write_tsv(path, rows):
    if not rows:
        raise ValueError("No output rows")
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, sort_keys=True) if isinstance(value, (list, dict)) else value for key, value in row.items()})


def summarize(cases):
    counts = defaultdict(Counter)
    for case in cases:
        for method, item in case["methods"].items():
            for stratum in (case["fixture_set"]+"/"+case["truth_stratum"], "ALL_EXPOSED_CASES"):
                row = counts[(stratum, method)]
                row["cases"] += 1
                row["truth_classes"] += len(item["truth_classes"])
                row["recovered_classes"] += len(item["recovered_classes"])
                row["missed_classes"] += len(item["missed_classes"])
                row["additional_false_classes"] += len(item["additional_false_classes"])
                row["positive_truth_fully_recovered_no_extra_cases"] += item["positive_truth_fully_recovered_no_extra"]
                row["no_class_call_cases"] += item["no_class_call"]
                row["negative_false_class_cases"] += item["negative_case_any_false_class"]
                row["ambiguous_class_count"] += len(item["ambiguous_classes"])
    return [dict(stratum=stratum, method=method, **dict(row)) for (stratum, method), row in sorted(counts.items())]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", required=True, type=Path)
    parser.add_argument("--panel-metadata", required=True, type=Path)
    parser.add_argument("--gate-config", required=True, type=Path)
    parser.add_argument("--results-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--history-json", type=Path, required=True,
                        help="Caller-supplied local archive; no remote access")
    args = parser.parse_args(argv)
    if args.output.exists():
        raise FileExistsError("Refuse to replace existing evaluation directory")
    inventory = json.loads(args.inventory.read_text(encoding="utf-8"))
    metadata = json.loads(args.panel_metadata.read_text(encoding="utf-8"))
    normalize_metadata(metadata)
    gate = normalize_gate(json.loads(args.gate_config.read_text(encoding="utf-8")))
    archived = json.loads(args.history_json.read_text(encoding="utf-8"))
    historical = {row["case_id"]: row for row in archived["cases"]}
    if len(historical) != len(archived["cases"]):
        raise ValueError("Duplicate historical case IDs")
    cases = inventory["cases"]
    if len({c["case_id"] for c in cases}) != len(cases) or set(historical) != {c["case_id"] for c in cases}:
        raise ValueError("Duplicate or missing case IDs")
    bindings = [file_record(p) for p in (args.inventory, args.panel_metadata, args.gate_config, Path(__file__))]
    if args.history_json:
        bindings.append(file_record(args.history_json))
    evaluated = []
    for case in cases:
        path = args.results_root / case["case_id"] / "RESULT.json"
        if not path.is_file():
            raise FileNotFoundError("Missing case retained as incomplete; no reduced denominator: " + str(path))
        bindings.append(file_record(path))
        evaluated.append(evaluate_case(case, json.loads(path.read_text(encoding="utf-8")), metadata, gate, historical[case["case_id"]]))
    if len(evaluated) != inventory["case_count"]:
        raise ValueError("Case denominator mismatch")
    positive = sum(c["truth_and_assembly"]["truth_alleles"] for c in evaluated)
    if positive != inventory["positive_alleles"]:
        raise ValueError("Positive allele denominator mismatch")
    summary = summarize(evaluated)
    output = {"schema": SCHEMA, "status": "COMPLETE", "case_count": len(evaluated), "positive_alleles": positive,
              "gate": gate, "scope": "Retrospective exposed-source development, conditional on shared QC/recruitment; descriptive architecture length classes",
              "assembly_exact_full": sum(c["truth_and_assembly"]["assembly_full_exact_alleles"] for c in evaluated),
              "assembly_exact_inner": sum(c["truth_and_assembly"]["assembly_inner_exact_alleles"] for c in evaluated),
              "cases": evaluated, "method_summary": summary}
    args.output.mkdir(parents=True)
    (args.output/"EVALUATION.json").write_text(json.dumps(output, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    (args.output/"HISTORY.json").write_text(json.dumps(archived, indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    rows, class_rows = [], []
    for case in evaluated:
        for method, item in case["methods"].items():
            rows.append({"case_id": case["case_id"], "fixture_set": case["fixture_set"], "truth_stratum": case["truth_stratum"], "method": method,
                         **{k: v for k, v in item.items() if k != "class_evidence"}})
            for group, evidence in item.get("class_evidence", {}).items():
                class_rows.append({"case_id": case["case_id"], "method": method, "class": group,
                                   "truth_present": group in item["truth_classes"], "bestfit_tv": item["bestfit_tv"],
                                   "accepted_fragments": item["accepted_fragments"], **evidence})
    write_tsv(args.output/"METHOD_CASES.tsv", rows)
    write_tsv(args.output/"CLASS_EVIDENCE.tsv", class_rows)
    write_tsv(args.output/"METHOD_SUMMARY.tsv", summary)
    (args.output/"SOURCE_MANIFEST.json").write_text(json.dumps({"inputs": bindings, "outputs": [file_record(p) for p in sorted(args.output.iterdir())]}, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": "COMPLETE", "cases": len(evaluated), "positive_alleles": positive, "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
