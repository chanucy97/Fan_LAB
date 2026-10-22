#!/usr/bin/env python3
"""Independently score v2 using immutable v1 observations and historical truth.

No optimization, alignment or model development is performed here. Observed
counts and no-frequency bounds are independently rebuilt from feature masks.
Model-dependent allocation minima remain fractional fragment quantities.
"""
from __future__ import annotations

import argparse
import copy
import json
import math
from pathlib import Path
from collections import Counter

import evaluate_benchmark as V1

GATE = {"min_accepted_fragments": 20, "max_joint_tv": .2,
        "min_group_lower_weight_exclusive": .02, "min_fractional_boundary_lower": 3,
        "min_informative_qnames": 3, "min_possible_boundary_qnames": 3,
        "baseline_breadth_cutoffs": [.5, .8, .95]}


def normalize_gate(raw):
    gate = raw.get("evaluation_gate", raw)
    if not set(GATE) <= set(gate):
        raise ValueError("Missing frozen v2 gate fields")
    for name in ("min_accepted_fragments", "min_informative_qnames", "min_possible_boundary_qnames"):
        if type(gate[name]) is not int or gate[name] < 1:
            raise ValueError("Invalid integer gate")
    V1.finite(gate["max_joint_tv"], "joint TV gate")
    V1.finite(gate["min_group_lower_weight_exclusive"], "weight gate")
    V1.finite(gate["min_fractional_boundary_lower"], "fractional boundary gate", high=1e9)
    cutoffs = gate["baseline_breadth_cutoffs"]
    if not isinstance(cutoffs, list) or not cutoffs or len(cutoffs) != len(set(cutoffs)):
        raise ValueError("Invalid cutoff grid")
    for value in cutoffs:
        V1.finite(value, "breadth cutoff")
    return gate


def parse_feature(feature, metadata):
    rows = json.loads(feature)
    if not isinstance(rows, list) or not rows:
        raise ValueError("Accepted feature requires at least one candidate")
    result = {}
    for row in rows:
        if not isinstance(row, list) or len(row) not in (1, 3):
            raise ValueError("Malformed feature entry")
        candidate = row[0]
        if candidate not in metadata or candidate in result:
            raise ValueError("Unknown or duplicate feature candidate")
        direct, bridge = row[1:] if len(row) == 3 else ([], [])
        for positions in (direct, bridge):
            if not isinstance(positions, list) or len(positions) != len(set(positions)) or not set(positions) <= set(metadata[candidate]["junctions"]):
                raise ValueError("Invalid observed boundary positions")
        result[candidate] = {"direct": set(direct), "bridge": set(bridge)}
    return result


def _positive(values, integral=False):
    if not isinstance(values, dict):
        raise ValueError("Feature count dictionary required")
    result = {}
    for key, value in values.items():
        if not isinstance(key, str) or not key or isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError("Invalid feature count")
        if integral and int(value) != value:
            raise ValueError("Observed QNAME counts must be integral")
        if value:
            result[key] = float(value)
    return result


def rebuild_masks(counts, profiles, metadata):
    """Closed-form row-wise allocation ranges without profile frequencies."""
    meta = V1.normalize_metadata(metadata)
    if set(profiles) != set(meta):
        raise ValueError("Profile/panel candidate universes differ")
    observed = _positive(counts, integral=True)
    unexplained = observed.pop("__unexplained__", 0)
    n = int(sum(observed.values()))
    normalized, discarded = {}, {}
    for candidate, raw in profiles.items():
        values = _positive(raw)
        discarded[candidate] = values.pop("__unexplained__", 0)
        total = sum(values.values())
        if total <= 0:
            raise ValueError("Empty accepted source profile")
        normalized[candidate] = {key: value/total for key, value in values.items()}
    features = sorted(set(observed).union(*(set(p) for p in normalized.values())))
    parsed = {key: parse_feature(key, meta) for key in features}
    groups = {}
    for candidate, item in meta.items():
        groups.setdefault(item["evaluation_class"], set()).add(candidate)
    profile_sources = {f: {c for c, p in normalized.items() if p.get(f, 0) > 0} for f in features}
    masks = {"observed_mask_no_profile": {f: set(parsed[f]) for f in features},
             "same_mask_no_profile": {f: set(parsed[f]) | profile_sources[f] for f in features}}
    retention = {c: sum(p for f, p in values.items() if c not in parsed[f]) for c, values in normalized.items()}
    support, no_frequency = {}, {name: {"group_ranges": {}, "boundary_bounds": {}} for name in masks}
    for group, members in groups.items():
        row = {"informative_qnames": sum(int(observed.get(f, 0)) for f in features if set(parsed[f]) <= members)}
        boundary_sets = {}
        for side in ("left", "right"):
            field = "inner_start" if side == "left" else "inner_end"
            direct, bridge, any_support = {}, {}, {}
            for feature in features:
                evidence = parsed[feature]
                direct[feature] = {c for c in members & set(evidence) if meta[c][field] in evidence[c]["direct"]}
                bridge[feature] = {c for c in members & set(evidence) if meta[c][field] in evidence[c]["bridge"]}
                any_support[feature] = direct[feature] | bridge[feature]
            boundary_sets[side] = {"": any_support, "direct_": direct}
            row[side] = {}
            for prefix, sets in (("", any_support), ("direct_", direct), ("bridged_", bridge)):
                row[side][prefix+"possible_qnames"] = sum(int(observed.get(f, 0)) for f in features if sets[f])
                row[side][prefix+"exclusive_qnames"] = sum(int(observed.get(f, 0)) for f in features if set(parsed[f]) <= sets[f])
        support[group] = row
        for name, mask in masks.items():
            if not n:
                continue
            lo = sum(observed.get(f, 0) for f in features if mask[f] <= members)/n
            hi = sum(observed.get(f, 0) for f in features if mask[f] & members)/n
            no_frequency[name]["group_ranges"][group] = [lo, hi]
            no_frequency[name]["boundary_bounds"][group] = {}
            for side in ("left", "right"):
                values = {}
                for prefix, supported in boundary_sets[side].items():
                    low = sum(observed.get(f, 0) for f in features if mask[f] <= supported[f])
                    high = sum(observed.get(f, 0) for f in features if mask[f] & supported[f])
                    values.update({prefix+"lower_mass": low/n, prefix+"upper_mass": high/n,
                                   prefix+"lower_fractional_fragments": float(low), prefix+"upper_fractional_fragments": float(high)})
                no_frequency[name]["boundary_bounds"][group][side] = values
    return {"accepted_fragments": n, "unexplained_count": int(unexplained), "observed_support": support,
            "no_frequency": no_frequency if n else {name: None for name in masks},
            "profile_only_probability_by_candidate": retention, "profile_unexplained_counts": discarded,
            "feature_count": len(features), "groups": {g: sorted(m) for g, m in groups.items()}}


def alias_map(fit, independent):
    configured = fit["config"]["candidate_groups"]
    aliases = {}
    for group, members in independent["groups"].items():
        names = [name for name, values in configured.items() if set(values) == set(members)]
        if len(names) != 1:
            raise ValueError("Fit must contain each full evaluation class once")
        aliases[names[0]] = group
    if set(aliases) != set(configured):
        raise ValueError("Unexpected overlapping or extra runtime group")
    return aliases


def assert_close(actual, expected, label, tolerance=1e-7):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            raise ValueError("Dictionary mismatch: " + label)
        return sum(assert_close(actual[k], v, label+"/"+str(k), tolerance) for k, v in expected.items())
    if isinstance(expected, (list, tuple)):
        if not isinstance(actual, (list, tuple)) or len(actual) != len(expected):
            raise ValueError("List mismatch: " + label)
        return sum(assert_close(a, b, label, tolerance) for a, b in zip(actual, expected))
    if expected is None or isinstance(expected, (str, bool)):
        if actual != expected:
            raise ValueError("Value mismatch: " + label)
    elif isinstance(actual, bool) or not isinstance(actual, (int, float)) or not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance):
        raise ValueError("Numeric mismatch: " + label)
    return 1


def verify_fit(fit, independent):
    n = independent["accepted_fragments"]
    aliases = alias_map(fit, independent)
    checked = assert_close(fit["accepted_fragment_count"], n, "accepted")
    checked += assert_close(fit["unexplained_local_fragment_count"], independent["unexplained_count"], "unexplained")
    checked += assert_close(fit["feature_count"], independent["feature_count"], "feature count")
    checked += assert_close(fit["profile_only_probability_by_candidate"], independent["profile_only_probability_by_candidate"], "profile retention")
    checked += assert_close(fit["profile_unexplained_counts"], independent["profile_unexplained_counts"], "profile unexplained")
    normalized = {"group_ranges": {}, "boundary_bounds": {}, "joint_tv": fit["joint_bestfit_tv"]}
    for alias, group in aliases.items():
        checked += assert_close(fit["observed_support"][alias], independent["observed_support"][group], "observed support/"+group)
    for mode, expected in independent["no_frequency"].items():
        actual = fit[mode]
        if not n:
            checked += assert_close(actual, None, mode)
            continue
        for alias, group in aliases.items():
            checked += assert_close(actual["group_ranges"][alias], expected["group_ranges"][group], mode+"/group/"+group)
            checked += assert_close(actual["boundary_bounds"][alias], expected["boundary_bounds"][group], mode+"/boundaries/"+group)
    if not n:
        if fit["status"] != "no_accepted_evidence" or any(fit[k] is not None for k in ("group_ranges", "boundary_bounds", "joint_bestfit_tv", "bestfit_weights")):
            raise ValueError("No-evidence fit contains numeric inference")
        return normalized, aliases, checked
    if fit["status"] != "fitted_descriptive":
        raise ValueError("Unexpected v2 fit state")
    D = V1.finite(fit["joint_bestfit_tv"], "joint TV")
    marginal = V1.finite(fit["marginal_tv_at_joint_bestfit"], "marginal TV at joint optimum")
    budget = V1.finite(fit["tv_budget"], "TV budget")
    if marginal > D+1e-7 or D > budget+1e-7:
        raise ValueError("Joint/marginal/budget inconsistency")
    validation = fit["solver_validation"]
    residual = max(validation[k] for k in ("max_equality_residual", "max_inequality_residual", "max_nonnegative_violation"))
    if residual > validation["tolerance"]:
        raise ValueError("Recorded solve residual exceeds declared tolerance")
    checked += assert_close(sum(fit["bestfit_weights"].values()), 1, "point weights sum", tolerance=max(1e-7, validation["tolerance"]))
    tol = max(1e-6, n*validation["tolerance"]*2)
    for alias, group in aliases.items():
        bounds = fit["group_ranges"][alias]
        lo, hi = [V1.finite(v, "group bound") for v in bounds]
        point = sum(fit["bestfit_weights"][c] for c in independent["groups"][group])
        nflo, nfhi = independent["no_frequency"]["same_mask_no_profile"]["group_ranges"][group]
        if not nflo-1e-7 <= lo <= point+1e-7 or not point-1e-7 <= hi <= nfhi+1e-7:
            raise ValueError("Model group bounds are inconsistent with point/no-frequency bounds")
        normalized["group_ranges"][group] = bounds
        normalized["boundary_bounds"][group] = fit["boundary_bounds"][alias]
        for side in ("left", "right"):
            row = fit["boundary_bounds"][alias][side]
            for prefix in ("", "direct_"):
                lower, point_b = row[prefix+"lower_fractional_fragments"], row[prefix+"bestfit_fractional_fragments"]
                possible = independent["observed_support"][group][side][prefix+"possible_qnames"]
                nf = independent["no_frequency"]["same_mask_no_profile"]["boundary_bounds"][group][side][prefix+"lower_fractional_fragments"]
                checked += assert_close(lower, row[prefix+"lower_mass"]*n, "boundary units", tolerance=1e-7)
                checked += assert_close(point_b, row[prefix+"bestfit_mass"]*n, "point boundary units", tolerance=1e-7)
                if not nf-tol <= lower <= point_b+tol or not point_b <= possible+tol or lower > hi*n+tol:
                    raise ValueError("Model boundary bound inconsistent with observed possible support")
            if row["direct_lower_fractional_fragments"] > row["lower_fractional_fragments"]+tol:
                raise ValueError("Direct-only lower bound exceeds all-boundary lower bound")
    return normalized, aliases, checked


def decision(bounds, independent, gate, require_tv=True):
    n = independent["accepted_fragments"]
    called, ambiguous, details = [], [], {}
    for group in independent["groups"]:
        pair = bounds["group_ranges"].get(group) if bounds else None
        boundary = bounds["boundary_bounds"].get(group) if bounds else None
        lo, hi = pair if pair else (None, None)
        failures = []
        if n < gate["min_accepted_fragments"]:
            failures.append("insufficient_accepted_fragments")
        if require_tv and (bounds["joint_tv"] is None or bounds["joint_tv"] > gate["max_joint_tv"]):
            failures.append("poor_or_unavailable_joint_fit")
        if lo is None or not lo > gate["min_group_lower_weight_exclusive"]:
            failures.append("group_lower_bound_not_above_threshold")
        values = {}
        for side in ("left", "right"):
            value = boundary[side]["lower_fractional_fragments"] if boundary else None
            values[side+"_lower_fractional_fragments"] = value
            if value is None or value < gate["min_fractional_boundary_lower"]:
                failures.append("insufficient_"+side+"_allocation_lower_bound")
        if not failures:
            called.append(group)
        elif hi is not None and hi > gate["min_group_lower_weight_exclusive"]:
            ambiguous.append(group)
        details[group] = {"lower_weight": lo, "upper_weight": hi, **values,
                          "supported": not failures, "failure_reasons": failures,
                          "observed_support": independent["observed_support"][group]}
    return {"supported_classes": sorted(called), "ambiguous_classes": sorted(ambiguous), "classes": details,
            "accepted_fragments": n, "joint_tv": bounds.get("joint_tv") if bounds else None,
            "fractional_unit": "model-dependent lower allocation; not observed exclusive QNAMEs"}


def possible_baselines(evidence, independent, gate):
    methods = {}
    for field in ("eligible_base_breadth_max", "near_best_base_breadth_max"):
        for cutoff in gate["baseline_breadth_cutoffs"]:
            support, details = [], {}
            for group, row in evidence.items():
                counts = independent["observed_support"][group]
                reasons = []
                if independent["accepted_fragments"] < gate["min_accepted_fragments"]:
                    reasons.append("insufficient_accepted_fragments")
                if counts["informative_qnames"] < gate["min_informative_qnames"]:
                    reasons.append("insufficient_informative_fragments")
                for side in ("left", "right"):
                    if counts[side]["possible_qnames"] < gate["min_possible_boundary_qnames"]:
                        reasons.append("insufficient_"+side+"_possible_support")
                if row[field] < cutoff:
                    reasons.append("insufficient_breadth")
                if not reasons:
                    support.append(group)
                details[group] = reasons
            key = "possible_boundary_"+field.removesuffix("_base_breadth_max")+"_"+format(cutoff, ".2f")
            methods[key] = {"supported_classes": sorted(support), "ambiguous_classes": [],
                            "field": field, "cutoff": cutoff, "class_failure_reasons": details}
    return methods


def crosscheck_calls(result, aliases, methods):
    def translated(names):
        if not isinstance(names, list) or len(names) != len(set(names)) or not set(names) <= set(aliases):
            raise ValueError("Invalid runtime v2 class list")
        return sorted(aliases[n] for n in names)
    raw = result["decision_v2"]
    expected = methods["coupled_assignment_model"]
    checks = assert_close(translated(raw["supported_classes"]), expected["supported_classes"], "v2 calls")
    if set(raw["classes"]) != set(aliases):
        raise ValueError("Runtime v2 classes missing")
    for alias, item in raw["classes"].items():
        target = expected["classes"][aliases[alias]]
        for key in ("lower_weight", "upper_weight", "left_lower_fractional_fragments", "right_lower_fractional_fragments", "supported"):
            checks += assert_close(item[key], target[key], "v2 decision class "+key)
    for name, values in methods.items():
        if name.startswith("possible_boundary_"):
            field = values["field"].removesuffix("_max")
            observed = result["reference_baseline_v2"]["possible_support"][field][str(values["cutoff"])]
            checks += assert_close(translated(observed), values["supported_classes"], name)
    for mode in ("same_mask_no_profile", "observed_mask_no_profile"):
        checks += assert_close(translated(result["no_profile_decisions"][mode]), methods[mode]["supported_classes"], mode)
    return checks


def evaluate_case(case, v1result, v1evaluation, v2result, profiles, metadata, gate):
    cid = case["case_id"]
    if any(row["case_id"] != cid for row in (v1result, v1evaluation, v2result)):
        raise ValueError("Mismatched case identity")
    independent = rebuild_masks(v1result["extraction"]["counts"], profiles, metadata)
    evidence = V1.class_evidence(v1result["extraction"], metadata)
    if independent["accepted_fragments"] != v1result["extraction"]["audit"]["accepted_fragments"]:
        raise ValueError("Counts and ledger denominators differ")
    for group, row in evidence.items():
        observed = independent["observed_support"][group]
        if observed["informative_qnames"] != row["class_exclusive_qnames"]:
            raise ValueError("Feature/ledger informative count mismatch")
        for side in ("left", "right"):
            if observed[side]["possible_qnames"] != row[side+"_any_qnames"] or observed[side]["exclusive_qnames"] != row[side+"_robust_qnames"]:
                raise ValueError("Feature/ledger boundary count mismatch")
    normalized, aliases, checked = verify_fit(v2result["assignment_fit"], independent)
    methods = {"coupled_assignment_model": decision(normalized, independent, gate)}
    methods.update(possible_baselines(evidence, independent, gate))
    for mode, values in independent["no_frequency"].items():
        methods[mode] = decision(values, independent, gate, require_tv=False)
    checked += crosscheck_calls(v2result, aliases, methods)
    truth = v1evaluation["truth_and_assembly"]
    for item in methods.values():
        item.update(V1.score_sets(item["supported_classes"], truth["truth_classes"]))
    for name, item in v1evaluation["methods"].items():
        methods["v1_"+name] = copy.deepcopy(item)
    main = methods["coupled_assignment_model"]
    comparisons = {}
    for name, item in methods.items():
        if name == "coupled_assignment_model":
            continue
        hits, other = set(main["recovered_classes"]), set(item["recovered_classes"])
        comparisons[name] = {"new_true_classes": sorted(hits-other), "lost_true_classes": sorted(other-hits),
                             "new_false_classes": sorted(set(main["additional_false_classes"])-set(item["additional_false_classes"])),
                             "removed_false_classes": sorted(set(item["additional_false_classes"])-set(main["additional_false_classes"])),
                             "abstention_without_new_recovery": main["no_class_call"] and not item["no_class_call"] and not hits}
    fit = v2result["assignment_fit"]
    return {"case_id": cid, "fixture_set": case["fixture_set"], "truth_stratum": v1evaluation["truth_stratum"],
            "source_coverages": v1evaluation["source_coverages"], "error_rate": case["error_rate"],
            "truth_and_assembly": copy.deepcopy(truth), "methods": methods, "new_method_vs": comparisons,
            "independent_checks": checked, "status": "PASS", "joint_bestfit_tv": fit["joint_bestfit_tv"],
            "marginal_tv_at_joint_bestfit": fit["marginal_tv_at_joint_bestfit"], "tv_budget": fit["tv_budget"],
            "bestfit_profile_only_allocation_mass": fit["bestfit_profile_only_allocation_mass"],
            "observed_support": independent["observed_support"], "no_frequency": independent["no_frequency"]}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("inventory", "panel-metadata", "profiles", "gate-config", "v1-evaluation", "v1-results-root", "v2-results-root", "output"):
        p.add_argument("--"+name, type=Path, required=True)
    args = p.parse_args(argv)
    if args.output.exists():
        raise FileExistsError("Refuse to overwrite evaluation output")
    def load(path):
        return json.loads(path.read_text(encoding="utf-8"))
    inventory, metadata, profiles = load(args.inventory), load(args.panel_metadata), load(args.profiles)
    gate = normalize_gate(load(args.gate_config))
    first = load(args.v1_evaluation)
    if first["status"] != "COMPLETE" or first["case_count"] != inventory["case_count"]:
        raise ValueError("Incomplete immutable v1 evaluation")
    v1cases = {row["case_id"]: row for row in first["cases"]}
    sources = [V1.file_record(path) for path in (args.inventory, args.panel_metadata, args.profiles, args.gate_config, args.v1_evaluation, Path(__file__), Path(V1.__file__))]
    cases = []
    for case in inventory["cases"]:
        cid = case["case_id"]
        old = args.v1_results_root/cid/"RESULT.json"
        new = args.v2_results_root/cid/"RESULT_V2.json"
        if not old.is_file() or not new.is_file():
            raise FileNotFoundError("Missing case; denominator must not shrink: "+cid)
        sources.extend([V1.file_record(old), V1.file_record(new)])
        cases.append(evaluate_case(case, load(old), v1cases[cid], load(new), profiles, metadata, gate))
    if len(cases) != inventory["case_count"] or sum(c["truth_and_assembly"]["truth_alleles"] for c in cases) != inventory["positive_alleles"]:
        raise ValueError("Fixed denominators differ")
    summary = V1.summarize(cases)
    result = {"schema": "independent-joint-assignment-evaluation.v2", "status": "COMPLETE", "case_count": len(cases),
              "positive_interval_alleles": inventory["positive_alleles"], "gate": gate, "cases": cases, "method_summary": summary,
              "independent_equal_values": sum(c["independent_checks"] for c in cases),
              "scope": "Retrospective v2 development after v1 outcome inspection; same immutable observations and candidate panel; no independent external validation"}
    args.output.mkdir(parents=True)
    (args.output/"EVALUATION_V2.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    rows, class_rows = [], []
    for case in cases:
        for method, values in case["methods"].items():
            rows.append({"case_id": case["case_id"], "fixture_set": case["fixture_set"], "truth_stratum": case["truth_stratum"], "method": method,
                         **{k: v for k, v in values.items() if k not in ("classes", "class_evidence")}})
            for group, value in values.get("classes", {}).items():
                class_rows.append({"case_id": case["case_id"], "method": method, "class": group,
                                   "truth_present": group in values["truth_classes"], **value})
    V1.write_tsv(args.output/"METHOD_CASES_V2.tsv", rows)
    V1.write_tsv(args.output/"CLASS_EVIDENCE_V2.tsv", class_rows)
    V1.write_tsv(args.output/"METHOD_SUMMARY_V2.tsv", summary)
    (args.output/"SOURCE_MANIFEST.json").write_text(json.dumps({"inputs": sources, "outputs": [V1.file_record(q) for q in sorted(args.output.iterdir())]}, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": "COMPLETE", "cases": len(cases), "independent_equal_values": result["independent_equal_values"], "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
