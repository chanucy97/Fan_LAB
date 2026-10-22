"""Experimental fragment-equivalence inference for fixed local candidates.

One QNAME is one observed sequencing fragment, not one cell or UMI molecule.
Fits and TV-budget ranges are descriptive, conditional on accepted local pairs.
They are neither confidence intervals nor biological presence/absence calls.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import groupby
import json
import math
from pathlib import Path
import re

SCHEMA = "glycan-evidence.structure-inference.research.v1"
UNEXPLAINED = "__unexplained__"
DEFAULTS = dict(near_best_delta=6, min_aligned_fraction=.9, min_map_identity=.9,
                min_anchor=15, max_insert=1000, mode="paired_junction",
                ledger_mode="accepted", max_ledger=None)


def _config(config):
    cfg = {**DEFAULTS, **(config or {})}
    if cfg["mode"] not in {"paired_junction", "paired"}:
        raise ValueError("mode must be paired_junction or paired")
    if cfg["ledger_mode"] not in {"accepted", "all", "none"}:
        raise ValueError("Invalid ledger_mode")
    for key in ("near_best_delta", "min_aligned_fraction", "min_map_identity"):
        value = cfg[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError("Invalid " + key)
    for key in ("min_aligned_fraction", "min_map_identity"):
        if cfg[key] > 1:
            raise ValueError("Invalid " + key)
    for key in ("min_anchor", "max_insert"):
        if type(cfg[key]) is not int or cfg[key] < 1:
            raise ValueError("Invalid " + key)
    if cfg["max_ledger"] is not None and (type(cfg["max_ledger"]) is not int or cfg["max_ledger"] < 0):
        raise ValueError("Invalid max_ledger")
    return cfg


def _metadata(metadata):
    if not isinstance(metadata, dict) or not metadata:
        raise ValueError("Nonempty candidate metadata required")
    result = {}
    for name, item in metadata.items():
        if not isinstance(name, str) or not name or any(c.isspace() for c in name) or name == UNEXPLAINED:
            raise ValueError("Invalid candidate ID")
        if not isinstance(item, dict) or type(item.get("length")) is not int or item["length"] < 1:
            raise ValueError("Invalid candidate length")
        junctions = item.get("junctions", [])
        if not isinstance(junctions, (list, tuple)) or any(type(p) is not int or not 0 < p < item["length"] for p in junctions):
            raise ValueError("Junctions must be zero-based interbase coordinates inside sequence")
        if len(set(junctions)) != len(junctions):
            raise ValueError("Duplicate junction")
        result[name] = dict(item, junctions=sorted(junctions))
    return result


def _cigar(text, start):
    parts = re.findall(r"(\d+)([MIDNSHP=X])", text)
    if not parts or "".join(n + op for n, op in parts) != text or any(int(n) < 1 for n, _ in parts):
        raise ValueError("Invalid CIGAR")
    rpos, qlength, aligned, denominator = start, 0, 0, 0
    blocks, block_start = [], None
    for raw, op in parts:
        length = int(raw)
        if op in "M=X":
            if block_start is None:
                block_start = rpos
            rpos += length
            aligned += length
            qlength += length
            denominator += length
        else:
            if block_start is not None:
                blocks.append((block_start, rpos))
                block_start = None
            if op in "IS":
                qlength += length
            if op == "I":
                aligned += length
                denominator += length
            if op in "DN":
                rpos += length
            if op == "D":
                denominator += length
    if block_start is not None:
        blocks.append((block_start, rpos))
    return blocks, qlength, aligned, denominator, rpos


def _record(line, metadata):
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 11:
        raise ValueError("Truncated SAM record")
    name, flag = fields[0], int(fields[1])
    if not name or name == "*" or any(c.isspace() for c in name) or not 0 <= flag <= 0xFFF:
        raise ValueError("Invalid QNAME/FLAG")
    if bool(flag & 64) == bool(flag & 128):
        raise ValueError("Every record must declare exactly one mate")
    if not flag & 1:
        raise ValueError("Paired records required")
    mapq = int(fields[4])
    if not 0 <= mapq <= 255:
        raise ValueError("Invalid MAPQ")
    result = dict(qname=name, flag=flag, mate=1 if flag & 64 else 2, mapq=mapq)
    if flag & 4:
        return result
    ref, start = fields[2], int(fields[3]) - 1
    if ref not in metadata or start < 0:
        raise ValueError("Unknown mapped reference/coordinate")
    blocks, qlength, aligned, denominator, end = _cigar(fields[5], start)
    if end > metadata[ref]["length"] or qlength <= 0 or denominator <= 0:
        raise ValueError("Invalid mapped CIGAR extent")
    if fields[9] != "*" and len(fields[9]) != qlength:
        raise ValueError("SEQ and CIGAR disagree")
    tags = {}
    for value in fields[11:]:
        bits = value.split(":", 2)
        if len(bits) != 3 or bits[0] in tags:
            raise ValueError("Malformed or duplicate SAM tag")
        tags[bits[0]] = bits[1:]
    for key in ("AS", "NM"):
        if key not in tags or tags[key][0] != "i":
            raise ValueError("Mapped records require integer AS and NM")
    score, nm = int(tags["AS"][1]), int(tags["NM"][1])
    if nm < 0 or nm > denominator:
        raise ValueError("Invalid NM")
    ys = int(tags["YS"][1]) if "YS" in tags and tags["YS"][0] == "i" else None
    result.update(ref=ref, start=start, end=end, blocks=blocks,
                  mate_ref=ref if fields[6] == "=" else fields[6], mate_start=int(fields[7])-1,
                  tlen=int(fields[8]), score=score, nm=nm, ys=ys,
                  query_length=qlength,
                  aligned_fraction=aligned / qlength, identity=1 - nm / denominator)
    return result


def _pair_key(record):
    r = record
    if r["mate"] == 1:
        return (r["ref"], r["start"], r["mate_start"], bool(r["flag"] & 16), bool(r["flag"] & 32))
    return (r["ref"], r["mate_start"], r["start"], bool(r["flag"] & 32), bool(r["flag"] & 16))


def _junctions(one, two, junctions, flank):
    direct, bridged = [], []
    for position in junctions:
        if any(lo <= position-flank and hi >= position+flank for r in (one, two) for lo, hi in r["blocks"]):
            direct.append(position)
        # Pair bridges are physical fragment evidence, not sequenced junction bases.
        left_one = any(hi <= position and hi-lo >= flank for lo, hi in one["blocks"])
        right_one = any(lo >= position and hi-lo >= flank for lo, hi in one["blocks"])
        left_two = any(hi <= position and hi-lo >= flank for lo, hi in two["blocks"])
        right_two = any(lo >= position and hi-lo >= flank for lo, hi in two["blocks"])
        if (left_one and right_two) or (left_two and right_one):
            bridged.append(position)
    return direct, bridged


def _eligible_pairs(records, metadata, cfg, audit):
    primary = Counter(r["mate"] for r in records if not r["flag"] & (256 | 2048))
    if any(n > 1 for n in primary.values()):
        raise ValueError("Duplicate primary mate records within QNAME")
    placements = defaultdict(lambda: {1: {}, 2: {}})
    for r in records:
        flag = r["flag"]
        if flag & (4 | 8 | 2048 | 512) or not flag & 2:
            audit["nonconcordant_or_filtered_records"] += 1
            continue
        if r["mate_ref"] != r["ref"] or r["mate_start"] < 0 or bool(flag & 16) == bool(flag & 32):
            raise ValueError("Proper pair has invalid mate geometry")
        key = _pair_key(r)
        variant = (r["score"], r["nm"], tuple(r["blocks"]), r["tlen"], r["ys"],
                   r["end"], r["query_length"], r["aligned_fraction"], r["identity"])
        if variant in placements[key][r["mate"]]:
            audit["duplicate_placement_records_collapsed"] += 1
        else:
            placements[key][r["mate"]][variant] = r
    accepted = []
    for pair in placements.values():
        if not pair[1] or not pair[2]:
            raise ValueError("Concordant placement lacks reciprocal mate record")
        matched_one, matched_two, invalid_reasons = set(), set(), set()
        # Same coordinates can support distinct gapped CIGARs. Keep all mate
        # combinations consistent with reciprocal coordinates, extent and tags.
        # This is conservative for junction intersection, unlike taking a union.
        for key_one, one in pair[1].items():
            for key_two, two in pair[2].items():
                if one["tlen"] != -two["tlen"] or abs(one["tlen"]) != max(one["end"], two["end"])-min(one["start"], two["start"]):
                    invalid_reasons.add("TLEN")
                    continue
                if any(r["ys"] is not None and r["ys"] != mate["score"] for r, mate in ((one, two), (two, one))):
                    invalid_reasons.add("AS/YS")
                    continue
                matched_one.add(key_one)
                matched_two.add(key_two)
                if (abs(one["tlen"]) > cfg["max_insert"] or
                        any(r["aligned_fraction"] < cfg["min_aligned_fraction"] or r["identity"] < cfg["min_map_identity"] for r in (one, two))):
                    audit["quality_rejected_pair_placements"] += 1
                    continue
                direct, bridged = _junctions(one, two, metadata[one["ref"]]["junctions"], cfg["min_anchor"])
                accepted.append(dict(candidate=one["ref"], score=one["score"]+two["score"],
                                     direct_junctions=direct, bridged_junctions=bridged,
                                     records=(one, two)))
        if len(matched_one) != len(pair[1]) or len(matched_two) != len(pair[2]):
            raise ValueError("Concordant alignment variant has no compatible reciprocal mate: " + "/".join(sorted(invalid_reasons)))
        if len(pair[1]) > 1 or len(pair[2]) > 1:
            audit["same_coordinate_variant_groups"] += 1
    return accepted


def extract_fragment_features(sam_path, candidate_metadata, config=None):
    """Read a complete QNAME-collated paired SAM, retaining near-best placements.

    SAM may omit unaligned pairs. Its residual denominator is then only QNAMEs
    present in this local alignment file, never total metagenomic read count.
    MAPQ 0/255 and secondary alignments are retained for equivalence analysis.
    """
    cfg, metadata = _config(config), _metadata(candidate_metadata)
    path = Path(sam_path)
    counts, audit, ledger, seen = Counter(), Counter(), [], set()
    coverage = {name: bytearray(item["length"]) for name, item in metadata.items()}
    near_coverage = {name: bytearray(item["length"]) for name, item in metadata.items()}
    candidate_fragments, near_fragments = Counter(), Counter()
    read_lengths, insert_lengths = Counter(), Counter()
    header = {}
    with path.open(encoding="utf-8") as handle:
        first = None
        for line in handle:
            if line.startswith("@SQ\t"):
                tags = dict(field.split(":", 1) for field in line.rstrip().split("\t")[1:])
                if tags.get("SN") in header or "SN" not in tags or "LN" not in tags:
                    raise ValueError("Invalid/duplicate SQ header")
                header[tags["SN"]] = int(tags["LN"])
            elif line.startswith("@"):
                continue
            elif line.strip():
                first = line
                break
        if header != {name: item["length"] for name, item in metadata.items()}:
            raise ValueError("SAM SQ IDs/lengths differ from candidate metadata")
        def rows():
            if first is not None:
                yield _record(first, metadata)
            for line in handle:
                if line.startswith("@"):
                    raise ValueError("SAM header after records")
                if line.strip():
                    yield _record(line, metadata)
        for name, group in groupby(rows(), key=lambda r: r["qname"]):
            if name in seen:
                raise ValueError("QNAME repeated in noncontiguous groups")
            seen.add(name)
            audit["qname_groups"] += 1
            records = list(group)
            audit["alignment_records"] += len(records)
            audit["secondary_records"] += sum(bool(r["flag"] & 256) for r in records)
            audit["mapq_zero_records_retained"] += sum(r["mapq"] == 0 and not r["flag"] & 4 for r in records)
            pairs = _eligible_pairs(records, metadata, cfg, audit)
            for candidate in {p["candidate"] for p in pairs}:
                candidate_fragments[candidate] += 1
            for pair in pairs:
                for r in pair["records"]:
                    for lo, hi in r["blocks"]:
                        coverage[pair["candidate"]][lo:hi] = b"\1" * (hi-lo)
            evidence, best = {}, None
            if not pairs:
                feature = UNEXPLAINED
                audit["unexplained_fragments"] += 1
            else:
                best = max(p["score"] for p in pairs)
                near = [p for p in pairs if p["score"] >= best-cfg["near_best_delta"]]
                representative = min((p for p in near if p["score"] == best), key=lambda p:
                    (p["candidate"], tuple((r["start"], r["end"], tuple(r["blocks"]), r["query_length"])
                                          for r in p["records"])))
                insert_lengths[abs(representative["records"][0]["tlen"])] += 1
                read_lengths.update(r["query_length"] for r in representative["records"])
                for candidate in sorted({p["candidate"] for p in near}):
                    subset = [p for p in near if p["candidate"] == candidate]
                    near_fragments[candidate] += 1
                    for pair in subset:
                        for r in pair["records"]:
                            for lo, hi in r["blocks"]:
                                near_coverage[candidate][lo:hi] = b"\1" * (hi-lo)
                    # Multimapping locations cannot jointly manufacture junction support.
                    direct = sorted(set.intersection(*(set(p["direct_junctions"]) for p in subset)))
                    bridge = sorted(set.intersection(*(set(p["bridged_junctions"]) for p in subset)))
                    evidence[candidate] = dict(score=max(p["score"] for p in subset),
                                               direct_junctions=direct, bridged_junctions=bridge,
                                               placements=len(subset))
                entries = [[candidate, evidence[candidate]["direct_junctions"], evidence[candidate]["bridged_junctions"]]
                           if cfg["mode"] == "paired_junction" else [candidate] for candidate in evidence]
                feature = json.dumps(entries, separators=(",", ":"), ensure_ascii=True)
                audit["accepted_fragments"] += 1
            counts[feature] += 1
            keep = cfg["ledger_mode"] == "all" or (cfg["ledger_mode"] == "accepted" and bool(evidence))
            if keep and (cfg["max_ledger"] is None or len(ledger) < cfg["max_ledger"]):
                ledger.append(dict(qname=name, feature=feature, candidates=list(evidence),
                                   best_pair_score=best, candidate_evidence=evidence))
            elif keep:
                audit["ledger_rows_omitted_by_cap"] += 1
    for key in ("qname_groups", "alignment_records", "accepted_fragments", "unexplained_fragments", "ledger_rows_omitted_by_cap"):
        audit.setdefault(key, 0)
    def distribution(histogram):
        total = sum(histogram.values())
        return dict(count=total, minimum=min(histogram) if total else None,
                    maximum=max(histogram) if total else None,
                    mean=sum(value*count for value, count in histogram.items())/total if total else None,
                    histogram={str(value): count for value, count in sorted(histogram.items())})
    audit_result = dict(audit)
    audit_result["technical_summary"] = dict(
        selection="one_deterministic_max_score_pair_placement_per_accepted_QNAME; ties_by_candidate_and_coordinates",
        scope="accepted_local_pairs_only; not_whole_library_or_independent_molecules",
        mate_query_lengths=distribution(read_lengths), absolute_template_lengths=distribution(insert_lengths))
    return dict(schema=SCHEMA, config=cfg, count_unit="one_QNAME_sequencing_fragment",
                scope="conditional_local_alignment_file; not_whole_community_denominator",
                counts=dict(sorted(counts.items())), audit=audit_result, fragment_ledger=ledger,
                candidate_coverage={name: dict(length=metadata[name]["length"],
                    eligible_fragments=candidate_fragments[name], near_best_fragments=near_fragments[name],
                    eligible_base_breadth=sum(coverage[name])/len(coverage[name]),
                    near_best_base_breadth=sum(near_coverage[name])/len(near_coverage[name])) for name in sorted(metadata)})


def _count_map(values, name):
    if not isinstance(values, dict):
        raise ValueError(name + " must be a feature-count dict")
    for key, value in values.items():
        if not isinstance(key, str) or not key or isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError("Invalid " + name)
    return {key: float(value) for key, value in values.items() if value}


def fit_identifiability(counts, profiles, config=None):
    """TV-optimal simplex fit and min/max ranges at optimal TV + fixed slack.

    profiles maps candidate ID to raw feature counts or nonnegative weights.
    Unknown accepted feature classes remain in the fit with A-row zero.
    Explicit unexplained local QNAMEs are reported separately and excluded.
    """
    import numpy as np
    from scipy import sparse
    from scipy.optimize import linprog

    cfg = {"tv_slack": .02, "numerical_tolerance": 1e-9, "candidate_order": None,
           "candidate_groups": {}, "diagnose_profiles": True, **(config or {})}
    if type(cfg["diagnose_profiles"]) is not bool:
        raise ValueError("diagnose_profiles must be boolean")
    for key in ("tv_slack", "numerical_tolerance"):
        val = cfg[key]
        if isinstance(val, bool) or not isinstance(val, (float, int)) or not math.isfinite(val) or not 0 <= val <= 1:
            raise ValueError("Invalid " + key)
    if not isinstance(profiles, dict) or not profiles:
        raise ValueError("Nonempty candidate profiles required")
    if any(not isinstance(candidate, str) or not candidate or candidate == UNEXPLAINED for candidate in profiles):
        raise ValueError("Invalid candidate profile ID")
    order = sorted(profiles) if cfg["candidate_order"] is None else cfg["candidate_order"]
    if not isinstance(order, (list, tuple)) or len(order) != len(set(order)) or set(order) != set(profiles):
        raise ValueError("candidate_order must list every profile exactly once")
    observed = _count_map(counts, "counts")
    if any(not value.is_integer() for value in observed.values()):
        raise ValueError("Observed counts must be integral QNAME fragment counts")
    normalized = {}
    discarded = {}
    for candidate in order:
        profile = _count_map(profiles[candidate], "profile")
        discarded[candidate] = profile.pop(UNEXPLAINED, 0)
        total = sum(profile.values())
        if total <= 0:
            raise ValueError("Profile has no accepted local fragments: " + candidate)
        normalized[candidate] = {key: value/total for key, value in profile.items()}
    groups = cfg["candidate_groups"]
    if not isinstance(groups, dict):
        raise ValueError("candidate_groups must be a mapping")
    for name, members in groups.items():
        if not isinstance(name, str) or not name or not isinstance(members, (list, tuple)) or not members or len(set(members)) != len(members) or not set(members) <= set(order):
            raise ValueError("Invalid candidate group")
    unexplained = observed.pop(UNEXPLAINED, 0)
    n = sum(observed.values())
    features = sorted(set(observed).union(*(set(p) for p in normalized.values())))
    matrix = np.array([[normalized[c].get(f, 0) for c in order] for f in features], dtype=float)
    tolerance = cfg["numerical_tolerance"]
    duplicate_groups, consumed = [], set()
    for j, candidate in enumerate(order):
        if candidate in consumed:
            continue
        group = [order[k] for k in range(j, len(order)) if np.max(np.abs(matrix[:, j]-matrix[:, k])) <= tolerance]
        consumed.update(group)
        if len(group) > 1:
            duplicate_groups.append(group)
    result = dict(schema=SCHEMA, status="no_accepted_evidence" if not n else "fitted_descriptive",
                  candidate_order=list(order), count_unit="one_QNAME_sequencing_fragment",
                  weight_unit="fraction_of_accepted_local_fragment_signal",
                  scope="fixed_candidate_conditional_local_fit_not_cell_or_strain_abundance",
                  calibration="none; TV-slack ranges are descriptive feasible ranges, not confidence intervals",
                  confidence_interval=None, p_value=None, biological_call=None,
                  accepted_fragment_count=n, unexplained_local_fragment_count=unexplained,
                  unexplained_local_fraction=unexplained/(n+unexplained) if n+unexplained else None,
                  profile_unexplained_counts=discarded, feature_count=len(features),
                  duplicate_profile_groups=duplicate_groups, config=cfg,
                  bestfit_tv=None, bestfit_weights=None, weight_ranges=None, group_ranges=None,
                  unmodeled_feature_count=sum(value for key, value in observed.items() if all(key not in p for p in normalized.values())),
                  unmodeled_feature_fraction=None, convex_hull_distance_tv={},
                  convex_hull_diagnostics_status="computed" if cfg["diagnose_profiles"] else "not_computed")
    def lp(target, design, objective=None, budget=None):
        rows, cols = design.shape
        design_sparse, negative_eye = sparse.csr_matrix(design), -sparse.eye(rows, format="csr")
        ub = sparse.vstack([sparse.hstack([design_sparse, negative_eye]),
                            sparse.hstack([-design_sparse, negative_eye])], format="csr")
        rhs = np.concatenate([target, -target])
        if budget is not None:
            budget_row = sparse.csr_matrix(np.concatenate([np.zeros(cols), np.ones(rows)*.5]).reshape(1, -1))
            ub = sparse.vstack([ub, budget_row], format="csr")
            rhs = np.append(rhs, budget)
        cost = np.concatenate([np.zeros(cols), np.ones(rows)*.5]) if objective is None else np.concatenate([objective, np.zeros(rows)])
        solution = linprog(cost, A_ub=ub, b_ub=rhs,
                           A_eq=np.array([np.concatenate([np.ones(cols), np.zeros(rows)])]),
                           b_eq=[1.], bounds=(0, None), method="highs")
        if not solution.success:
            raise RuntimeError("Identifiability LP failed: " + solution.message)
        return solution, cols
    if cfg["diagnose_profiles"]:
        for j, candidate in enumerate(order):
            if len(order) == 1:
                result["convex_hull_distance_tv"][candidate] = None
            else:
                other = np.delete(matrix, j, axis=1)
                solution, _ = lp(matrix[:, j], other)
                result["convex_hull_distance_tv"][candidate] = float(max(0, solution.fun))
    if not n:
        return result
    y = np.array([observed.get(f, 0)/n for f in features])
    optimum, columns = lp(y, matrix)
    best = float(max(0, optimum.fun))
    budget = min(1., best+cfg["tv_slack"]+tolerance)
    weights = np.maximum(0., optimum.x[:columns])
    weights /= weights.sum()
    result.update(bestfit_tv=best, bestfit_weights={c: float(w) for c, w in zip(order, weights)},
                  tv_budget=budget, unmodeled_feature_fraction=result["unmodeled_feature_count"]/n)
    def interval(members):
        objective = np.array([1. if c in members else 0. for c in order])
        lower, _ = lp(y, matrix, objective, budget)
        upper, _ = lp(y, matrix, -objective, budget)
        return [float(max(0., min(1., lower.fun))), float(max(0., min(1., -upper.fun)))]
    result["weight_ranges"] = {c: interval([c]) for c in order}
    result["group_ranges"] = {name: interval(members) for name, members in groups.items()}
    return result
