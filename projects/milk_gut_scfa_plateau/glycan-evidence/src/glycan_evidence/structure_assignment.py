"""Experimental joint feature/source assignment with explicit ambiguity.

The LP couples empirical fragment allocations to fixed reference profiles.
Allocation bounds are model-dependent fractional fragments, not observed
class-exclusive QNAMEs, confidence intervals, or physical long-range phasing.
"""
from __future__ import annotations

import json
import math

SCHEMA = "glycan-evidence.structure-assignment.research.v2"
UNEXPLAINED = "__unexplained__"


def _values(values, label, integral=False):
    if not isinstance(values, dict):
        raise ValueError(label + " must be a feature-count dictionary")
    result = {}
    for feature, count in values.items():
        if (not isinstance(feature, str) or not feature or isinstance(count, bool) or
                not isinstance(count, (float, int)) or not math.isfinite(count) or count < 0 or
                (integral and not float(count).is_integer())):
            raise ValueError("Invalid " + label)
        if count:
            result[feature] = float(count)
    return result


def _parse_feature(feature, metadata):
    try:
        rows = json.loads(feature)
    except (TypeError, ValueError) as exc:
        raise ValueError("Feature is not JSON") from exc
    if not isinstance(rows, list) or not rows:
        raise ValueError("Accepted feature must contain known candidate alignments")
    evidence = {}
    for row in rows:
        if not isinstance(row, list) or len(row) not in {1, 3}:
            raise ValueError("Invalid feature entry")
        candidate = row[0]
        if not isinstance(candidate, str) or candidate not in metadata or candidate in evidence:
            raise ValueError("Unknown or duplicate feature candidate")
        direct, bridge = ([], []) if len(row) == 1 else row[1:]
        for positions in (direct, bridge):
            if (not isinstance(positions, list) or any(type(x) is not int for x in positions) or
                    len(set(positions)) != len(positions) or
                    not set(positions) <= set(metadata[candidate]["junctions"])):
                raise ValueError("Invalid candidate-relative feature junction")
        evidence[candidate] = (set(direct), set(bridge))
    return evidence


def fit_assignment(counts, profiles, metadata, config=None):
    """Fit joint-TV allocation and descriptive group/boundary lower bounds.

    counts uses the full JSON feature keys from v1 extract_fragment_features.
    profiles maps each metadata candidate to accepted-feature counts/weights.
    The primary edge mask is observed-candidate OR positive profile probability.
    Literal __unexplained__ is reported separately and excluded from this fit.
    """
    import numpy as np
    from scipy import sparse
    from scipy.optimize import linprog

    cfg = {"tv_slack": .02, "numerical_tolerance": 1e-9, "candidate_order": None,
           "candidate_groups": None, **(config or {})}
    for name in ("tv_slack", "numerical_tolerance"):
        value = cfg[name]
        if isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("Invalid " + name)
    if not isinstance(metadata, dict) or not metadata or not isinstance(profiles, dict) or set(metadata) != set(profiles):
        raise ValueError("Nonempty identical candidate metadata/profile universes required")
    for candidate, item in metadata.items():
        if not isinstance(candidate, str) or not candidate or candidate == UNEXPLAINED or not isinstance(item, dict):
            raise ValueError("Invalid candidate metadata")
        length, junctions = item.get("length"), item.get("junctions")
        if (type(length) is not int or length < 1 or not isinstance(junctions, (list, tuple)) or
                len(junctions) != 2 or any(type(p) is not int for p in junctions) or
                not 0 < junctions[0] < junctions[1] < length):
            raise ValueError("Two ordered candidate-relative inner boundaries required")
    order = sorted(metadata) if cfg["candidate_order"] is None else cfg["candidate_order"]
    if not isinstance(order, (tuple, list)) or len(order) != len(set(order)) or set(order) != set(metadata):
        raise ValueError("candidate_order must list every candidate exactly once")
    groups = cfg["candidate_groups"]
    if groups is None:
        groups = {}
        for candidate in order:
            name = metadata[candidate].get("class_id", candidate)
            groups.setdefault(name, []).append(candidate)
    if not isinstance(groups, dict) or not groups:
        raise ValueError("Nonempty candidate groups required")
    for name, members in groups.items():
        if (not isinstance(name, str) or not name or not isinstance(members, (tuple, list)) or
                not members or len(members) != len(set(members)) or not set(members) <= set(order)):
            raise ValueError("Invalid candidate group")
    observed = _values(counts, "observed integral fragment counts", True)
    unexplained = observed.pop(UNEXPLAINED, 0)
    normalized, discarded = {}, {}
    for candidate in order:
        values = _values(profiles[candidate], "candidate profile")
        discarded[candidate] = values.pop(UNEXPLAINED, 0)
        total = math.fsum(values.values())
        if not total:
            raise ValueError("Candidate profile has no accepted local evidence")
        normalized[candidate] = {feature: value/total for feature, value in values.items()}
    features = sorted(set(observed).union(*(set(p) for p in normalized.values())))
    parsed = [_parse_feature(feature, metadata) for feature in features]
    n = math.fsum(observed.values())
    A = np.array([[normalized[candidate].get(feature, 0.) for candidate in order] for feature in features])
    O = np.array([[candidate in feature for candidate in order] for feature in parsed], dtype=bool)
    E = O | (A > 0)
    counts_vector = np.array([observed.get(feature, 0.) for feature in features])
    y = counts_vector/n if n else np.zeros(len(features))
    group_members = {name: set(members) for name, members in groups.items()}
    group_indices = {name: [i for i, candidate in enumerate(order) if candidate in members]
                     for name, members in group_members.items()}
    boundary_masks = {}
    observed_support = {}
    for group, members in group_members.items():
        support = {"informative_qnames": 0}
        for side_index, side in enumerate(("left", "right")):
            mask, direct_mask = np.zeros_like(E), np.zeros_like(E)
            possible = exclusive = direct_possible = direct_exclusive = bridge_possible = bridge_exclusive = 0
            for f, evidence in enumerate(parsed):
                supported, direct_supported, bridge_supported = set(), set(), set()
                for c, candidate in enumerate(order):
                    if candidate not in members or candidate not in evidence:
                        continue
                    boundary = metadata[candidate]["junctions"][side_index]
                    direct, bridge = evidence[candidate]
                    if boundary in direct:
                        direct_supported.add(candidate)
                        direct_mask[f, c] = True
                    if boundary in bridge:
                        bridge_supported.add(candidate)
                    if boundary in direct or boundary in bridge:
                        supported.add(candidate)
                        mask[f, c] = True
                count = int(counts_vector[f])
                possible += count * bool(supported)
                direct_possible += count * bool(direct_supported)
                bridge_possible += count * bool(bridge_supported)
                exclusive += count * bool(set(evidence) <= supported)
                direct_exclusive += count * bool(set(evidence) <= direct_supported)
                bridge_exclusive += count * bool(set(evidence) <= bridge_supported)
            boundary_masks[(group, side)] = mask
            boundary_masks[(group, side, "direct")] = direct_mask
            support[side] = dict(possible_qnames=possible, exclusive_qnames=exclusive,
                                 direct_possible_qnames=direct_possible, direct_exclusive_qnames=direct_exclusive,
                                 bridged_possible_qnames=bridge_possible, bridged_exclusive_qnames=bridge_exclusive)
        support["informative_qnames"] = sum(int(counts_vector[f]) for f, feature in enumerate(parsed) if set(feature) <= members)
        observed_support[group] = support
    retention = []
    for f, feature in enumerate(features):
        for c, candidate in enumerate(order):
            if A[f, c] > 0 and not O[f, c]:
                retention.append(dict(feature=feature, candidate=candidate, probability=float(A[f, c])))
    result = dict(schema=SCHEMA, status="no_accepted_evidence" if not n else "fitted_descriptive",
                  config={**cfg, "candidate_order": list(order), "candidate_groups": groups},
                  allocation_mask="observed_near_best_candidate_OR_positive_generating_profile",
                  count_unit="one_QNAME_sequencing_fragment", accepted_fragment_count=int(n),
                  unexplained_local_fragment_count=int(unexplained),
                  unexplained_local_fraction=unexplained/(n+unexplained) if n+unexplained else None,
                  profile_unexplained_counts=discarded, candidate_order=list(order),
                  feature_count=len(features), observed_edge_count=int(O.sum()), allowed_edge_count=int(E.sum()),
                  profile_only_edges=retention, profile_only_probability_by_candidate={candidate: float(A[~O[:, c], c].sum()) for c, candidate in enumerate(order)},
                  observed_support=observed_support, group_ranges=None, boundary_bounds=None,
                  bestfit_weights=None, joint_bestfit_tv=None, marginal_tv_at_joint_bestfit=None,
                  marginal_tv_definition="evaluated_at_joint_optimal_weights; not_independently_minimized_marginal_TV",
                  bestfit_profile_only_allocation_mass=None, tv_budget=None,
                  same_mask_no_profile=None, observed_mask_no_profile=None,
                  confidence_interval=None, p_value=None, biological_call=None,
                  scope="conditional_model_source_allocation; not_exclusive_QNAME_counts_or_physical_phasing_or_organism_abundance",
                  calibration="none; fixed joint-TV slack yields descriptive feasible bounds",
                  boundary_count_unit="model_dependent_fractional_fragment_allocation; never_round_to_observed_exclusive_QNAMEs")
    if not n:
        return result

    def allocation_only(mask):
        answer = {"group_ranges": {}, "boundary_bounds": {},
                  "mask": "union_observed_or_profile" if mask is E else "observed_near_best_only"}
        for group, indices in group_indices.items():
            member_columns = np.zeros(len(order), dtype=bool)
            member_columns[indices] = True
            any_member = (mask & member_columns).any(axis=1)
            all_member = ~(mask & ~member_columns).any(axis=1)
            answer["group_ranges"][group] = [float(y[all_member].sum()), float(y[any_member].sum())]
            answer["boundary_bounds"][group] = {}
            for side in ("left", "right"):
                values = {}
                for suffix in (None, "direct"):
                    B = boundary_masks[(group, side) if suffix is None else (group, side, suffix)]
                    forced = ~(mask & ~B).any(axis=1)
                    possible = (mask & B).any(axis=1)
                    prefix = "" if suffix is None else "direct_"
                    lo, hi = float(y[forced].sum()), float(y[possible].sum())
                    values.update({prefix+"lower_mass": lo, prefix+"upper_mass": hi,
                                   prefix+"lower_fractional_fragments": lo*n, prefix+"upper_fractional_fragments": hi*n})
                answer["boundary_bounds"][group][side] = values
        return answer
    result["same_mask_no_profile"] = allocation_only(E)
    result["observed_mask_no_profile"] = allocation_only(O)

    # z_e, w_c, t_e; forbidden edges have both A=0 and O=0 and add zero L1 mass.
    edge_f, edge_c = np.nonzero(E)
    edge_A = A[edge_f, edge_c]
    edges, candidates, rows = len(edge_f), len(order), len(features)
    width = 2*edges+candidates
    z_index = np.arange(edges)
    w_index = edges+np.arange(candidates)
    t_index = edges+candidates+np.arange(edges)
    eq_rows = np.concatenate([edge_f, rows+edge_c, rows+np.arange(candidates), np.repeat(rows+candidates, candidates)])
    eq_cols = np.concatenate([z_index, z_index, w_index, w_index])
    eq_data = np.concatenate([np.ones(edges), np.ones(edges), -np.ones(candidates), np.ones(candidates)])
    equalities = sparse.coo_matrix((eq_data, (eq_rows, eq_cols)), shape=(rows+candidates+1, width)).tocsr()
    equality_rhs = np.concatenate([y, np.zeros(candidates), [1.]])
    ineq_rows = np.concatenate([z_index, z_index, z_index, edges+z_index, edges+z_index, edges+z_index])
    ineq_cols = np.concatenate([z_index, edges+edge_c, t_index, z_index, edges+edge_c, t_index])
    ineq_data = np.concatenate([np.ones(edges), -edge_A, -np.ones(edges), -np.ones(edges), edge_A, -np.ones(edges)])
    inequalities = sparse.coo_matrix((ineq_data, (ineq_rows, ineq_cols)), shape=(2*edges, width)).tocsr()
    objective_D = np.zeros(width); objective_D[t_index] = .5
    validation_tolerance = max(1e-8, 10*cfg["numerical_tolerance"])
    solver_tolerance = min(1e-7, max(1e-10, cfg["numerical_tolerance"]))
    validations = []
    def solve(objective, budget=None):
        constraints, rhs = inequalities, np.zeros(2*edges)
        if budget is not None:
            constraints = sparse.vstack([inequalities, sparse.csr_matrix(objective_D.reshape(1, -1))], format="csr")
            rhs = np.append(rhs, budget)
        solution = linprog(objective, A_ub=constraints, b_ub=rhs, A_eq=equalities, b_eq=equality_rhs,
                           bounds=(0, None), method="highs", options={"primal_feasibility_tolerance": solver_tolerance,
                                                                      "dual_feasibility_tolerance": solver_tolerance})
        if not solution.success:
            raise RuntimeError("Joint assignment LP failed: " + solution.message)
        x = solution.x
        z, w, t = x[z_index], x[w_index], x[t_index]
        equality_error = float(np.max(np.abs(equalities@x-equality_rhs)))
        inequality_error = float(max(0., np.max(constraints@x-rhs)))
        nonnegative_error = float(max(0., -np.min(x)))
        actual_D = float(.5*np.abs(z-edge_A*w[edge_c]).sum())
        if (max(equality_error, inequality_error, nonnegative_error) > validation_tolerance or
                actual_D > float(.5*t.sum())+validation_tolerance or
                (budget is not None and actual_D > budget+validation_tolerance) or
                abs(float(objective@x)-float(solution.fun)) > validation_tolerance):
            raise RuntimeError("Joint assignment LP solution failed residual validation")
        validations.append((equality_error, inequality_error, nonnegative_error))
        return solution, actual_D

    optimal, best = solve(objective_D)
    if abs(optimal.fun-best) > validation_tolerance:
        raise RuntimeError("Joint objective epigraph does not match optimal discrepancy")
    budget = min(1., best+cfg["tv_slack"]+cfg["numerical_tolerance"])
    point = optimal.x
    point_z, point_w = point[z_index], point[w_index]
    result.update(joint_bestfit_tv=best, tv_budget=budget,
                  marginal_tv_at_joint_bestfit=float(.5*np.abs(y-A@point_w).sum()),
                  bestfit_weights={candidate: float(max(0., point_w[c])) for c, candidate in enumerate(order)},
                  bestfit_profile_only_allocation_mass=float(point_z[~O[edge_f, edge_c]].sum()),
                  group_ranges={}, boundary_bounds={})
    if result["marginal_tv_at_joint_bestfit"] > best+validation_tolerance:
        raise RuntimeError("Marginal discrepancy exceeds joint discrepancy")
    bound_cache = {}
    def bound(objective, maximize=False):
        cache_key = (tuple(np.flatnonzero(objective)), maximize)
        if cache_key in bound_cache:
            return bound_cache[cache_key]
        solution, _ = solve(-objective if maximize else objective, budget)
        value = -solution.fun if maximize else solution.fun
        bound_cache[cache_key] = float(max(0., min(1., value)))
        return bound_cache[cache_key]
    for group, indices in group_indices.items():
        objective = np.zeros(width); objective[w_index[indices]] = 1.
        result["group_ranges"][group] = [bound(objective), bound(objective, True)]
        result["boundary_bounds"][group] = {}
        for side in ("left", "right"):
            values = {}
            for suffix in (None, "direct"):
                B = boundary_masks[(group, side) if suffix is None else (group, side, suffix)]
                objective = np.zeros(width); objective[z_index] = B[edge_f, edge_c]
                prefix = "" if suffix is None else "direct_"
                possible = observed_support[group][side]["possible_qnames" if suffix is None else "direct_possible_qnames"]
                lower = bound(objective) if possible else 0.
                point_value = float(objective@point)
                values.update({prefix+"lower_mass": lower, prefix+"lower_fractional_fragments": lower*n,
                               prefix+"bestfit_mass": point_value, prefix+"bestfit_fractional_fragments": point_value*n})
            result["boundary_bounds"][group][side] = values
    result["solver_validation"] = dict(solves=len(validations), unique_bound_objectives=len(bound_cache), tolerance=validation_tolerance,
        solver_feasibility_tolerance=solver_tolerance,
        max_equality_residual=max(v[0] for v in validations), max_inequality_residual=max(v[1] for v in validations),
        max_nonnegative_violation=max(v[2] for v in validations))
    return result
