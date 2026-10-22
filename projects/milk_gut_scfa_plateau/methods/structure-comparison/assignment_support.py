"""Frozen v2 decisions and strong support-only comparators; no truth input."""


def summarize_assignment(extraction, fit, metadata, gate):
    n = fit['accepted_fragment_count']
    assert n == extraction['audit']['accepted_fragments']
    tv = fit['joint_bestfit_tv']
    global_pass = n >= gate['min_accepted_fragments'] and tv is not None and tv <= gate['max_joint_tv']
    groups = {}
    for candidate, item in metadata.items():
        groups.setdefault(item['class_id'], []).append(candidate)
    rows = {}
    for group, members in groups.items():
        lo, hi = (fit.get('group_ranges') or {}).get(group, [None, None])
        bounds = (fit.get('boundary_bounds') or {}).get(group, {})
        left = bounds.get('left', {}).get('lower_fractional_fragments')
        right = bounds.get('right', {}).get('lower_fractional_fragments')
        called = bool(global_pass and lo is not None and lo > gate['min_group_lower_weight_exclusive']
                      and left is not None and right is not None
                      and min(left, right) >= gate['min_fractional_boundary_lower'])
        support = fit['observed_support'][group]
        row = dict(members=members, lower_weight=lo, upper_weight=hi,
                   left_lower_fractional_fragments=left, right_lower_fractional_fragments=right,
                   supported=called, observed_support=support)
        for field in ('eligible_base_breadth', 'near_best_base_breadth'):
            row[field] = max(extraction['candidate_coverage'][c][field] for c in members)
        rows[group] = row
    decision = dict(status='supported_classes' if any(x['supported'] for x in rows.values()) else 'unresolved',
                    supported_classes=sorted(g for g, x in rows.items() if x['supported']),
                    accepted_fragments=n, joint_bestfit_tv=tv, global_gate_pass=global_pass, classes=rows,
                    scope='conditional_reference_architecture_assignment_bounds')
    possible = {}
    for field in ('eligible_base_breadth', 'near_best_base_breadth'):
        possible[field] = {}
        for cutoff in gate['baseline_breadth_cutoffs']:
            possible[field][str(cutoff)] = sorted(g for g, row in rows.items()
                if n >= gate['min_accepted_fragments'] and row[field] >= cutoff
                and row['observed_support']['informative_qnames'] >= gate['min_informative_qnames']
                and all(row['observed_support'][side]['possible_qnames'] >= gate['min_possible_boundary_qnames']
                        for side in ('left', 'right')))
    ablations = {}
    for name in ('same_mask_no_profile', 'observed_mask_no_profile'):
        obj = fit.get(name) or {}
        ablations[name] = sorted(group for group, interval in obj.get('group_ranges', {}).items()
            if n >= gate['min_accepted_fragments'] and interval[0] > gate['min_group_lower_weight_exclusive']
            and all(obj['boundary_bounds'][group][side]['lower_fractional_fragments']
                    >= gate['min_fractional_boundary_lower'] for side in ('left', 'right')))
    return decision, {'possible_support': possible}, ablations
