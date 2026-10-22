"""Prespecified support summaries, independent of benchmark source truth."""
from collections import defaultdict


def summarize(extraction, fit, metadata, gate):
    groups = defaultdict(list)
    for candidate, item in metadata.items():
        groups[item['class_id']].append(candidate)
    support = {name:{'members':members,'left':set(),'right':set(),'left_direct':set(),
        'right_direct':set()} for name,members in groups.items()}
    for fragment in extraction['fragment_ledger']:
        names = fragment['candidates']
        classes = {metadata[c]['class_id'] for c in names}
        if len(classes) != 1:
            continue
        group = next(iter(classes))
        for i, side in enumerate(('left','right')):
            evidence = [fragment['candidate_evidence'][c] for c in names]
            boundaries = [metadata[c]['junctions'][i] for c in names]
            if all(b in e['direct_junctions'] or b in e['bridged_junctions'] for b,e in zip(boundaries,evidence)):
                support[group][side].add(fragment['qname'])
            if all(b in e['direct_junctions'] for b,e in zip(boundaries,evidence)):
                support[group][side+'_direct'].add(fragment['qname'])
    n = extraction['audit']['accepted_fragments']
    tv = fit['bestfit_tv']
    global_pass = n >= gate['min_accepted_fragments'] and tv is not None and tv <= gate['max_bestfit_tv']
    called, ambiguous, rows = [], [], {}
    for group, values in support.items():
        interval = (fit.get('group_ranges') or {}).get(group)
        lower, upper = interval if interval else (None,None)
        supported = (global_pass and lower is not None and lower > gate['min_group_lower_weight_exclusive']
            and all(len(values[side]) >= gate['min_class_exclusive_boundary_qnames'] for side in ('left','right')))
        if supported:
            called.append(group)
        elif upper is not None and upper > gate['min_group_lower_weight_exclusive']:
            ambiguous.append(group)
        row = {k:len(v) if isinstance(v,set) else v for k,v in values.items()}
        row.update(lower_weight=lower,upper_weight=upper,supported=supported)
        for field in ('eligible_base_breadth','near_best_base_breadth'):
            row[field] = max(extraction['candidate_coverage'][c][field] for c in values['members'])
        rows[group] = row
    baselines = {}
    for field in ('eligible_base_breadth','near_best_base_breadth'):
        baselines[field] = {str(cutoff):sorted(g for g,row in rows.items() if row[field]>=cutoff)
                           for cutoff in gate['baseline_breadth_cutoffs']}
    baselines['matched_support']={}
    for field in ('eligible_base_breadth','near_best_base_breadth'):
        baselines['matched_support'][field]={str(cutoff):sorted(g for g,row in rows.items()
            if n >= gate['min_accepted_fragments'] and row[field]>=cutoff and
            all(row[side]>=gate['min_class_exclusive_boundary_qnames'] for side in ('left','right')))
            for cutoff in gate['baseline_breadth_cutoffs']}
    decision = dict(status='supported_classes' if called else 'unresolved',supported_classes=sorted(called),
        ambiguous_classes=sorted(ambiguous),global_gate_pass=global_pass,accepted_fragments=n,
        bestfit_tv=tv,classes=rows,scope='descriptive_reference_architecture_support')
    return decision,baselines
