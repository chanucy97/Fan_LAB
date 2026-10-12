"""Evaluate explicit known-truth cases without treating technical repeats as strains."""
import csv
import math
from collections import defaultdict
from pathlib import Path
from .integrity import digest, write_json

TRUTH_FIELDS = {'sample', 'target', 'present', 'genome', 'cluster', 'split', 'stratum'}


def read_tsv(path, required):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames) or not required <= set(reader.fieldnames):
            raise ValueError(f'Missing or duplicate TSV columns in {path}: {sorted(required)}')
        rows = list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError(f'Malformed TSV row in {path}')
    return rows


def summarize(rows):
    counts = dict(TP=0, FN=0, FP=0, TN=0)
    per_cluster = defaultdict(list)
    for r in rows:
        predicted = r['identity_state'] == 'detected'
        positive = r['present'] == '1'
        category = ('TP' if predicted else 'FN') if positive else ('FP' if predicted else 'TN')
        counts[category] += 1
        per_cluster[r['cluster']].append((positive, predicted))
    positive_rates, negative_rates = [], []
    for values in per_cluster.values():
        truth = {x[0] for x in values}
        if len(truth) != 1:
            raise ValueError('One target/independent cluster has contradictory binary truth labels')
        rate = sum(x[1] for x in values) / len(values)
        (positive_rates if next(iter(truth)) else negative_rates).append(rate)
    return dict(**counts, technical_conditions=len(rows),
                unique_genomes=len({r['genome'] for r in rows}), independent_clusters=len(per_cluster),
                condition_sensitivity=counts['TP']/(counts['TP']+counts['FN']) if counts['TP']+counts['FN'] else None,
                condition_false_positive_rate=counts['FP']/(counts['FP']+counts['TN']) if counts['FP']+counts['TN'] else None,
                cluster_macro_sensitivity=sum(positive_rates)/len(positive_rates) if positive_rates else None,
                cluster_macro_false_positive_rate=sum(negative_rates)/len(negative_rates) if negative_rates else None,
                positive_clusters=len(positive_rates), negative_clusters=len(negative_rates))


def evaluate(calls, truth, out, margin=0, split='test'):
    if not isinstance(margin, (int, float)) or not math.isfinite(margin) or margin not in {0, 6, 12, 18}:
        raise ValueError('Margin must be one of 0, 6, 12, 18')
    predictions = read_tsv(calls, {'sample', 'target', 'identity_state', 'margin'})
    labels = read_tsv(truth, TRUTH_FIELDS)
    if not labels:
        raise ValueError('Truth table is empty')
    lookup, seen_truth, cluster_splits, genome_clusters = {}, set(), defaultdict(set), defaultdict(set)
    for r in predictions:
        if r['identity_state'] not in {'detected', 'insufficient_evidence'}:
            raise ValueError('Unsupported identity state')
        try:
            m = float(r['margin'])
        except ValueError as e:
            raise ValueError('Invalid prediction margin') from e
        if not math.isfinite(m) or m not in {0, 6, 12, 18}:
            raise ValueError('Invalid prediction margin')
        if m != margin:
            continue
        key = (r['sample'], r['target'])
        if key in lookup:
            raise ValueError('Duplicate prediction for sample/target/margin')
        lookup[key] = r
    selected, excluded, expected = [], [], set()
    for r in labels:
        if any(not r[k].strip() for k in TRUTH_FIELDS):
            raise ValueError('Truth metadata fields must be explicit and nonempty')
        if r['present'] not in {'0', '1', 'unresolved'}:
            raise ValueError('Truth present must be 0, 1, or unresolved')
        key = (r['sample'], r['target'])
        if key in seen_truth:
            raise ValueError('Duplicate truth sample/target')
        seen_truth.add(key)
        cluster_splits[r['cluster']].add(r['split'])
        genome_clusters[r['genome']].add(r['cluster'])
        if r['split'] != split:
            continue
        expected.add(key)
        if key not in lookup:
            raise ValueError(f'Missing prediction for declared case: {key}')
        if r['present'] == 'unresolved':
            excluded.append(dict(sample=r['sample'], target=r['target'], reason='truth_unresolved'))
        else:
            selected.append({**r, 'identity_state': lookup[key]['identity_state']})
    if any(len(s)>1 for s in cluster_splits.values()):
        raise ValueError('Independent cluster appears in multiple splits')
    if any(len(s)>1 for s in genome_clusters.values()):
        raise ValueError('Genome assigned to multiple independent clusters')
    if not expected or not selected:
        raise ValueError('No evaluable known-truth cases in requested split')
    extra = sorted(set(lookup) - seen_truth)
    if extra:
        raise ValueError('Predictions contain undeclared cases; supply a complete truth ledger')
    groups=defaultdict(list)
    totals=defaultdict(list)
    for r in selected:
        groups[(r['target'],r['stratum'])].append(r)
        totals[r['target']].append(r)
    result=dict(schema_version=1,status='COMPLETE',split=split,margin=margin,
                inputs=dict(calls_sha256=digest(calls),truth_sha256=digest(truth)),
                targets={t:summarize(v) for t,v in sorted(totals.items())},
                strata=[dict(target=t,stratum=s,**summarize(v)) for (t,s),v in sorted(groups.items())],
                declared_sample_target_cases=len(expected),
                declared_unique_samples=len({s for s,t in expected}),truth_unresolved=excluded,
                scope='Known-truth detection evaluation. Insufficient evidence is a nondetection outcome, not a biological absence assertion. Technical conditions are not independent genomes. No technical-replicate confidence interval is reported.')
    out=Path(out)
    out.mkdir(parents=True,exist_ok=False)
    write_json(out/'benchmark.json',result)
    (out/'COMPLETE').write_text('Known-truth report complete\n',encoding='utf-8')
    return result
