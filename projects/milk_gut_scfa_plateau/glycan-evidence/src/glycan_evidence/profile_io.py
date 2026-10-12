"""Portable file interface for descriptive full-marker diagnostics."""
import csv
from dataclasses import fields
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from uuid import uuid4

from . import __version__, profile
from .identity import verify_completed_metrics
from .integrity import write_json


def _read(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f'Duplicate JSON key: {key}')
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding='utf-8-sig'), object_pairs_hook=unique)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_from_dict(value):
    if not isinstance(value, dict) or set(value) != {f.name for f in fields(profile.Source)}:
        raise ValueError('Source requires source_id, cluster_id, positive_targets and exposure')
    if not isinstance(value['positive_targets'], list):
        raise ValueError('positive_targets must be a JSON array')
    return profile.Source(**value)


def load_expected(path):
    value = _read(path)
    required = {f.name for f in fields(profile.ExpectedProfile)} | {
        'schema', 'panel_sha256', 'count_unit', 'scope', 'estimator', 'smoothing'}
    if not isinstance(value, dict) or set(value) != required:
        raise ValueError('Unexpected expected-profile schema')
    constants = dict(schema=profile.SCHEMA, count_unit=profile.COUNT_UNIT, scope=profile.SCOPE,
                     estimator='equal_source_mean_of_within_source_pooled_count_fractions',
                     smoothing='none; unsupported reference probability remains zero')
    if any(value[k] != expected for k, expected in constants.items()):
        raise ValueError('Expected-profile interpretation contract changed')
    if not isinstance(value['panel'], dict) or set(value['panel']) != {f.name for f in fields(profile.MarkerPanel)}:
        raise ValueError('Invalid model panel schema')
    vectors = ['probabilities', 'training_sources', 'training_samples', 'training_input_sha256', 'source_probabilities']
    if any(not isinstance(value[k], list) for k in vectors):
        raise ValueError('Model vectors must be JSON arrays')
    if not isinstance(value['panel']['markers'], list) or any(
            not isinstance(row, list) or len(row) != 2 for row in value['panel']['markers']):
        raise ValueError('Model markers must be ID/length pairs')
    if any(not isinstance(row, list) for row in value['source_probabilities']):
        raise ValueError('Source probabilities must be arrays')
    if any(isinstance(p, bool) or not isinstance(p, (int, float))
           for row in [value['probabilities'], *value['source_probabilities']] for p in row):
        raise ValueError('Probabilities must be numeric')
    panel = profile.MarkerPanel(**value['panel'])
    if panel.sha256 != value['panel_sha256']:
        raise ValueError('Model panel SHA mismatch')
    expected = profile.ExpectedProfile(panel, tuple(value['probabilities']),
                                       tuple(source_from_dict(s) for s in value['training_sources']),
                                       tuple(value['training_samples']), tuple(value['training_input_sha256']),
                                       tuple(tuple(p) for p in value['source_probabilities']))
    # A serialized model must implement the estimator it declares.
    for index, p in enumerate(expected.probabilities):
        mean = sum(row[index] for row in expected.source_probabilities) / len(expected.source_probabilities)
        if not math.isclose(mean, p, rel_tol=1e-10, abs_tol=1e-12):
            raise ValueError('Expected probabilities differ from declared equal-source mean')
    return expected


def run_profile(metrics, model, source, out):
    """Diagnose complete registered-v2 metrics with an explicit source ledger.

    Models are user supplied, hash recorded and source declared. This command
    does not certify their calibration, truth labels or external validity.
    """
    path = Path(metrics).resolve()
    if path.name != 'gene_metrics.tsv':
        raise ValueError('Expected complete-run gene_metrics.tsv')
    sid = path.parent.name
    verified = verify_completed_metrics(path.parent.parent, sid)
    expected = load_expected(model)
    origin = source_from_dict(_read(source))
    with path.open(newline='', encoding='utf-8-sig') as handle:
        rows = [r for r in csv.DictReader(handle, delimiter='\t')
                if r['channel'] == 'identity' and r['taxon_group'] == expected.panel.target
                and r['eligible_identity'] == '1' and r['pseudo'] == '0'
                and int(r['margin']) == expected.panel.margin]
    actual_panel = profile.MarkerPanel(expected.panel.target,
                                       tuple((r['feature_id'], int(r['length'])) for r in rows),
                                       verified['rules']['mapping_contract_sha256'], expected.panel.margin)
    if actual_panel != expected.panel:
        raise ValueError('Model does not match the complete target marker panel/mapping/margin')
    observation = profile.observation_from_rows(actual_panel, origin, sid,
                    [{k: r[k] for k in profile.ROW_FIELDS} for r in rows], input_sha256=_sha(path))
    diagnostic = profile.diagnose(observation, expected)
    destination = Path(out).resolve()
    if destination.exists():
        raise FileExistsError(f'Output path already exists: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = destination.with_name(f'.{destination.name}.running-{uuid4().hex}')
    stage.mkdir()
    receipt = dict(schema='glycan-evidence.profile-run.v1', package_version=__version__,
                   status='RUNNING', sample=sid, scope=profile.SCOPE,
                   inputs={name: {'path': str(Path(p).resolve()), 'sha256': _sha(p)}
                           for name,p in [('metrics', path), ('model', model), ('source', source)]},
                   model_validation='schema, complete panel, estimator and declared source separation; no calibration certification',
                   historical_raw_input_scope=verified['provenance']['raw_input_verification'])
    try:
        write_json(stage/'diagnostic.json', diagnostic)
        write_json(stage/'expected_profile.json', expected.to_dict())
        receipt.update(status='COMPLETE', completed_utc=datetime.now(timezone.utc).isoformat(),
                       diagnostic_sha256=_sha(stage/'diagnostic.json'))
        write_json(stage/'receipt.json', receipt)
        (stage/'COMPLETE').write_text(receipt['completed_utc']+'\n', encoding='utf-8')
        stage.rename(destination)
    except BaseException as error:
        receipt.update(status='FAILED', error_type=type(error).__name__, error=str(error))
        (stage/'COMPLETE').unlink(missing_ok=True)
        write_json(stage/'FAILED.json', receipt)
        raise
    return receipt
