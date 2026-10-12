"""Sealed-metrics bridge to the separately frozen local-support candidate."""
from datetime import datetime, timezone
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

from . import __version__, compatibility_io as baseline, identity, integrity
from .integrity import digest, write_json

VERSION = 'local-reference-support.research.candidate20261011'
PAYLOAD = frozenset({'CANDIDATE_MODEL.json', 'LOCAL_REFERENCE_ORDER.json',
                     'local_support_api.py', 'local_candidate_api.py'})


def _verified_candidate(bundle, expected):
    """Read data only; never import an unverified candidate payload."""
    bundle = Path(bundle).resolve()
    baseline._require(baseline._sha(expected), 'External candidate manifest SHA256 required')
    baseline._require(digest(bundle/'LOCAL_MANIFEST.json') == expected, 'Candidate manifest changed')
    manifest = identity._json(bundle/'LOCAL_MANIFEST.json')
    baseline._require(manifest.get('version') == VERSION, 'Unsupported local-support candidate')
    files = manifest.get('files_sha256')
    allowed = PAYLOAD | {'baseline/'+name for name in baseline.PAYLOAD | {'MANIFEST.json'}}
    baseline._require(type(files) is dict and set(files) == allowed, 'Incomplete or unrecognized candidate payload')
    for name, expected_file in files.items():
        path = bundle/name
        baseline._require(baseline._sha(expected_file) and path.is_file() and not path.is_symlink()
                          and path.resolve().is_relative_to(bundle) and digest(path) == expected_file,
                          'Candidate payload changed: '+name)
    base, base_manifest, markers, mapping = baseline._verified_bundle(
        bundle/'baseline', manifest.get('baseline_manifest_sha256'))
    return bundle, manifest, base, base_manifest, markers, mapping


def run_local_support(metrics, bundle, expected_manifest_sha256, out):
    """Derive query N/counts/hashes from complete metrics, without user overrides."""
    destination = Path(out).resolve()
    if destination.exists():
        raise FileExistsError(str(destination))
    bundle, manifest, base, base_manifest, markers, mapping = _verified_candidate(bundle, expected_manifest_sha256)
    query, inputs, preprocessing = baseline._measured_query(metrics, base, base_manifest, markers, mapping)
    destination.mkdir(parents=True, exist_ok=False)
    receipt = dict(schema='glycan-evidence.local-support-run.v1', status='RUNNING',
                   package_version=__version__, started_utc=datetime.now(timezone.utc).isoformat(),
                   inputs=inputs, preprocessing=preprocessing, code=baseline._record(__file__),
                   helper_code={name:baseline._record(module.__file__) for name,module in
                                [('compatibility_io',baseline),('identity',identity),('integrity',integrity)]},
                   candidate_manifest_sha256=expected_manifest_sha256, commands=[], identity_call=None,
                   scope='frozen development candidate; conditional local reference support; not a taxonomic classifier')
    write_json(destination/'receipt.json', receipt)
    try:
        snapshot = destination/'verified_candidate'
        snapshot.mkdir()
        for name in sorted(set(manifest['files_sha256']) | {'LOCAL_MANIFEST.json'}):
            target = snapshot/name
            target.parent.mkdir(exist_ok=True, parents=True)
            shutil.copyfile(bundle/name, target)
        _verified_candidate(snapshot, expected_manifest_sha256)
        write_json(destination/'query.json', query)
        argv = [sys.executable, '-I', '-B', '-X', 'utf8', str(snapshot/'local_candidate_api.py'),
                '--bundle', str(snapshot), '--expected-manifest-sha256', expected_manifest_sha256,
                '--query', str(destination/'query.json'), '--out', str(destination/'result.json')]
        with (destination/'api.stdout.log').open('x', encoding='utf-8') as stdout, \
             (destination/'api.stderr.log').open('x', encoding='utf-8') as stderr:
            completed = subprocess.run(argv, cwd=destination, stdout=stdout, stderr=stderr, check=False)
        receipt['commands'].append({'argv':argv, 'returncode':completed.returncode})
        if completed.returncode:
            raise RuntimeError(f'Local-support API failed ({completed.returncode}); inspect api.stderr.log')
        result = identity._json(destination/'result.json')
        baseline._require(result.get('version') == VERSION and
                          result.get('candidate_manifest_sha256') == expected_manifest_sha256 and
                          result.get('baseline_manifest_sha256') == manifest['baseline_manifest_sha256'] and
                          result.get('input_binding_sha256') == hashlib.sha256(baseline._canonical(query).encode()).hexdigest(),
                          'Prediction does not bind the verified query and candidate')
        baseline._require(result.get('input_pairs') == query['input_pairs'] and
                          result.get('composite_support') == sum(query['counts']), 'Prediction denominator mismatch')
        model = identity._json(snapshot/'CANDIDATE_MODEL.json')
        baseline._require(result.get('local_support') == sum(query['counts'][j] for j in model['Z_indices']),
                          'Prediction local-support mismatch')
        observation = result.get('observation', {})
        baseline._require(observation.get('state') in {'conditional_target_local_support','insufficient_evidence'}
                          and type(observation.get('reasons')) is list, 'Invalid observation state')
        baseline._require(observation['state'] != 'conditional_target_local_support' or not observation['reasons'],
                          'Conditional support contradicts abstention reasons')
        baseline._require(all(k in obj and obj[k] is None for obj in [result, observation]
                          for k in ['identity_call','p_value','confidence_interval','organism_abundance']),
                          'Prediction exceeded research scope')
        _verified_candidate(snapshot, expected_manifest_sha256)
        for source in inputs.values():
            baseline._require(baseline._record(source['path']) == source, 'Input changed during execution')
        for source in [receipt['code'], *receipt['helper_code'].values()]:
            baseline._require(baseline._record(source['path']) == source, 'Bridge helper changed during execution')
        receipt.update(status='COMPLETE', completed_utc=datetime.now(timezone.utc).isoformat(),
                       state=observation['state'], input_binding_sha256=result['input_binding_sha256'],
                       files={p.relative_to(destination).as_posix():baseline._record(p)
                              for p in sorted(destination.rglob('*')) if p.is_file() and p.name != 'receipt.json'})
        write_json(destination/'receipt.json', receipt)
        (destination/'COMPLETE').write_text(digest(destination/'receipt.json')+'\n', encoding='ascii')
    except BaseException as error:
        receipt.update(status='FAILED', error_type=type(error).__name__, error=str(error))
        (destination/'COMPLETE').unlink(missing_ok=True)
        write_json(destination/'receipt.json', receipt)
        raise
    return receipt
