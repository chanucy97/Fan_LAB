"""Replay one case from a SHA-bound, truth-free normalized evidence roster."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

from . import __version__
from .structure_recovery import recover_structure, validate_file_binding


def record(path):
    p=Path(path).resolve()
    with p.open('rb') as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(p),bytes=p.stat().st_size,sha256=h)


def run_structure_recovery(roster,case_id,expected_roster_sha256,out):
    if not isinstance(expected_roster_sha256,str) or not re.fullmatch('[a-f0-9]{64}',expected_roster_sha256):
        raise ValueError('Expected lowercase roster SHA256 required')
    source=record(roster)
    if source['sha256']!=expected_roster_sha256:
        raise ValueError('Roster SHA256 mismatch')
    data=json.loads(Path(roster).read_text(encoding='utf-8'))
    if data.get('schema')!='truth-free-existing-assembly-roster.v1':
        raise ValueError('Unsupported recovery roster schema')
    rows=data.get('cases')
    if not isinstance(rows,list) or not rows:
        raise ValueError('Nonempty roster cases required')
    ids=[r.get('case_id') for r in rows]
    if any(not isinstance(x,str) or not x for x in ids) or len(ids)!=len(set(ids)):
        raise ValueError('Unique nonempty case IDs required')
    selected=[r for r in rows if r['case_id']==case_id]
    if len(selected)!=1:raise ValueError('Requested case absent from roster')
    item=selected[0]
    bindings=[data['source_inventory'],*data['panel_sources'],item['v2_source']]
    for candidate in item['assembly_candidates']:
        bindings.extend(candidate['source_bindings'].values())
    by_path={}
    for binding in bindings:
        validate_file_binding(binding)
        path=str(Path(binding['path']).resolve())
        if path in by_path and by_path[path]!=binding:
            raise ValueError('Conflicting bindings for one source file')
        by_path[path]=binding
    original=json.loads(Path(item['v2_source']['path']).read_text(encoding='utf-8'))
    if (original.get('case_id')!=case_id or original.get('status')!='COMPLETE' or
            sorted(original.get('decision_v2',{}).get('supported_classes',[]))!=sorted(item['reference_calls'])):
        raise ValueError('Reference call source does not match the requested case and calls')
    # This is a reproducible normalized-evidence replay, not a fresh mapper or
    # raw-table adapter. Exact declared source artifacts are rehashed here.
    result=recover_structure(case_id,item['reference_calls'],data['panel'],item['assembly_candidates'])
    destination=Path(out).resolve();destination.mkdir(parents=True,exist_ok=False)
    def put(name,value):
        (destination/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    receipt=dict(status='RUNNING',schema='glycan-evidence.structure-recovery-io/1',package_version=__version__,
                 case_id=case_id,roster=source,source_files=list(by_path.values()),
                 started_utc=datetime.now(timezone.utc).isoformat(),
                 scope='SHA_bound_normalized_evidence_replay; adapter_transformation_not_recomputed; no_new_sequence_reconstruction')
    put('receipt.json',receipt)
    try:
        put('recovery.json',result)
        if record(roster)!=source:raise ValueError('Roster changed during execution')
        for binding in by_path.values():validate_file_binding(binding)
        receipt.update(status='COMPLETE',completed_utc=datetime.now(timezone.utc).isoformat(),
                       output=record(destination/'recovery.json'))
        put('receipt.json',receipt)
        (destination/'COMPLETE').write_text(record(destination/'receipt.json')['sha256']+'\n',encoding='ascii')
    except BaseException as error:
        (destination/'COMPLETE').unlink(missing_ok=True)
        receipt.update(status='FAILED',error_type=type(error).__name__,error=str(error));put('receipt.json',receipt)
        raise
    return receipt
