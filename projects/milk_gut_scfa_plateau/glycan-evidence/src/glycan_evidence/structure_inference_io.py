"""File interface for experimental local fragment inference."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from . import __version__
from .structure_inference import extract_fragment_features, fit_identifiability


def _record(path):
    p=Path(path).resolve()
    with p.open('rb') as stream:
        digest=hashlib.file_digest(stream,'sha256').hexdigest()
    return dict(path=str(p),bytes=p.stat().st_size,sha256=digest)


def run_structure_inference(sam,candidates,profiles,config,out):
    return _run_local_fit(sam,candidates,profiles,config,out,assignment=False)


def run_structure_assignment(sam,candidates,profiles,config,out):
    return _run_local_fit(sam,candidates,profiles,config,out,assignment=True)


def _run_local_fit(sam,candidates,profiles,config,out,*,assignment):
    paths=dict(sam=sam,candidates=candidates,profiles=profiles,config=config)
    inputs={k:_record(v) for k,v in paths.items()}
    def load(path):
        return json.loads(Path(path).read_text(encoding='utf-8'))
    metadata,expected,cfg=load(candidates),load(profiles),load(config)
    if set(metadata)!=set(expected):
        raise ValueError('Candidate and expected-profile IDs must match exactly')
    destination=Path(out).resolve()
    destination.mkdir(parents=True,exist_ok=False)
    def put(name,obj):
        (destination/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    receipt=dict(status='RUNNING',schema='structure-assignment-io.v2' if assignment else 'structure-inference-io.v1',package_version=__version__,
        inputs=inputs,started_utc=datetime.now(timezone.utc).isoformat(),
        scope='conditional_source_assignment_and_boundary_lower_bounds' if assignment else 'fixed_candidate_descriptive_local_fragment_fit',biological_call=None)
    put('receipt.json',receipt)
    try:
        extraction=extract_fragment_features(sam,metadata,cfg.get('extraction',{}))
        if assignment:
            from .structure_assignment import fit_assignment
            fit=fit_assignment(extraction['counts'],expected,metadata,cfg.get('fit',{}))
        else:
            fit=fit_identifiability(extraction['counts'],expected,cfg.get('fit',{}))
        put('extraction.json',extraction);put('fit.json',fit)
        for key,path in paths.items():
            if _record(path)!=inputs[key]:
                raise ValueError('Input changed during execution: '+key)
        receipt.update(status='COMPLETE',completed_utc=datetime.now(timezone.utc).isoformat(),
            outputs={name:_record(destination/name) for name in ['extraction.json','fit.json']})
        put('receipt.json',receipt)
        (destination/'COMPLETE').write_text(_record(destination/'receipt.json')['sha256']+'\n',encoding='ascii')
    except BaseException as exc:
        receipt.update(status='FAILED',error_type=type(exc).__name__,error=str(exc))
        put('receipt.json',receipt)
        raise
    return receipt
