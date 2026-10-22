"""One pre-specified, hash-bound candidate readback; not a variant-discovery tool.

Public adapter v1: sample identity is supplied in reviewed metadata. Importing
performs no IO; a run requires frozen inputs and explicit authorization.
"""
from __future__ import annotations
import argparse
from collections import Counter
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import pickle
import platform
import sys
import traceback
import core


def utc(): return datetime.now(timezone.utc).isoformat()


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()


def write_json(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def canonical(value): return json.dumps(value,sort_keys=True,separators=(",",":"))


def bound_json(root,rec):
    # Deployed immutable metadata live beside the rules, never at arbitrary relative paths.
    if Path(rec["file"]).name != rec["file"]: raise ValueError("Metadata name must be a basename")
    p=root/rec["file"]
    if digest(p)!=rec["sha256"]: raise ValueError("Metadata SHA mismatch: "+str(p))
    return json.loads(p.read_text(encoding="utf-8-sig"))


def verify_authorization(rules,actual_sha,expected_sha,authorized):
    if not authorized or rules.get("status")!="FROZEN_FOR_SINGLE_CANDIDATE_READBACK" or actual_sha!=expected_sha:
        raise ValueError("Requires root-frozen rules, exact supplied SHA and explicit authorized-run flag")
    if rules.get("schema")!="single-candidate-readback-rules.public-adapter.v1" or rules.get("candidate_count")!=1 or not isinstance(rules.get("sample_id"),str) or not rules["sample_id"].strip():
        raise ValueError("Requires exactly one pre-specified candidate and a nonempty sample ID")


def stat_bound(rec):
    p=Path(rec["path"]);s=p.stat()
    if s.st_size!=rec["historic_bytes"] or s.st_mtime_ns!=rec["historic_mtime_ns"]:
        raise ValueError("Original source stat changed: "+str(p))
    return {"path":str(p),"bytes":s.st_size,"mtime_ns":s.st_mtime_ns,"historic_sha256":rec["sha256"],"current_full_sha256":None}


def verify_hashed_record(rec):
    p=Path(rec["path"])
    if p.stat().st_size!=rec["bytes"] or digest(p)!=rec["sha256"]:raise ValueError("Source SHA mismatch: "+str(p))


class PrimitiveUnpickler(pickle.Unpickler):
    def find_class(self,module,name):
        raise pickle.UnpicklingError("Rdic must contain only builtin primitive containers, not globals")


class HashingReader:
    def __init__(self,handle):self.handle=handle;self.hash=hashlib.sha256()
    def read(self,size=-1):
        data=self.handle.read(size);self.hash.update(data);return data
    def readline(self,size=-1):
        data=self.handle.readline(size);self.hash.update(data);return data
    def readinto(self,b):
        n=self.handle.readinto(b);self.hash.update(memoryview(b)[:n]);return n


def load_rdic(path,expected_sha,references):
    # Pickle is not region indexed. Load it once, hashing the very same byte stream.
    with Path(path).open("rb") as f:
        wrapped=HashingReader(f);data=PrimitiveUnpickler(wrapped).load()
        if wrapped.read(1):raise ValueError("Trailing bytes after Rdic pickle")
        if wrapped.hash.hexdigest()!=expected_sha:raise ValueError("Rdic SHA mismatch")
    if type(data) is not dict:raise ValueError("Unexpected Rdic container")
    selected={}
    for ref in references:
        if ref not in data or type(data[ref]) is not set or any(type(q) is not str for q in data[ref]):
            raise ValueError("Original skip-mm Rdic set unavailable/invalid: "+ref)
        selected[ref]=data[ref]
    return selected


def read_cov_window(handle,reference,start,end):
    import numpy as np
    keys=[k for k in handle.keys() if k.rsplit("::",1)[0]==reference]
    if not keys:return None
    if keys != [reference+"::0"]:raise ValueError("Unexpected covT mm buckets")
    ds=handle[keys[0]]
    if ds.ndim!=2 or ds.shape[0]!=2:raise ValueError("Unexpected covT shape")
    # inStrain's special HDF5 stores [values, positions], as gene-v2 consumes it.
    positions=ds[1,:]
    if not np.all(np.isfinite(positions)) or not np.all(positions==np.floor(positions)):
        raise ValueError("Noninteger/nonfinite covT positions")
    indices=np.flatnonzero((positions>=start)&(positions<end))
    values=ds[0,indices] if len(indices) else []
    result={}
    for i,depth in zip(indices,values):
        pos=int(positions[i])
        if pos in result or not np.isfinite(depth) or depth!=int(depth) or depth<=0:
            raise ValueError("Invalid sparse covT window")
        result[pos]=int(depth)
    return result


def parity_report(bam,cov,accepted,windows):
    results=[]
    for w in windows:
        ref=w["reference"];a,b=w["comparison_window0"]
        observed=[]
        for start,end in w["original_splits0_inclusive"]:
            # Preserve original split origin/cap/overlap state, retain only the frozen window.
            for pos,depth in core.original_split_coverage(bam,ref,start,end,accepted[ref]):
                if a<=pos<b:observed.append((pos,depth))
        expected=read_cov_window(cov,ref,a,b)
        outcome=core.compare_covt_depth(observed,expected,a,b)
        outcome.update(w)
        outcome["observed_depth_by_position0"]={str(p):d for p,d in observed}
        outcome["expected_sparse_depth_by_position0"]=None if expected is None else {str(p):d for p,d in expected.items()}
        results.append(outcome)
    return {"all_windows_match":all(r["status"]=="MATCH" for r in results),"windows":results,
            "mismatch_possible_causes_not_assumed":["split/boundary or max_depth context","runtime overlap/quality semantics","source or accepted-set drift"]}


def perform_event_fetch(bam,assays,accepted):
    calls=[]
    for a in assays:
        # Indexed fetch, never whole BAM traversal and never a new alignment.
        for r in bam.fetch(a.reference,a.start0,a.end0):
            if r.cigarstring is None or r.query_sequence is None:raise ValueError("Malformed source alignment in frozen event window")
            raw=core.Read(a.sample_id,r.query_name,a.reference,r.reference_start,r.cigarstring,
                          r.query_sequence,None if r.query_qualities is None else tuple(r.query_qualities),r.flag,r.mapping_quality)
            calls.append(core.classify_read(raw,a,accepted[a.reference]))
    return calls


def interpret_calls(calls,assays,technical,parity):
    if not parity["all_windows_match"]:
        return {"status":"BLOCKED_COVT_PARITY", "technical_label":None,"event_counts_generated":False},[],[]
    main=core.collapse_qnames(calls,"main");diagnostic=core.collapse_qnames(calls,"diagnostic")
    summary=core.summarize_fragments(main,core.TechnicalRules(**technical),[k for k,_ in assays[0].templates])
    summary.update(status="COMPLETE_PARITY_MATCH_READBACK_REPORTED",event_counts_generated=True,
                   diagnostic_qnames=len(diagnostic),diagnostic_record_rdic_membership=dict(Counter(str(c["accepted_original_rdic"]) for c in calls if not c["main_record_eligible"])),
                   diagnostic_counts_excluded_from_all_main_denominators=True,
                   mapping_uniqueness_asserted=False,original_primary_mapping_scope_only=True)
    return summary,main,diagnostic


def validate_metadata(rules,candidate,bindings,selection):
    if bindings["sample_id"]!=rules["sample_id"] or candidate["candidate_id"]!=rules["candidate_id"] or candidate["sample_id"]!=rules["sample_id"]:
        raise ValueError("Frozen candidate/sample identity mismatch")
    if len(selection)!=1 or selection[0]["sample_id"]!=rules["sample_id"] or selection[0]["contig_id"]!=candidate["sample_contig"] or selection[0]["readback_selected"] is not True:
        raise ValueError("Frozen complete-candidate selection mismatch")
    a=core.assays_from_candidate(candidate,**rules["single_read"])
    if canonical([asdict(x) for x in a])!=canonical(rules["assays"]):raise ValueError("Frozen assay differs from sequence receipt")
    lengths={x["reference_scaffold_id"]:x["reference_length"] for x in bindings["reference_panel"]["two_event_background_references"]}
    windows=[]
    for ref in candidate["reference_comparisons"]:
        start,end=ref["reference_window0"];name=ref["subject_id"]
        splits=[list(s) for s in core.original_splits(lengths[name],bindings["split_definition"]["window_length_default"]) if s[0]<end and s[1]>=start]
        windows.append({"reference":name,"reference_length":lengths[name],"comparison_window0":[start,end],"original_splits0_inclusive":splits})
    if canonical(windows)!=canonical(rules["parity_windows"]):raise ValueError("Frozen original split/window definition mismatch")
    return a


def run(args):
    rp=Path(args.rules).resolve();rules=json.loads(rp.read_text(encoding="utf-8-sig"))
    verify_authorization(rules,digest(rp),args.rules_sha256,args.root_authorized_run)
    root=rp.parent
    for name,sha in rules["implementation_sha256"].items():
        if Path(name).name!=name or digest(Path(__file__).parent/name)!=sha:raise ValueError("Implementation drift: "+name)
    candidate=bound_json(root,rules["inputs"]["candidate"])
    bindings=bound_json(root,rules["inputs"]["bindings"])
    sources=bound_json(root,rules["inputs"]["filter_sources"])
    manifest=bound_json(root,rules["inputs"]["frozen_candidate_manifest"])
    selection=bound_json(root,rules["inputs"]["frozen_selected_cases"])
    if rules["inputs"]["frozen_selected_cases"]["sha256"]!=next(x["sha256"] for x in manifest["files"] if x["file"]=="SELECTED_CASES.json"):
        raise ValueError("Selection not bound by frozen candidate manifest")
    assays=validate_metadata(rules,candidate,bindings,selection)
    out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=False)
    run_record={"schema":"single-candidate-readback-run.public-adapter.v1","status":"RUNNING","started_utc":utc(),"rules_sha256":args.rules_sha256,"rules":rules,"command":sys.argv,"inputs_verified":[],"source_profile_files_modified":False,"FA_read":False}
    write_json(out/"RUN.json",run_record)
    def log(message):
        with (out/"execution.log").open("a",encoding="utf-8") as f:f.write(utc()+" "+message+"\n")
        print(message,flush=True)
    try:
        log("Validating frozen source metadata/runtime; no event counts yet")
        for rec in sources["remote_source_files"]+sources["current_pysam_runtime_files"]+bindings["split_definition"]["current_installed_definition_sources"]:verify_hashed_record(rec)
        for rec in bindings["current_stat_check"]["small_source_files_current_sha256"]:verify_hashed_record(rec)
        active_keys=("profile_bam","profile_bai","accepted_Rdic","covT","attributes","scaffold_lengths")
        verified={k:stat_bound(bindings["file_bindings"][k]) for k in active_keys}
        for k in ("profile_bai","covT","attributes","scaffold_lengths"):
            verified[k]["current_full_sha256"]=digest(verified[k]["path"])
            if verified[k]["current_full_sha256"]!=verified[k]["historic_sha256"]:raise ValueError("Input SHA mismatch: "+k)
        import pysam
        import pysam.version
        import h5py
        import inStrain
        versions={"inStrain":inStrain.__version__,"pysam":pysam.__version__,"htslib_reported_by_pysam":pysam.version.__htslib_version__,"samtools_reported_by_pysam":pysam.__samtools_version__}
        if versions!=sources["versions"]:raise ValueError("Runtime version drift")
        run_record["runtime"]={**versions,"python":platform.python_version(),"h5py":h5py.__version__}
        log("Loading original skip-mm Rdic once, verifying SHA on the same stream")
        accepted=load_rdic(verified["accepted_Rdic"]["path"],verified["accepted_Rdic"]["historic_sha256"],[a.reference for a in assays])
        verified["accepted_Rdic"]["current_full_sha256"]=verified["accepted_Rdic"]["historic_sha256"]
        run_record["inputs_verified"]=verified
        write_json(out/"RUN.json",run_record)
        log("Replaying original pileup only in the two original event-containing splits")
        with pysam.AlignmentFile(verified["profile_bam"]["path"],"rb",index_filename=verified["profile_bai"]["path"]) as bam, h5py.File(verified["covT"]["path"],"r") as cov:
            for w in rules["parity_windows"]:
                if bam.get_reference_length(w["reference"])!=w["reference_length"]:raise ValueError("BAM reference length mismatch")
            parity=parity_report(bam,cov,accepted,rules["parity_windows"])
            write_json(out/"PARITY.json",parity)
            calls=[]
            if parity["all_windows_match"]:
                log("Both frozen windows match covT; fetching event-overlapping original reads")
                calls=perform_event_fetch(bam,assays,accepted)
            else:log("Parity did not match; event classification is blocked")
        summary,main,diagnostic=interpret_calls(calls,assays,rules["technical_support"],parity)
        summary.update(candidate_id=rules["candidate_id"],sample_id=rules["sample_id"],parity_all_windows_match=parity["all_windows_match"])
        write_json(out/"SUMMARY.json",summary)
        for filename,rows in (("READ_CALLS.jsonl",calls),("MAIN_QNAMES.jsonl",main),("DIAGNOSTIC_QNAMES.jsonl",diagnostic)):
            with (out/filename).open("x",encoding="utf-8") as f:
                for row in rows:f.write(json.dumps(row,ensure_ascii=False,sort_keys=True)+"\n")
        for k in active_keys:stat_bound(bindings["file_bindings"][k])
        run_record.update(status=summary["status"],finished_utc=utc())
        log("Finished bounded readback: "+summary["status"])
        write_json(out/"RUN.json",run_record)
        outputs=[{"file":p.name,"bytes":p.stat().st_size,"sha256":digest(p)} for p in sorted(out.iterdir()) if p.is_file()]
        receipt={"schema":"single-candidate-readback-completion.public-adapter.v1","status":summary["status"],"rules_sha256":args.rules_sha256,"outputs":outputs}
        write_json(out/"MANIFEST.json",receipt)
        (out/"COMPLETE").write_text(digest(out/"MANIFEST.json")+"\n",encoding="ascii")
        return 0 if parity["all_windows_match"] else 2
    except BaseException as e:
        run_record.update(status="FAILED",finished_utc=utc(),error=repr(e),traceback=traceback.format_exc())
        write_json(out/"RUN.json",run_record)
        log("FAILED; partial output retained without COMPLETE: "+repr(e))
        raise


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rules",required=True)
    parser.add_argument("--rules-sha256",required=True)
    parser.add_argument("--out",required=True)
    parser.add_argument("--root-authorized-run",action="store_true")
    sys.exit(run(parser.parse_args()))
