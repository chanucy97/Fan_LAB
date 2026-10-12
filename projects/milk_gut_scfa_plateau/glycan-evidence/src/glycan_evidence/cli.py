"""Command-line interfaces; expensive dependencies are imported only when used."""
import argparse
import csv
import json
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

from . import __version__
from .integrity import write_json, seal_outputs, verify_outputs


def sample_list(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        reader=csv.DictReader(stream, delimiter='\t')
        if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames) or 'sample' not in reader.fieldnames:
            raise ValueError('Sample TSV requires a unique sample column')
        rows=list(reader)
    if any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError('Malformed sample TSV row')
    samples=[r['sample'] for r in rows]
    if not samples or any(not s or s != s.strip() or any(c in s for c in '/\\:\t\r\n') or s in {'.','..'} for s in samples):
        raise ValueError('Sample IDs must be nonempty safe directory names')
    if len(set(samples))!=len(samples):
        raise ValueError('Duplicate sample IDs')
    return samples


def doctor(bin_dir=None, structural=False, require=False):
    requested=['bowtie2','fastp'] + (['blastn','bowtie2-build','samtools'] if structural else [])
    tools={}
    for name in requested:
        p=Path(bin_dir)/name if bin_dir else None
        if p and not p.is_file() and p.with_suffix('.exe').is_file():
            p=p.with_suffix('.exe')
        executable=str(p) if p and p.is_file() else (shutil.which(name) if not bin_dir else None)
        entry=dict(available=bool(executable),path=executable)
        if executable:
            flag='-version' if name in {'blastn','makeblastdb'} else '--version'
            try:
                run=subprocess.run([executable,flag],capture_output=True,text=True,timeout=20)
                entry.update(returncode=run.returncode,version=(run.stdout+'\n'+run.stderr).strip()[:1200])
                entry['available']=run.returncode==0
            except (OSError,subprocess.TimeoutExpired) as error:
                entry.update(available=False,error=str(error))
        tools[name]=entry
    complete=all(x['available'] for x in tools.values())
    return dict(status='PASS' if complete else 'MISSING_DEPENDENCIES',package_version=__version__,
                python=sys.version,tools=tools,scope='Dependency inspection, not benchmark or scientific validation'), (0 if complete or not require else 1)


def parser():
    p=argparse.ArgumentParser(description='Separate identity-marker and glycan-locus sequence evidence (development release).')
    p.add_argument('--version',action='version',version=__version__)
    s=p.add_subparsers(dest='command',required=True)
    d=s.add_parser('doctor',help='Inspect external tools without changing the environment')
    d.add_argument('--bin-dir'); d.add_argument('--structural',action='store_true'); d.add_argument('--require',action='store_true')
    r=s.add_parser('run',help='QC, screen, compete and apply the frozen v2 identity rule')
    for k in ['r1','r2','sample','panel','rules','out']:
        r.add_argument('--'+k,required=True)
    r.add_argument('--threads',type=int,default=1); r.add_argument('--bin-dir'); r.add_argument('--skip-qc',action='store_true',help='Only for explicitly pre-QC or simulated reads')
    c=s.add_parser('call',help='Apply frozen rules to complete, provenance-verified metrics')
    for k in ['metrics-root','samples','rules','out']:
        c.add_argument('--'+k,required=True)
    t=s.add_parser('structure',help='Inspect configured anchor geometry and optional readback')
    for k in ['contigs','anchors','out']:
        t.add_argument('--'+k,required=True)
    for k in ['r1','r2','bam','config','bin-dir']:
        t.add_argument('--'+k)
    t.add_argument('--threads',type=int,default=1)
    a=s.add_parser('assemble',help='QC and assemble reads recruited against an explicit reference index')
    for k in ['r1','r2','sample','recruit-index','reference-fasta','recruitment-manifest','out','megahit']:
        a.add_argument('--'+k,required=True)
    a.add_argument('--threads',type=int,default=12); a.add_argument('--bin-dir')
    a.add_argument('--expected-qc-pairs',type=int)
    a.add_argument('--memory-bytes',type=int,default=24000000000)
    a.add_argument('--min-contig-len',type=int,default=500)
    a.add_argument('--validation-tmpdir',help='Existing scratch directory for the exact QNAME uniqueness check')
    b=s.add_parser('benchmark',help='Evaluate an explicit known-truth ledger; retain unknown labels')
    for k in ['calls','truth','out']:
        b.add_argument('--'+k,required=True)
    b.add_argument('--margin',type=int,choices=[0,6,12,18],default=0); b.add_argument('--split',default='test')
    f=s.add_parser('profile',help='Describe full-marker shape with explicit source-separated expectations; no new identity call')
    for k in ['metrics','model','source','out']:
        f.add_argument('--'+k,required=True)
    h=s.add_parser('compatibility',help='Verify completed same-input metrics and run a fixed reference-family research bundle')
    for k in ['metrics','bundle','expected-manifest-sha256','out']:
        h.add_argument('--'+k,required=True)
    l=s.add_parser('local-support',help='Run the frozen local-support candidate from complete same-input metrics')
    for k in ['metrics','bundle','expected-manifest-sha256','out']:
        l.add_argument('--'+k,required=True)
    dv=s.add_parser('discrepancy-verify',help='Replay sealed complete-QC discrepancy certificates without optimization or new calls')
    for k in ['root','sample','expected-freeze-sha256','expected-results-sha256','out']:
        dv.add_argument('--'+k,required=True)
    q=s.add_parser('bases',help='Audit every candidate base against the bound readback BAM, collapsing consistent mates by QNAME')
    for k in ['readback-receipt','structures','out']:
        q.add_argument('--'+k,required=True)
    cr=s.add_parser('cds-replay',help='Replay coding-sequence checks from bound candidate sequences and saved HSPs')
    for k in ['input-manifest','expected-manifest-sha256','out']:
        cr.add_argument('--'+k,required=True)
    v=s.add_parser('verify',help='Verify the small-output integrity manifest')
    v.add_argument('--directory',required=True)
    return p


def execute(args):
    if args.command=='doctor':
        return doctor(args.bin_dir,args.structural,args.require)
    if args.command=='verify':
        return verify_outputs(args.directory),0
    if args.command=='benchmark':
        from .benchmark import evaluate
        receipt=evaluate(args.calls,args.truth,args.out,args.margin,args.split)
    elif args.command=='call':
        from .identity import apply_rules
        receipt=apply_rules(args.metrics_root,sample_list(args.samples),args.rules,args.out)
    elif args.command=='structure':
        from .structure import run_structure
        receipt=run_structure(args.contigs,args.anchors,args.out,r1=args.r1,r2=args.r2,bam=args.bam,
                              config=args.config,bin_dir=args.bin_dir,threads=args.threads)
    elif args.command=='assemble':
        from .assembly import run_assembly
        receipt=run_assembly(args.r1,args.r2,args.sample,args.recruit_index,args.reference_fasta,
                             args.recruitment_manifest,args.out,args.megahit,threads=args.threads,
                             bin_dir=args.bin_dir,expected_qc_pairs=args.expected_qc_pairs,
                             memory_bytes=args.memory_bytes,min_contig_len=args.min_contig_len,
                             validation_tmpdir=args.validation_tmpdir)
    elif args.command=='profile':
        from .profile_io import run_profile
        receipt=run_profile(args.metrics,args.model,args.source,args.out)
    elif args.command=='compatibility':
        from .compatibility_io import run_compatibility
        receipt=run_compatibility(args.metrics,args.bundle,args.expected_manifest_sha256,args.out)
    elif args.command=='local-support':
        from .local_support_io import run_local_support
        receipt=run_local_support(args.metrics,args.bundle,args.expected_manifest_sha256,args.out)
    elif args.command=='discrepancy-verify':
        from .discrepancy_verify import run_discrepancy_verify
        receipt=run_discrepancy_verify(args.root,args.sample,args.expected_freeze_sha256,
                                       args.expected_results_sha256,args.out)
    elif args.command=='bases':
        from .base_concordance import run_base_concordance
        receipt=run_base_concordance(args.readback_receipt,args.structures,args.out)
    elif args.command=='cds-replay':
        from .cds import run_cds_replay
        receipt=run_cds_replay(args.input_manifest,args.out,
                               expected_manifest_sha256=args.expected_manifest_sha256)
    elif args.command=='run':
        from .identity import run_sample,apply_rules,load_rules
        # Validate a single sample with the same constraints as manifest samples.
        if not args.sample or args.sample != args.sample.strip() or args.sample in {'.','..'} or any(c in args.sample for c in '/\\:\t\r\n'):
            raise ValueError('Unsafe sample ID')
        if args.threads<1:
            raise ValueError('Threads must be positive')
        load_rules(args.rules)
        out=Path(args.out).resolve()
        out.mkdir(parents=True,exist_ok=False)
        receipt=dict(schema_version=1,package_version=__version__,status='RUNNING',sample=args.sample,
                     started_utc=datetime.now(timezone.utc).isoformat(),scope='Frozen v2 identity and continuous reference-locus evidence; not absolute abundance or function')
        write_json(out/'receipt.json',receipt)
        try:
            run_sample(args.r1,args.r2,args.sample,args.panel,out/'metrics'/args.sample,
                       threads=args.threads,skip_qc=args.skip_qc,bin_dir=args.bin_dir)
            apply_rules(out/'metrics',[args.sample],args.rules,out/'calls')
            receipt.update(status='COMPLETE',completed_utc=datetime.now(timezone.utc).isoformat())
            write_json(out/'receipt.json',receipt)
            (out/'COMPLETE').write_text(receipt['completed_utc']+'\n',encoding='utf-8')
        except BaseException as error:
            receipt.update(status='FAILED',error_type=type(error).__name__,error=str(error))
            write_json(out/'receipt.json',receipt)
            raise
    else:
        raise ValueError('Unknown command')
    try:
        integrity=seal_outputs(args.out)
    except BaseException as error:
        # A completed scientific substage is not a completed top-level command.
        # Keep its reports for diagnosis, but revoke the marker consumed by batch gates.
        out=Path(args.out).resolve()
        try:
            (out/'COMPLETE').unlink(missing_ok=True)
        except OSError as cleanup_error:
            raise RuntimeError(f'Output sealing failed: {error}; could not revoke COMPLETE: {cleanup_error}') from error
        write_json(out/'CLI_FAILED.json',dict(status='FAILED',stage='seal_outputs',
                                             error_type=type(error).__name__,error=str(error),
                                             completed_scientific_substages_preserved=True))
        raise
    return dict(status='COMPLETE',command=args.command,out=str(Path(args.out).resolve()),
                package_version=__version__,integrity=integrity),0


def main(argv=None):
    args=parser().parse_args(argv)
    try:
        result,status=execute(args)
        print(json.dumps(result,ensure_ascii=False,allow_nan=False))
        return status
    except (OSError,ValueError,RuntimeError,KeyError,ImportError,csv.Error,subprocess.SubprocessError) as error:
        print(json.dumps(dict(status='FAILED',error_type=type(error).__name__,error=str(error)),ensure_ascii=False),file=sys.stderr)
        return 2
