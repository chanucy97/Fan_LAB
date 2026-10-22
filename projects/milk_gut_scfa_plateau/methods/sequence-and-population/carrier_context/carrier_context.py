"""Frozen descriptive H5/background contrast; no phenotype or carrier classification."""
import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

FIELDS = 'qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore gaps gapopen btop'.split()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()


def rec(path):
    p=Path(path)
    return {'path':str(p.resolve()),'bytes':p.stat().st_size,'sha256':sha(p)}


def verify(r):
    p=Path(r['path'])
    if p.stat().st_size!=r['bytes'] or sha(p)!=r['sha256']:raise ValueError('source changed: '+str(p))


def fasta(path):
    name=None;parts=[]
    with Path(path).open(encoding='ascii') as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(parts).upper()
                name=line[1:].split()[0];parts=[]
            elif line.strip():
                if name is None:raise ValueError('sequence before header')
                parts.append(line.strip())
    if name is not None:yield name,''.join(parts).upper()


def table(path, rows, fields=None):
    fields=fields or list(rows[0])
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)


def regions(case, config):
    n=case['contig_length'];a=case['interval_start0'];b=case['interval_end0']
    if not 0<=a<b<=n:raise ValueError('invalid local interval')
    gap=config['outer_interval_buffer'];width=config['flank_length']
    left_stop=max(0,a-gap);right_start=min(n,b+gap)
    def reg(name,parts,requested):
        parts=[(x,y) for x,y in parts if y>x]
        length=sum(y-x for x,y in parts)
        return {'region':name,'segments':parts,'length':length,'requested_length':requested,
                'full_requested_length_available':length==requested}
    return [reg('whole_contig',[(0,n)],n),reg('local_interval',[(a,b)],b-a),
            reg('masked_background',[(0,left_stop),(right_start,n)],left_stop+n-right_start),
            reg('left_flank',[(max(0,left_stop-width),left_stop)],width),
            reg('right_flank',[(right_start,min(n,right_start+width))],width)]


def positions(h):
    """Yield (query0, subject0 or None, exact), including query deletion bases."""
    qs,qe,ss,se=(int(h[k]) for k in ('qstart','qend','sstart','send'))
    qstep=1 if qe>=qs else -1;sstep=1 if se>=ss else -1
    q=qs-1;s=ss-1
    tokens=re.findall(r'\d+|[^\d]{2}',h['btop'])
    if ''.join(tokens)!=h['btop']:raise ValueError('malformed BTOP')
    for token in tokens:
        if token.isdigit():
            for _ in range(int(token)):
                yield q,s,True;q+=qstep;s+=sstep
        else:
            qb,sb=token
            if qb=='-' and sb=='-':raise ValueError('double gap')
            if qb!='-':
                yield q,None if sb=='-' else s,False;q+=qstep
            if sb!='-':s+=sstep
    if q!=qe-1+qstep or s!=se-1+sstep:raise ValueError('BTOP coordinate span disagreement')


def masks(hits,n):
    """Highest HSP bitscore per query position; exact status ANDed across score ties."""
    from array import array
    covered=bytearray(n);matched=bytearray(n);tie=bytearray(n)
    score=array('d',[-math.inf])*n
    for h in sorted(hits,key=lambda x:(-float(x['bitscore']),x['sseqid'],int(x['qstart']),int(x['sstart']),x['btop'])):
        bits=float(h['bitscore'])
        if not math.isfinite(bits):raise ValueError('nonfinite bitscore')
        for q,s,match in positions(h):
            if not 0<=q<n:raise ValueError('query coordinate outside contig')
            if not covered[q]:covered[q]=1;matched[q]=match;score[q]=bits
            elif score[q]==bits and matched[q]!=int(match):tie[q]=1;matched[q]=0
    return covered,matched,tie


def clipped_hsp(h,start,end,subject_lengths):
    lo=min(int(h['qstart']),int(h['qend']))-1;hi=max(int(h['qstart']),int(h['qend']))
    if hi<=start or lo>=end:return None
    selected=sorted((q,s,match) for q,s,match in positions(h) if start<=q<end)
    # Terminal query bases opposite a subject gap have no reference coordinate;
    # trim only those bases, retaining internal gaps as covered nonmatches.
    while selected and selected[0][1] is None:selected.pop(0)
    while selected and selected[-1][1] is None:selected.pop()
    if not selected:return None
    qmin,slow,_=selected[0];qmax,shigh,_=selected[-1]
    covered=len(selected);matched=sum(int(p[2]) for p in selected)
    strand=1 if (int(h['qend'])-int(h['qstart']))*(int(h['send'])-int(h['sstart']))>=0 else -1
    sid=h['sseqid'];length=subject_lengths[sid]
    if not 0<=min(slow,shigh)<=max(slow,shigh)<length:raise ValueError('subject bounds')
    return {'subject':sid,'strand':strand,'qstart0':qmin,'qend0':qmax+1,
            'subject_at_query_start0':slow,'subject_at_query_end0':shigh,
            'tstart':strand*slow,'tend':strand*shigh+1,
            'covered':covered,'matched':matched,'bitscore':float(h['bitscore'])}


def chain(hits,start,end,subject_lengths,config):
    blocks=[b for h in hits if (b:=clipped_hsp(h,start,end,subject_lengths)) is not None]
    blocks.sort(key=lambda b:(b['subject'],b['strand'],b['qend0'],b['qstart0'],b['tstart'],-b['matched']))
    dp=[];paths=[];possible_origin=False
    for i,b in enumerate(blocks):
        best=(b['matched'],b['covered']);path=[i]
        for j,a in enumerate(blocks[:i]):
            if a['subject']!=b['subject'] or a['strand']!=b['strand']:continue
            qgap=b['qstart0']-a['qend0'];sgap=b['tstart']-a['tend']
            if qgap<0 or qgap>config['chain_max_query_gap']:continue
            if sgap<0:
                wrap=sgap+subject_lengths[b['subject']]
                if 0<=wrap<=config['chain_max_subject_gap'] and abs(wrap-qgap)<=config['chain_max_gap_discrepancy']:
                    possible_origin=True
                continue
            if sgap>config['chain_max_subject_gap'] or abs(sgap-qgap)>config['chain_max_gap_discrepancy']:continue
            score=(dp[j][0]+b['matched'],dp[j][1]+b['covered'])
            if score>best:best=score;path=paths[j]+[i]
        dp.append(best);paths.append(path)
    if not blocks:return {'covered':0,'matched':0,'blocks':[],'possible_circular_origin':False,'near_equal_other_locus':False}
    winner=max(range(len(dp)),key=lambda i:dp[i]);chosen=[blocks[i] for i in paths[winner]]
    first,last=chosen[0],chosen[-1]
    alt=False;threshold=dp[winner][0]-config['unique_margin']*(end-start)
    for i,score in enumerate(dp):
        if i==winner or score[0]<threshold:continue
        other=[blocks[j] for j in paths[i]];a,z=other[0],other[-1]
        same=(a['subject']==first['subject'] and a['strand']==first['strand'])
        disjoint=not same or z['tend']<=first['tstart'] or a['tstart']>=last['tend']
        if disjoint:alt=True
    return {'covered':dp[winner][1],'matched':dp[winner][0],'subject':first['subject'],'strand':first['strand'],
            'qstart0':first['qstart0'],'qend0':last['qend0'],
            'subject_start0':min(first['subject_at_query_start0'],last['subject_at_query_end0']),
            'subject_end0':max(first['subject_at_query_start0'],last['subject_at_query_end0'])+1,
            'tstart':first['tstart'],'tend':last['tend'],
            'blocks':chosen,'possible_circular_origin':possible_origin,'near_equal_other_locus':alt}


def fractions(covered,matched,length):
    return {'covered_bases':covered,'exact_matching_bases':matched,'coverage':covered/length if length else None,
            'covered_identity':matched/covered if covered else None,'whole_region_matching_fraction':matched/length if length else None}


def score_case(case,hits,manifest,config):
    refs=[r['accession'] for r in manifest['references']]
    subject_lengths={r['subject_id']:int(r['length']) for r in manifest['subjects']}
    subject_accession={r['subject_id']:r['accession'] for r in manifest['subjects']}
    grouped=defaultdict(list)
    for h in hits:
        if h['qseqid']!=case['contig'] or int(h['qlen'])!=case['contig_length']:raise ValueError('unexpected query')
        if not 1<=min(int(h['sstart']),int(h['send']))<=max(int(h['sstart']),int(h['send']))<=subject_lengths[h['sseqid']]:raise ValueError('subject coordinate outside reference')
        grouped[subject_accession[h['sseqid']]].append(h)
    rr=regions(case,config);scores=[];chains={}
    for acc in refs:
        cov,mat,ties=masks(grouped[acc],case['contig_length'])
        for region in rr:
            c=sum(sum(cov[a:b]) for a,b in region['segments']);m=sum(sum(mat[a:b]) for a,b in region['segments']);t=sum(sum(ties[a:b]) for a,b in region['segments'])
            base={'case_id':case['case_id'],'sample_id':case['sample_id'],'scope':case['scope'],'reference':acc,
                  'region':region['region'],'method':'query_position_union','region_length':region['length'],
                  'requested_length':region['requested_length'],'full_requested_length_available':region['full_requested_length_available']}
            scores.append({**base,**fractions(c,m,region['length']),'equal_bitscore_conflicting_positions':t})
            if region['region'] in ('left_flank','right_flank'):
                ch=chain(grouped[acc],*region['segments'][0],subject_lengths,config) if region['segments'] else {'covered':0,'matched':0,'blocks':[],'possible_circular_origin':False,'near_equal_other_locus':False}
                chains[region['region']+'|'+acc]=ch
                scores.append({**base,'method':'collinear_chain',**fractions(ch['covered'],ch['matched'],region['length']),'equal_bitscore_conflicting_positions':None})
    ranked={}
    for name,method in [('whole_contig','query_position_union'),('local_interval','query_position_union'),('masked_background','query_position_union'),('left_flank','collinear_chain'),('right_flank','collinear_chain')]:
        relevant=sorted([r for r in scores if r['region']==name and r['method']==method],key=lambda r:(-(r['whole_region_matching_fraction'] or 0),r['reference']))
        top,second=relevant[:2];value=top['whole_region_matching_fraction'] or 0;margin=value-(second['whole_region_matching_fraction'] or 0)
        available=top['region_length']>0 and top['full_requested_length_available']
        qualified=available and top['coverage']>=config['minimum_region_coverage'] and (top['covered_identity'] or 0)>=config['minimum_covered_identity']
        unique=qualified and margin+1e-12>=config['unique_margin']
        ranked[name]={'best_reference':top['reference'],'runner_up_reference':second['reference'],'best_score':top['whole_region_matching_fraction'],
                      'runner_up_score':second['whole_region_matching_fraction'],'margin':margin,'coverage_qualified':qualified,'unique_reference_affinity':unique,
                      'region_status':'insufficient_flank_length' if not available else 'unique_affinity' if unique else 'ambiguous_affinity' if qualified else 'insufficient_alignment_support',
                      'near_best_references':[r['reference'] for r in relevant if value-(r['whole_region_matching_fraction'] or 0)<config['unique_margin']-1e-12],
                      'full_requested_length_available':available}
    reasons=[];required=['masked_background','left_flank','right_flank']
    if case['scope']!='primary_complete_H5':reasons.append('secondary_case_not_primary_complete_H5')
    for name in required:
        if not ranked[name]['unique_reference_affinity']:reasons.append(name+':'+ranked[name]['region_status'])
    winners={ranked[name]['best_reference'] for name in required}
    if len(winners)!=1:reasons.append('different_best_references_between_regions')
    left=chains['left_flank|'+ranked['left_flank']['best_reference']];right=chains['right_flank|'+ranked['right_flank']['best_reference']]
    for name,ch in [('left_flank',left),('right_flank',right)]:
        if ch['near_equal_other_locus']:reasons.append(name+':multiple_reference_locations')
        if ch['possible_circular_origin']:reasons.append('REVIEW_CIRCULAR_ORIGIN:'+name)
    geometry={'subject_order_consistent':False,'possible_circular_origin':False}
    if left['blocks'] and right['blocks']:
        same=left['subject']==right['subject'] and left['strand']==right['strand']
        geometry.update(same_subject_and_direction=same,left_subject=left['subject'],right_subject=right['subject'],left_strand=left['strand'],right_strand=right['strand'],left_subject_start0=left['subject_start0'],left_subject_end0=left['subject_end0'],right_subject_start0=right['subject_start0'],right_subject_end0=right['subject_end0'])
        if not same:reasons.append('different_subject_or_direction_between_flanks')
        else:
            qgap=right['qstart0']-left['qend0'];sgap=right['tstart']-left['tend'];geometry.update(query_gap=qgap,subject_oriented_gap=sgap)
            if sgap<0:
                wrapped=sgap+subject_lengths[left['subject']]
                possible=wrapped>=0 and abs(wrapped-qgap)<=config['cross_flank_max_gap_discrepancy']
                geometry.update(possible_circular_origin=possible,circular_adjusted_gap=wrapped)
                reasons.append('REVIEW_CIRCULAR_ORIGIN:between_flanks' if possible else 'noncollinear_subject_order')
            elif abs(sgap-qgap)>config['cross_flank_max_gap_discrepancy']:reasons.append('distant_reference_locations_between_flanks')
            else:geometry['subject_order_consistent']=True
    else:reasons.append('insufficient_collinear_chain_geometry')
    decision={'case_id':case['case_id'],'sample_id':case['sample_id'],'scope':case['scope'],'region_rankings':ranked,'cross_flank_geometry':geometry,
              'status':'COHERENT_PANEL_REFERENCE_AFFINITY' if not reasons else 'REVIEW_OR_UNRESOLVED','reasons':reasons,
              'affinity_reference':next(iter(winners)) if not reasons else None,'read_connection_verified':False,
              'carrier_taxon_or_strain_asserted':False,'HGT_asserted':False,'same_nutrient_enzyme_set_asserted':False}
    return scores,chains,decision,rr


def run_cmd(argv,dest,commands):
    r=subprocess.run([str(x) for x in argv],capture_output=True,text=True)
    (dest/('command_%02d.stdout.txt'%len(commands))).write_text(r.stdout,encoding='utf-8')
    (dest/('command_%02d.stderr.txt'%len(commands))).write_text(r.stderr,encoding='utf-8')
    commands.append({'argv':[str(x) for x in argv],'returncode':r.returncode})
    if r.returncode:raise RuntimeError('command failed: '+str(argv[0]))


def prepare(inputs,config_path,out):
    m=read(inputs);config=read(config_path);verify(m['reference_fasta'])
    if out.exists():raise FileExistsError('fresh prepared directory required')
    out.mkdir(parents=True);commands=[]
    expected={r['subject_id']:int(r['length']) for r in m['subjects']};observed={}
    for name,seq in fasta(m['reference_fasta']['path']):
        if name in observed:raise ValueError('duplicate reference subject')
        observed[name]=len(seq)
    if observed!=expected:raise ValueError('reference FASTA does not match frozen subject roster')
    run_cmd([config['makeblastdb'],'-version'],out,commands)
    run_cmd([config['blastn'],'-version'],out,commands)
    run_cmd([config['makeblastdb'],'-in',m['reference_fasta']['path'],'-dbtype','nucl','-parse_seqids','-out',out/'reference21'],out,commands)
    dbfiles=[rec(p) for p in sorted(out.glob('reference21.*'))]
    if not dbfiles:raise ValueError('empty database')
    state={'status':'COMPLETE','inputs':rec(inputs),'config':rec(config_path),'code':rec(__file__),'reference_fasta':m['reference_fasta'],
           'database_prefix':str((out/'reference21').resolve()),'database_files':dbfiles,'commands':commands,'subjects_verified':len(expected),
           'python_version':sys.version,'command_logs':[rec(f) for f in sorted(out.glob('command_*.txt'))]}
    write(out/'PREPARED.json',state);(out/'COMPLETE').write_text(sha(out/'PREPARED.json')+'\n')
    return state


def execute(inputs,config_path,prepared,outroot,index):
    m=read(inputs);config=read(config_path);p=read(prepared)
    if index is None or not 0<=index<len(m['cases']):raise ValueError('invalid case index')
    case=m['cases'][index]
    if p['inputs']!=rec(inputs) or p['config']!=rec(config_path) or p['code']!=rec(__file__):raise ValueError('preparation freeze mismatch')
    if (prepared.parent/'COMPLETE').read_text().strip()!=sha(prepared):raise ValueError('preparation incomplete')
    for r in p['database_files']:verify(r)
    if not re.fullmatch('[A-Za-z0-9_.-]+',case['case_id']):raise ValueError('unsafe case id')
    out=outroot/case['case_id']
    if out.exists():
        if (out/'COMPLETE').exists():
            previous=read(out/'RUN.json')
            if (out/'COMPLETE').read_text().strip()!=sha(out/'RUN.json') or previous['inputs']!=rec(inputs) or previous['code']!=rec(__file__) or previous['config']!=rec(config_path) or previous['prepared']!=rec(prepared) or previous['case_id']!=case['case_id']:raise ValueError('completed run mismatch')
            for r in previous['outputs']:verify(r)
            return {'status':'ALREADY_COMPLETE_VERIFIED','case_id':case['case_id']}
        raise FileExistsError('failed or partial directory retained; use a new output root')
    out.mkdir(parents=True);commands=[]
    state={'status':'RUNNING','case_id':case['case_id'],'inputs':rec(inputs),'config':rec(config_path),'code':rec(__file__),'prepared':rec(prepared),'source':case,'commands':commands}
    try:
        verify(case['assembly']);selected=[s for name,s in fasta(case['assembly']['path']) if name==case['contig']]
        if len(selected)!=1:raise ValueError('contig absent or duplicated')
        seq=selected[0]
        if len(seq)!=case['contig_length'] or hashlib.sha256(seq.encode()).hexdigest()!=case['contig_sha256']:raise ValueError('contig sequence changed')
        local=seq[case['interval_start0']:case['interval_end0']]
        if case['local_sequence_hash_orientation']=='-':local=local.translate(str.maketrans('ACGT','TGCA'))[::-1]
        if hashlib.sha256(local.encode()).hexdigest()!=case['local_sequence_sha256']:raise ValueError('interval sequence/coordinates changed')
        query=out/'selected_contig.fna';query.write_text('>'+case['contig']+'\n'+seq+'\n',encoding='ascii')
        hp=out/'all_HSPs_btop.tsv'
        run_cmd([config['blastn'],'-task',config['blast_task'],'-word_size',config['blast_word_size'],'-query',query,'-db',p['database_prefix'],'-evalue',config['blast_evalue'],'-num_threads',config['threads'],'-max_target_seqs',config['blast_max_target_seqs'],'-outfmt','6 '+' '.join(FIELDS),'-out',hp],out,commands)
        with hp.open(newline='') as f:hits=list(csv.DictReader(f,delimiter='\t',fieldnames=FIELDS))
        scores,chains,decision,rr=score_case(case,hits,m,config)
        table(out/'REGION_REFERENCE_SCORES.tsv',scores);write(out/'FLANK_CHAINS.json',chains);write(out/'DECISION.json',decision);write(out/'REGIONS.json',rr)
        state.update(status='COMPLETE',completed_utc=datetime.now(timezone.utc).isoformat(),HSPs=len(hits),score_rows=len(scores),phenotype_inputs=[],BAM_read=False,
                     outputs=[rec(f) for f in sorted(out.iterdir()) if f.is_file()])
        write(out/'RUN.json',state);(out/'COMPLETE').write_text(sha(out/'RUN.json')+'\n')
        return {'status':'COMPLETE','case_id':case['case_id'],'decision':decision['status']}
    except BaseException as e:
        state.update(status='FAILED',error=repr(e));write(out/'FAILED.json',state);raise


def same_content(receipt,path):
    current=rec(path)
    return receipt['bytes']==current['bytes'] and receipt['sha256']==current['sha256']


def collect(inputs,config_path,root,out):
    """Aggregate frozen results; allow downloaded copies by unambiguous basename/hash."""
    m=read(inputs)
    if out.exists():raise FileExistsError('fresh collection directory required')
    decisions=[];all_scores=[];all_regions=[];bindings=[]
    for case in m['cases']:
        folder=root/case['case_id'];runpath=folder/'RUN.json';run=read(runpath)
        if (folder/'COMPLETE').read_text().strip()!=sha(runpath) or run['status']!='COMPLETE' or run['case_id']!=case['case_id']:raise ValueError('incomplete case '+case['case_id'])
        if run['source']!=case or not same_content(run['inputs'],inputs) or not same_content(run['config'],config_path) or not same_content(run['code'],__file__):raise ValueError('case source/freeze mismatch')
        names=set()
        for r in run['outputs']:
            name=Path(r['path']).name
            if name in names:raise ValueError('ambiguous output basename')
            names.add(name)
            if not same_content(r,folder/name):raise ValueError('changed downloaded result '+name)
        if not {'REGION_REFERENCE_SCORES.tsv','DECISION.json','REGIONS.json','FLANK_CHAINS.json','all_HSPs_btop.tsv','selected_contig.fna'}<=names:raise ValueError('missing case products')
        d=read(folder/'DECISION.json')
        if d['case_id']!=case['case_id'] or d['scope']!=case['scope']:raise ValueError('decision source mismatch')
        decisions.append(d);bindings.append(rec(runpath))
        with (folder/'REGION_REFERENCE_SCORES.tsv').open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
        if len(rows)!=len(m['references'])*7 or any(r['case_id']!=case['case_id'] for r in rows):raise ValueError('score denominator mismatch')
        all_scores.extend(rows)
        for region in read(folder/'REGIONS.json'):
            all_regions.append({'case_id':case['case_id'],'scope':case['scope'],**region})
    summary=[]
    for d in decisions:
        row={k:d[k] for k in ['case_id','sample_id','scope','status','affinity_reference']}
        row['reasons']=';'.join(d['reasons'])
        for name,v in d['region_rankings'].items():
            for k in ['best_reference','runner_up_reference','best_score','runner_up_score','margin','region_status']:
                row[name+'__'+k]=v[k]
        summary.append(row)
    out.mkdir(parents=True)
    table(out/'ALL_CASES.tsv',summary);table(out/'ALL_REFERENCE_SCORES.tsv',all_scores)
    write(out/'ALL_DECISIONS.json',decisions);write(out/'ALL_REGIONS.json',all_regions)
    report={'status':'COMPLETE','schema':'H5-context-collection.v1','cases':len(decisions),'primary_complete_H5':sum(d['scope']=='primary_complete_H5' for d in decisions),
            'primary_coherent_affinity':sum(d['scope']=='primary_complete_H5' and d['status']=='COHERENT_PANEL_REFERENCE_AFFINITY' for d in decisions),
            'reference_score_rows':len(all_scores),'cohort_denominator':m['cohort_denominator'],'inputs':rec(inputs),'config':rec(config_path),'code':rec(__file__),
            'case_runs':bindings,'read_connection_verified':False,'interpretation':'Descriptive reference-panel affinity only; incomplete flanks, mapping ambiguity, origin crossings and noncollinearity remain review states.',
            'outputs':[rec(p) for p in sorted(out.iterdir()) if p.is_file()]}
    write(out/'COLLECTION.json',report);(out/'COMPLETE').write_text(sha(out/'COLLECTION.json')+'\n')
    return report


def main():
    a=argparse.ArgumentParser();a.add_argument('mode',choices=['prepare','run','collect']);a.add_argument('--inputs',type=Path,required=True);a.add_argument('--config',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--prepared',type=Path);a.add_argument('--index',type=int);a.add_argument('--results-root',type=Path)
    args=a.parse_args()
    if args.mode=='prepare':result=prepare(args.inputs,args.config,args.out)
    elif args.mode=='run':result=execute(args.inputs,args.config,args.prepared,args.out,args.index)
    else:result=collect(args.inputs,args.config,args.results_root,args.out)
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()
