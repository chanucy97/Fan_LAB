"""Bounded, descriptive enzyme-contig reference affinity; no population/FFA calls."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

DEPENDENCY_SHA = '99995a09554e42563a67ad345dedef4d0be03532bc015168887041ec036706d3'
CONFIG_SHA = 'ebabc782561dd2d7456a1c2ed1d9ee0d7c1b2e2b2d72e7837c1f140e300c473d'
POSITIVE = {'AfcA_Bbif_JCM1254', 'AfcB_Bbif_JCM1254', 'LnbB_Bbif_JCM1254', 'NanH2_Binf_ATCC15697', 'Bga42A_Binf_ATCC15697', 'Bga2A_Binf_ATCC15697'}
COMPARATORS = {'NanH1_Binf_ATCC15697', 'Bga42B_Binf_ATCC15697', 'Bga42C_Binf_ATCC15697'}
REGION_NAMES = ['whole_contig', 'local_interval', 'masked_background', 'left_flank', 'right_flank']
SCORE_FIELDS = ['case_id','sample_id','gene_id','trait_id','reference','region','method','region_length','requested_length','raw_window_length','removed_by_masks_nt','segments','full_requested_length_available','mask_interrupted','covered_bases','exact_matching_bases','coverage','covered_identity','whole_region_matching_fraction','equal_bitscore_conflicting_positions','possible_circular_origin','near_equal_other_locus','single_chain_asserted']
POINTER_FIELDS = ['case_id','sample_id','gene_id','trait_id','reference','subject_id','hsp_row_number','query_start0','query_end0','subject_start0','subject_end0','relative_alignment_strand','sample_gene_strand','subject_coding_strand','covered_query_nt','exact_matching_query_nt','overlapping_known_reference_homologue_ids','overlapping_same_trait_reference_homologue_ids','orthology_asserted']


def now(): return datetime.now(timezone.utc).isoformat()
def read(path): return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def write(path,obj): Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda:handle.read(4194304),b''): h.update(block)
    return h.hexdigest()
def seqsha(sequence): return hashlib.sha256(sequence.encode('ascii')).hexdigest()
def rec(path):
    p=Path(path).resolve();before=p.stat();digest=sha(p);after=p.stat()
    if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns): raise ValueError('file changed while hashing')
    return {'path':str(p),'bytes':after.st_size,'sha256':digest}
def same(a,b): return a['bytes']==b['bytes'] and a['sha256']==b['sha256']
def verify(record):
    if not same(record,rec(record['path'])): raise ValueError('source changed: '+record['path'])
def table(path,rows,fields):
    with Path(path).open('w',encoding='utf-8',newline='') as handle:
        writer=csv.DictWriter(handle,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
def read_table(path):
    with Path(path).open(encoding='utf-8',newline='') as handle:return list(csv.DictReader(handle,delimiter='\t'))
def fresh(out):
    out=Path(out).resolve();out.mkdir(parents=True,exist_ok=False);return out


def dependency(path):
    if sha(path)!=DEPENDENCY_SHA: raise ValueError('unexpected carrier dependency bytes')
    spec=importlib.util.spec_from_file_location('enzyme_context_carrier_dependency',path)
    module=importlib.util.module_from_spec(spec);old=sys.dont_write_bytecode;sys.dont_write_bytecode=True
    try:spec.loader.exec_module(module)
    finally:sys.dont_write_bytecode=old
    return module


def merge_intervals(intervals):
    out=[]
    for start,end in sorted(intervals):
        if start>=end:continue
        if out and start<=out[-1][1]:out[-1][1]=max(out[-1][1],end)
        else:out.append([start,end])
    return out


def subtract(intervals,masks):
    result=[]
    for start,end in intervals:
        cursor=start
        for left,right in masks:
            if right<=cursor:continue
            if left>=end:break
            if left>cursor:result.append([cursor,min(left,end)])
            cursor=max(cursor,right)
            if cursor>=end:break
        if cursor<end:result.append([cursor,end])
    return result


def region_set(case,gene,config):
    n=case['contig_length'];gap=config['outer_interval_buffer'];width=config['flank_length']
    masks=merge_intervals((max(0,g['start0']-gap),min(n,g['end0']+gap)) for g in case['mask_orfs'])
    left=max(0,gene['start0']-gap);right=min(n,gene['end0']+gap)
    raw={'whole_contig':[[0,n]],'local_interval':[[gene['start0'],gene['end0']]],'masked_background':[[0,n]],
         'left_flank':[[max(0,left-width),left]],'right_flank':[[right,min(n,right+width)]]}
    regions=[]
    for name in REGION_NAMES:
        original=[p for p in raw[name] if p[0]<p[1]]
        parts=subtract(original,masks) if name in ('masked_background','left_flank','right_flank') else original
        length=sum(b-a for a,b in parts);raw_length=sum(b-a for a,b in original)
        requested=width if name.endswith('_flank') else length
        regions.append({'region':name,'segments':parts,'length':length,'requested_length':requested,
                        'raw_window_length':raw_length,'removed_by_masks_nt':raw_length-length,
                        'full_requested_length_available':length==requested,
                        'mask_interrupted':name.endswith('_flank') and length<raw_length})
    return masks,regions


def validate_manifest(m):
    if m['schema']!='enzyme-context.inputs.v1' or m['cohort_denominator']!=50:raise ValueError('unexpected input schema/cohort')
    if m['dependency']['sha256']!=DEPENDENCY_SHA or m['config_sha256']!=CONFIG_SHA:raise ValueError('dependency/config mismatch')
    if len(m['cases'])!=m['contig_count'] or sum(len(c['genes']) for c in m['cases'])!=m['candidate_count']:raise ValueError('candidate denominator mismatch')
    if sum(len(c['mask_orfs']) for c in m['cases'])!=m['mask_orf_count']:raise ValueError('mask denominator mismatch')
    ids=set();keys=set();genes_seen=set()
    for case in m['cases']:
        if not re.fullmatch(r'[A-Za-z0-9_.-]+',case['case_id']) or case['case_id'] in ids:raise ValueError('unsafe/duplicate case id')
        ids.add(case['case_id']);key=(case['sample_id'],case['contig'])
        if key in keys:raise ValueError('duplicate sample/contig; BLAST must run once per contig')
        keys.add(key);mask={g['gene_id']:g for g in case['mask_orfs']}
        if len(mask)!=len(case['mask_orfs']) or not case['genes']:raise ValueError('duplicate mask or empty case')
        for g in case['mask_orfs']:
            if not 0<=g['start0']<g['end0']<=case['contig_length'] or g['strand'] not in (-1,1):raise ValueError('invalid ORF coordinates')
            expected='positive_enzyme_seed' if g['trait_id'] in POSITIVE else 'specificity_comparator' if g['trait_id'] in COMPARATORS else None
            if g['trait_role']!=expected or expected is None or g['mapping_status'] not in ('high_similarity_candidate','exact_seed_peptide'):raise ValueError('unaccepted target mask')
            if g.get('competitor_ambiguity') not in (None,False,'False'):raise ValueError('ambiguous target mask')
            if g['partial'] not in ('00','01','10','11') or type(g['complete_predicted_orf']) is not bool:raise ValueError('invalid ORF flags')
            for field in ('cds_sha256','protein_sha256'):
                if not re.fullmatch('[0-9a-f]{64}',g[field]):raise ValueError('invalid sequence hash')
        for gene in case['genes']:
            key=(case['sample_id'],gene['gene_id'])
            if key in genes_seen or gene!=mask.get(gene['gene_id']) or gene['trait_role']!='positive_enzyme_seed':raise ValueError('reported gene missing from masks or duplicated')
            genes_seen.add(key)
    references=[r['accession'] for r in m['references']]
    if len(references)<2 or len(set(references))!=len(references):raise ValueError('reference accession aliases must be deduplicated')
    subjects={s['subject_id']:s for s in m['subjects']}
    if len(subjects)!=len(m['subjects']) or any(s['accession'] not in references or s['length']<1 for s in subjects.values()):raise ValueError('invalid subjects')
    for item in m.get('reference_orthologs',[]):
        if item['subject_id'] not in subjects or not 0<=int(item['start0'])<int(item['end0'])<=subjects[item['subject_id']]['length']:raise ValueError('invalid reference homologue pointer')


def inputs_config(inputs,config_path):
    m=read(inputs);config=read(config_path);validate_manifest(m)
    if sha(config_path)!=CONFIG_SHA or config['phenotype_inputs_allowed'] or config['readback_linkage_done']:raise ValueError('frozen config changed')
    return m,config


def selected_contig(case,dep):
    verify(case['assembly']);matches=[seq for name,seq in dep.fasta(case['assembly']['path']) if name==case['contig']]
    if len(matches)!=1:raise ValueError('selected contig absent/duplicated')
    seq=matches[0]
    if len(seq)!=case['contig_length'] or seqsha(seq)!=case['contig_sha256']:raise ValueError('contig source mismatch')
    for gene in case['mask_orfs']:
        cds=seq[gene['start0']:gene['end0']]
        if gene['strand']==-1:cds=cds.translate(str.maketrans('ACGTRYMKSWBDHVN','TGCAYRKMSWVHDBN'))[::-1]
        if seqsha(cds)!=gene['cds_sha256']:raise ValueError('target CDS source mismatch')
    return seq


def parse_hits(path,case,m,dep):
    subjects={s['subject_id']:s for s in m['subjects']};out=[]
    with Path(path).open(encoding='ascii') as handle:
        for line_no,line in enumerate(handle,1):
            values=line.rstrip('\n').split('\t')
            if len(values)!=len(dep.FIELDS):raise ValueError('malformed BLAST row')
            h=dict(zip(dep.FIELDS,values))
            for k in ('length','qlen','qstart','qend','sstart','send','gaps','gapopen'):h[k]=int(h[k])
            for k in ('pident','evalue','bitscore'):
                h[k]=float(h[k])
                if not math.isfinite(h[k]):raise ValueError('nonfinite HSP metric')
            if h['qseqid']!=case['contig'] or h['qlen']!=case['contig_length'] or h['sseqid'] not in subjects:raise ValueError('unbound HSP sequence')
            if not 1<=min(h['qstart'],h['qend'])<=max(h['qstart'],h['qend'])<=h['qlen']:raise ValueError('query bounds')
            if not 1<=min(h['sstart'],h['send'])<=max(h['sstart'],h['send'])<=subjects[h['sseqid']]['length']:raise ValueError('subject bounds')
            if not (0<=h['pident']<=100 and h['bitscore']>0 and h['evalue']>=0 and h['length']>0 and h['gaps']>=0 and h['gapopen']>=0):raise ValueError('invalid HSP metrics')
            # Exhausting the frozen iterator validates BTOP endpoints, including gaps.
            for q,s,match in dep.positions(h):
                if not 0<=q<h['qlen'] or s is not None and not 0<=s<subjects[h['sseqid']]['length']:raise ValueError('BTOP position outside sequence')
            h['hsp_row_number']=line_no;out.append(h)
    return out


def rank(rows,region,config):
    ordered=sorted(rows,key=lambda r:(-(r['whole_region_matching_fraction'] or 0),r['reference']))
    top,second=ordered[:2];value=top['whole_region_matching_fraction'] or 0;margin=value-(second['whole_region_matching_fraction'] or 0)
    available=top['region_length']>0 and top['full_requested_length_available']
    qualified=available and top['coverage']>=config['minimum_region_coverage'] and (top['covered_identity'] or 0)>=config['minimum_covered_identity']
    geometric=not(top['possible_circular_origin'] or top['near_equal_other_locus'])
    unique=qualified and geometric and margin+1e-12>=config['unique_margin']
    if top['region_length']==0:status='no_available_background_or_flank' if region!='local_interval' else 'empty_local_interval'
    elif top['mask_interrupted']:status='masked_flank_unresolved'
    elif not available:status='insufficient_flank_length'
    elif not qualified:status='insufficient_alignment_support'
    elif not geometric:status='reference_locus_geometry_unresolved'
    elif unique:status='unique_panel_affinity'
    else:status='near_tied_panel_affinity'
    return {'best_reference':top['reference'],'runner_up_reference':second['reference'],'best_score':top['whole_region_matching_fraction'],
            'runner_up_score':second['whole_region_matching_fraction'],'margin':margin,'region_length':top['region_length'],
            'requested_length':top['requested_length'],'coverage':top['coverage'],'covered_identity':top['covered_identity'],
            'coverage_qualified':qualified,'unique_panel_affinity':unique,'region_status':status,
            'near_best_references':[r['reference'] for r in ordered if value-(r['whole_region_matching_fraction'] or 0)<config['unique_margin']-1e-12],
            'full_requested_length_available':available,'mask_interrupted':top['mask_interrupted']}


def score(case,hits,m,config,dep):
    subjects={s['subject_id']:s for s in m['subjects']};lengths={k:v['length'] for k,v in subjects.items()}
    grouped=defaultdict(list)
    for h in hits:grouped[subjects[h['sseqid']]['accession']].append(h)
    scores=[];chains={};pointers=[];regions_out={}
    # Union masks are made once per reference/contig, not once per gene.
    for ref in m['references']:
        acc=ref['accession'];refhits=grouped[acc];cov,matched,ties=dep.masks(refhits,case['contig_length'])
        for gene in case['genes']:
            masks,regions=region_set(case,gene,config);regions_out[gene['gene_id']]={'target_masks':masks,'regions':regions,'gene':gene}
            for region in regions:
                parts=region['segments'];isflank=region['region'].endswith('_flank')
                common={'case_id':case['case_id'],'sample_id':case['sample_id'],'gene_id':gene['gene_id'],'trait_id':gene['trait_id'],
                        'reference':acc,'region':region['region'],'region_length':region['length'],'requested_length':region['requested_length'],
                        'raw_window_length':region['raw_window_length'],'removed_by_masks_nt':region['removed_by_masks_nt'],
                        'segments':';'.join(f'{a}:{b}' for a,b in parts),'full_requested_length_available':region['full_requested_length_available'],
                        'mask_interrupted':region['mask_interrupted']}
                c=sum(sum(cov[a:b]) for a,b in parts);mcount=sum(sum(matched[a:b]) for a,b in parts);t=sum(sum(ties[a:b]) for a,b in parts)
                scores.append({**common,'method':'query_position_union',**dep.fractions(c,mcount,region['length']),
                               'equal_bitscore_conflicting_positions':t,'possible_circular_origin':False,'near_equal_other_locus':False,'single_chain_asserted':False})
                if isflank:
                    pieces=[]
                    for start,end in parts:
                        ch=dep.chain(refhits,start,end,lengths,config)
                        # Original block records lack raw-row IDs; recover stable pointers without HSP copies.
                        pointers_by_block=defaultdict(list)
                        for h in refhits:
                            clipped=dep.clipped_hsp(h,start,end,lengths)
                            if clipped is not None:pointers_by_block[tuple(sorted(clipped.items()))].append(h['hsp_row_number'])
                        for block in ch['blocks']:
                            block['hsp_row_numbers']=pointers_by_block[tuple(sorted(block.items()))]
                        pieces.append({'start0':start,'end0':end,'chain':ch})
                    chcovered=sum(p['chain']['covered'] for p in pieces);chmatched=sum(p['chain']['matched'] for p in pieces)
                    origin=any(p['chain']['possible_circular_origin'] for p in pieces);other=any(p['chain']['near_equal_other_locus'] for p in pieces)
                    key=gene['gene_id']+'|'+acc+'|'+region['region']
                    chains[key]={'segments':pieces,'cross_mask_chain_asserted':False,'mask_interrupted':region['mask_interrupted']}
                    scores.append({**common,'method':'independent_segment_chains',**dep.fractions(chcovered,chmatched,region['length']),
                                   'equal_bitscore_conflicting_positions':None,'possible_circular_origin':origin,'near_equal_other_locus':other,
                                   'single_chain_asserted':len(pieces)==1 and bool(pieces[0]['chain']['blocks'])})
            for h in refhits:
                block=dep.clipped_hsp(h,gene['start0'],gene['end0'],lengths)
                if block is None:continue
                start=min(block['subject_at_query_start0'],block['subject_at_query_end0']);end=max(block['subject_at_query_start0'],block['subject_at_query_end0'])+1
                overlaps=[r for r in m.get('reference_orthologs',[]) if r['subject_id']==h['sseqid'] and int(r['start0'])<end and int(r['end0'])>start]
                pointers.append({'case_id':case['case_id'],'sample_id':case['sample_id'],'gene_id':gene['gene_id'],'trait_id':gene['trait_id'],
                                 'reference':acc,'subject_id':h['sseqid'],'hsp_row_number':h['hsp_row_number'],
                                 'query_start0':block['qstart0'],'query_end0':block['qend0'],'subject_start0':start,'subject_end0':end,
                                 'relative_alignment_strand':block['strand'],'sample_gene_strand':gene['strand'],'subject_coding_strand':gene['strand']*block['strand'],
                                 'covered_query_nt':block['covered'],'exact_matching_query_nt':block['matched'],
                                 'overlapping_known_reference_homologue_ids':';'.join(r.get('gene_id',r.get('region_id','')) for r in overlaps),
                                 'overlapping_same_trait_reference_homologue_ids':';'.join(r.get('gene_id',r.get('region_id','')) for r in overlaps if r['trait_id']==gene['trait_id']),
                                 'orthology_asserted':False})
    summaries=[]
    for gene in case['genes']:
        rankings={}
        for region in REGION_NAMES:
            method='independent_segment_chains' if region.endswith('_flank') else 'query_position_union'
            relevant=[r for r in scores if r['gene_id']==gene['gene_id'] and r['region']==region and r['method']==method]
            rankings[region]=rank(relevant,region,config)
        summaries.append({'case_id':case['case_id'],'sample_id':case['sample_id'],'contig':case['contig'],'gene':gene,'region_rankings':rankings,
                          'local_reference_pointer_count':sum(p['gene_id']==gene['gene_id'] for p in pointers),
                          'interpretation':'continuous reference-panel affinity; no sample carrier, strain, orthology, HGT or activity assignment',
                          'population_inputs':[],'phenotype_inputs':[],'read_connection_verified':False})
    return scores,chains,pointers,regions_out,summaries


def output_records(out,exclude):
    records=[]
    for p in sorted(out.iterdir()):
        if p.is_symlink() or p.is_dir():raise ValueError('unexpected output directory/symlink')
        if p.name not in exclude:r=rec(p);r['path']=p.name;records.append(r)
    return records


def finish(out,filename,state):
    if (out/'COMPLETE').exists() or (out/'FAILED.json').exists():raise ValueError('cannot rewrite completed/failed output')
    state.update(status='COMPLETE',finished_utc=now(),outputs=output_records(out,{filename,'COMPLETE','FAILED.json'}))
    write(out/filename,state);(out/'COMPLETE').write_text(sha(out/filename)+'\n',encoding='ascii');return state


def bundle(folder,filename,check_outputs=True):
    folder=Path(folder);path=folder/filename;state=read(path)
    if state['status']!='COMPLETE' or (folder/'COMPLETE').read_text().strip()!=sha(path) or (folder/'FAILED.json').exists():raise ValueError('incomplete result')
    seen=set()
    for record in state['outputs']:
        name=record['path']
        if Path(name).name!=name or '/' in name or '\\' in name or name in seen:raise ValueError('unsafe output record')
        seen.add(name)
        if check_outputs:
            path=folder/name
            if path.is_symlink() or not same(record,rec(path)):raise ValueError('output hash mismatch')
    if check_outputs and seen!={p.name for p in folder.iterdir() if p.name not in {filename,'COMPLETE'}}:raise ValueError('unbound output')
    return state


def command(argv,out,commands):
    ident=len(commands);stdout=out/f'command_{ident:02d}.stdout.txt';stderr=out/f'command_{ident:02d}.stderr.txt'
    entry={'argv':[str(x) for x in argv],'started_utc':now(),'stdout':stdout.name,'stderr':stderr.name};commands.append(entry);write(out/'COMMANDS.json',commands)
    with stdout.open('wb') as so,stderr.open('wb') as se:result=subprocess.run(entry['argv'],stdout=so,stderr=se,check=False,cwd=out)
    entry.update(returncode=result.returncode,finished_utc=now());write(out/'COMMANDS.json',commands)
    if result.returncode:raise RuntimeError('command failed; retained logs: '+str(argv[0]))


def freeze_bindings(inputs,config_path,m,config):
    return {'inputs':rec(inputs),'config':rec(config_path),'core':rec(__file__),'dependency':m['dependency'],
            'tools':{k:rec(shutil.which(config[k]) or config[k]) for k in ('blastn','makeblastdb')}}


def compatible_bindings(saved,inputs,config_path,m):
    for key,path in [('inputs',inputs),('config',config_path),('core',__file__)]:
        if not same(saved[key],rec(path)):raise ValueError('frozen '+key+' changed')
    if saved['dependency']!=m['dependency']:raise ValueError('frozen dependency changed')


def prepare(inputs,config_path,out):
    m,config=inputs_config(inputs,config_path);dep=dependency(m['dependency']['path']);out=fresh(out);commands=[]
    try:
        verify(m['dependency']);verify(m['reference_fasta'])
        for record in m['source_receipts']:verify(record)
        observed={}
        for name,seq in dep.fasta(m['reference_fasta']['path']):
            if name in observed:raise ValueError('duplicate FASTA subject')
            observed[name]=(len(seq),seqsha(seq))
        if observed!={s['subject_id']:(s['length'],s['sequence_sha256']) for s in m['subjects']}:raise ValueError('subject sequence roster mismatch')
        # Bound all selected contigs and all masks, including comparators, before DB creation.
        for case in m['cases']:selected_contig(case,dep)
        bindings=freeze_bindings(inputs,config_path,m,config)
        command([bindings['tools']['makeblastdb']['path'],'-version'],out,commands);command([bindings['tools']['blastn']['path'],'-version'],out,commands)
        command([bindings['tools']['makeblastdb']['path'],'-in',m['reference_fasta']['path'],'-dbtype','nucl','-parse_seqids','-out',out/'reference'],out,commands)
        if not list(out.glob('reference.*')):raise ValueError('database was not created')
        if bindings!=freeze_bindings(inputs,config_path,m,config):raise ValueError('source/tool changed during prepare')
        verify(m['reference_fasta'])
        return finish(out,'PREPARED.json',{'schema':'enzyme-context.prepared.v1','created_utc':now(),'bindings':bindings,
                      'reference_fasta':m['reference_fasta'],'source_receipts':m['source_receipts'],'parameters':config,
                      'database_prefix':str(out/'reference'),'case_count':len(m['cases']),'candidate_count':m['candidate_count'],
                      'alignment_run':False,'phenotype_inputs':[],'BAM_read':False})
    except Exception as exc:write(out/'FAILED.json',{'status':'FAILED','error':repr(exc),'commands':commands});raise


def validate_prepared(p,inputs,config_path,m,config):
    if p['schema']!='enzyme-context.prepared.v1' or p['status']!='COMPLETE':raise ValueError('wrong preparation')
    compatible_bindings(p['bindings'],inputs,config_path,m)
    if p['parameters']!=config or p['reference_fasta']!=m['reference_fasta'] or p['source_receipts']!=m['source_receipts']:raise ValueError('preparation sources mismatch')
    if p['case_count']!=len(m['cases']) or p['candidate_count']!=m['candidate_count']:raise ValueError('preparation denominator mismatch')


def execute(inputs,config_path,prepared,outroot,index):
    inputs=Path(inputs);config_path=Path(config_path);prepared=Path(prepared)
    m,config=inputs_config(inputs,config_path)
    if index is None or not 0<=index<len(m['cases']):raise ValueError('invalid case index')
    p=bundle(prepared.parent,prepared.name);validate_prepared(p,inputs,config_path,m,config)
    for tool in p['bindings']['tools'].values():verify(tool)
    dep=dependency(m['dependency']['path']);case=m['cases'][index];out=fresh(Path(outroot)/case['case_id']);commands=[]
    state={'schema':'enzyme-context.case.v1','status':'RUNNING','created_utc':now(),'case_index':index,'case_id':case['case_id'],
           'source':case,'bindings':p['bindings'],'prepared':rec(prepared),'phenotype_inputs':[],'population_inputs':[],'BAM_read':False}
    write(out/'RUN.json',state)
    try:
        shutil.copyfile(prepared,out/'PREPARED_SNAPSHOT.json')
        sequence=selected_contig(case,dep);query=out/'selected_contig.fna';query.write_text('>'+case['contig']+'\n'+sequence+'\n',encoding='ascii')
        path=out/'all_HSPs_btop.tsv'
        command([p['bindings']['tools']['blastn']['path'],'-task',config['blast_task'],'-word_size',config['blast_word_size'],'-query',query,'-db',p['database_prefix'],
                 '-evalue',config['blast_evalue'],'-num_threads',config['threads'],'-max_target_seqs',config['blast_max_target_seqs'],
                 '-outfmt','6 '+' '.join(dep.FIELDS),'-out',path],out,commands)
        hits=parse_hits(path,case,m,dep);scores,chains,pointers,regions,summaries=score(case,hits,m,config,dep)
        table(out/'REGION_REFERENCE_SCORES.tsv',scores,SCORE_FIELDS);table(out/'LOCAL_REFERENCE_POINTERS.tsv',pointers,POINTER_FIELDS)
        write(out/'SEGMENT_CHAINS.json',chains);write(out/'REGIONS.json',regions);write(out/'GENE_SUMMARIES.json',summaries)
        compatible_bindings(state['bindings'],inputs,config_path,m);verify(m['dependency']);verify(case['assembly'])
        for tool in state['bindings']['tools'].values():verify(tool)
        if rec(prepared)!=state['prepared'] or not same(state['prepared'],rec(out/'PREPARED_SNAPSHOT.json')):raise ValueError('preparation changed during run')
        # Database and command-log bytes remain bound to the same verified preparation.
        bundle(prepared.parent,prepared.name)
        return finish(out,'RUN.json',{**state,'HSPs':len(hits),'gene_count':len(case['genes']),'score_rows':len(scores),'local_pointer_rows':len(pointers)})
    except Exception as exc:write(out/'FAILED.json',{**state,'status':'FAILED','error':repr(exc),'commands':commands});raise


def collect(inputs,config_path,results_root,out):
    m,config=inputs_config(inputs,config_path);out=fresh(out)
    try:
        all_scores=[];all_pointers=[];summaries=[];runs=[];prepared_hashes=set()
        expected_regions={(name,'independent_segment_chains' if name.endswith('_flank') else 'query_position_union') for name in REGION_NAMES}
        expected_regions|={(name,'query_position_union') for name in ('left_flank','right_flank')}
        for index,case in enumerate(m['cases']):
            folder=Path(results_root)/case['case_id'];run=bundle(folder,'RUN.json')
            if run['schema']!='enzyme-context.case.v1' or run['case_id']!=case['case_id'] or run['case_index']!=index or run['source']!=case:raise ValueError('case identity mismatch')
            compatible_bindings(run['bindings'],inputs,config_path,m)
            snapshot=folder/'PREPARED_SNAPSHOT.json'
            if not same(run['prepared'],rec(snapshot)):raise ValueError('preparation snapshot mismatch')
            p=read(snapshot);validate_prepared(p,inputs,config_path,m,config)
            if run['bindings']!=p['bindings']:raise ValueError('case/preparation binding mismatch')
            prepared_hashes.add(sha(snapshot));runs.append(rec(folder/'RUN.json'))
            case_scores=read_table(folder/'REGION_REFERENCE_SCORES.tsv')
            expected={(g['gene_id'],r['accession'],name,method) for g in case['genes'] for r in m['references'] for name,method in expected_regions}
            observed={(r['gene_id'],r['reference'],r['region'],r['method']) for r in case_scores}
            if len(case_scores)!=len(expected) or observed!=expected or any(r['case_id']!=case['case_id'] or r['sample_id']!=case['sample_id'] for r in case_scores):raise ValueError('score denominator mismatch')
            current=read(folder/'GENE_SUMMARIES.json')
            if len(current)!=len(case['genes']) or {x['gene']['gene_id'] for x in current}!={g['gene_id'] for g in case['genes']}:raise ValueError('summary gene denominator mismatch')
            for row in current:
                gene=next(g for g in case['genes'] if g['gene_id']==row['gene']['gene_id'])
                if row['gene']!=gene or row['case_id']!=case['case_id'] or row['sample_id']!=case['sample_id'] or set(row['region_rankings'])!=set(REGION_NAMES):raise ValueError('summary source mismatch')
            pointers=read_table(folder/'LOCAL_REFERENCE_POINTERS.tsv')
            if any(r['case_id']!=case['case_id'] or r['gene_id'] not in {g['gene_id'] for g in case['genes']} for r in pointers):raise ValueError('pointer source mismatch')
            all_scores.extend(case_scores);all_pointers.extend(pointers);summaries.extend(current)
        if len(prepared_hashes)!=1:raise ValueError('mixed preparation generations')
        table(out/'ALL_REFERENCE_SCORES.tsv',all_scores,SCORE_FIELDS);table(out/'ALL_LOCAL_REFERENCE_POINTERS.tsv',all_pointers,POINTER_FIELDS)
        write(out/'ALL_GENE_SUMMARIES.json',summaries)
        flattened=[]
        for row in summaries:
            item={k:row[k] for k in ('case_id','sample_id','contig')};item.update(row['gene'])
            for region,ranking in row['region_rankings'].items():
                for key,value in ranking.items():item[region+'__'+key]=';'.join(value) if isinstance(value,list) else value
            flattened.append(item)
        table(out/'ALL_GENES.tsv',flattened,list(flattened[0]))
        return finish(out,'COLLECTION.json',{'schema':'enzyme-context.collection.v1','created_utc':now(),'bindings':{'inputs':rec(inputs),'config':rec(config_path),'core':rec(__file__),'dependency':m['dependency']},
                      'case_runs':runs,'prepared_sha256':next(iter(prepared_hashes)),'cohort_denominator':50,'contigs':len(m['cases']),
                      'genes':len(summaries),'score_rows':len(all_scores),'reference_count':len(m['references']),
                      'phenotype_inputs':[],'population_inputs':[],'BAM_read':False,'carrier_or_strain_asserted':False,'HGT_asserted':False})
    except Exception as exc:write(out/'FAILED.json',{'status':'FAILED','error':repr(exc)});raise


def main():
    parser=argparse.ArgumentParser(description=__doc__);sub=parser.add_subparsers(dest='mode',required=True)
    for mode in ('prepare','run','collect'):
        p=sub.add_parser(mode);p.add_argument('--inputs',type=Path,required=True);p.add_argument('--config',type=Path,required=True)
        if mode=='run':
            p.add_argument('--prepared',type=Path,required=True);p.add_argument('--outroot',type=Path,required=True);p.add_argument('--index',type=int,required=True)
        else:p.add_argument('--out',type=Path,required=True)
        if mode=='collect':p.add_argument('--results-root',type=Path,required=True)
    args=parser.parse_args()
    if args.mode=='prepare':result=prepare(args.inputs,args.config,args.out)
    elif args.mode=='run':result=execute(args.inputs,args.config,args.prepared,args.outroot,args.index)
    else:result=collect(args.inputs,args.config,args.results_root,args.out)
    print(json.dumps({'status':result['status'],'schema':result['schema']}))


if __name__=='__main__':main()
