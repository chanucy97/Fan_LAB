"""Frozen inStrain 1.11 comparison and metadata-only descriptive summaries.

prepare requires every one of the 50 profile jobs to be COMPLETE, including
explicit zero-analysable-pair outcomes. No nutrient/phenotype values are used.
run executes one reference genome. summarize retains every possible sample pair.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
from itertools import combinations
import json
import math
import os
from pathlib import Path
import statistics
import subprocess

GATE = {'breadth_minCov': .5, 'min_cov': 5, 'min_freq': .05, 'fdr': 1e-6,
        'joint_whole_genome_fraction': .5, 'popANI': .99999}
PROFILE_KEYS = ['min_mapq', 'min_read_ani', 'min_cov', 'min_freq', 'snv_fdr',
                'min_insert', 'max_insert_relative', 'bowtie_max_insert', 'seed', 'skip_mm_profiling']
PAIR_COLUMNS = ['sample1','sample2','dyad1','dyad2','window1','window2',
                'nominal_day1','nominal_day2','pair_type','consecutive_observed',
                'nominal_day_gap','actual_day_gap_known']


def now(): return datetime.now(timezone.utc).isoformat()


def record(path):
    p=Path(path).resolve(strict=True);h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return {'path':str(p),'bytes':p.stat().st_size,'sha256':h.hexdigest()}


def verify(rec):
    actual=record(rec['path'])
    if actual['bytes']!=rec['bytes'] or actual['sha256']!=rec['sha256']:
        raise ValueError('Frozen file changed: '+rec['path'])


def read_json(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def write_json(path,data):
    Path(path).write_text(json.dumps(data,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')


def read_tsv(path,delimiter='\t'):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))


def write_tsv(path,rows,columns=None):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=columns or list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


def number(value,nullable=False):
    if value in (None,'','nan','NaN') and nullable:return None
    x=float(value)
    if not math.isfinite(x):
        if nullable:return None
        raise ValueError('Non-finite numeric field')
    return x


def true(value):return value is True or str(value).lower()=='true'


def metadata(samples,within_path=None,enforce_cohort=True):
    """Only identity and nominal timing fields enter the comparison metadata."""
    cleaned=[{k:r[k] for k in ['sample_id','dyad_id','window','nominal_day']} for r in samples]
    ids=[r['sample_id'] for r in cleaned]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate sample IDs')
    groups=defaultdict(list)
    for r in cleaned:
        r['dyad_id']=str(r['dyad_id']);r['nominal_day']=int(r['nominal_day'])
        groups[r['dyad_id']].append(r)
    consecutive=set()
    for rs in groups.values():
        ordered=sorted(rs,key=lambda r:(r['nominal_day'],r['sample_id']))
        if len({r['nominal_day'] for r in ordered})!=len(ordered):raise ValueError('Duplicated nominal visit in infant')
        consecutive.update(tuple(sorted([a['sample_id'],b['sample_id']])) for a,b in zip(ordered,ordered[1:]))
    pairs=[]
    for a,b in combinations(sorted(cleaned,key=lambda r:r['sample_id']),2):
        same=a['dyad_id']==b['dyad_id'];key=(a['sample_id'],b['sample_id'])
        pairs.append(dict(zip(PAIR_COLUMNS,[*key,a['dyad_id'],b['dyad_id'],a['window'],b['window'],
            a['nominal_day'],b['nominal_day'],'within_infant' if same else 'between_infants',
            key in consecutive,abs(a['nominal_day']-b['nominal_day']),False])))
    within={tuple(sorted([r['sample1'],r['sample2']])):r for r in pairs if r['pair_type']=='within_infant'}
    if within_path:
        # Deliberately ignore all phenotype columns in the historical pair table.
        selected=[{k:r[k] for k in ['dyad_id','sample_before','sample_after','window_before','window_after',
                                   'nominal_day_gap','actual_day_gap_known','consecutive_observed']}
                  for r in read_tsv(within_path)]
        if len(selected)!=len(within):raise ValueError('Historical within-pair row count mismatch')
        seen=set()
        lookup={r['sample_id']:r for r in cleaned}
        for r in selected:
            key=tuple(sorted([r['sample_before'],r['sample_after']]))
            if key not in within or key in seen:raise ValueError('Historical within-pair membership mismatch')
            seen.add(key);calculated=within[key]
            if str(r['dyad_id'])!=calculated['dyad1']:raise ValueError('Historical dyad identity mismatch')
            if int(r['nominal_day_gap'])!=calculated['nominal_day_gap']:raise ValueError('Nominal gap mismatch')
            if true(r['consecutive_observed'])!=calculated['consecutive_observed']:raise ValueError('Consecutive observed mismatch')
            if true(r['actual_day_gap_known']):raise ValueError('Unexpected actual-date metadata; review before changing protocol')
            if r['window_before']!=lookup[r['sample_before']]['window'] or r['window_after']!=lookup[r['sample_after']]['window']:
                raise ValueError('Historical visit window mismatch')
    shape={'samples':len(cleaned),'infants':len(groups),'repeated_infants':sum(len(v)>1 for v in groups.values()),
           'all_sample_pairs':len(pairs),'within_pairs':len(within),'consecutive_within_pairs':len(consecutive),
           'between_pairs':len(pairs)-len(within)}
    expected={'samples':50,'infants':30,'repeated_infants':16,'all_sample_pairs':1225,
              'within_pairs':25,'consecutive_within_pairs':20,'between_pairs':1200}
    if enforce_cohort and shape!=expected:raise ValueError('Fixed cohort shape differs: '+str(shape))
    return cleaned,pairs,shape


def latest_genome_rows(rows):
    grouped=defaultdict(list)
    for r in rows:grouped[r['genome']].append(r)
    result={}
    for genome,items in grouped.items():
        ranks=[number(r.get('mm',0)) for r in items]
        top=max(ranks);chosen=[r for r,x in zip(items,ranks) if x==top]
        if len(chosen)!=1:raise ValueError('Duplicated maximum-mm genome rows: '+genome)
        result[genome]=chosen[0]
    return result


def reference_lengths(fasta,stb_path):
    stb={}
    for line in Path(stb_path).read_text(encoding='utf-8').splitlines():
        scaffold,genome=line.split('\t')
        if scaffold in stb:raise ValueError('Duplicate scaffold in STB')
        stb[scaffold]=genome
    lengths={};name=None
    with Path(fasta).open(encoding='ascii') as f:
        for line in f:
            if line.startswith('>'):
                name=line[1:].split()[0]
                if name in lengths:raise ValueError('Duplicate scaffold in FASTA')
                lengths[name]=0
            elif line.strip():
                if name is None:raise ValueError('FASTA sequence before header')
                lengths[name]+=len(line.strip())
    if set(lengths)!=set(stb):raise ValueError('FASTA and STB scaffold sets differ')
    genomes=Counter()
    for scaffold,length in lengths.items():genomes[stb[scaffold]]+=length
    return dict(genomes)


def profile_metadata(root,sample,expected_refs,genomes):
    out=Path(root)/sample['sample_id'];receipt=out/'RUN.json';run=read_json(receipt)
    if run['status']!='COMPLETE' or not (out/'COMPLETE').is_file():raise ValueError('Profile incomplete: '+str(out))
    if run['sample']['sample_id']!=sample['sample_id']:raise ValueError('Profile identity mismatch')
    actual_refs={r['sha256'] for r in run['references']}
    if actual_refs!=set(expected_refs.values()):raise ValueError('Profiles do not share frozen reference files')
    p={'sample_id':sample['sample_id'],'profile_dir':str((out/'instrain').resolve()),
       'status':run['profile_status'],'run_receipt':record(receipt),'compact_bindings':[],
       'profile_parameters':{k:run['config'][k] for k in PROFILE_KEYS},
       'profile_script_sha256':run['bindings']['script']['sha256'],
       'base_quality_source_sha256':run['instrain_base_quality_source']['sha256'],
       'retained_pairs':run['stage_counts']['pair_filter']['retained_pairs']}
    missing={g:{'breadth_minCov':0.,'coverage':0.,'nucl_diversity':None,'nucl_diversity_rarefied':None,
                'profile_row_present':False,'sample_eligible':False}
             for g in genomes}
    if run['profile_status']=='NO_ANALYSABLE_PAIRS_NOT_BIOLOGICAL_ABSENCE':
        if p['retained_pairs']!=0:raise ValueError('No-analysable status contradicts retained pairs')
        p.update(bam_name=None,genomes=missing);return p
    profile=Path(p['profile_dir']);attrs_path=profile/'raw_data/attributes.tsv'
    attrs={r['name']:r for r in read_tsv(attrs_path)}
    if attrs['version']['value']!='1.11.0':raise ValueError('Wrong profile version')
    if Path(attrs['location']['value']).resolve()!=profile.resolve():raise ValueError('Relocated inStrain profile requires explicit review')
    name=Path(attrs['bam_loc']['value']).name
    if name!=sample['sample_id']+'.bam':raise ValueError('inStrain comparison BAM basename must be unique sample_id.bam')
    raw_info=profile/'raw_data'/Path(attrs['genome_level_info']['value']).name
    # SNVprofile._store_pandas stores comma-separated indexed CSV, despite TSV exports.
    raw=latest_genome_rows(read_tsv(raw_info,delimiter=','))
    outputs=list((profile/'output').glob('*genome_info.tsv'))
    if len(outputs)!=1:raise ValueError('Exactly one compact genome_info output expected')
    compact=latest_genome_rows(read_tsv(outputs[0]))
    known_outputs={str(Path(r['path']).resolve()):r for r in run['outputs']}
    for f in [attrs_path,raw_info,outputs[0]]:
        r=record(f);saved=known_outputs.get(r['path'])
        if not saved or r['sha256']!=saved['sha256']:raise ValueError('Profile output changed after completion: '+str(f))
        p['compact_bindings'].append(r)
    if set(raw)-set(genomes):raise ValueError('Unknown genome in profile')
    for genome,row in raw.items():
        breadth=number(row['breadth_minCov']);coverage=number(row.get('coverage'),True)
        if not 0<=breadth<=1:raise ValueError('Invalid breadth')
        if genome not in compact or not math.isclose(breadth,number(compact[genome]['breadth_minCov']),abs_tol=1e-10):
            raise ValueError('Compact and stored genome-level breadth disagree')
        missing[genome]={'breadth_minCov':breadth,'coverage':coverage,'nucl_diversity':number(row.get('nucl_diversity'),True),
                         'nucl_diversity_rarefied':number(row.get('nucl_diversity_rarefied'),True),'profile_row_present':True,
                         'sample_eligible':breadth>=GATE['breadth_minCov']}
    p.update(bam_name=name,genomes=missing)
    return p


def prepare(config_path):
    config=read_json(config_path);dest=Path(config['prepared_dir'])
    samples,pairs,shape=metadata(read_json(config['sample_manifest']),config['within_pairs'])
    gs=read_tsv(config['genome_metadata'])
    genomes={r['assembly_accession']:{'length':int(r['bases']),'taxon_group':r['taxon_group']} for r in gs}
    if len(genomes)!=11 or len(gs)!=11:raise ValueError('Expected frozen eleven-genome panel')
    if reference_lengths(config['reference_fasta'],config['reference_stb'])!={g:r['length'] for g,r in genomes.items()}:
        raise ValueError('Genome metadata lengths differ from entire frozen reference FASTA')
    bindings={k:record(config[k]) for k in ['sample_manifest','within_pairs','genome_metadata','reference_fasta','reference_stb','reference_genes']}
    refs={k:bindings[k]['sha256'] for k in ['reference_fasta','reference_stb','reference_genes']}
    profiles=[profile_metadata(config['profiles_root'],r,refs,genomes) for r in samples]
    for key in ['profile_parameters','profile_script_sha256','base_quality_source_sha256']:
        if len({json.dumps(p[key],sort_keys=True) for p in profiles})!=1:raise ValueError('Heterogeneous profile protocol: '+key)
    if profiles[0]['profile_parameters']['min_cov']!=GATE['min_cov']:raise ValueError('Profile depth differs from comparison depth')
    names=[p['bam_name'] for p in profiles if p['bam_name']]
    if len(names)!=len(set(names)):raise ValueError('Duplicate inStrain BAM basename')
    result={'schema':'strain-comparison-prepared.v1','created_utc':now(),'gate':GATE,'config':config,
            'config_binding':record(config_path),'script_binding':record(__file__),'source_bindings':bindings,
            'samples':samples,'pairs':pairs,'cohort_shape':shape,'genomes':genomes,'profiles':profiles,
            'phenotype_values_used':False}
    dest.mkdir(parents=True,exist_ok=False)
    write_json(dest/'PREPARED.json',result)
    write_tsv(dest/'ALL_PAIR_METADATA.tsv',pairs,PAIR_COLUMNS)
    eligible=[{'sample_id':p['sample_id'],'genome':g,'profile_status':p['status'],**r} for p in profiles for g,r in p['genomes'].items()]
    write_tsv(dest/'SAMPLE_GENOME_ELIGIBILITY.tsv',eligible)
    write_json(dest/'PREPARATION_COMPLETE.json',{'status':'COMPLETE','shape':shape,'bindings':[record(dest/n) for n in ['PREPARED.json','ALL_PAIR_METADATA.tsv','SAMPLE_GENOME_ELIGIBILITY.tsv']]})
    return {'status':'COMPLETE','prepared':str(dest/'PREPARED.json'),'shape':shape}


def compare_argv(executable,inputs,out,stb,genome,threads):
    return [executable,'compare','-i',str(inputs),'-o',str(out),'-p',str(threads),'-s',str(stb),
            '--genome',genome,'--database_mode','--breadth','0.5','-c','5','-f','0.05','-fdr','1e-6',
            '-cov','0.5','-ani','0.99999','--skip_plot_generation','--store_coverage_overlap']


def run_genome(prepared_path,index):
    prepared_path=Path(prepared_path);p=read_json(prepared_path);config=p['config']
    genomes=list(p['genomes'])
    if not 0<=index<len(genomes):raise ValueError('Genome array index outside frozen panel')
    if p['gate']!=GATE:raise ValueError('Comparison gate changed')
    verify(p['script_binding']);verify(p['config_binding'])
    for k in ['reference_fasta','reference_stb','reference_genes','genome_metadata']:verify(p['source_bindings'][k])
    genome=genomes[index];out=Path(config['output_root'])/genome
    selected=[x for x in p['profiles'] if x['genomes'][genome]['sample_eligible']]
    receipt={'schema':'strain-comparison-run.v1','status':'RUNNING','started_utc':now(),'genome':genome,
             'reference':p['genomes'][genome],'gate':GATE,'prepared_binding':record(prepared_path),
             'eligible_samples':[x['sample_id'] for x in selected],'eligible_sample_count':len(selected),
             'expected_pairs':len(selected)*(len(selected)-1)//2,'command':None}
    if out.exists():
        saved=read_json(out/'RUN.json')
        if saved.get('status')!='COMPLETE' or saved['prepared_binding']!=receipt['prepared_binding']:
            raise FileExistsError('Preserve partial or differently configured comparison; use a new output root')
        for item in saved['outputs']:verify(item)
        return {'status':'ALREADY_COMPLETE_HASH_VERIFIED','genome':genome}
    out.mkdir(parents=True,exist_ok=False);write_json(out/'RUN.json',receipt)
    try:
        for profile in selected:
            verify(profile['run_receipt'])
            for item in profile['compact_bindings']:verify(item)
        if len(selected)<2:
            receipt['comparison_status']='NO_ELIGIBLE_PROFILE' if not selected else 'SINGLE_ELIGIBLE_PROFILE'
        else:
            inputs=out/'INPUT_PROFILES.txt'
            inputs.write_text(''.join(x['profile_dir']+'\n' for x in selected),encoding='utf-8')
            environment=dict(os.environ)
            bins=[str(Path(config['instrain']).parent)]
            if config.get('samtools'):bins.append(str(Path(config['samtools']).parent))
            environment['PATH']=os.pathsep.join(bins+[environment.get('PATH','')])
            receipt['child_path_prefixes']=bins
            version=subprocess.run([config['instrain'],'compare','--version'],capture_output=True,text=True,check=True,env=environment)
            (out/'VERSION.txt').write_text(version.stdout+version.stderr,encoding='utf-8')
            if '1.11.0' not in version.stdout+version.stderr:raise ValueError('Unexpected inStrain version')
            argv=compare_argv(config['instrain'],inputs,out/'instrain',config['reference_stb'],genome,config['threads'])
            receipt['command']=argv;write_json(out/'RUN.json',receipt)
            with (out/'COMPARE.stdout.log').open('wb') as stdout,(out/'COMPARE.stderr.log').open('wb') as stderr:
                proc=subprocess.run(argv,stdout=stdout,stderr=stderr,check=False,env=environment)
            receipt['returncode']=proc.returncode
            if proc.returncode:raise RuntimeError('inStrain compare failed; logs and output preserved')
            tables=list((out/'instrain/output').glob('*genomeWide_compare.tsv'))
            if len(tables)!=1:raise RuntimeError('Expected one genome-wide comparison table')
            receipt['comparison_table']=str(tables[0]);receipt['comparison_status']='COMPARISON_GENERATED'
        receipt['outputs']=[record(f) for f in sorted(out.rglob('*')) if f.is_file() and f.name not in ['RUN.json','COMPLETE']]
        receipt.update(status='COMPLETE',completed_utc=now());write_json(out/'RUN.json',receipt)
        (out/'COMPLETE').write_text(receipt['completed_utc']+'\n',encoding='ascii')
        return {'status':'COMPLETE','genome':genome,'eligible_samples':len(selected),'comparison_status':receipt['comparison_status']}
    except BaseException as exc:
        receipt.update(status='FAILED',error_type=type(exc).__name__,error=str(exc),failed_utc=now())
        write_json(out/'RUN.json',receipt);raise


def normalized_comparisons(rows,genome,length,name_to_sample):
    result={}
    for row in rows:
        if row['genome']!=genome:raise ValueError('Unexpected genome in comparison result')
        names=[row['name1'],row['name2']]
        if any(n not in name_to_sample for n in names):raise ValueError('Unknown sample name in comparison result')
        key=tuple(sorted(name_to_sample[n] for n in names))
        if key[0]==key[1] or key in result:raise ValueError('Self or duplicated unordered pair in result')
        compared=number(row['compared_bases_count'])
        if compared!=int(compared) or not 0<=compared<=length:raise ValueError('Invalid compared base count')
        compared=int(compared);ani=number(row.get('popANI'),True)
        snps=number(row.get('population_SNPs'),True)
        con_ani=number(row.get('conANI'),True);con_snps=number(row.get('consensus_SNPs'),True)
        if snps is not None and (snps!=int(snps) or not 0<=snps<=compared):raise ValueError('Invalid population SNP count')
        if ani is not None and not 0<=ani<=1:raise ValueError('Invalid popANI')
        if compared>0 and ani is None:raise ValueError('Missing popANI despite compared positions')
        if compared>0 and snps is not None and not math.isclose(ani,1-snps/compared,abs_tol=1e-10,rel_tol=0):
            raise ValueError('popANI inconsistent with population SNP count')
        if con_snps is not None and (con_snps!=int(con_snps) or not 0<=con_snps<=compared):raise ValueError('Invalid consensus SNP count')
        if con_ani is not None and not 0<=con_ani<=1:raise ValueError('Invalid conANI')
        if compared>0 and con_ani is not None and con_snps is not None and not math.isclose(con_ani,1-con_snps/compared,abs_tol=1e-10,rel_tol=0):
            raise ValueError('conANI inconsistent with consensus SNP count')
        fraction=compared/length
        reported=number(row.get('percent_compared'),True)
        result[key]={'compared_positions':compared,'full_reference_length':length,
                     'joint_whole_genome_fraction':fraction,'popANI':ani,'population_SNPs':snps,
                     'conANI':con_ani,'consensus_SNPs':con_snps,
                     'inStrain_percent_compared':reported,
                     'reported_denominator_matches_full_reference':reported is not None and math.isclose(fraction,reported,abs_tol=1e-8,rel_tol=0),
                     'coverage_overlap':number(row.get('coverage_overlap'),True)}
    return result


def pair_evidence(pair,genome,profiles,comparison,length):
    a=profiles[pair['sample1']];b=profiles[pair['sample2']]
    ea=a['genomes'][genome];eb=b['genomes'][genome]
    both=ea['sample_eligible'] and eb['sample_eligible']
    values={'compared_positions':None,'full_reference_length':length,'joint_whole_genome_fraction':None,
            'popANI':None,'population_SNPs':None,'conANI':None,'consensus_SNPs':None,'inStrain_percent_compared':None,
            'reported_denominator_matches_full_reference':None,'coverage_overlap':None}
    if comparison is not None and not both:raise ValueError('Compared a profile excluded by frozen breadth gate')
    if not both:
        status='NO_ANALYSABLE_PAIRS' if not a['retained_pairs'] or not b['retained_pairs'] else 'SAMPLE_BREADTH_INELIGIBLE'
    elif comparison is None:status='EXPECTED_COMPARISON_MISSING'
    else:
        values.update(comparison)
        if comparison['joint_whole_genome_fraction']<GATE['joint_whole_genome_fraction']:status='JOINT_COVERAGE_INELIGIBLE'
        elif comparison['popANI'] is None:raise ValueError('Eligible joint genome without ANI')
        elif comparison['popANI']>=GATE['popANI']:status='ELIGIBLE_HIGH_SIMILARITY'
        else:status='ELIGIBLE_BELOW_SIMILARITY_THRESHOLD'
    return {'genome':genome,**pair,'breadth_minCov1':ea['breadth_minCov'],'breadth_minCov2':eb['breadth_minCov'],
            'coverage1':ea.get('coverage'),'coverage2':eb.get('coverage'),
            'nucl_diversity1':ea.get('nucl_diversity'),'nucl_diversity2':eb.get('nucl_diversity'),
            'nucl_diversity_rarefied1':ea.get('nucl_diversity_rarefied'),'nucl_diversity_rarefied2':eb.get('nucl_diversity_rarefied'),
            'both_sample_breadth_eligible':both,**values,'comparison_status':status,
            'jointly_evaluable':status.startswith('ELIGIBLE_'),'high_similarity':status=='ELIGIBLE_HIGH_SIMILARITY'}


def describe(rows):
    eligible=[r for r in rows if r['jointly_evaluable']]
    high=[r for r in eligible if r['high_similarity']]
    anis=[r['popANI'] for r in eligible]
    infants={str(r[k]) for r in eligible for k in ['dyad1','dyad2']}
    infant_pairs={tuple(sorted([str(r['dyad1']),str(r['dyad2'])])) for r in eligible}
    high_infants={str(r[k]) for r in high for k in ['dyad1','dyad2']}
    high_infant_pairs={tuple(sorted([str(r['dyad1']),str(r['dyad2'])])) for r in high}
    return {'possible_visit_pairs':len(rows),'sample_breadth_eligible_pairs':sum(r['both_sample_breadth_eligible'] for r in rows),
            'jointly_evaluable_pairs':len(eligible),'high_similarity_pairs':len(high),
            'eligible_below_similarity_pairs':len(eligible)-len(high),'unresolved_pairs':len(rows)-len(eligible),
            'evaluable_infants':len(infants),'evaluable_infant_pairs':len(infant_pairs),
            'high_similarity_infants':len(high_infants),'high_similarity_infant_pairs':len(high_infant_pairs),
            'high_similarity_fraction_among_evaluable':len(high)/len(eligible) if eligible else None,
            'popANI_median':statistics.median(anis) if anis else None,'popANI_min':min(anis) if anis else None,
            'popANI_max':max(anis) if anis else None,
            'conANI_median':statistics.median([r['conANI'] for r in eligible if r['conANI'] is not None]) if any(r['conANI'] is not None for r in eligible) else None,
            'joint_fraction_median':statistics.median(r['joint_whole_genome_fraction'] for r in eligible) if eligible else None}


def between_shared_context(evidence):
    """Direct within-infant compatibility only; never construct transitive clones."""
    within=defaultdict(list)
    for r in evidence:
        if r['pair_type']=='within_infant' and r['high_similarity']:
            within[(r['genome'],r['sample1'])].append((r['sample2'],r['nominal_day_gap']))
            within[(r['genome'],r['sample2'])].append((r['sample1'],r['nominal_day_gap']))
    result=[]
    for r in evidence:
        if r['pair_type']!='between_infants' or not r['high_similarity']:continue
        row=dict(r)
        for side in [1,2]:
            related=within[(r['genome'],r['sample'+str(side)])]
            row['related_within_samples'+str(side)]=';'.join(x[0] for x in sorted(related))
            row['related_within_pairs'+str(side)]=len(related)
            row['related_within_min_nominal_gap'+str(side)]=min(x[1] for x in related) if related else None
            row['related_within_max_nominal_gap'+str(side)]=max(x[1] for x in related) if related else None
        result.append(row)
    return result


def summarize(prepared_path,destination):
    p=read_json(prepared_path);dest=Path(destination);profiles={r['sample_id']:r for r in p['profiles']}
    names={r['bam_name']:r['sample_id'] for r in p['profiles'] if r['bam_name']}
    evidence=[];run_bindings=[];genome_status=[]
    for genome,meta in p['genomes'].items():
        runpath=Path(p['config']['output_root'])/genome/'RUN.json';run=read_json(runpath)
        if run['status']!='COMPLETE' or run['prepared_binding']!=record(prepared_path):raise ValueError('Comparison task incomplete or mismatched')
        for r in run['outputs']:verify(r)
        raw=read_tsv(run['comparison_table']) if run.get('comparison_table') else []
        normalized=normalized_comparisons(raw,genome,meta['length'],names)
        genome_status.append({'genome':genome,'taxon_group':meta['taxon_group'],
                              'eligible_samples':run['eligible_sample_count'],'expected_comparison_pairs':run['expected_pairs'],
                              'observed_comparison_pairs':len(normalized),'task_status':run['comparison_status']})
        evidence.extend(pair_evidence(pair,genome,profiles,normalized.get((pair['sample1'],pair['sample2'])),meta['length']) for pair in p['pairs'])
        run_bindings.append(record(runpath))
    grouped=defaultdict(list);infant_grouped=defaultdict(list)
    for row in evidence:
        grouped[(row['genome'],row['pair_type'])].append(row)
        if row['pair_type']=='within_infant' and row['consecutive_observed']:
            grouped[(row['genome'],'within_infant_consecutive')].append(row)
        if row['pair_type']=='within_infant':
            infant_grouped[(row['genome'],row['dyad1'])].append(row)
    descriptives=[{'genome':g,'pair_group':kind,**describe(rows)} for (g,kind),rows in sorted(grouped.items())]
    infant_rows=[{'genome':g,'dyad_id':infant,**describe(rows)} for (g,infant),rows in sorted(infant_grouped.items())]
    statuses=Counter(r['comparison_status'] for r in evidence)
    shared=between_shared_context(evidence)
    result={'schema':'strain-comparison-summary.v1','status':'PASS' if not statuses['EXPECTED_COMPARISON_MISSING'] else 'REVIEW_MISSING_COMPARISONS',
            'created_utc':now(),'gate':GATE,'cohort_shape':p['cohort_shape'],'genome_count':len(p['genomes']),
            'sample_pair_genome_rows':len(evidence),'comparison_status_counts':dict(statuses),
            'shared_compatible_between_visit_pairs':len(shared),
            'shared_compatible_between_infant_pairs_across_genomes':len({tuple(sorted([r['dyad1'],r['dyad2']])) for r in shared}),
            'genome_status':genome_status,'descriptives':descriptives,'prepared_binding':record(prepared_path),
            'run_bindings':run_bindings,'script_binding':record(__file__),'phenotype_values_used':False,
            'interpretation':'High similarity is evaluated only where >=50% of the entire frozen reference is jointly compared at depth >=5. Below-threshold similarity is not itself evidence of replacement. Low coverage or missing comparison is unresolved. Repeated visit pairs are dependent; summaries are descriptive and have no P values.'}
    dest.mkdir(parents=True,exist_ok=False)
    write_tsv(dest/'ALL_PAIR_GENOME_EVIDENCE.tsv',evidence)
    write_tsv(dest/'WITHIN_INFANT_EVIDENCE.tsv',[r for r in evidence if r['pair_type']=='within_infant'])
    shared_columns=list(evidence[0])+[base+str(side) for side in [1,2] for base in ['related_within_samples','related_within_pairs','related_within_min_nominal_gap','related_within_max_nominal_gap']]
    write_tsv(dest/'shared_compatible_between_infant_pairs.tsv',shared,shared_columns)
    write_tsv(dest/'WITHIN_VS_BETWEEN_DESCRIPTIVES.tsv',descriptives)
    write_tsv(dest/'WITHIN_INFANT_SUMMARIES.tsv',infant_rows)
    write_tsv(dest/'GENOME_COMPARISON_ELIGIBILITY.tsv',genome_status)
    write_json(dest/'SUMMARY.json',result)
    text=['# 菌株比较描述性汇总','',f"50 份文库、30 位婴儿；预先固定 25 个同婴儿访问对（16 位婴儿），其中 20 个相邻观测对；另有 1,200 个跨婴儿访问对。",'',
          '每个样本先要求参考基因组至少 50% 的位置达到深度 5；样本对还必须有至少整个参考基因组 50% 的位置共同可比较。只有满足这一资格后，popANI ≥0.99999 才记为高度相似。',
          '',f"本次状态：{result['status']}。共保留 {len(evidence):,} 个访问对×参考基因组条目。",'',
          '|参考组|访问对组|可比较/全部|高度相似|有资格但低于阈值|未决|',
          '|---|---|---:|---:|---:|---:|']
    for r in descriptives:
        text.append(f"|{p['genomes'][r['genome']]['taxon_group']}|{r['pair_group']}|{r['jointly_evaluable_pairs']}/{r['possible_visit_pairs']}|{r['high_similarity_pairs']}|{r['eligible_below_similarity_pairs']}|{r['unresolved_pairs']}|")
    text+=['','低相似性本身不等于菌株替换，低覆盖不等于缺失。高度相似可用于描述群体序列连续性，但不能单独证明单一菌株、定植途径或传递关系。访问对存在同一婴儿重复观测的依赖，因此本表不提供独立样本假设下的 P 值。营养指标未参与筛选或比较。',
           '',f"跨婴儿固定标准下有 {len(shared)} 个访问对×参考基因组的兼容命中，完整保留于 shared_compatible_between_infant_pairs.tsv；附 conANI、两端多样性和各端直接相关的同婴儿相似访问间隔。该表不建立传递链或克隆聚类。",
           '','完整条目见 ALL_PAIR_GENOME_EVIDENCE.tsv；同婴儿条目保留名义日龄和相邻观测标记，实际采样日间隔仍记为未知。']
    (dest/'REPORT_ZH.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    write_json(dest/'MANIFEST.json',{'status':result['status'],'files':[record(f) for f in sorted(dest.iterdir()) if f.is_file()]})
    return {'status':result['status'],'rows':len(evidence),'out':str(dest)}


def main():
    parser=argparse.ArgumentParser(description=__doc__);subs=parser.add_subparsers(dest='command',required=True)
    a=subs.add_parser('prepare');a.add_argument('--config',required=True)
    b=subs.add_parser('run');b.add_argument('--prepared',required=True);b.add_argument('--index',required=True,type=int)
    c=subs.add_parser('summarize');c.add_argument('--prepared',required=True);c.add_argument('--output',required=True)
    args=parser.parse_args()
    if args.command=='prepare':result=prepare(args.config)
    elif args.command=='run':result=run_genome(args.prepared,args.index)
    else:result=summarize(args.prepared,args.output)
    print(json.dumps(result,ensure_ascii=False))


if __name__=='__main__':main()
