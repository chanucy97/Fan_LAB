"""Isolated synthetic evidence checks. Never execute production alignments."""
import copy
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import core

HERE=Path(__file__).resolve().parent
CARRIER=HERE.parents[1]/'carrier_context'


class ContextTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dep=core.dependency(CARRIER/'carrier_context.py')
        cls.config=core.read(CARRIER/'config.json')

    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)

    def gene(self,ident='g1',start=11000,end=12000,trait='Bga42A_Binf_ATCC15697',strand=1):
        return {'gene_id':ident,'trait_id':trait,'trait_role':'positive_enzyme_seed' if trait in core.POSITIVE else 'specificity_comparator',
                'mapping_status':'high_similarity_candidate','start0':start,'end0':end,'strand':strand,
                'cds_sha256':core.seqsha(('A' if strand==1 else 'T')*(end-start)),
                'protein_sha256':'a'*64,'partial':'00','complete_predicted_orf':True}

    def case(self,n=30000,genes=None,masks=None):
        genes=genes or [self.gene()]
        return {'case_id':'synthetic__contig','sample_id':'synthetic','contig':'ctg','contig_length':n,
                'contig_sha256':core.seqsha('A'*n),'genes':genes,'mask_orfs':masks if masks is not None else copy.deepcopy(genes)}

    def manifest(self,case):
        return {'schema':'enzyme-context.inputs.v1','cohort_denominator':50,'candidate_count':len(case['genes']),
                'contig_count':1,'mask_orf_count':len(case['mask_orfs']),'dependency':core.rec(CARRIER/'carrier_context.py'),
                'config_sha256':core.CONFIG_SHA,'cases':[case],'references':[{'accession':'A'},{'accession':'B'}],
                'subjects':[{'subject_id':k,'accession':k,'length':40000,'sequence_sha256':'a'*64} for k in ('A','B')],
                'reference_orthologs':[],'source_receipts':[]}

    def hit(self,case,start,end,subject='A',subject_start=None,reverse=False,row=1):
        length=end-start;subject_start=start if subject_start is None else subject_start
        return {'qseqid':'ctg','sseqid':subject,'pident':100.,'length':length,'qlen':case['contig_length'],
                'qstart':start+1,'qend':end,'sstart':subject_start+length if reverse else subject_start+1,
                'send':subject_start+1 if reverse else subject_start+length,'evalue':1e-30,'bitscore':length*2.,
                'gaps':0,'gapopen':0,'btop':str(length),'hsp_row_number':row}

    def test_multiple_masks_merge_and_subtract(self):
        self.assertEqual(core.merge_intervals([(2,5),(4,8),(8,10),(14,17)]),[[2,10],[14,17]])
        self.assertEqual(core.subtract([[0,20]],[[2,10],[14,17]]),[[0,2],[10,14],[17,20]])
        other=self.gene('comparator',5000,6000,'Bga42B_Binf_ATCC15697')
        case=self.case(masks=[self.gene(),other]);masks,regions=core.region_set(case,case['genes'][0],self.config)
        self.assertEqual(masks,[[4000,7000],[10000,13000]])
        left=next(r for r in regions if r['region']=='left_flank')
        self.assertEqual(left['segments'],[[0,4000],[7000,10000]])
        self.assertEqual(left['length'],7000);self.assertTrue(left['mask_interrupted'])
        self.assertEqual(next(r for r in regions if r['region']=='masked_background')['length'],24000)

    def test_background_excludes_all_targets_and_buffers(self):
        genes=[self.gene(),self.gene('g2',17000,18000,'Bga2A_Binf_ATCC15697')]
        case=self.case(genes=genes);m=self.manifest(case)
        hits=[self.hit(case,10000,13000,row=1),self.hit(case,16000,19000,row=2),self.hit(case,0,10000,'B',row=3),
              self.hit(case,13000,16000,'B',row=4),self.hit(case,19000,30000,'B',row=5)]
        scores,_,_,_,summaries=core.score(case,hits,m,self.config,self.dep)
        for gene in genes:
            a=next(r for r in scores if r['gene_id']==gene['gene_id'] and r['reference']=='A' and r['region']=='masked_background')
            self.assertEqual(a['covered_bases'],0)
            b=next(r for r in scores if r['gene_id']==gene['gene_id'] and r['reference']=='B' and r['region']=='masked_background')
            self.assertEqual(b['whole_region_matching_fraction'],1)
        self.assertEqual(summaries[0]['region_rankings']['local_interval']['best_reference'],'A')
        self.assertEqual(summaries[0]['region_rankings']['masked_background']['best_reference'],'B')

    def test_segment_chains_never_cross_other_target_mask(self):
        other=self.gene('comp',5000,6000,'Bga42B_Binf_ATCC15697')
        case=self.case(masks=[self.gene(),other]);m=self.manifest(case)
        scores,chains,_,_,summary=core.score(case,[self.hit(case,0,30000)],m,self.config,self.dep)
        ch=chains['g1|A|left_flank'];self.assertEqual(len(ch['segments']),2)
        self.assertFalse(ch['cross_mask_chain_asserted'])
        row=next(r for r in scores if r['gene_id']=='g1' and r['reference']=='A' and r['region']=='left_flank' and r['method']=='independent_segment_chains')
        self.assertEqual(row['covered_bases'],7000);self.assertFalse(row['single_chain_asserted'])
        self.assertEqual(summary[0]['region_rankings']['left_flank']['region_status'],'masked_flank_unresolved')
        self.assertFalse(summary[0]['region_rankings']['left_flank']['unique_panel_affinity'])
        for segment in ch['segments']:
            for block in segment['chain']['blocks']:
                self.assertEqual(block['hsp_row_numbers'],[1])
                self.assertTrue(block['qend0']<=4000 or block['qstart0']>=7000)

    def test_negative_coordinate_reference_pointers(self):
        gene=self.gene(strand=-1);case=self.case(genes=[gene]);m=self.manifest(case)
        m['reference_orthologs']=[{'subject_id':'A','start0':18000,'end0':19000,'gene_id':'ref_g','trait_id':gene['trait_id']}]
        _,_,pointers,_,_=core.score(case,[self.hit(case,0,30000,reverse=True)],m,self.config,self.dep)
        pointer=pointers[0]
        self.assertEqual((pointer['subject_start0'],pointer['subject_end0']),(18000,19000))
        self.assertEqual(pointer['relative_alignment_strand'],-1);self.assertEqual(pointer['subject_coding_strand'],1)
        self.assertEqual(pointer['overlapping_same_trait_reference_homologue_ids'],'ref_g')
        self.assertFalse(pointer['orthology_asserted'])

    def test_empty_and_short_flanks_remain_unqualified(self):
        case=self.case(5000,[self.gene(start=10,end=100)]);m=self.manifest(case)
        scores,_,_,_,summaries=core.score(case,[self.hit(case,0,5000)],m,self.config,self.dep)
        ranked=summaries[0]['region_rankings']
        self.assertEqual(ranked['left_flank']['region_length'],0)
        self.assertEqual(ranked['right_flank']['region_length'],3900)
        self.assertFalse(ranked['right_flank']['coverage_qualified'])
        self.assertEqual(ranked['right_flank']['region_status'],'insufficient_flank_length')

    def test_zero_background_retains_all_references(self):
        case=self.case(1000,[self.gene(start=100,end=900)]);m=self.manifest(case)
        scores,_,_,_,summaries=core.score(case,[],m,self.config,self.dep)
        self.assertEqual(len(scores),14)
        self.assertIsNone(summaries[0]['region_rankings']['masked_background']['best_score'])
        self.assertEqual(summaries[0]['region_rankings']['masked_background']['near_best_references'],['A','B'])

    def test_tie_and_exact_margin_boundary(self):
        base={'region_length':10000,'requested_length':10000,'full_requested_length_available':True,'mask_interrupted':False,
              'coverage':1.,'covered_identity':1.,'possible_circular_origin':False,'near_equal_other_locus':False}
        rows=[{**base,'reference':'A','whole_region_matching_fraction':1.},{**base,'reference':'B','whole_region_matching_fraction':.99}]
        result=core.rank(rows,'masked_background',self.config)
        self.assertEqual(result['near_best_references'],['A','B']);self.assertFalse(result['unique_panel_affinity'])
        rows[1]['whole_region_matching_fraction']=.98
        result=core.rank(rows,'masked_background',self.config)
        self.assertTrue(result['unique_panel_affinity']);self.assertEqual(result['near_best_references'],['A'])

    def test_duplicate_contig_and_unmasked_gene_rejected(self):
        case=self.case();m=self.manifest(case);core.validate_manifest(m)
        dupe=copy.deepcopy(case);dupe['case_id']='different_name';m['cases'].append(dupe);m.update(contig_count=2,candidate_count=2,mask_orf_count=2)
        with self.assertRaisesRegex(ValueError,'duplicate sample/contig'):core.validate_manifest(m)
        m=self.manifest(case);m['cases'][0]['mask_orfs']=[];m['mask_orf_count']=0
        with self.assertRaisesRegex(ValueError,'missing from masks'):core.validate_manifest(m)

    def test_btop_malformed_and_reverse_coordinates(self):
        case=self.case();m=self.manifest(case);h=self.hit(case,0,100,reverse=True)
        p=self.root/'hits.tsv';p.write_text('\t'.join(str(h[k]) for k in self.dep.FIELDS)+'\n')
        self.assertEqual(core.parse_hits(p,case,m,self.dep)[0]['hsp_row_number'],1)
        h['btop']='99';p.write_text('\t'.join(str(h[k]) for k in self.dep.FIELDS)+'\n')
        with self.assertRaisesRegex(ValueError,'BTOP coordinate'):core.parse_hits(p,case,m,self.dep)

    def test_case_sequence_negative_strand_hash(self):
        case=self.case(1000,[self.gene(start=100,end=900,strand=-1)])
        p=self.root/'assembly.fa';p.write_text('>ctg\n'+'A'*1000+'\n');case['assembly']=core.rec(p)
        self.assertEqual(len(core.selected_contig(case,self.dep)),1000)
        case['mask_orfs'][0]['cds_sha256']='0'*64
        with self.assertRaisesRegex(ValueError,'CDS source'):core.selected_contig(case,self.dep)

    def test_fresh_outputs_and_mutated_receipts_fail(self):
        out=core.fresh(self.root/'out');(out/'data').write_text('synthetic')
        core.finish(out,'RUN.json',{'schema':'synthetic'})
        self.assertEqual(core.bundle(out,'RUN.json')['status'],'COMPLETE')
        with self.assertRaises(FileExistsError):core.fresh(out)
        (out/'data').write_text('changed')
        with self.assertRaisesRegex(ValueError,'hash'):core.bundle(out,'RUN.json')

    def test_public_path_tools_are_bound_and_missing_tool_fails(self):
        inputs=self.root/'inputs.json';inputs.write_text('{}')
        tool=self.root/'synthetic_tool';tool.write_text('synthetic binary content')
        m={'dependency':core.rec(CARRIER/'carrier_context.py')}
        with patch.object(core.shutil,'which',return_value=str(tool)):
            bound=core.freeze_bindings(inputs,CARRIER/'config.json',m,self.config)
        self.assertEqual(bound['tools']['blastn'],core.rec(tool))
        self.assertEqual(bound['tools']['makeblastdb'],core.rec(tool))
        missing={**self.config,'blastn':str(self.root/'missing_tool')}
        with patch.object(core.shutil,'which',return_value=None):
            with self.assertRaises(FileNotFoundError):
                core.freeze_bindings(inputs,CARRIER/'config.json',m,missing)

    def test_one_contig_one_command_and_original_profile_untouched(self):
        genes=[self.gene(),self.gene('g2',17000,18000,'Bga2A_Binf_ATCC15697')]
        case=self.case(genes=genes);assembly=self.root/'assembly.fa';assembly.write_text('>ctg\n'+'A'*30000+'\n');case['assembly']=core.rec(assembly)
        m=self.manifest(case);reference=self.root/'ref.fa';reference.write_text('>A\nAAAA\n>B\nAAAA\n');m['reference_fasta']=core.rec(reference)
        inputs=self.root/'INPUTS.json';core.write(inputs,m);config_path=CARRIER/'config.json'
        sentinel=self.root/'original_profile';sentinel.mkdir();(sentinel/'RUN.json').write_text('{"status":"COMPLETE"}')
        sentinel_before=core.rec(sentinel/'RUN.json')
        tool=self.root/'synthetic_tool';tool.write_text('not executable: test double')
        bindings={'inputs':core.rec(inputs),'config':core.rec(config_path),'core':core.rec(core.__file__),'dependency':m['dependency'],
                  'tools':{'blastn':core.rec(tool),'makeblastdb':core.rec(tool)}}
        prepared=core.fresh(self.root/'prepared');(prepared/'reference.nsq').write_text('synthetic database')
        core.finish(prepared,'PREPARED.json',{'schema':'enzyme-context.prepared.v1','bindings':bindings,'reference_fasta':m['reference_fasta'],
                    'source_receipts':[],'parameters':self.config,'case_count':1,'candidate_count':2,'database_prefix':str(prepared/'reference')})
        calls=[]
        def fake_command(argv,out,commands):
            calls.append(argv);self.assertIn('-query',argv)
            self.assertEqual(argv[0],bindings['tools']['blastn']['path'])
            h=self.hit(case,0,30000);path=Path(argv[argv.index('-out')+1]);path.write_text('\t'.join(str(h[k]) for k in self.dep.FIELDS)+'\n')
            commands.append({'argv':[str(x) for x in argv],'returncode':0,'synthetic':True});core.write(out/'COMMANDS.json',commands)
        with patch.object(core,'command',side_effect=fake_command):
            result=core.execute(inputs,config_path,prepared/'PREPARED.json',self.root/'results',0)
        self.assertEqual(len(calls),1);self.assertEqual(result['gene_count'],2)
        collected=core.collect(inputs,config_path,self.root/'results',self.root/'collection')
        self.assertEqual(collected['genes'],2);self.assertEqual(collected['score_rows'],28)
        self.assertEqual(core.rec(sentinel/'RUN.json'),sentinel_before)
        self.assertEqual(core.rec(assembly),case['assembly'])
        # A forged preparation status cannot pass collection even with a refreshed outer receipt.
        folder=self.root/'results'/case['case_id'];snapshot=folder/'PREPARED_SNAPSHOT.json';p=core.read(snapshot);p['status']='FAILED';core.write(snapshot,p)
        run=core.read(folder/'RUN.json');run['prepared'].update(bytes=snapshot.stat().st_size,sha256=core.sha(snapshot))
        for item in run['outputs']:
            if item['path']==snapshot.name:item.update(bytes=snapshot.stat().st_size,sha256=core.sha(snapshot))
        core.write(folder/'RUN.json',run);(folder/'COMPLETE').write_text(core.sha(folder/'RUN.json')+'\n')
        with self.assertRaisesRegex(ValueError,'wrong preparation'):core.collect(inputs,config_path,self.root/'results',self.root/'bad_collection')


if __name__=='__main__':unittest.main()
