import importlib.util
from pathlib import Path
import tempfile
import unittest
import os
import subprocess

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('compare_strains',HERE/'compare_strains.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def base_pair():
    return dict(zip(m.PAIR_COLUMNS,['S1','S2','1','1','D05','D14',5,14,'within_infant',True,9,False]))


def profiles(breadths=(.5,.5),retained=(10,10)):
    return {s:{'retained_pairs':n,'genomes':{'G':{'sample_eligible':b>=.5,'breadth_minCov':b}}}
            for s,b,n in zip(['S1','S2'],breadths,retained)}


def normalized(compared=100000,popani=.99999,snps=1,length=200000,reported=.5):
    return m.normalized_comparisons([{'genome':'G','name1':'S1.bam','name2':'S2.bam',
        'compared_bases_count':str(compared),'popANI':str(popani),'population_SNPs':str(snps),
        'percent_compared':str(reported),'coverage_overlap':'.8'}],'G',length,{'S1.bam':'S1','S2.bam':'S2'})[('S1','S2')]


class CompareEvidence(unittest.TestCase):
    def test_joint_half_and_ani_exact_boundary_inclusive(self):
        row=m.pair_evidence(base_pair(),'G',profiles(),normalized(),200000)
        self.assertTrue(row['jointly_evaluable']);self.assertTrue(row['high_similarity'])

    def test_perfect_ani_insufficient_whole_genome_not_evaluable(self):
        row=m.pair_evidence(base_pair(),'G',profiles(),normalized(compared=99999,popani=1,snps=0,reported=.499995),200000)
        self.assertEqual(row['comparison_status'],'JOINT_COVERAGE_INELIGIBLE')
        self.assertFalse(row['high_similarity'])

    def test_use_whole_genome_denominator_not_reported_coverage(self):
        row=normalized(compared=99000,popani=1,snps=0,reported=.99)
        self.assertEqual(row['joint_whole_genome_fraction'],.495)
        self.assertFalse(row['reported_denominator_matches_full_reference'])

    def test_lower_ani_not_called_replacement(self):
        row=m.pair_evidence(base_pair(),'G',profiles(),normalized(popani=.99998,snps=2),200000)
        self.assertEqual(row['comparison_status'],'ELIGIBLE_BELOW_SIMILARITY_THRESHOLD')
        self.assertTrue(row['jointly_evaluable']);self.assertFalse(row['high_similarity'])

    def test_missing_compare_retained_not_negative(self):
        row=m.pair_evidence(base_pair(),'G',profiles(),None,200000)
        self.assertEqual(row['comparison_status'],'EXPECTED_COMPARISON_MISSING')
        self.assertIsNone(row['popANI'])

    def test_breadth_failure_and_no_analysable_distinct(self):
        self.assertEqual(m.pair_evidence(base_pair(),'G',profiles((.49,.8)),None,200000)['comparison_status'],'SAMPLE_BREADTH_INELIGIBLE')
        self.assertEqual(m.pair_evidence(base_pair(),'G',profiles((0,.8),(0,10)),None,200000)['comparison_status'],'NO_ANALYSABLE_PAIRS')

    def test_ineligible_input_comparison_rejected(self):
        with self.assertRaises(ValueError):m.pair_evidence(base_pair(),'G',profiles((.1,.8)),normalized(),200000)

    def test_ani_consistency_checked(self):
        with self.assertRaises(ValueError):normalized(popani=1,snps=1)

    def test_unknown_name_and_duplicate_pair_rejected(self):
        row={'genome':'G','name1':'S1.bam','name2':'S2.bam','compared_bases_count':'10','popANI':'1','population_SNPs':'0'}
        with self.assertRaises(ValueError):m.normalized_comparisons([row],'G',100,{'S1.bam':'S1'})
        with self.assertRaises(ValueError):m.normalized_comparisons([row,row],'G',100,{'S1.bam':'S1','S2.bam':'S2'})

    def test_zero_overlap_nan_ani_explicit(self):
        row=normalized(compared=0,popani='nan',snps=0,reported=0)
        result=m.pair_evidence(base_pair(),'G',profiles(),row,200000)
        self.assertEqual(result['comparison_status'],'JOINT_COVERAGE_INELIGIBLE')
        self.assertIsNone(result['popANI'])

    def test_latest_mm_selected_not_average(self):
        rows=[{'genome':'G','mm':'1','breadth_minCov':'.1'},{'genome':'G','mm':'2','breadth_minCov':'.7'}]
        self.assertEqual(m.latest_genome_rows(rows)['G']['breadth_minCov'],'.7')
        with self.assertRaises(ValueError):m.latest_genome_rows(rows+rows)

    def test_metadata_only_and_consecutive_not_allpairs(self):
        samples=[{'sample_id':'A','dyad_id':1,'window':'D05','nominal_day':5,'FFA':999},
                 {'sample_id':'B','dyad_id':1,'window':'D14','nominal_day':14,'FFA':999},
                 {'sample_id':'C','dyad_id':1,'window':'D90','nominal_day':90,'FFA':999},
                 {'sample_id':'D','dyad_id':2,'window':'D05','nominal_day':5,'FFA':999}]
        clean,pairs,shape=m.metadata(samples,enforce_cohort=False)
        self.assertEqual(shape['within_pairs'],3);self.assertEqual(shape['consecutive_within_pairs'],2)
        self.assertEqual(shape['between_pairs'],3)
        self.assertNotIn('FFA',str(clean));self.assertNotIn('FFA',str(pairs))

    def test_descriptives_keep_unresolved_denominator(self):
        one=m.pair_evidence(base_pair(),'G',profiles(),normalized(),200000)
        two=m.pair_evidence(base_pair(),'G',profiles(),None,200000)
        r=m.describe([one,two]);self.assertEqual(r['possible_visit_pairs'],2)
        self.assertEqual(r['jointly_evaluable_pairs'],1);self.assertEqual(r['unresolved_pairs'],1)
        self.assertEqual(r['high_similarity_fraction_among_evaluable'],1)
        self.assertEqual(r['evaluable_infants'],1)

    def test_conani_retained_without_changing_popani_gate(self):
        raw={'genome':'G','name1':'S1.bam','name2':'S2.bam','compared_bases_count':'100000',
             'popANI':'0.99999','population_SNPs':'1','conANI':'0.999','consensus_SNPs':'100'}
        row=m.normalized_comparisons([raw],'G',200000,{'S1.bam':'S1','S2.bam':'S2'})[('S1','S2')]
        final=m.pair_evidence(base_pair(),'G',profiles(),row,200000)
        self.assertTrue(final['high_similarity']);self.assertEqual(final['conANI'],.999)

    def test_cross_infant_shared_context_is_direct_only(self):
        within=m.pair_evidence(base_pair(),'G',profiles(),normalized(),200000)
        cross=dict(within,sample2='S3',dyad2='2',pair_type='between_infants')
        unrelated=dict(within,sample1='S2',sample2='S4',nominal_day_gap=30)
        shared=m.between_shared_context([within,cross,unrelated])
        self.assertEqual(len(shared),1)
        self.assertEqual(shared[0]['related_within_samples1'],'S2')
        self.assertEqual(shared[0]['related_within_min_nominal_gap1'],9)
        self.assertNotIn('S4',shared[0]['related_within_samples1'])
        self.assertIsNone(shared[0]['related_within_min_nominal_gap2'])

    def test_command_matches_frozen_cli(self):
        argv=m.compare_argv('/bin/inStrain','list','out','stb','G',16)
        for option,value in [('--breadth','0.5'),('-c','5'),('-f','0.05'),('-fdr','1e-6'),('-cov','0.5'),('-ani','0.99999'),('--genome','G')]:
            self.assertEqual(argv[argv.index(option)+1],value)
        self.assertIn('--database_mode',argv);self.assertIn('--store_coverage_overlap',argv)
        self.assertNotIn('--include_self_comparisons',argv)

if __name__=='__main__':unittest.main(verbosity=2)
