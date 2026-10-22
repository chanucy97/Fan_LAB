import json
from pathlib import Path
import tempfile
import unittest

from glycan_evidence.structure_conditioning import extract_locus_conditioned
from glycan_evidence.structure_inference import extract_fragment_features


def pair(name="q", ref="A", start=20, mate_start=70, length=20, score=0,
         secondary=False, cigar=None, nm=0, span=None):
    cigar = cigar or f"{length}M"
    span = mate_start+length-start if span is None else span
    extra = 256 if secondary else 0
    return ["\t".join(map(str, [name, flag+extra, ref, pos+1, 0, cigar, "=", other+1,
        tlen, "A"*length, "I"*length, f"AS:i:{score}", f"YS:i:{score}", f"NM:i:{nm}"]))
        for flag, pos, other, tlen in ((99,start,mate_start,span),(147,mate_start,start,-span))]


class ConditioningTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.metadata={"A":{"length":1000,"junctions":[200,700],"class_id":"a"},
                       "B":{"length":1000,"junctions":[200,700],"class_id":"b"}}

    def run_case(self, lines, config=None):
        path=Path(self.tmp.name)/"x.sam"
        path.write_text("\n".join(["@HD\tVN:1.6",*[f"@SQ\tSN:{c}\tLN:{m['length']}" for c,m in self.metadata.items()],*lines])+"\n")
        original=extract_fragment_features(path,self.metadata,config)
        result=extract_locus_conditioned(path,self.metadata,config,expected_full_counts=original["counts"])
        self.assertEqual(result['original_candidate_coverage'],original['candidate_coverage'])
        self.assertEqual(result['full_counts'],original['counts'])
        json.dumps(result,allow_nan=False)
        return result,original,path

    def test_anchor_only_removed_and_partitioned(self):
        r,old,_=self.run_case(pair())
        self.assertEqual(r['counts'],{})
        self.assertEqual(r['anchor_only_counts'],old['counts'])
        self.assertEqual(r['selection_audit']['anchor_only_fragments'],1)

    def test_internal_contiguous_block_at_least_min_anchor(self):
        r,_,_=self.run_case(pair(start=185,mate_start=300,length=30))
        self.assertEqual(sum(r['counts'].values()),1)
        self.assertTrue(r['fragment_ledger'][0]['internal_evidence'])

    def test_exact_internal_overlap_threshold(self):
        # Both mates end inside the left boundary; only the second enters it.
        yes,_,_=self.run_case(pair(start=100,mate_start=195,length=20))
        no,_,_=self.run_case(pair(start=100,mate_start=194,length=20))
        self.assertEqual(sum(yes['counts'].values()),1)
        self.assertEqual(no['counts'],{})

    def test_direct_boundary_crossing_qualifies(self):
        r,_,_=self.run_case(pair(start=180,mate_start=300,length=40))
        self.assertTrue(r['fragment_ledger'][0]['boundary_evidence'])
        self.assertIn(200,r['fragment_ledger'][0]['candidate_evidence']['A']['direct_junctions'])

    def test_pure_anchor_pair_bridge_qualifies_without_internal_bases(self):
        r,_,_=self.run_case(pair(start=150,mate_start=750,length=30))
        item=r['fragment_ledger'][0]
        self.assertTrue(item['boundary_evidence']);self.assertFalse(item['internal_evidence'])
        self.assertEqual(r['selection_audit']['qualification_reasons'],{'boundary_only':1})
        self.assertEqual(r['candidate_coverage']['A']['eligible_inner_breadth'],0)

    def test_any_near_candidate_keeps_original_ambiguity(self):
        r,old,_=self.run_case(pair()+pair(ref='B',start=250,mate_start=300,secondary=True))
        self.assertEqual(r['counts'],old['counts'])
        self.assertEqual(r['fragment_ledger'][0]['candidates'],['A','B'])
        self.assertEqual(r['fragment_ledger'][0]['qualifying_candidates'],['B'])
        self.assertEqual(r['candidate_coverage']['A']['near_best_inner_breadth'],0)
        self.assertGreater(r['candidate_coverage']['B']['near_best_inner_breadth'],0)

    def test_eligible_but_not_near_cannot_select(self):
        r,_,_=self.run_case(pair()+pair(ref='B',start=250,mate_start=300,score=-4,secondary=True))
        self.assertEqual(r['counts'],{})
        self.assertGreater(r['original_candidate_coverage']['B']['eligible_base_breadth'],0)

    def test_near_score_boundary_inclusive(self):
        r,_,_=self.run_case(pair()+pair(ref='B',start=250,mate_start=300,score=-3,secondary=True))
        self.assertEqual(sum(r['counts'].values()),1)

    def test_same_candidate_multimap_any_selection_preserves_intersection(self):
        r,old,_=self.run_case(pair(length=40)+pair(start=180,mate_start=300,length=40,secondary=True))
        self.assertEqual(r['counts'],old['counts'])
        item=r['fragment_ledger'][0]
        self.assertTrue(item['boundary_evidence'])
        self.assertEqual(item['candidate_evidence']['A']['direct_junctions'],[])

    def test_one_qname_not_double_counted_with_two_qualifiers(self):
        r,_,_=self.run_case(pair(start=250,mate_start=300)+pair(ref='B',start=400,mate_start=450,secondary=True))
        self.assertEqual(sum(r['counts'].values()),1)
        self.assertEqual(r['audit']['technical_summary']['mate_query_lengths']['count'],2)

    def test_gapped_short_blocks_not_joined_for_internal_threshold(self):
        # 10M1I10M blocks enter the interval, each below 15; no boundary bridge.
        r,_,_=self.run_case(pair(start=250,mate_start=300,length=21,cigar='10M1I10M',nm=1,span=70))
        self.assertEqual(r['counts'],{})

    def test_coverage_inner_is_union_of_clipped_blocks(self):
        r,_,_=self.run_case(pair(start=195,mate_start=230,length=20))
        self.assertAlmostEqual(r['candidate_coverage']['A']['eligible_inner_breadth'],35/500)
        self.assertAlmostEqual(r['candidate_coverage']['A']['near_best_inner_breadth'],35/500)

    def test_unexplained_separate_from_selected_and_anchor_only(self):
        r,old,_=self.run_case(pair(name='a')+pair(name='b',start=250,mate_start=300)+pair(name='bad',nm=4))
        self.assertEqual(r['unexplained_fragments'],1)
        self.assertEqual(sum(r['counts'].values()),1)
        self.assertEqual(sum(r['anchor_only_counts'].values()),1)
        self.assertNotIn('__unexplained__',r['counts'])
        self.assertEqual(sum(r['full_counts'].values()),3)

    def test_expected_full_counts_mismatch_fails(self):
        _,_,p=self.run_case(pair())
        with self.assertRaisesRegex(ValueError,'historical'):
            extract_locus_conditioned(p,self.metadata,expected_full_counts={})

    def test_empty_header_only_sam_is_explicit_empty(self):
        r,_,_=self.run_case([])
        self.assertEqual(r['counts'],{})
        self.assertEqual(r['selection_audit']['selected_fragments'],0)

    def test_ledger_cap_does_not_change_conditioned_counts(self):
        r,_,_=self.run_case(pair(name='q1',start=250,mate_start=300)+pair(name='q2',start=300,mate_start=350),{'max_ledger':1})
        self.assertEqual(sum(r['counts'].values()),2)
        self.assertEqual(len(r['fragment_ledger']),1)
        self.assertEqual(r['audit']['ledger_rows_omitted_by_cap'],1)


if __name__=='__main__':unittest.main()
