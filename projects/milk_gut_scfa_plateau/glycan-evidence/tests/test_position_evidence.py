import unittest
from glycan_evidence.position_evidence import assess_position


def row(assembly='A', values=None):
    result={'assembled_genomic_base':assembly}
    for base in 'ACGT':
        f,r,b=(values or {}).get(base,(0,0,0))
        for key,value in zip(('F','R','both','total'),(f,r,b,f+r+b)):
            result[base+'_qname_'+key]=value
    result['acgt_qname_total']=sum(result[b+'_qname_total'] for b in 'ACGT')
    return result


class PositionEvidenceTests(unittest.TestCase):
    def test_tie_is_reported_without_inventing_correction(self):
        out=assess_position(row(values={'A':(1,0,0),'G':(0,1,0)}))
        self.assertTrue(out['assembly_tied_for_top'])
        self.assertFalse(out['nonassembly_exceeds_assembly'])
        self.assertFalse(out['strong_nonassembly_discordance'])
        self.assertEqual(out['assembly_fraction'],.5)

    def test_both_orientations_count_once(self):
        out=assess_position(row(values={'A':(1,0,0),'G':(1,1,3)}))
        self.assertEqual(out['acgt_qnames'],6)
        self.assertEqual(out['strong_nonassembly_bases'],['G'])

    def test_single_orientation_is_not_strong(self):
        out=assess_position(row(values={'G':(10,0,0)}))
        self.assertTrue(out['nonassembly_exceeds_assembly'])
        self.assertFalse(out['strong_nonassembly_discordance'])

    def test_empty_is_no_evidence(self):
        out=assess_position(row())
        self.assertTrue(out['no_acgt_evidence'])
        self.assertIsNone(out['assembly_fraction'])
        self.assertEqual(out['top_bases'],[])

    def test_unknown_assembly_base_is_not_confident_discordance(self):
        out=assess_position(row(assembly='N',values={'G':(5,5,0)}))
        self.assertTrue(out['assembly_base_unknown'])
        self.assertIsNone(out['assembly_fraction'])
        self.assertFalse(out['strong_nonassembly_discordance'])
        self.assertFalse(out['nonassembly_exceeds_assembly'])

    def test_bad_denominator_rejected(self):
        item=row(values={'A':(2,3,1)});item['acgt_qname_total']=7
        with self.assertRaises(ValueError):assess_position(item)

    def test_overlapping_orientation_tallies_rejected(self):
        item=row(values={'A':(2,3,1)});item['A_qname_total']=5
        with self.assertRaises(ValueError):assess_position(item)


if __name__=='__main__':unittest.main()
