import json
from pathlib import Path
import tempfile
import unittest

from glycan_evidence.structure_inference import extract_fragment_features
from glycan_evidence.structure_inference_io import run_structure_inference


class InterfaceTest(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.metadata={'A':dict(length=500,junctions=[100,175])}
        self.sam=self.root/'reads.sam'
        lines=['@SQ\tSN:A\tLN:500']
        for flag,pos,pnext,tlen in [(99,51,201,250),(147,201,51,-250)]:
            lines.append('\t'.join(map(str,['opaque',flag,'A',pos,0,'100M','=',pnext,tlen,
                'A'*100,'I'*100,'AS:i:0','YS:i:0','NM:i:0'])))
        self.sam.write_text('\n'.join(lines)+'\n')
        counts=extract_fragment_features(self.sam,self.metadata)['counts']
        for name,obj in [('candidates',self.metadata),('profiles',{'A':counts}),
                         ('config',{'fit':{'candidate_groups':{'G':['A']}}})]:
            (self.root/(name+'.json')).write_text(json.dumps(obj))

    def run_io(self):
        return run_structure_inference(self.sam,self.root/'candidates.json',self.root/'profiles.json',
            self.root/'config.json',self.root/'out')

    def test_complete_binds_inputs_and_outputs(self):
        receipt=self.run_io()
        self.assertEqual(receipt['status'],'COMPLETE')
        self.assertEqual(len(receipt['inputs']),4)
        result=json.loads((self.root/'out/fit.json').read_text())
        self.assertAlmostEqual(result['bestfit_weights']['A'],1)
        self.assertTrue((self.root/'out/COMPLETE').is_file())

    def test_existing_destination_refused(self):
        self.run_io()
        with self.assertRaises(FileExistsError):self.run_io()

    def test_different_candidate_universe_refused(self):
        (self.root/'profiles.json').write_text(json.dumps({'B':{'x':1}}))
        with self.assertRaisesRegex(ValueError,'IDs must match'):self.run_io()
        self.assertFalse((self.root/'out').exists())

    def test_malformed_sam_records_failure_without_complete(self):
        self.sam.write_text('@SQ\tSN:B\tLN:500\n')
        with self.assertRaises(ValueError):self.run_io()
        receipt=json.loads((self.root/'out/receipt.json').read_text())
        self.assertEqual(receipt['status'],'FAILED')
        self.assertFalse((self.root/'out/COMPLETE').exists())


if __name__=='__main__':unittest.main()
