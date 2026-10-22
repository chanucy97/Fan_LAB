import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from glycan_evidence.structure_recovery_io import record,run_structure_recovery
from glycan_evidence.cli import main
from glycan_evidence.integrity import verify_outputs
from test_structure_recovery import fixture


class RecoveryIOTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        source=self.root/'source.txt';source.write_text('fixed source\n')
        self.binding=record(source)
        v2=self.root/'v2.json'
        v2.write_text(json.dumps(dict(case_id='sample',status='COMPLETE',decision_v2=dict(supported_classes=['short_824']))))
        panel,candidate=fixture()
        candidate['source_bindings']={k:copy.deepcopy(self.binding) for k in candidate['source_bindings']}
        self.data=dict(schema='truth-free-existing-assembly-roster.v1',source_inventory=self.binding,
             panel_sources=[self.binding],panel=panel,cases=[dict(case_id='sample',reference_calls=['short_824'],
             v2_source=record(v2),assembly_candidates=[candidate])])
        self.roster=self.root/'roster.json';self.save()

    def save(self):
        self.roster.write_text(json.dumps(self.data),encoding='utf-8')
        self.sha=record(self.roster)['sha256']

    def run_case(self):
        return run_structure_recovery(self.roster,'sample',self.sha,self.root/'out')

    def test_actual_cli_and_integrity(self):
        self.assertEqual(main(['structure-recover','--roster',str(self.roster),'--case-id','sample',
                  '--expected-roster-sha256',self.sha,'--out',str(self.root/'out')]),0)
        self.assertEqual(verify_outputs(self.root/'out')['status'],'PASS')
        result=json.loads((self.root/'out/recovery.json').read_text())
        self.assertEqual(result['concordant_classes'],['short_824'])

    def test_roster_digest_mismatch(self):
        self.roster.write_text(self.roster.read_text()+' ')
        with self.assertRaises(ValueError):self.run_case()
        self.assertFalse((self.root/'out').exists())

    def test_bound_file_change_is_rejected(self):
        (self.root/'source.txt').write_text('changed source\n')
        with self.assertRaises(ValueError):self.run_case()
        self.assertFalse((self.root/'out').exists())

    def test_duplicate_case_is_rejected(self):
        self.data['cases'].append(copy.deepcopy(self.data['cases'][0]));self.save()
        with self.assertRaises(ValueError):self.run_case()

    def test_cross_library_candidate_is_rejected(self):
        self.data['cases'][0]['assembly_candidates'][0]['case_id']='different';self.save()
        with self.assertRaises(ValueError):self.run_case()

    def test_wrong_reference_source_case_is_rejected(self):
        path=self.root/'v2.json'
        path.write_text(json.dumps(dict(case_id='different',status='COMPLETE',decision_v2=dict(supported_classes=['short_824']))))
        self.data['cases'][0]['v2_source']=record(path);self.save()
        with self.assertRaises(ValueError):self.run_case()

    def test_declared_calls_must_equal_source(self):
        self.data['cases'][0]['reference_calls']=[];self.save()
        with self.assertRaises(ValueError):self.run_case()

    def test_existing_output_preserved(self):
        self.run_case();old=(self.root/'out/recovery.json').read_bytes()
        with self.assertRaises(FileExistsError):self.run_case()
        self.assertEqual((self.root/'out/recovery.json').read_bytes(),old)

    def test_partial_complete_marker_removed_on_write_failure(self):
        original=Path.write_text
        def fail_marker(path,text,*args,**kwargs):
            if path.name=='COMPLETE':
                original(path,'partial',*args,**kwargs)
                raise OSError('simulated marker write failure')
            return original(path,text,*args,**kwargs)
        with patch.object(Path,'write_text',fail_marker):
            with self.assertRaises(OSError):self.run_case()
        self.assertFalse((self.root/'out/COMPLETE').exists())
        self.assertEqual(json.loads((self.root/'out/receipt.json').read_text())['status'],'FAILED')


if __name__=='__main__':unittest.main()
