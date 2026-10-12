"""Input/output integrity of the new bridge; synthetic fixtures, no performance claims."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import unittest
from unittest.mock import patch

from glycan_evidence import compatibility_io as base, local_support_io as bridge
from glycan_evidence.integrity import digest
import test_compatibility_io


class LocalSupportBridgeTests(unittest.TestCase):
    def setUp(self):
        self.fixture = test_compatibility_io.CompatibilityIOTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.doCleanups)
        self.bundle = self.fixture.root/'local_bundle'
        self.bundle.mkdir()
        shutil.copytree(self.fixture.bundle, self.bundle/'baseline')
        self.fixture.save(self.bundle/'CANDIDATE_MODEL.json', {'Z_indices':[0,1]})
        self.fixture.save(self.bundle/'LOCAL_REFERENCE_ORDER.json', ['synthetic_source'])
        for name in ['local_candidate_api.py','local_support_api.py']:
            (self.bundle/name).write_text('synthetic unexecuted fixture\n')
        manifest = {'version':bridge.VERSION, 'baseline_manifest_sha256':self.fixture.expected,
                    'files_sha256':{p.relative_to(self.bundle).as_posix():digest(p)
                                    for p in self.bundle.rglob('*') if p.is_file()}}
        self.fixture.save(self.bundle/'LOCAL_MANIFEST.json',manifest)
        self.expected = digest(self.bundle/'LOCAL_MANIFEST.json')
        self.out = self.fixture.root/'local_result'
        self.fixture.reseal_run()

    def fake(self, argv, **kwargs):
        self.assertEqual(argv[1:5],['-I','-B','-X','utf8'])
        snapshot = Path(argv[argv.index('--bundle')+1])
        self.assertNotEqual(snapshot,self.bundle)
        self.assertFalse((snapshot/'numpy.py').exists())
        query = self.fixture.load(Path(argv[argv.index('--query')+1]))
        nulls = dict(identity_call=None,p_value=None,confidence_interval=None,organism_abundance=None)
        result = dict(version=bridge.VERSION,candidate_manifest_sha256=self.expected,
                      baseline_manifest_sha256=self.fixture.expected,
                      input_binding_sha256=hashlib.sha256(base._canonical(query).encode()).hexdigest(),
                      input_pairs=query['input_pairs'],composite_support=sum(query['counts']),
                      local_support=sum(query['counts'][:2]),
                      observation=dict(state='insufficient_evidence',reasons=['synthetic_fixture'],**nulls),**nulls)
        self.fixture.save(Path(argv[argv.index('--out')+1]),result)
        return subprocess.CompletedProcess(argv,0)

    def run_bridge(self, callback=None):
        with patch.object(bridge.subprocess,'run',side_effect=callback or self.fake):
            return bridge.run_local_support(self.fixture.metrics,self.bundle,self.expected,self.out)

    def test_complete_measured_query_and_abstention_preserved(self):
        result = self.run_bridge()
        query = self.fixture.load(self.out/'query.json')
        self.assertEqual(query['input_pairs'],50)
        self.assertEqual(len(query['counts']),220)
        self.assertEqual(query['read_sha256'],['a'*64,'b'*64])
        self.assertEqual(result['state'],'insufficient_evidence')
        self.assertEqual(set(result['helper_code']),{'compatibility_io','identity','integrity'})
        self.assertEqual((self.out/'COMPLETE').read_text().strip(),digest(self.out/'receipt.json'))

    def test_inconsistent_N_rejected_before_execution(self):
        self.fixture.change_adapter(lambda x:x['raw_input_validation'].update(pairs_n=49))
        with patch.object(bridge.subprocess,'run') as run:
            with self.assertRaisesRegex(ValueError,'differs from the no-QC'):
                bridge.run_local_support(self.fixture.metrics,self.bundle,self.expected,self.out)
            run.assert_not_called()
        self.assertFalse(self.out.exists())

    def test_changed_candidate_code_rejected_before_execution(self):
        (self.bundle/'local_support_api.py').write_text('changed')
        with patch.object(bridge.subprocess,'run') as run:
            with self.assertRaisesRegex(ValueError,'Candidate payload changed'):
                bridge.run_local_support(self.fixture.metrics,self.bundle,self.expected,self.out)
            run.assert_not_called()

    def test_post_QC_counts_do_not_borrow_raw_N(self):
        path = self.fixture.folder/'run_manifest.json'
        run = self.fixture.load(path); run['skip_qc'] = False
        self.fixture.save(path,run); self.fixture.reseal_run()
        with self.assertRaisesRegex(ValueError,'skip_qc=true'):
            self.run_bridge()
        self.assertFalse(self.out.exists())

    def test_wrong_output_binding_fails_without_COMPLETE(self):
        def altered(argv,**kwargs):
            completed = self.fake(argv,**kwargs)
            path = Path(argv[argv.index('--out')+1]); data = self.fixture.load(path)
            data['input_binding_sha256'] = 'f'*64; self.fixture.save(path,data)
            return completed
        with self.assertRaisesRegex(ValueError,'does not bind'):
            self.run_bridge(altered)
        self.assertFalse((self.out/'COMPLETE').exists())
        self.assertEqual(self.fixture.load(self.out/'receipt.json')['status'],'FAILED')

    def test_failure_keeps_logs_without_complete_marker(self):
        with self.assertRaisesRegex(RuntimeError,'API failed'):
            self.run_bridge(lambda argv,**kwargs:subprocess.CompletedProcess(argv,2))
        self.assertFalse((self.out/'COMPLETE').exists())
        self.assertTrue((self.out/'api.stderr.log').exists())

    def test_conditional_support_cannot_keep_abstention_reasons(self):
        def altered(argv,**kwargs):
            completed = self.fake(argv,**kwargs)
            path = Path(argv[argv.index('--out')+1]); data = self.fixture.load(path)
            data['observation']['state'] = 'conditional_target_local_support'; self.fixture.save(path,data)
            return completed
        with self.assertRaisesRegex(ValueError,'contradicts abstention'):
            self.run_bridge(altered)

    def test_input_change_during_execution_revokes_complete(self):
        def altered(argv,**kwargs):
            completed = self.fake(argv,**kwargs)
            with self.fixture.metrics.open('a') as stream: stream.write('\n')
            return completed
        with self.assertRaisesRegex(ValueError,'Input changed'):
            self.run_bridge(altered)
        self.assertFalse((self.out/'COMPLETE').exists())


if __name__ == '__main__':
    unittest.main()
