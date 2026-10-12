"""CLI wiring and completion boundaries for independent discrepancy replay."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from glycan_evidence import cli
from glycan_evidence.integrity import write_json,verify_outputs


class DiscrepancyCliTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.out=Path(self.temp.name)/'verified'
        self.args=['discrepancy-verify','--root','input-root','--sample','toy',
                   '--expected-freeze-sha256','a'*64,'--expected-results-sha256','b'*64,
                   '--out',str(self.out)]

    def fake_run(self,root,sample,freeze_sha,results_sha,out):
        self.assertEqual((root,sample,freeze_sha,results_sha),('input-root','toy','a'*64,'b'*64))
        path=Path(out)
        path.mkdir(exist_ok=False)
        receipt={'status':'COMPLETE','optimizer_calls':0}
        write_json(path/'receipt.json',receipt)
        (path/'COMPLETE').write_text('verified\n',encoding='utf-8')
        return receipt

    def test_success_dispatch_and_integrity_sealing(self):
        with patch('glycan_evidence.discrepancy_verify.run_discrepancy_verify',side_effect=self.fake_run) as run,contextlib.redirect_stdout(io.StringIO()):
            code=cli.main(self.args)
        self.assertEqual(code,0)
        run.assert_called_once()
        self.assertEqual(verify_outputs(self.out)['status'],'PASS')

    def test_rejected_input_is_not_complete(self):
        with patch('glycan_evidence.discrepancy_verify.run_discrepancy_verify',side_effect=ValueError('input hash mismatch')),contextlib.redirect_stderr(io.StringIO()):
            code=cli.main(self.args)
        self.assertEqual(code,2)
        self.assertFalse(self.out.exists())

    def test_sealing_failure_revokes_completion(self):
        with patch('glycan_evidence.discrepancy_verify.run_discrepancy_verify',side_effect=self.fake_run),patch('glycan_evidence.cli.seal_outputs',side_effect=OSError('disk full')),contextlib.redirect_stderr(io.StringIO()):
            code=cli.main(self.args)
        self.assertEqual(code,2)
        self.assertFalse((self.out/'COMPLETE').exists())
        self.assertTrue((self.out/'receipt.json').is_file())
        self.assertEqual(json.loads((self.out/'CLI_FAILED.json').read_text())['stage'],'seal_outputs')


if __name__=='__main__':
    unittest.main()
