"""CLI behavior at failure boundaries; no sequencing data in these fixtures."""
import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from glycan_evidence import cli
from glycan_evidence.identity import DEFAULT_RULES


class CliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_sample_manifest_errors(self):
        path = self.root / 'samples.tsv'
        for data in ['sample\nx\nx\n', 'sample\n../x\n', 'sample\tbatch\nx\n', 'sample\nx\textra\n', 'sample\n']:
            with self.subTest(data=data):
                path.write_text(data)
                with self.assertRaises(ValueError):
                    cli.sample_list(path)
        path.write_text('sample\tbatch\nx\t1\ny\t2\n')
        self.assertEqual(cli.sample_list(path), ['x', 'y'])

    def test_doctor_missing_dependency_is_actionable(self):
        info, status = cli.doctor(str(self.root), structural=True, require=True)
        self.assertEqual(status, 1)
        self.assertEqual(info['status'], 'MISSING_DEPENDENCIES')
        self.assertEqual(set(info['tools']), {'bowtie2', 'fastp', 'blastn', 'bowtie2-build', 'samtools'})

    def test_failed_run_never_has_complete(self):
        out = self.root / 'output'
        stderr = io.StringIO()
        with patch('glycan_evidence.identity.run_sample', side_effect=ValueError('bad paired input')), contextlib.redirect_stderr(stderr):
            code = cli.main(['run', '--r1', 'r1', '--r2', 'r2', '--sample', 'toy', '--panel', 'panel', '--rules', str(DEFAULT_RULES), '--out', str(out)])
        self.assertEqual(code, 2)
        self.assertFalse((out / 'COMPLETE').exists())
        self.assertFalse((out / 'artifact_manifest.json').exists())
        self.assertEqual(json.loads((out / 'receipt.json').read_text())['status'], 'FAILED')
        self.assertIn('bad paired input', stderr.getvalue())

    def test_unsafe_sample_rejected_before_output_creation(self):
        out = self.root / 'output'
        with contextlib.redirect_stderr(io.StringIO()):
            code = cli.main(['run', '--r1', 'r1', '--r2', 'r2', '--sample', '../toy', '--panel', 'panel', '--rules', 'rules', '--out', str(out)])
        self.assertEqual(code, 2)
        self.assertFalse(out.exists())

    def test_seal_failure_revokes_top_level_completion(self):
        calls=self.root/'calls.tsv'
        truth=self.root/'truth.tsv'
        calls.write_text('sample\ttarget\tidentity_state\tmargin\na\tinfantis\tdetected\t0\n')
        truth.write_text('sample\ttarget\tpresent\tgenome\tcluster\tsplit\tstratum\na\tinfantis\t1\tg\tg\ttest\thigh\n')
        out=self.root/'report'
        with patch('glycan_evidence.cli.seal_outputs',side_effect=OSError('disk failure')), contextlib.redirect_stderr(io.StringIO()):
            code=cli.main(['benchmark','--calls',str(calls),'--truth',str(truth),'--out',str(out)])
        self.assertEqual(code,2)
        self.assertFalse((out/'COMPLETE').exists())
        self.assertEqual(json.loads((out/'CLI_FAILED.json').read_text())['stage'],'seal_outputs')
        self.assertTrue((out/'benchmark.json').is_file())

    def test_invalid_rules_rejected_before_expensive_work(self):
        out=self.root/'bad_rules_output'
        with patch('glycan_evidence.identity.run_sample') as run, contextlib.redirect_stderr(io.StringIO()):
            code=cli.main(['run','--r1','r1','--r2','r2','--sample','toy','--panel','panel',
                           '--rules',str(self.root/'missing_rules.json'),'--out',str(out)])
        self.assertEqual(code,2)
        run.assert_not_called()
        self.assertFalse(out.exists())


if __name__ == '__main__':
    unittest.main()
