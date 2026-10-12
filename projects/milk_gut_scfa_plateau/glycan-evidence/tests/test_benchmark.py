"""Keep biological denominators, uncertain truth, and leakage visible."""
import csv
import json
from pathlib import Path
import tempfile
import unittest

from glycan_evidence.benchmark import evaluate


class BenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.truth = [
            dict(sample='a1', target='infantis', present='1', genome='A', cluster='A', split='test', stratum='low'),
            dict(sample='a2', target='infantis', present='1', genome='A', cluster='A', split='test', stratum='high'),
            dict(sample='b1', target='infantis', present='1', genome='B', cluster='B', split='test', stratum='high'),
            dict(sample='c1', target='infantis', present='0', genome='C', cluster='C', split='test', stratum='high'),
            dict(sample='u1', target='infantis', present='unresolved', genome='U', cluster='U', split='test', stratum='high')]
        self.calls = [dict(sample=r['sample'], target=r['target'], margin='0', identity_state=s)
                      for r, s in zip(self.truth, ['insufficient_evidence', 'detected', 'detected', 'detected', 'detected'])]

    def write(self, filename, rows):
        path = self.root / filename
        with path.open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        return path

    def run_report(self):
        return evaluate(self.write('calls.tsv', self.calls), self.write('truth.tsv', self.truth), self.root / 'out')

    def test_condition_and_cluster_denominators_differ(self):
        result = self.run_report()
        t = result['targets']['infantis']
        self.assertEqual((t['TP'], t['FN'], t['FP'], t['TN']), (2, 1, 1, 0))
        self.assertAlmostEqual(t['condition_sensitivity'], 2 / 3)
        self.assertEqual(t['cluster_macro_sensitivity'], .75)
        self.assertEqual((t['technical_conditions'], t['unique_genomes'], t['independent_clusters']), (4, 3, 3))
        self.assertEqual(result['declared_unique_samples'], 5)
        self.assertEqual(result['truth_unresolved'], [dict(sample='u1', target='infantis', reason='truth_unresolved')])
        self.assertTrue((self.root / 'out/COMPLETE').is_file())

    def test_missing_predictions_not_dropped(self):
        self.calls.pop()
        with self.assertRaisesRegex(ValueError, 'Missing prediction'):
            self.run_report()
        self.assertFalse((self.root / 'out').exists())

    def test_truth_uncertainty_not_scored_as_negative(self):
        result = self.run_report()
        self.assertEqual(result['targets']['infantis']['negative_clusters'], 1)
        self.assertEqual(result['targets']['infantis']['FP'], 1)

    def test_cluster_leakage_rejected_even_outside_selected_split(self):
        self.truth.append(dict(self.truth[0], sample='training', split='development'))
        with self.assertRaisesRegex(ValueError, 'multiple splits'):
            self.run_report()

    def test_contradictory_cluster_truth_rejected(self):
        self.truth[1]['present'] = '0'
        with self.assertRaisesRegex(ValueError, 'contradictory'):
            self.run_report()

    def test_undeclared_and_duplicate_predictions_rejected(self):
        self.calls.append(dict(self.calls[0], sample='hidden'))
        with self.assertRaisesRegex(ValueError, 'undeclared'):
            self.run_report()
        self.calls[-1] = self.calls[0].copy()
        with self.assertRaisesRegex(ValueError, 'Duplicate prediction'):
            self.run_report()

    def test_nonfinite_and_unknown_margins_rejected(self):
        for value in ['nan', 'inf', '-inf', '3', 'error']:
            with self.subTest(value=value):
                self.calls[0]['margin'] = value
                with self.assertRaisesRegex(ValueError, 'margin'):
                    self.run_report()

    def test_duplicate_truth_rejected(self):
        self.truth.append(self.truth[0].copy())
        with self.assertRaisesRegex(ValueError, 'Duplicate truth'):
            self.run_report()

    def test_report_never_overwrites_existing_output(self):
        self.run_report()
        before = (self.root / 'out/benchmark.json').read_bytes()
        with self.assertRaises(FileExistsError):
            self.run_report()
        self.assertEqual((self.root / 'out/benchmark.json').read_bytes(), before)


if __name__ == '__main__':
    unittest.main()
