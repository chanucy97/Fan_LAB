import copy
import json
from pathlib import Path
import tempfile
import unittest

from glycan_evidence import profile
from glycan_evidence.profile_io import load_expected, source_from_dict


class ProfileIOTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / 'model.json'
        panel = profile.MarkerPanel('longum', (('a', 100), ('b', 100)), 'f'*64)
        source = profile.Source('reference', 'reference_cluster', ('longum',))
        obs = profile.Observation(panel, source, 'development_sample', (2, 3), (50, 75), (75, 100), 'a'*64)
        self.model = profile.fit_expected([obs]).to_dict()

    def save(self, value):
        self.path.write_text(json.dumps(value))
        return self.path

    def test_model_roundtrip_preserves_diagnostic(self):
        expected = load_expected(self.save(self.model))
        heldout = profile.Observation(expected.panel, profile.Source('heldout', 'different_cluster', ()),
                                      'heldout_sample', (2, 0), (60, 0), (80, 0), 'b'*64)
        diagnostic = profile.diagnose(heldout, expected)
        self.assertEqual(diagnostic['training_sources'], ['reference'])
        self.assertAlmostEqual(diagnostic['zero_opportunity_mass'], .6)
        self.assertIsNone(diagnostic['identity_call'])

    def test_inconsistent_estimator_and_changed_panel_rejected(self):
        for field in ['probabilities', 'panel_sha256']:
            with self.subTest(field=field):
                value = copy.deepcopy(self.model)
                value[field] = [.5, .5] if field == 'probabilities' else 'b'*64
                with self.assertRaises(ValueError):
                    load_expected(self.save(value))

    def test_malformed_and_duplicate_fields_rejected(self):
        for value in [dict(self.model, probabilities=None), dict(self.model, source_probabilities=[None]),
                      dict(self.model, probabilities=[True, False])]:
            with self.assertRaises(ValueError):
                load_expected(self.save(value))
        self.path.write_text('{"schema": "x", "schema": "y"}')
        with self.assertRaises(ValueError):
            load_expected(self.path)

    def test_source_exposure_and_explicit_unknown_truth(self):
        value = dict(source_id='sample', cluster_id='subject', positive_targets=[], exposure='unlabelled_application')
        self.assertEqual(source_from_dict(value).positive_targets, ())
        for change in [{'positive_targets': None}, {'positive_targets': ['longum']}, {'exposure': 'unknown'}]:
            with self.assertRaises(ValueError):
                source_from_dict({**value, **change})


if __name__ == '__main__':
    unittest.main()
