import json
import unittest

from glycan_evidence import cli
from glycan_evidence.structure_inference_io import run_structure_assignment
import test_structure_inference_io as fixture


class AssignmentInterfaceTest(fixture.InterfaceTest):
    def run_io(self):
        return run_structure_assignment(self.sam,self.root/'candidates.json',self.root/'profiles.json',
            self.root/'config.json',self.root/'out')

    def test_public_cli_seals_joint_output(self):
        argv = ['structure-assign']
        for name, path in [('sam', self.sam), ('candidates', self.root/'candidates.json'),
                           ('profiles', self.root/'profiles.json'), ('config', self.root/'config.json'),
                           ('out', self.root/'cli_out')]:
            argv += ['--'+name, str(path)]
        self.assertEqual(cli.main(argv), 0)
        self.assertEqual(cli.main(['verify','--directory',str(self.root/'cli_out')]), 0)
        fit = json.loads((self.root/'cli_out/fit.json').read_text())
        self.assertEqual(fit['schema'], 'glycan-evidence.structure-assignment.research.v2')
        self.assertAlmostEqual(fit['joint_bestfit_tv'], 0)
        self.assertEqual(fit['observed_support']['G']['informative_qnames'], 1)


if __name__ == '__main__':
    unittest.main()
