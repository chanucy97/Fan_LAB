"""Detect changed reports and hostile manifest paths without touching raw reads."""
import json
from pathlib import Path
import tempfile
import unittest

from glycan_evidence.integrity import seal_outputs, verify_outputs


class IntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        (self.root / 'COMPLETE').write_text('done\n')
        (self.root / 'metrics.tsv').write_text('sample\tvalue\na\t1\n')
        (self.root / 'raw.fastq').write_text('not a report\n')

    def test_small_output_scope_and_mutation(self):
        self.assertEqual(seal_outputs(self.root)['report_files'], 2)
        (self.root / 'raw.fastq').write_text('raw deliberately out of scope\n')
        self.assertEqual(verify_outputs(self.root)['status'], 'PASS')
        (self.root / 'metrics.tsv').write_text('changed\n')
        with self.assertRaisesRegex(ValueError, 'integrity'):
            verify_outputs(self.root)

    def test_missing_complete_rejected(self):
        (self.root / 'COMPLETE').unlink()
        with self.assertRaisesRegex(ValueError, 'completed'):
            seal_outputs(self.root)

    def test_compressed_base_report_is_bound_without_including_raw_fastq(self):
        import gzip
        positions = self.root / 'candidate1.positions.tsv.gz'
        positions.write_bytes(gzip.compress(b'position\tA_qname_total\n0\t2\n'))
        (self.root / 'reads.fastq.gz').write_bytes(gzip.compress(b'raw input remains outside report scope'))
        self.assertEqual(seal_outputs(self.root)['report_files'], 3)
        positions.write_bytes(gzip.compress(b'position\tA_qname_total\n0\t9\n'))
        with self.assertRaisesRegex(ValueError, 'integrity'):
            verify_outputs(self.root)

    def test_reseal_does_not_bless_modified_reports(self):
        seal_outputs(self.root)
        (self.root / 'metrics.tsv').write_text('changed\n')
        with self.assertRaises(FileExistsError):
            seal_outputs(self.root)

    def test_paths_and_duplicate_entries_rejected(self):
        seal_outputs(self.root)
        path = self.root / 'artifact_manifest.json'
        original = json.loads(path.read_text())
        for bad in ['../secret.txt', 'C:/secret.txt', '/secret.txt', '..\\secret.txt']:
            with self.subTest(path=bad):
                data = json.loads(json.dumps(original))
                data['files'][0]['path'] = bad
                path.write_text(json.dumps(data))
                with self.assertRaisesRegex(ValueError, 'Unsafe'):
                    verify_outputs(self.root)
        original['files'].append(original['files'][0])
        path.write_text(json.dumps(original))
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            verify_outputs(self.root)


if __name__ == '__main__':
    unittest.main()
