"""Small synthetic boundary tests; real fixed-query replay has separate receipts."""
import csv
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from glycan_evidence import cds
from glycan_evidence.integrity import seal_outputs, verify_outputs


def blob(value):
    return (json.dumps(value, indent=2)+'\n').encode()


def table(rows, fields):
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode()


def sha(value):
    return hashlib.sha256(value).hexdigest()


def hsp(q='q', sid='c', **changes):
    values = dict(zip(cds.FIELDS, [q, sid, '100', '9', '0', '0', '1', '9',
        '1', '9', '1e-20', '40', '9', '9', '1', '9']))
    values.update(changes)
    return ('\t'.join(values[name] for name in cds.FIELDS)+'\n').encode()


def fixture(directory):
    """Synthetic query hashes are patched only by the test caller."""
    sequence = 'ATGAAATAA'
    names = cds.ROLES + ['context_'+str(i) for i in range(26)]
    nt = ''.join('>'+name+'\n'+sequence+'\n' for name in names+[cds.EXCLUDED]).encode()
    aa = ''.join('>'+name+'\nMK\n' for name in names).encode()
    row = dict(candidate_id='c', sample='sample', full_sequence_sha256=sha(sequence.encode()),
               interval_sequence_sha256=sha(sequence.encode()), anchored_interval_length_bp='9',
               candidate_evidence='unresolved_partial_or_below_screen_anchors',
               interval_start_0based='0', interval_end_0based_exclusive='9')
    files = {'all_candidate_full.fasta': ('>c\n'+sequence+'\n').encode(),
             'CANDIDATE_SEQUENCE_LEDGER.tsv': table([row], list(row)),
             'query_nt.fasta': nt, 'query_aa.fasta': aa,
             'blastn_all_hsp.tsv': hsp(q=cds.ROLES[0]), 'tblastn_all_hsp.tsv': b''}
    record = lambda name: dict(bytes=len(files[name]), sha256=sha(files[name]))
    freeze = dict(mode='SYNTHETIC_TEST', execution_freeze_sha256='a'*64,
                  files={name: record(name) for name in ['all_candidate_full.fasta', 'CANDIDATE_SEQUENCE_LEDGER.tsv']},
                  libraries=[dict(sample='sample', status='COMPLETE', candidate_count=1),
                             dict(sample='zero', status='COMPLETE', candidate_count=0)],
                  candidate_count=1, candidate_evidence_counts={'unresolved_partial_or_below_screen_anchors': 1})
    files['INPUT_FREEZE.json'] = blob(freeze)
    files['SUMMARY.json'] = blob(dict(status='COMPLETE', mode=freeze['mode'],
        input_freeze_sha256=sha(files['INPUT_FREEZE.json']), execution_freeze_sha256='a'*64,
        libraries=2, candidate_count=1, candidate_evidence_counts=freeze['candidate_evidence_counts'],
        nt_queries=33, aa_queries=32, six_role_cds_rows=6,
        cds_status_counts={'COMPLETE_OPEN_CODING_SEQUENCE': 1, 'NO_FULL_QUERY_SPAN_HSP': 5}))
    commands = []
    for tool, query, extra in [('blastn', 'query_nt.fasta', ['-task', 'blastn', '-word_size', '7', '-dust', 'no']),
                               ('tblastn', 'query_aa.fasta', ['-seg', 'no'])]:
        commands.append(dict(returncode=0, argv=[tool, '-query', query, '-subject', 'all_candidate_full.fasta',
            '-evalue', '1e-5', '-outfmt', '6 '+' '.join(cds.FIELDS), '-num_threads', '2', '-out', tool+'_all_hsp.tsv', *extra]))
    files['COMMANDS.json'] = blob(commands)
    files['OUTPUT_MANIFEST.json'] = blob(dict(files={name: record(name) for name in files if name != 'INPUT_FREEZE.json'}))
    manifest = dict(schema=cds.SCHEMA, protocol=cds.PROTOCOL,
                    files={name: dict(path=name, **record(name)) for name in cds.INPUT_NAMES})
    for name, data in files.items():
        (directory/name).write_bytes(data)
    path = directory/'replay.json'
    path.write_bytes(blob(manifest))
    return path, manifest, files


class CodingBoundaryTests(unittest.TestCase):
    def run_coding(self, sequence, **changes):
        hit = dict(status='HSP_PRESENT', query_span_fraction=1, sstart=1, send=len(sequence))
        hit.update(changes)
        return cds.coding('candidate', sequence, 'reference', 'ATGAAATAA', 'MK', hit)

    def test_complete_reverse_and_alternative_start(self):
        row, seq, protein = self.run_coding('TTATTTCAT', sstart=9, send=1)
        self.assertEqual((row['status'], seq, protein, row['strand']),
                         ('COMPLETE_OPEN_CODING_SEQUENCE', 'ATGAAATAA', 'MK', '-'))
        row, _, protein = self.run_coding('GTGAAATAA')
        self.assertTrue(row['same_reference_encoded_protein'])
        self.assertFalse(row['same_reference_nt'])
        self.assertEqual(protein, 'MK')

    def test_stop_truncation_unknown_and_internal_stop_fail(self):
        for seq in ('ATGAAAAAA', 'ATGAAATA', 'ATGTAAAAATAA', 'ATGNNNTAA'):
            with self.subTest(sequence=seq):
                row, _, protein = self.run_coding(seq)
                self.assertEqual(row['status'], 'FULL_SPAN_HSP_CODING_CHECK_FAILED')
                self.assertIsNone(protein)

    def test_partial_full_length_sequence_is_not_promoted(self):
        row, seq, _ = self.run_coding('ATGAAATAA', query_span_fraction=0.9)
        self.assertEqual(row['status'], 'NO_FULL_QUERY_SPAN_HSP')
        self.assertIsNone(seq)

    def test_coordinate_bounds_rejected(self):
        for bounds in [(0, 9), (1, 10), (-1, 9)]:
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                self.run_coding('ATGAAATAA', sstart=bounds[0], send=bounds[1])

    def test_best_hit_tie_order_and_full_matrix(self):
        raw = hsp(q='q1', pident='90') + hsp(q='q1', pident='90', btop='8AG')
        rows = cds.best_hits(raw, {'q1':'ATGAAATAA', 'q2':'ATGAAATAA'}, {'c':'ATGAAATAA', 'c2':'ATGAAATAA'})
        self.assertEqual(len(rows), 4)
        self.assertEqual((rows[0]['hsps_n'], rows[0]['btop']), (2, '9'))
        self.assertEqual(rows[-1]['status'], 'NO_HSP_AT_PARAMETERS')

    def test_malformed_hsp_rejected(self):
        for changes in ({'sstart':'0'}, {'send':'10'}, {'qlen':'8'}, {'qstart':'10'},
                        {'bitscore':'nan'}, {'pident':'101'}, {'sframe':'0'}, {'qseqid':'other'}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                cds.best_hits(hsp(**changes), {'q':'ATGAAATAA'}, {'c':'ATGAAATAA'})

    def test_duplicate_and_invalid_fasta_rejected(self):
        for data in (b'>c\nATG\n>c\nATG\n', b'>c\n', b'ATG\n', b'>c\nAZG\n'):
            with self.subTest(data=data), self.assertRaises(ValueError):
                cds._fasta(data)


class ReplayIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path, self.manifest, self.files = fixture(self.root)

    def load(self):
        return cds._load_bound_inputs(self.path, sha(self.path.read_bytes()))

    def rewrite_manifest(self):
        self.path.write_bytes(blob(self.manifest))

    def test_manifest_anchor_and_version_rejected(self):
        with self.assertRaisesRegex(ValueError, 'SHA256 mismatch'):
            cds._load_bound_inputs(self.path, 'a'*64)
        self.manifest['protocol'] = 'changed'
        self.rewrite_manifest()
        with self.assertRaisesRegex(ValueError, 'protocol'):
            self.load()

    def test_absolute_traversal_and_duplicate_paths_rejected(self):
        original = self.manifest['files']['all_candidate_full.fasta']['path']
        for path in ('../other', '/absolute', 'C:/absolute', '..\\other', './file', 'a//file',
                     self.manifest['files']['query_nt.fasta']['path']):
            with self.subTest(path=path):
                self.manifest['files']['all_candidate_full.fasta']['path'] = path
                self.rewrite_manifest()
                with self.assertRaises(ValueError):
                    self.load()
        self.manifest['files']['all_candidate_full.fasta']['path'] = original

    def test_all_hashes_checked_before_source_interpretation(self):
        # Earlier source bytes are malformed but bound; a later file has a wrong hash.
        (self.root/'INPUT_FREEZE.json').write_bytes(b'invalid JSON')
        self.manifest['files']['INPUT_FREEZE.json'].update(bytes=12, sha256=sha(b'invalid JSON'))
        self.manifest['files']['COMMANDS.json']['sha256'] = 'b'*64
        self.rewrite_manifest()
        with patch.object(cds, '_source_bindings', side_effect=AssertionError('parsed too soon')):
            with self.assertRaisesRegex(ValueError, 'Input integrity mismatch: COMMANDS.json'):
                cds.run_cds_replay(self.path, self.root/'out', expected_manifest_sha256=sha(self.path.read_bytes()))
        self.assertFalse((self.root/'out'/'COMPLETE').exists())

    def test_source_cross_bundle_mismatch_rejected(self):
        data = dict(self.files)
        summary = json.loads(data['SUMMARY.json'])
        summary['input_freeze_sha256'] = 'b'*64
        data['SUMMARY.json'] = blob(summary)
        with self.assertRaisesRegex(ValueError, 'input-freeze binding'):
            cds._source_bindings(data)

    def test_duplicate_json_key_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate JSON'):
            cds._json(b'{"a":1,"a":2}')

    def test_missing_duplicate_and_count_mismatched_samples_rejected(self):
        freeze = json.loads(self.files['INPUT_FREEZE.json'])
        seq = cds._fasta(self.files['all_candidate_full.fasta'])
        ledger = cds._table(self.files['CANDIDATE_SEQUENCE_LEDGER.tsv'])
        for mutation in ('missing', 'duplicate', 'count', 'candidate_duplicate', 'geometry'):
            f, rows = json.loads(json.dumps(freeze)), [dict(row) for row in ledger]
            if mutation == 'missing': f['libraries'][0]['sample'] = 'unknown'
            if mutation == 'duplicate': f['libraries'].append(dict(f['libraries'][0]))
            if mutation == 'count': f['libraries'][0]['candidate_count'] = 0
            if mutation == 'candidate_duplicate': rows.append(dict(rows[0]))
            if mutation == 'geometry': rows[0]['interval_end_0based_exclusive'] = '8'
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                cds._validate_candidates(seq, rows, f)

    def test_fixed_queries_fail_without_test_patch(self):
        with self.assertRaisesRegex(ValueError, 'Frozen reference query bytes'):
            cds.run_cds_replay(self.path, self.root/'untrusted', expected_manifest_sha256=sha(self.path.read_bytes()))

    def test_synthetic_replay_retains_partial_and_zero_denominators(self):
        with patch.object(cds, 'NT_SHA', sha(self.files['query_nt.fasta'])), patch.object(cds, 'AA_SHA', sha(self.files['query_aa.fasta'])):
            receipt = cds.run_cds_replay(self.path, self.root/'out', expected_manifest_sha256=sha(self.path.read_bytes()))
        self.assertEqual((receipt['candidate_count'], receipt['libraries'], receipt['new_alignment_calls']), (1, 2, 0))
        rows = cds._table((self.root/'out'/'ALL_LIBRARY_CDS_LEDGER.tsv').read_bytes())
        self.assertEqual((rows[0]['partial_inner_candidate_count'], rows[0]['complete_inner_candidate_count'],
                          rows[0]['complete_cds_query_roles']), ('1', '0', '1'))
        self.assertEqual(rows[1]['cds_state'], 'NOT_APPLICABLE_NO_CANDIDATES')
        seal_outputs(self.root/'out')
        self.assertEqual(verify_outputs(self.root/'out')['status'], 'PASS')

    def test_failed_write_has_no_complete_marker(self):
        with patch.object(cds, 'NT_SHA', sha(self.files['query_nt.fasta'])), patch.object(cds, 'AA_SHA', sha(self.files['query_aa.fasta'])), patch.object(cds, '_write_table', side_effect=OSError('simulated write failure')):
            with self.assertRaises(OSError):
                cds.run_cds_replay(self.path, self.root/'out', expected_manifest_sha256=sha(self.path.read_bytes()))
        self.assertFalse((self.root/'out'/'COMPLETE').exists())


if __name__ == '__main__':
    unittest.main()
