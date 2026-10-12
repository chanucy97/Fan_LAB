"""Replay frozen reference-localized CDS checks from saved BLAST HSPs.

The v1 replay contract preserves Step6 selection/translation rules. It executes
no aligner and makes no change to the upstream paired-anchor evidence state.
"""
from collections import Counter, defaultdict
import csv
import hashlib
import io
import itertools
import json
import math
from pathlib import Path, PurePosixPath
import re

from . import __version__
from .integrity import write_json


SCHEMA = 'glycan-evidence.cds-replay-input/1'
PROTOCOL = 'h5-cds-step6-v3'
FIELDS = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen sframe btop'.split()
ROLES = ['Blon_2175', 'Blon_2176', 'Blon_2177',
         'BT1_RY67_RS01230_Blon_2175_alternative',
         'BT1_RY67_RS01235_Blon_2176_alternative',
         'BT1_RY67_RS01240_Blon_2177_alternative']
EXCLUDED = 'GCF_001281305.1__RY67_RS13345'
NT_SHA = '6947da8f94af1b3e940aff65fc6f4edc8273472eea3f54ca09c0f2ae428f4ceb'
AA_SHA = '9398220e80751c3d6e7db271130181f0def28b9014ca1dc35fb0b8665f4ea8be'
INPUT_NAMES = ('all_candidate_full.fasta', 'CANDIDATE_SEQUENCE_LEDGER.tsv',
               'query_nt.fasta', 'query_aa.fasta', 'blastn_all_hsp.tsv',
               'tblastn_all_hsp.tsv', 'INPUT_FREEZE.json', 'OUTPUT_MANIFEST.json',
               'SUMMARY.json', 'COMMANDS.json')
CDS_FIELDS = ['candidate_id', 'reference_query', 'status', 'start_0based', 'end_0based',
              'strand', 'nt_length', 'in_frame_length', 'start_codon', 'stop_codon',
              'internal_stop_n', 'unknown_aa_n', 'protein_length_aa', 'same_reference_nt',
              'same_reference_encoded_protein', 'nt_sha256', 'encoded_protein_sha256']
BEST_FIELDS = ['query_id', 'candidate_id', 'status', 'hsps_n', *FIELDS[2:],
               'query_span_fraction', 'region_interpretation']
CODONS = dict(zip((''.join(c) for c in itertools.product('TCAG', repeat=3)),
                 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
EVIDENCE = {'complete_inner_edges', 'unresolved_partial_or_below_screen_anchors'}
SCOPE = ('Reference-localized coding consistency from saved HSPs; no new alignment, '
         'readback, expression, function, strain or paired-anchor inference')


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _object(pairs):
    result = {}
    for key, value in pairs:
        _require(key not in result, 'Duplicate JSON key: ' + key)
        result[key] = value
    return result


def _json(data):
    return json.loads(data.decode('utf-8'), object_pairs_hook=_object,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError('Nonfinite JSON value')))


def _is_sha(value):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{64}', value) is not None


def _relative_file(root, name):
    _require(isinstance(name, str) and name, 'Missing input path')
    rel = PurePosixPath(name)
    _require(not rel.is_absolute() and not any(p in {'', '.', '..'} for p in name.split('/'))
             and '\\' not in name and ':' not in name, 'Input path must be safe and relative')
    path = root.joinpath(*rel.parts)
    _require(path.resolve().is_relative_to(root), 'Input path escapes manifest directory')
    _require(not any(p.is_symlink() for p in [path, *path.parents] if p != root.parent),
             'Symbolic links are not valid replay input paths')
    _require(path.is_file(), 'Input file missing: ' + name)
    return path


def _load_bound_inputs(input_manifest, expected_manifest_sha256):
    """Hash every declared file before interpreting any source file bytes."""
    _require(_is_sha(expected_manifest_sha256), 'Expected manifest SHA256 is required')
    path = Path(input_manifest)
    _require(path.is_file() and not path.is_symlink(), 'Input manifest must be a regular file')
    manifest_bytes = path.read_bytes()
    _require(_sha(manifest_bytes) == expected_manifest_sha256, 'Input manifest SHA256 mismatch')
    manifest = _json(manifest_bytes)
    _require(manifest.get('schema') == SCHEMA and manifest.get('protocol') == PROTOCOL,
             'Unsupported CDS replay schema or protocol')
    files = manifest.get('files')
    _require(isinstance(files, dict) and set(files) == set(INPUT_NAMES), 'Incomplete replay file contract')
    root = path.resolve().parent
    data, seen = {}, set()
    for name in INPUT_NAMES:
        record = files[name]
        _require(isinstance(record, dict), 'Invalid input record')
        p = _relative_file(root, record.get('path'))
        canonical = str(p.resolve()).casefold()
        _require(canonical not in seen, 'Duplicate input file path')
        seen.add(canonical)
        _require(type(record.get('bytes')) is int and record['bytes'] >= 0 and _is_sha(record.get('sha256')),
                 'Invalid input size or SHA256')
        blob = p.read_bytes()
        _require(len(blob) == record['bytes'] and _sha(blob) == record['sha256'], 'Input integrity mismatch: ' + name)
        data[name] = blob
    return manifest, data


def _table(data):
    reader = csv.DictReader(io.StringIO(data.decode('utf-8')), delimiter='\t')
    fields = reader.fieldnames
    _require(fields and len(fields) == len(set(fields)) and all(fields), 'Invalid TSV header')
    rows = list(reader)
    _require(all(None not in row and all(v is not None for v in row.values()) for row in rows), 'Malformed TSV row')
    return rows


def _fasta(data, *, protein=False):
    result, key = {}, None
    for line in data.decode('ascii').splitlines():
        if line.startswith('>'):
            _require(bool(line[1:].split()), 'Empty FASTA identifier')
            key = line[1:].split()[0]
            _require(key not in result, 'Duplicate FASTA identifier: ' + key)
            result[key] = ''
        elif line.strip():
            _require(key is not None, 'Sequence before FASTA header')
            result[key] += line.strip().upper()
    alphabet = set('ACDEFGHIKLMNPQRSTVWYBXZJUO*') if protein else set('ACGTRYMKBDHVNSW')
    _require(all(seq and set(seq) <= alphabet for seq in result.values()), 'Empty or invalid FASTA sequence')
    return result


def translate(sequence):
    if len(sequence) % 3:
        return None
    aa = ''.join(CODONS.get(sequence[i:i+3], 'X') for i in range(0, len(sequence), 3))
    return 'M' + aa[1:] if sequence[:3] in {'ATG', 'GTG', 'TTG'} else aa


def rc(sequence):
    return sequence.translate(str.maketrans('ACGTRYMKBDHVNSW', 'TGCAYRKMVHDBNSW'))[::-1]


def coding(sid, sequence, query, nt, aa, hit):
    """Apply the frozen full-query-span and open-coding checks, without promotion."""
    row = dict.fromkeys(CDS_FIELDS, '')
    row.update(candidate_id=sid, reference_query=query, status='NO_FULL_QUERY_SPAN_HSP',
               same_reference_nt=False, same_reference_encoded_protein=False)
    cds = protein = None
    if hit['status'] == 'HSP_PRESENT' and float(hit['query_span_fraction']) == 1:
        a, b = int(hit['sstart']), int(hit['send'])
        x, y = min(a, b)-1, max(a, b)
        _require(0 <= x < y <= len(sequence), 'CDS slice outside candidate')
        strand = '+' if a < b else '-'
        cds = sequence[x:y]
        if strand == '-':
            cds = rc(cds)
        translated = translate(cds)
        valid = (translated is not None and cds[:3] in {'ATG', 'GTG', 'TTG'} and
                 translated.endswith('*') and '*' not in translated[:-1] and 'X' not in translated)
        protein = translated[:-1] if valid else None
        row.update(status='COMPLETE_OPEN_CODING_SEQUENCE' if valid else 'FULL_SPAN_HSP_CODING_CHECK_FAILED',
                   start_0based=x, end_0based=y, strand=strand, nt_length=len(cds),
                   in_frame_length=len(cds)%3 == 0, start_codon=cds[:3], stop_codon=cds[-3:],
                   internal_stop_n=translated[:-1].count('*') if translated else '',
                   unknown_aa_n=translated.count('X') if translated else '',
                   protein_length_aa=len(translated)-1 if valid else '', same_reference_nt=cds == nt,
                   same_reference_encoded_protein=protein == aa if valid else False,
                   nt_sha256=_sha(cds.encode('ascii')), encoded_protein_sha256=_sha(protein.encode('ascii')) if valid else '')
    return row, cds, protein


def best_hits(data, queries, candidates):
    """Retain query x candidate denominator; break exact score ties by input order."""
    groups = defaultdict(list)
    for values in csv.reader(io.StringIO(data.decode('ascii')), delimiter='\t'):
        _require(len(values) == len(FIELDS), 'Malformed BLAST row')
        row = dict(zip(FIELDS, values))
        q, sid = row['qseqid'], row['sseqid']
        _require(q in queries and sid in candidates, 'Unknown BLAST identifier')
        for field in ('length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'qlen', 'slen', 'sframe'):
            _require(re.fullmatch('-?[0-9]+', row[field]) is not None, 'Invalid integer in BLAST row')
        for field in ('pident', 'evalue', 'bitscore'):
            _require(math.isfinite(float(row[field])) and float(row[field]) >= 0, 'Invalid numeric BLAST field')
        _require(float(row['pident']) <= 100 and int(row['length']) > 0 and int(row['mismatch']) >= 0 and
                 int(row['gapopen']) >= 0 and abs(int(row['sframe'])) in {1, 2, 3}, 'Invalid BLAST measurement')
        _require(int(row['qlen']) == len(queries[q]) and int(row['slen']) == len(candidates[sid]), 'BLAST length mismatch')
        _require(1 <= int(row['qstart']) <= int(row['qend']) <= int(row['qlen']), 'Query coordinates invalid')
        _require(1 <= min(int(row['sstart']), int(row['send'])) <= max(int(row['sstart']), int(row['send'])) <= int(row['slen']), 'Subject coordinates invalid')
        row['query_span_fraction'] = (int(row['qend'])-int(row['qstart'])+1)/int(row['qlen'])
        groups[q, sid].append(row)
    rows = []
    for q in queries:
        for sid in candidates:
            hits = groups[q, sid]
            row = dict(query_id=q, candidate_id=sid, status='NO_HSP_AT_PARAMETERS', hsps_n=0,
                       **dict.fromkeys(FIELDS[2:], ''), query_span_fraction='',
                       region_interpretation='EXCLUDED_CORRECTED_ANNOTATION' if q == EXCLUDED else 'RETAINED')
            if hits:
                h = max(hits, key=lambda x: (float(x['bitscore']), x['query_span_fraction'], float(x['pident'])))
                row.update(status='HSP_PRESENT', hsps_n=len(hits), **{k: h[k] for k in FIELDS[2:]},
                           query_span_fraction=h['query_span_fraction'])
            rows.append(row)
    return rows


def _count(value):
    _require(type(value) is int or isinstance(value, str) and re.fullmatch('[0-9]+', value), 'Invalid denominator')
    value = int(value)
    _require(value >= 0, 'Negative denominator')
    return value


def _validate_candidates(sequences, ledger, freeze):
    libraries = freeze.get('libraries')
    _require(isinstance(libraries, list) and libraries, 'Library denominator missing')
    _require(all(isinstance(row, dict) and isinstance(row.get('sample'), str) and row['sample'] and
                 isinstance(row.get('status'), str) and row['status'] for row in libraries), 'Invalid library ledger')
    samples = [row['sample'] for row in libraries]
    _require(len(samples) == len(set(samples)), 'Duplicate library sample')
    _require(all(set(row) == set(libraries[0]) for row in libraries), 'Inconsistent library columns')
    _require(all(row['status'] not in {'NOT_STARTED', 'RUNNING', 'QUEUED', 'PENDING'} for row in libraries), 'Nonterminal library input')
    required = {'candidate_id', 'sample', 'full_sequence_sha256', 'interval_sequence_sha256',
                'anchored_interval_length_bp', 'candidate_evidence'}
    _require(all(required <= set(row) for row in ledger), 'Incomplete candidate ledger')
    ids = [row['candidate_id'] for row in ledger]
    _require(len(ids) == len(set(ids)) and set(ids) == set(sequences), 'Candidate denominator mismatch or duplicate')
    _require(_count(freeze['candidate_count']) == len(ids), 'Frozen candidate count mismatch')
    _require(dict(Counter(row['candidate_evidence'] for row in ledger)) == freeze.get('candidate_evidence_counts'), 'Frozen geometry state counts mismatch')
    for row in ledger:
        seq = sequences[row['candidate_id']]
        _require(row['sample'] in samples, 'Candidate sample absent from library denominator')
        _require(_sha(seq.encode('ascii')) == row['full_sequence_sha256'] == row['interval_sequence_sha256'], 'Candidate sequence SHA mismatch')
        _require(len(seq) == _count(row['anchored_interval_length_bp']), 'Candidate geometry length mismatch')
        _require(row['candidate_evidence'] in EVIDENCE, 'Unknown geometry evidence state')
        if 'interval_start_0based' in row and 'interval_end_0based_exclusive' in row:
            start, end = _count(row['interval_start_0based']), _count(row['interval_end_0based_exclusive'])
            _require(end > start and end-start == len(seq), 'Invalid geometry interval coordinates')
    for row in libraries:
        n = sum(candidate['sample'] == row['sample'] for candidate in ledger)
        _require(row.get('candidate_count') in {'', None} and n == 0 or
                 row.get('candidate_count') not in {'', None} and _count(row['candidate_count']) == n,
                 'Library candidate denominator mismatch')
    return libraries


def _source_bindings(data):
    freeze, summary, output = (_json(data[name]) for name in ('INPUT_FREEZE.json', 'SUMMARY.json', 'OUTPUT_MANIFEST.json'))
    _require(summary.get('status') == 'COMPLETE', 'Source Step6 is incomplete')
    _require(summary.get('input_freeze_sha256') == _sha(data['INPUT_FREEZE.json']), 'Source input-freeze binding mismatch')
    _require(_is_sha(freeze.get('execution_freeze_sha256')) and
             summary.get('execution_freeze_sha256') == freeze['execution_freeze_sha256'], 'Source execution binding mismatch')
    _require(summary.get('mode') == freeze.get('mode'), 'Source mode mismatch')
    for name in INPUT_NAMES:
        if name in {'INPUT_FREEZE.json', 'OUTPUT_MANIFEST.json'}:
            continue
        record = output.get('files', {}).get(name, {})
        _require(record.get('sha256') == _sha(data[name]) and record.get('bytes') == len(data[name]), 'Source output binding mismatch: ' + name)
    for name in ('all_candidate_full.fasta', 'CANDIDATE_SEQUENCE_LEDGER.tsv'):
        record = freeze.get('files', {}).get(name, {})
        _require(record.get('sha256') == _sha(data[name]) and record.get('bytes') == len(data[name]), 'Source candidate binding mismatch: ' + name)
    _require(_sha(data['query_nt.fasta']) == NT_SHA and _sha(data['query_aa.fasta']) == AA_SHA, 'Frozen reference query bytes changed')
    commands = _json(data['COMMANDS.json'])
    _require(isinstance(commands, list), 'Invalid source commands')
    tools = []
    for command in commands:
        argv = command.get('argv')
        _require(command.get('returncode') == 0 and isinstance(argv, list) and argv and
                 all(isinstance(value, str) for value in argv), 'Source alignment unsuccessful')
        tool = argv[0].replace('\\', '/').split('/')[-1]
        _require(tool in {'blastn', 'tblastn'} and len(argv[1:]) % 2 == 0, 'Unsupported source command')
        pairs = list(zip(argv[1::2], argv[2::2]))
        options = dict(pairs)
        _require(len(pairs) == len(options), 'Duplicate source command option')
        expected = {'-query': 'query_nt.fasta' if tool == 'blastn' else 'query_aa.fasta',
                    '-subject': 'all_candidate_full.fasta', '-evalue': '1e-5',
                    '-outfmt': '6 ' + ' '.join(FIELDS), '-num_threads': '2', '-out': tool+'_all_hsp.tsv'}
        expected.update({'-task': 'blastn', '-word_size': '7', '-dust': 'no'} if tool == 'blastn' else {'-seg': 'no'})
        _require(options == expected, 'Source alignment parameters differ from frozen protocol')
        tools.append(tool)
    _require(sorted(tools) == ['blastn', 'tblastn'] if _count(freeze['candidate_count']) else not tools,
             'Incomplete or duplicate source alignment commands')
    return freeze, summary


def _write_table(path, rows, fields):
    with path.open('x', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def run_cds_replay(input_manifest, out, *, expected_manifest_sha256):
    """Validate a portable frozen bundle and replay CDS tables into a new directory."""
    manifest, data = _load_bound_inputs(input_manifest, expected_manifest_sha256)
    freeze, source_summary = _source_bindings(data)
    sequences = _fasta(data['all_candidate_full.fasta'])
    ledger = _table(data['CANDIDATE_SEQUENCE_LEDGER.tsv'])
    libraries = _validate_candidates(sequences, ledger, freeze)
    nt, aa = _fasta(data['query_nt.fasta']), _fasta(data['query_aa.fasta'], protein=True)
    _require(len(nt) == 33 and len(aa) == 32 and EXCLUDED in nt and EXCLUDED not in aa, 'Reference query set mismatch')
    for query in ROLES:
        _require(query in nt and query in aa and translate(nt[query]) == aa[query]+'*', 'Reference coding translation mismatch')
    all_best = {tool: best_hits(data[tool+'_all_hsp.tsv'], queries, sequences)
                for tool, queries in [('blastn', nt), ('tblastn', aa)]}
    lookup = {(row['query_id'], row['candidate_id']): row for row in all_best['blastn']}
    rows, coding_sequences, proteins = [], {}, {}
    for sid, seq in sequences.items():
        for query in ROLES:
            row, cds, protein = coding(sid, seq, query, nt[query], aa[query], lookup[query, sid])
            rows.append(row)
            key = sid+'__'+query
            if cds is not None:
                coding_sequences[key] = cds
            if protein is not None:
                proteins[key] = protein
    library_rows = []
    for library in libraries:
        own = [row for row in ledger if row['sample'] == library['sample']]
        ids = {row['candidate_id'] for row in own}
        library_rows.append({**library, 'cds_candidate_count': len(ids),
            'complete_inner_candidate_count': sum(row['candidate_evidence'] == 'complete_inner_edges' for row in own),
            'partial_inner_candidate_count': sum(row['candidate_evidence'] == 'unresolved_partial_or_below_screen_anchors' for row in own),
            'complete_cds_query_roles': sum(row['candidate_id'] in ids and row['status'] == 'COMPLETE_OPEN_CODING_SEQUENCE' for row in rows),
            'cds_state': 'CHARACTERIZED' if own else ('NOT_APPLICABLE_NO_CANDIDATES' if library['status'] == 'COMPLETE' else 'NO_VALID_CANDIDATE_GEOMETRY_RETAIN_PARENT_STATE')})
    measured = {'libraries': len(libraries), 'candidate_count': len(sequences),
                'candidate_evidence_counts': dict(Counter(row['candidate_evidence'] for row in ledger)),
                'nt_queries': len(nt), 'aa_queries': len(aa), 'six_role_cds_rows': len(rows),
                'cds_status_counts': dict(Counter(row['status'] for row in rows))}
    _require(all(source_summary.get(key) == value for key, value in measured.items()), 'Replayed source summary mismatch')
    receipt = {'schema': 'glycan-evidence.cds-replay/1', 'status': 'COMPLETE', 'package_version': __version__,
               'protocol': PROTOCOL, 'input_manifest_sha256': expected_manifest_sha256,
               'inputs': manifest['files'], 'source_mode': freeze['mode'], 'source_execution_sha256': freeze['execution_freeze_sha256'],
               'new_alignment_calls': 0, 'query_nt_sha256': NT_SHA, 'query_aa_sha256': AA_SHA,
               'coverage_unit': 'reference query role x retained candidate; separate complete library denominator',
               'scope': SCOPE, **measured}
    destination = Path(out)
    _require(not any(path.is_symlink() for path in [destination, *destination.parents]), 'Output path contains symbolic links')
    destination.mkdir(parents=True, exist_ok=False)
    try:
        for tool, best in all_best.items():
            _write_table(destination/(tool+'_best_by_query_candidate.tsv'), best, BEST_FIELDS)
        _write_table(destination/'CANDIDATE_CDS.tsv', rows, CDS_FIELDS)
        _write_table(destination/'CANDIDATE_SEQUENCE_LEDGER.tsv', ledger,
                     list(ledger[0]) if ledger else list(csv.reader(io.StringIO(data['CANDIDATE_SEQUENCE_LEDGER.tsv'].decode()), delimiter='\t'))[0])
        _write_table(destination/'ALL_LIBRARY_CDS_LEDGER.tsv', library_rows, list(library_rows[0]))
        _write_table(destination/'QUERY_ROLES.tsv', [{'query_id': query, 'core_cds_role': query in ROLES,
            'region_interpretation': 'EXCLUDED_CORRECTED_ANNOTATION' if query == EXCLUDED else 'RETAINED',
            'has_protein_query': query in aa} for query in nt],
            ['query_id', 'core_cds_role', 'region_interpretation', 'has_protein_query'])
        for name, records in [('localized_coding_sequences.fna', coding_sequences), ('localized_encoded_proteins.faa', proteins)]:
            (destination/name).write_text(''.join('>'+key+'\n'+value+'\n' for key, value in records.items()), encoding='ascii', newline='\n')
        receipt['outputs'] = [{'path': path.name, 'bytes': path.stat().st_size, 'sha256': _sha(path.read_bytes())}
                              for path in sorted(destination.iterdir()) if path.is_file()]
        write_json(destination/'receipt.json', receipt)
        (destination/'COMPLETE').write_text(_sha((destination/'receipt.json').read_bytes())+'\n', encoding='ascii')
    except BaseException:
        (destination/'COMPLETE').unlink(missing_ok=True)
        raise
    return receipt
