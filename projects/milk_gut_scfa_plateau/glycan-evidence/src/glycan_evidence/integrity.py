"""Small-output integrity manifests, separate from scientific validation."""
from pathlib import Path, PurePosixPath
import hashlib
import json


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n', encoding='utf-8')
    temporary.replace(path)


def seal_outputs(directory):
    directory = Path(directory).resolve()
    if not (directory / 'COMPLETE').is_file():
        raise ValueError('Only completed output directories can be sealed')
    destination = directory / 'artifact_manifest.json'
    if destination.exists():
        raise FileExistsError(str(destination))
    items = []
    for p in sorted(directory.rglob('*')):
        if not p.is_file() or p == destination or p.suffix == '.tmp':
            continue
        if p.name != 'COMPLETE' and p.suffix not in {'.json', '.tsv', '.log', '.txt', '.fna', '.faa'} and not p.name.endswith('.positions.tsv.gz'):
            continue
        if p.is_symlink() or not p.resolve().is_relative_to(directory):
            raise ValueError('Output member is a symbolic link or escapes the output directory')
        items.append(dict(path=p.relative_to(directory).as_posix(), bytes=p.stat().st_size, sha256=digest(p)))
    write_json(destination, dict(schema_version=1, status='COMPLETE', files=items,
                                scope='Integrity of listed small reports; not independent scientific validation or raw-read checksum verification'))
    return verify_outputs(directory)


def verify_outputs(directory):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / 'artifact_manifest.json').read_text(encoding='utf-8'))
    if manifest.get('schema_version') != 1 or manifest.get('status') != 'COMPLETE':
        raise ValueError('Invalid completed output manifest')
    files = manifest.get('files')
    if not isinstance(files, list) or not files:
        raise ValueError('Manifest has no report files')
    seen = set()
    for item in files:
        name = item.get('path', '')
        rel = PurePosixPath(name)
        if not name or rel.is_absolute() or '..' in rel.parts or '\\' in name or ':' in name or name in seen:
            raise ValueError('Unsafe or duplicate output path')
        seen.add(name)
        p = directory.joinpath(*rel.parts)
        if p.is_symlink() or not p.resolve().is_relative_to(directory):
            raise ValueError('Output path escapes the report directory')
        if type(item.get('bytes')) is not int or item['bytes'] < 0:
            raise ValueError('Invalid artifact byte size')
        if not p.is_file() or p.stat().st_size != item['bytes'] or digest(p) != item.get('sha256'):
            raise ValueError(f'Artifact failed integrity check: {name}')
    if 'COMPLETE' not in seen:
        raise ValueError('Manifest does not bind the COMPLETE marker')
    return dict(status='PASS', report_files=len(files), scope=manifest.get('scope'))
