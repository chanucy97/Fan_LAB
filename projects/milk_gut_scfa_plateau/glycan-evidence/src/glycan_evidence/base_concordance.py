"""SHA-bound position-level assembly/read concordance, not variant calling.

Counting and inspection flags preserve the 2026-10-09 research audit. Input
validation and receipts are the package interface; no remapping is performed.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import gzip
import hashlib
from itertools import groupby
import json
from pathlib import Path
import re
import time

from . import __version__
from .integrity import digest, write_json
from .structure import file_record, read_fasta, utc_now


SCHEMA_VERSION = "glycan-evidence.base-concordance/1"
BASES = "ACGTN"
COMP = str.maketrans("ACGTN", "TGCAN")
REQUIRED_STRUCTURE_FIELDS = (
    "structure_id", "contig", "contig_sequence_sha256", "interval_start_0based",
    "interval_end_0based_exclusive", "orientation_to_reference", "interval_sequence_sha256",
)
POSITION_FIELDS = [
    "structure_id", "contig", "genomic_position_0based", "full_position_0based",
    "orientation_to_reference", "assembled_genomic_base", "assembled_oriented_base",
]
for _base in BASES:
    POSITION_FIELDS.extend([f"{_base}_read_F", f"{_base}_read_R"] +
                           [f"{_base}_qname_{s}" for s in ("F", "R", "both")] +
                           [f"{_base}_qname_total"])
POSITION_FIELDS.extend([
    "conflict_qname_F", "conflict_qname_R", "conflict_qname_both",
    "DEL_OR_SKIP_read_total", "LOW_BQ_read_total", "acgt_qname_total",
    "assembly_base_qname_n", "assembly_fraction", "zero_acgt_support",
    "assembly_base_no_qname_support", "alternative_exceeds_assembly",
    "alternate_ge2_qnames_both_orientations",
])


def _check_file(record, label, supplied_path=None):
    if not isinstance(record, dict) or not isinstance(record.get("path"), str) or not record["path"]:
        raise ValueError(f"Missing file binding: {label}")
    if type(record.get("bytes")) is not int or record["bytes"] < 0:
        raise ValueError(f"Invalid byte size: {label}")
    expected_sha = record.get("sha256", "")
    if not isinstance(expected_sha, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected_sha):
        raise ValueError(f"Invalid SHA-256: {label}")
    path = Path(supplied_path if supplied_path is not None else record["path"]).resolve()
    if not path.is_file() or path.stat().st_size != record["bytes"] or digest(path) != expected_sha.lower():
        raise ValueError(f"File SHA/size mismatch: {label}: {path}")
    return path


def _read_candidates(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if (not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames) or
                not set(REQUIRED_STRUCTURE_FIELDS).issubset(reader.fieldnames)):
            raise ValueError("Structures TSV requires unique columns and all base-concordance fields")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("Malformed structures TSV row")
    seen = set()
    for row in rows:
        sid = row["structure_id"]
        if (not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sid) or sid in seen or
                sid.upper().split(".")[0] in {"CON", "PRN", "AUX", "NUL", *[f"COM{i}" for i in range(1, 10)],
                                             *[f"LPT{i}" for i in range(1, 10)]}):
            raise ValueError("Unsafe or duplicate structure_id")
        seen.add(sid)
    return rows


def _require_readback_validation(receipt):
    rb = receipt.get("readback")
    if not isinstance(rb, dict) or rb.get("status") != "COMPLETE":
        raise ValueError("A completed source readback is required when an alignment is supplied")
    required = {
        "qname_validation": "PASS_contiguous_groups_no_duplicate_primary_mates",
        "boundary_pair_validation": "both_primary_mates_present_with_reciprocal_coordinates_and_strands",
        "sqlite_integrity_check": "ok",
    }
    for field, expected in required.items():
        if rb.get(field) != expected:
            raise ValueError(f"Source readback validation missing or failed: {field}")
    header = rb.get("header_validation")
    if not isinstance(header, dict) or header.get("reference_dictionary") != "PASS_all_contig_names_and_lengths":
        raise ValueError("Source readback header validation missing or failed")
    return rb


def collapse_group(observations):
    """One call per QNAME/position; strand belongs to contributing records."""
    bases = {base for base, strand in observations}
    strands = {strand for base, strand in observations}
    return (next(iter(bases)) if len(bases) == 1 else "CONFLICT",
            "both" if len(strands) == 2 else next(iter(strands)))


def _aligned_reference_pairs(cigar, start, query_length, reference_length):
    """Yield SAM query/reference coordinates without pysam's padding offset.

    Validate the entire CIGAR before yielding any evidence. P/H consume neither
    sequence; I/S consume query only, D/N reference only, M/=/X both. Padding
    must not shift either the nucleotide or its base-quality observation.
    """
    if (not isinstance(cigar, (list, tuple)) or not cigar or
            type(start) is not int or type(query_length) is not int or
            type(reference_length) is not int or start < 0 or query_length < 0 or
            reference_length < 0 or start > reference_length):
        raise ValueError("Invalid mapped CIGAR or query/reference coordinates")
    for index, item in enumerate(cigar):
        if (not isinstance(item, (list, tuple)) or len(item) != 2 or
                type(item[0]) is not int or type(item[1]) is not int or
                item[0] not in range(9) or item[1] <= 0):
            raise ValueError("CIGAR requires supported operators and positive integer lengths")
    operations = [op for op, length in cigar]
    for index, op in enumerate(operations):
        if op == 5 and index not in (0, len(operations) - 1):
            raise ValueError("Hard clipping is only permitted at CIGAR ends")
        if op == 4 and not (all(x == 5 for x in operations[:index]) or
                            all(x == 5 for x in operations[index + 1:])):
            raise ValueError("Soft clipping is only permitted at CIGAR ends, inside optional hard clipping")
    consumed_query = sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
    consumed_reference = sum(length for op, length in cigar if op in (0, 2, 3, 7, 8))
    if consumed_query != query_length:
        raise ValueError("Alignment SEQ length disagrees with CIGAR query consumption")
    if start + consumed_reference > reference_length:
        raise ValueError("CIGAR extends beyond reference length")
    qpos, rpos = 0, start
    for op, length in cigar:
        if op in (0, 7, 8):
            for offset in range(length):
                yield qpos + offset, rpos + offset
            qpos += length
            rpos += length
        elif op in (1, 4):
            qpos += length
        elif op in (2, 3):
            for offset in range(length):
                yield None, rpos + offset
            rpos += length
        # H and P consume neither query nor reference.


def count_bam(path, targets):
    """Preserve the frozen audit's primary/MAPQ20/BQ20 observation rules."""
    try:
        import pysam
    except ImportError as error:
        raise ImportError("Base concordance requires the structure extra (pysam) on the compute host") from error
    counts = {key: Counter() for key in targets}
    per_contig = defaultdict(set)
    for contig, position in targets:
        per_contig[contig].add(position)
    span = {c: (min(p), max(p) + 1) for c, p in per_contig.items()}
    ledger = Counter()
    seen = set()
    with pysam.AlignmentFile(str(path), "rb", threads=2) as alignment:
        if not alignment.is_bam:
            raise ValueError("Base concordance requires a BAM alignment, not SAM/CRAM")
        for name, records in groupby(alignment.fetch(until_eof=True), lambda a: a.query_name):
            observed = defaultdict(list)
            retained_group = False
            for record in records:
                ledger["alignment_records"] += 1
                if record.is_unmapped or record.reference_name not in span:
                    continue
                contig = record.reference_name
                low, high = span[contig]
                if record.reference_end is None or record.reference_start >= high or record.reference_end <= low:
                    continue
                if record.is_secondary or record.is_supplementary:
                    ledger["nonprimary_records_excluded"] += 1
                    continue
                if record.mapping_quality == 255 or record.mapping_quality < 20:
                    ledger["mapq_records_excluded"] += 1
                    continue
                seq, quality = record.query_sequence, record.query_qualities
                if seq is None or quality is None:
                    ledger["missing_sequence_or_quality"] += 1
                    continue
                if len(quality) != len(seq):
                    raise ValueError("Alignment quality length disagrees with SEQ")
                retained_group = True
                ledger["quality_primary_overlap_records"] += 1
                strand = "R" if record.is_reverse else "F"
                for qi, ri in _aligned_reference_pairs(record.cigartuples, record.reference_start,
                                                       len(seq), alignment.get_reference_length(contig)):
                    if ri not in per_contig[contig]:
                        continue
                    key = (contig, ri)
                    if qi is None:
                        counts[key][("read", "DEL_OR_SKIP", strand)] += 1
                    elif quality[qi] < 20:
                        counts[key][("read", "LOW_BQ", strand)] += 1
                    else:
                        base = seq[qi].upper()
                        if base not in BASES:
                            base = "N"
                        counts[key][("read", base, strand)] += 1
                        observed[key].append((base, strand))
            if retained_group:
                if name in seen:
                    raise ValueError(f"Retained QNAME groups are not contiguous: {name!r}")
                seen.add(name)
            for key, values in observed.items():
                base, strand = collapse_group(values)
                counts[key][("qname", base, strand)] += 1
    ledger["retained_unique_qnames"] = len(seen)
    return counts, dict(ledger)


def position_row(candidate, position, sequence, counts):
    reverse = candidate["orientation_to_reference"] == "-"
    assembly = sequence[position]
    start, end = int(candidate["interval_start_0based"]), int(candidate["interval_end_0based_exclusive"])
    row = {"structure_id": candidate["structure_id"], "contig": candidate["contig"],
           "genomic_position_0based": position, "full_position_0based": end - 1 - position if reverse else position - start,
           "orientation_to_reference": candidate["orientation_to_reference"], "assembled_genomic_base": assembly,
           "assembled_oriented_base": assembly.translate(COMP) if reverse else assembly}
    q = {}
    both_orientation_alternates = []
    for base in BASES:
        q[base] = sum(counts[("qname", base, s)] for s in ("F", "R", "both"))
        for strand in ("F", "R"):
            row[f"{base}_read_{strand}"] = counts[("read", base, strand)]
        for strand in ("F", "R", "both"):
            row[f"{base}_qname_{strand}"] = counts[("qname", base, strand)]
        row[f"{base}_qname_total"] = q[base]
        if base in "ACGT" and base != assembly and q[base] >= 2 and (
                counts[("qname", base, "F")] + counts[("qname", base, "both")] > 0) and (
                counts[("qname", base, "R")] + counts[("qname", base, "both")] > 0):
            both_orientation_alternates.append(base)
    for strand in ("F", "R", "both"):
        row[f"conflict_qname_{strand}"] = counts[("qname", "CONFLICT", strand)]
    for label in ("DEL_OR_SKIP", "LOW_BQ"):
        row[label + "_read_total"] = sum(counts[("read", label, s)] for s in ("F", "R"))
    denominator = sum(q[b] for b in "ACGT")
    assembly_n = q.get(assembly, 0)
    alternative_max = max(q[b] for b in "ACGT" if b != assembly)
    row.update(acgt_qname_total=denominator, assembly_base_qname_n=assembly_n,
               assembly_fraction=assembly_n / denominator if denominator and assembly in "ACGT" else "",
               zero_acgt_support=denominator == 0,
               assembly_base_no_qname_support=assembly_n == 0,
               alternative_exceeds_assembly=alternative_max > assembly_n,
               alternate_ge2_qnames_both_orientations=",".join(both_orientation_alternates))
    return row


def run_base_concordance(readback_receipt, structures, out):
    """Audit an explicitly bound completed structure run into a new directory.

    Inputs are paths. Only structure receipt schema /1 is accepted. The supplied
    structures table must have the exact bytes bound by that receipt. Original
    FASTA/BAM paths in the receipt must remain accessible on the compute host.
    """
    start_time = time.monotonic()
    out = Path(out).resolve()
    if out.exists():
        raise FileExistsError(f"Output directory already exists: {out}")
    receipt_path = Path(readback_receipt).resolve()
    source_record = file_record(receipt_path)
    marker = receipt_path.parent / "COMPLETE"
    if not marker.is_file() or not marker.read_text(encoding="utf-8").strip():
        raise ValueError("Source structure COMPLETE marker is missing or empty")
    marker_record = file_record(marker)
    source = json.loads(receipt_path.read_text(encoding="utf-8-sig"))
    if (not isinstance(source, dict) or source.get("schema_version") != "glycan-evidence.structure/1" or
            source.get("status") != "COMPLETE"):
        raise ValueError("A COMPLETE glycan-evidence.structure/1 receipt is required")
    outputs = source.get("outputs")
    if not isinstance(outputs, list):
        raise ValueError("Source receipt has no output file bindings")
    bound_tables = [item for item in outputs if isinstance(item, dict) and
                    isinstance(item.get("path"), str) and Path(item["path"]).name == "structures.tsv"]
    if len(bound_tables) != 1:
        raise ValueError("Source receipt must bind exactly one structures.tsv")
    structures_path = _check_file(bound_tables[0], "structures", structures)
    candidates = _read_candidates(structures_path)
    if type(source.get("structure_count")) is not int or source["structure_count"] != len(candidates):
        raise ValueError("Source structure_count differs from the bound structures table")
    contigs_path = _check_file(source.get("contigs"), "contigs")
    seqs = read_fasta(contigs_path)
    targets = set()
    for candidate in candidates:
        if candidate["contig"] not in seqs:
            raise ValueError("Candidate contig is absent from the bound FASTA")
        seq = seqs[candidate["contig"]]
        if hashlib.sha256(seq.encode("ascii")).hexdigest() != candidate["contig_sequence_sha256"]:
            raise ValueError("Candidate contig sequence SHA mismatch")
        if candidate["orientation_to_reference"] not in ("+", "-"):
            raise ValueError("Candidate orientation must be '+' or '-'")
        lo, hi = int(candidate["interval_start_0based"]), int(candidate["interval_end_0based_exclusive"])
        if not 0 <= lo < hi <= len(seq):
            raise ValueError("Invalid candidate interval coordinates")
        interval = seq[lo:hi]
        # The frozen audit defines its orientation transform over A/C/G/T/N.
        # Refuse unsupported ambiguity instead of emitting an incorrect base.
        if set(interval) - set(BASES):
            raise ValueError("Base-concordance intervals require A/C/G/T/N assembly bases")
        if candidate["orientation_to_reference"] == "-":
            interval = interval.translate(COMP)[::-1]
        if hashlib.sha256(interval.encode("ascii")).hexdigest() != candidate["interval_sequence_sha256"]:
            raise ValueError("Candidate interval sequence SHA mismatch")
        targets.update((candidate["contig"], p) for p in range(lo, hi))
    bam = None
    if source.get("alignment_input") is not None:
        _require_readback_validation(source)
        bam = _check_file(source["alignment_input"], "alignment_input")
    else:
        rb = source.get("readback")
        if candidates or not isinstance(rb, dict) or rb.get("status") != "NOT_REQUESTED":
            raise ValueError("Nonzero candidates require a bound, completed BAM readback; only zero candidates allow NOT_REQUESTED")
    bound_paths = [receipt_path, marker, structures_path, contigs_path] + ([bam] if bam else [])
    before = {p: (p.stat().st_size, p.stat().st_mtime_ns) for p in bound_paths}
    out.mkdir(parents=True, exist_ok=False)
    result = {
        "schema_version": SCHEMA_VERSION, "package_version": __version__, "status": "RUNNING",
        "started_utc": utc_now(), "readback_receipt": source_record, "source_complete_marker": marker_record,
        "structures": file_record(structures_path), "source_structures_binding": bound_tables[0],
        "contigs": source["contigs"], "bam": source.get("alignment_input"), "candidate_count": len(candidates),
        "bam_sha_reverified": bam is not None, "alignment_scan_performed": bool(candidates),
        "source_readback_validation": source.get("readback"), "code_sha256": digest(__file__),
        "parameters": {"min_mapq": 20, "min_baseq": 20, "mapq_255": "excluded_unknown"},
        "count_unit": "one_consensus_call_per_QNAME_per_position_not_original_molecules",
        "coordinate_orientation": "counts_and_alternate_letters_genomic_BAM; full_position_and_assembled_oriented_base_reference_oriented",
        "scope": "descriptive base concordance; no calibrated variant, identity, abundance, function or absence call; QNAME not molecule",
    }
    try:
        write_json(out / "receipt.json", result)
        counts, ledger = count_bam(bam, targets) if candidates else ({}, {})
        if any((p.stat().st_size, p.stat().st_mtime_ns) != state for p, state in before.items()):
            raise ValueError("Bound input changed during base-concordance analysis")
        summaries = []
        for candidate in candidates:
            rows = [position_row(candidate, pos, seqs[candidate["contig"]], counts[candidate["contig"], pos])
                    for pos in range(int(candidate["interval_start_0based"]), int(candidate["interval_end_0based_exclusive"]))]
            with gzip.open(out / (candidate["structure_id"] + ".positions.tsv.gz"), "wt", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, POSITION_FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            summaries.append({"structure_id": candidate["structure_id"], "positions": len(rows),
                              "zero_acgt_support_positions": sum(r["zero_acgt_support"] for r in rows),
                              "assembly_base_no_qname_support_positions": sum(r["assembly_base_no_qname_support"] for r in rows),
                              "alternative_exceeds_assembly_positions": sum(r["alternative_exceeds_assembly"] for r in rows),
                              "alternate_ge2_qnames_both_orientations_positions": sum(bool(r["alternate_ge2_qnames_both_orientations"]) for r in rows)})
        result.update(status="COMPLETE", completed_utc=utc_now(), ledger=ledger, candidates=summaries,
                      evidence_status="descriptive_concordance" if candidates else "not_applicable_no_candidates",
                      interpretation="assembly_conditional_observation" if candidates else "unresolved_no_geometrically_paired_anchors",
                      bam_unchanged_during_read=True if candidates else None,
                      alignment_revalidation="SHA_verified_source_readback_validations_bound" if bam else "not_revalidated_no_alignment_requested",
                      seconds=time.monotonic() - start_time,
                      files={p.name: {"bytes": p.stat().st_size, "sha256": digest(p)}
                             for p in sorted(out.glob("*.positions.tsv.gz"))})
        write_json(out / "receipt.json", result)
        (out / "COMPLETE").write_text(digest(out / "receipt.json") + "\n", encoding="ascii")
        return result
    except BaseException as error:
        (out / "COMPLETE").unlink(missing_ok=True)
        result.update(status="FAILED", failed_utc=utc_now(), error_type=type(error).__name__, error=str(error))
        write_json(out / "receipt.json", result)
        raise
