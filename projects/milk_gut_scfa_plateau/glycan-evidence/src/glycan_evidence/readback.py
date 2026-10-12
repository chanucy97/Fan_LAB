"""Strict alignment input and descriptive read depth / boundary observations.

Each accepted primary mate alignment contributes read depth, including both
overlapping mates. Boundary fragments are unique QNAMEs per candidate/side;
they are not deduplicated original DNA molecules. MAPQ is conditional on the
alignment reference and is never a species-uniqueness claim.
"""
from __future__ import annotations

import gzip
import itertools
import json
import os
import re
import sqlite3
import statistics
import subprocess
import tempfile
import time
from array import array
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

COVERAGE_FIELDS = ["structure_id", "contig", "segment", "filter", "start_0based",
                   "end_0based_exclusive", "length_bp", "breadth_ge_1", "breadth_ge_5",
                   "min_read_depth", "median_read_depth", "mean_read_depth", "coverage_unit"]
BOUNDARY_FIELDS = ["structure_id", "contig", "side", "anchor_id", "boundary_0based",
                   "complete_inner_anchor_edge", "unique_qnames", "distinct_read_start_orientation",
                   "distinct_pair_start_coordinates", "forward_primary_mate_alignments",
                   "reverse_primary_mate_alignments", "nm_histogram", "flank_bp_each_side",
                   "min_mapq", "min_baseq", "count_unit", "boundary_type"]
WINDOW_FIELDS = ["structure_id", "contig", "start_0based", "end_0based_exclusive",
                 "min_read_depth", "median_read_depth", "mean_read_depth", "breadth_ge_1", "breadth_ge_5"]
CIGAR_CODES = {c: i for i, c in enumerate("MIDNSHP=X")}


@contextmanager
def _unique_names(work_dir):
    """An exact disk-backed UNIQUE index avoids an unbounded Python name set."""
    handle = tempfile.NamedTemporaryFile(prefix="qname_", suffix=".sqlite", dir=work_dir, delete=False)
    path = Path(handle.name)
    handle.close()
    connection = sqlite3.connect(path)
    try:
        connection.execute("PRAGMA journal_mode=OFF")
        connection.execute("PRAGMA synchronous=OFF")
        connection.execute("CREATE TABLE names (name TEXT PRIMARY KEY) WITHOUT ROWID")
        yield connection
    finally:
        connection.close()
        path.unlink(missing_ok=True)


def _register_name(connection, name):
    try:
        connection.execute("INSERT INTO names(name) VALUES (?)", (name,))
    except sqlite3.IntegrityError as exc:
        raise ValueError(f"QNAME repeats in a non-contiguous group (or duplicate FASTQ fragment): {name!r}. "
                         "BAM/SAM must be name-collated, not coordinate-sorted; collate the complete alignment first.") from exc


def validate_header(header, contigs):
    """Require the complete sequence dictionary, not just selected targets."""
    sq = header.get("SQ", [])
    names = [r.get("SN") for r in sq]
    if len(names) != len(set(names)):
        raise ValueError("Duplicate @SQ sequence name")
    expected = {name: len(seq) for name, seq in contigs.items()}
    try:
        actual = {r["SN"]: int(r["LN"]) for r in sq}
    except (KeyError, ValueError, TypeError) as exc:
        raise ValueError("Malformed @SQ dictionary") from exc
    if actual != expected:
        raise ValueError("BAM/SAM @SQ names/lengths differ from the complete contig FASTA")
    return {"reference_dictionary": "PASS_all_contig_names_and_lengths",
            "declared_sort_order": header.get("HD", {}).get("SO", "unknown"),
            "contig_count": len(expected)}


def _parse_cigar(text):
    if text == "*":
        return []
    tokens = re.findall(r"([0-9]+)([MIDNSHP=X])", text)
    if not tokens or "".join(n+c for n, c in tokens) != text or any(int(n) <= 0 for n, _ in tokens):
        raise ValueError(f"Malformed CIGAR: {text}")
    return [(CIGAR_CODES[c], int(n)) for n, c in tokens]


def _sam_record(line):
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 11:
        raise ValueError("SAM alignment has fewer than 11 columns")
    qname, flag, ref, pos, mapq, cigar, mate_ref, mate_pos, tlen, seq, qual = fields[:11]
    tags = {}
    for item in fields[11:]:
        parts = item.split(":", 2)
        if len(parts) != 3:
            raise ValueError("Malformed SAM auxiliary tag")
        tags[parts[0]] = int(parts[2]) if parts[1] == "i" else parts[2]
    return {"qname": qname, "flag": int(flag), "contig": ref, "start": int(pos)-1,
            "mapq": int(mapq), "cigar": _parse_cigar(cigar),
            "mate_contig": ref if mate_ref == "=" else mate_ref, "mate_start": int(mate_pos)-1,
            "template_length": int(tlen), "seq": None if seq == "*" else seq,
            "qualities": None if qual == "*" else bytes(ord(q)-33 for q in qual), "nm": tags.get("NM", -1)}


@contextmanager
def alignment_records(path):
    """Yield (header, record iterator, reader metadata); BAM uses lazy pysam."""
    path = Path(path)
    if path.suffix.lower() == ".sam":
        with path.open(encoding="ascii") as handle:
            header = {"SQ": []}
            first = None
            for line in handle:
                if not line.strip():
                    continue
                if not line.startswith("@"):
                    first = line
                    break
                fields = line.rstrip("\r\n").split("\t")
                kind = fields[0][1:]
                if kind in ("SQ", "HD"):
                    data = dict(token.split(":", 1) for token in fields[1:])
                    if kind == "SQ":
                        header["SQ"].append(data)
                    else:
                        header["HD"] = data
            def records():
                for line in itertools.chain([first] if first else [], handle):
                    if line.strip():
                        if line.startswith("@"):
                            raise ValueError("SAM header appears after alignment records")
                        yield _sam_record(line)
            yield header, records(), {"reader": "stdlib_strict_sam"}
    else:
        try:
            import pysam
        except ImportError as exc:
            raise ImportError("BAM readback requires pysam on the Linux compute host; .sam fixtures use the standard library") from exc
        with pysam.AlignmentFile(str(path), "rb", check_sq=True, ignore_truncation=False) as handle:
            header = handle.header.to_dict()
            def records():
                for a in handle.fetch(until_eof=True):
                    yield {"qname": a.query_name, "flag": a.flag, "contig": a.reference_name,
                           "start": a.reference_start, "mapq": a.mapping_quality,
                           "cigar": a.cigartuples or [], "mate_contig": a.next_reference_name,
                           "mate_start": a.next_reference_start, "template_length": a.template_length,
                           "seq": a.query_sequence, "qualities": a.query_qualities,
                           "nm": a.get_tag("NM") if a.has_tag("NM") else -1}
            yield header, records(), {"reader": "pysam", "pysam_version": pysam.__version__}


def _record_valid(record, contigs):
    name, flag = record["qname"], record["flag"]
    if not name or name == "*" or any(c.isspace() for c in name) or flag < 0 or flag > 0xFFF:
        raise ValueError("Invalid alignment QNAME/FLAG")
    if not 0 <= record["mapq"] <= 255:
        raise ValueError("Invalid alignment MAPQ")
    if flag & 0x40 and flag & 0x80:
        raise ValueError(f"Both read1 and read2 flags set for {name}")
    if flag & 0x4:
        return
    if record["contig"] not in contigs or record["start"] < 0:
        raise ValueError("Mapped record has unknown reference or negative coordinate")
    cigar = record["cigar"]
    if not cigar or any(op not in range(9) or length <= 0 for op, length in cigar):
        raise ValueError("Mapped record has missing or invalid CIGAR")
    ref_length = sum(length for op, length in cigar if op in (0, 2, 3, 7, 8))
    query_length = sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
    if record["start"] + ref_length > len(contigs[record["contig"]]):
        raise ValueError("Alignment extends beyond reference length")
    if record["seq"] is not None and len(record["seq"]) != query_length:
        raise ValueError("Alignment SEQ length disagrees with CIGAR")
    if not (flag & (0x100 | 0x800)):
        if record["seq"] is None or record["qualities"] is None or len(record["qualities"]) != query_length:
            raise ValueError("Primary mapped alignment requires complete sequence/base qualities")
        if max(record["qualities"], default=0) > 93:
            raise ValueError("Unsupported base quality outside Phred 0..93")


def _blocks(record):
    """Merge adjacent aligned-base operators; every other operator breaks a run.

    M and =/X are alternate encodings of the same reference/query-contiguous
    alignment. Their operator boundaries are not sequence discontinuities.
    """
    rpos, qpos = record["start"], 0
    block = None
    for op, length in record["cigar"]:
        if op in (0, 7, 8):
            if block is None:
                block = (rpos, qpos)
            rpos += length
            qpos += length
        else:
            if block is not None:
                yield block[0], rpos, block[1]
                block = None
            if op in (1, 4):
                qpos += length
            elif op in (2, 3):
                rpos += length
    if block is not None:
        yield block[0], rpos, block[1]


def _quality_runs(lo, hi, query_start, qualities, minimum):
    values = qualities[query_start:query_start+hi-lo]
    if min(values) >= minimum:
        yield lo, hi
        return
    begin = None
    for offset, value in enumerate(values):
        if value >= minimum and begin is None:
            begin = lo + offset
        elif value < minimum and begin is not None:
            yield begin, lo + offset
            begin = None
    if begin is not None:
        yield begin, hi


def _depths_from_delta(delta):
    try:
        import numpy as np
    except ImportError:
        running, result = 0, array("q")
        for value in delta[:-1]:
            running += value
            result.append(running)
        if running + delta[-1] != 0:
            raise AssertionError("Unbalanced depth deltas")
        return result
    values = np.cumsum(np.frombuffer(delta, dtype=np.int64))
    if values[-1] != 0:
        raise AssertionError("Unbalanced depth deltas")
    return values[:-1]


def _stats(depth):
    if not len(depth):
        return {"breadth_ge_1": "", "breadth_ge_5": "", "min_read_depth": "", "median_read_depth": "", "mean_read_depth": ""}
    try:
        import numpy as np
    except ImportError:
        return {"breadth_ge_1": sum(x >= 1 for x in depth)/len(depth), "breadth_ge_5": sum(x >= 5 for x in depth)/len(depth),
                "min_read_depth": int(min(depth)), "median_read_depth": float(statistics.median(depth)), "mean_read_depth": sum(depth)/len(depth)}
    depth = np.asarray(depth)
    return {"breadth_ge_1": float(np.mean(depth >= 1)), "breadth_ge_5": float(np.mean(depth >= 5)),
            "min_read_depth": int(depth.min()), "median_read_depth": float(np.median(depth)), "mean_read_depth": float(depth.mean())}


def analyze_records(records, header, contigs, structures, config, work_dir=None):
    """Pure core for SAM/BAM and tiny test fixtures; consume every input record."""
    header_result = validate_header(header, contigs)
    selected = {s["contig"] for s in structures}
    for s in structures:
        if s["contig"] not in contigs:
            raise ValueError("Structure contig not in FASTA")
        bounds = [s["interval_start_0based"], s["left_anchor_boundary_0based"],
                  s["right_anchor_boundary_0based"], s["interval_end_0based_exclusive"]]
        if not 0 <= bounds[0] <= bounds[1] <= bounds[2] <= bounds[3] <= len(contigs[s["contig"]]):
            raise ValueError("Invalid structure interval coordinates")
    deltas = {label: {name: array("q", [0]) * (len(contigs[name])+1) for name in selected}
              for label in ("quality_filtered", "all_primary_base_quality_filtered")}
    boundary_info = {}
    by_contig = {}
    for s in structures:
        for side in ("left", "right"):
            key = (s["structure_id"], side)
            item = {"structure": s, "side": side, "position": s[f"{side}_anchor_boundary_0based"],
                    "qnames": 0, "starts": set(), "pair_starts": set(), "forward": 0, "reverse": 0, "nm": Counter()}
            boundary_info[key] = item
            by_contig.setdefault(s["contig"], []).append((key, item))
    counts, sql_seconds = Counter(), 0.0
    with _unique_names(work_dir) as names:
        for name, group in itertools.groupby(records, key=lambda r: r["qname"]):
            before = time.perf_counter()
            _register_name(names, name)
            sql_seconds += time.perf_counter() - before
            counts["qname_groups"] += 1
            primary_mates, supported_boundaries, primary_records = set(), set(), {}
            for record in group:
                counts["alignment_records"] += 1
                _record_valid(record, contigs)
                flag = record["flag"]
                if flag & (0x100 | 0x800):
                    counts["secondary_or_supplementary_records"] += 1
                    continue
                # Distinct mates in one group are allowed, duplicate primary
                # mate records are not silently counted as separate fragments.
                mate = 1 if flag & 0x40 else 2 if flag & 0x80 else 0
                if mate in primary_mates or (mate == 0 and primary_mates) or (mate != 0 and 0 in primary_mates):
                    raise ValueError(f"Duplicate or incompatible primary mate records for QNAME {name!r}")
                primary_mates.add(mate)
                primary_records[mate] = record
                if flag & 0x4:
                    counts["unmapped_primary_records"] += 1
                    continue
                counts["mapped_primary_mate_alignments"] += 1
                if record["mapq"] == 255:
                    counts["mapped_primary_mapq_unknown255"] += 1
                if flag & 0x400:
                    counts["duplicate_flag_primary_alignments_included"] += 1
                c = record["contig"]
                if c not in selected:
                    continue
                counts["selected_contig_primary_mate_alignments"] += 1
                if record["mapq"] == 255:
                    counts["selected_contig_primary_mapq_unknown255"] += 1
                quality_pass = record["mapq"] >= config["min_mapq"] and record["mapq"] != 255
                blocks = list(_blocks(record))
                for lo, hi, qo in blocks:
                    for start, end in _quality_runs(lo, hi, qo, record["qualities"], config["min_baseq"]):
                        d = deltas["all_primary_base_quality_filtered"][c]
                        d[start] += 1
                        d[end] -= 1
                        if quality_pass:
                            d = deltas["quality_filtered"][c]
                            d[start] += 1
                            d[end] -= 1
                proper = (flag & 0x1 and flag & 0x2 and not flag & 0x8 and
                          record["mate_contig"] == c and record["mate_start"] >= 0 and
                          bool(flag & 0x10) != bool(flag & 0x20))
                if not quality_pass or not proper:
                    continue
                flank = config["boundary_flank_bp"]
                for key, info in by_contig.get(c, []):
                    boundary = info["position"]
                    for lo, hi, qo in blocks:
                        if lo <= boundary-flank and hi >= boundary+flank:
                            qstart = qo + boundary-flank-lo
                            if min(record["qualities"][qstart:qstart+2*flank]) >= config["min_baseq"]:
                                supported_boundaries.add(key)
                                reverse = bool(flag & 0x10)
                                info["starts"].add((record["start"], reverse))
                                info["pair_starts"].add(tuple(sorted((record["start"], record["mate_start"]))))
                                info["reverse" if reverse else "forward"] += 1
                                info["nm"][record["nm"]] += 1
                            break
            if supported_boundaries:
                # A proper-pair flag alone is not evidence that the partner is
                # actually present in an externally supplied alignment file.
                if set(primary_records) != {1, 2}:
                    raise ValueError(f"Boundary-supporting QNAME {name!r} lacks its primary mate; supply the complete name-collated BAM/SAM")
                one, two = primary_records[1], primary_records[2]
                if (any(not r["flag"] & 0x1 or not r["flag"] & 0x2 or r["flag"] & (0x4 | 0x8) for r in (one, two)) or
                    one["mate_contig"] != two["contig"] or two["mate_contig"] != one["contig"] or
                    one["mate_start"] != two["start"] or two["mate_start"] != one["start"] or
                    bool(one["flag"] & 0x20) != bool(two["flag"] & 0x10) or
                    bool(two["flag"] & 0x20) != bool(one["flag"] & 0x10)):
                    raise ValueError(f"Boundary-supporting QNAME {name!r} has inconsistent primary mate coordinates/strands")
            for key in supported_boundaries:
                boundary_info[key]["qnames"] += 1
            if counts["qname_groups"] % 50000 == 0:
                before = time.perf_counter()
                names.commit()
                sql_seconds += time.perf_counter() - before
        before = time.perf_counter()
        names.commit()
        integrity = names.execute("PRAGMA integrity_check").fetchone()[0]
        sql_seconds += time.perf_counter() - before
        if integrity != "ok":
            raise RuntimeError("QNAME SQLite integrity_check failed")
    depths = {label: {c: _depths_from_delta(d) for c, d in values.items()} for label, values in deltas.items()}
    coverage, windows, boundaries = [], [], []
    for s in structures:
        c = s["contig"]
        for segment, start, end in (("full_anchored_interval", s["interval_start_0based"], s["interval_end_0based_exclusive"]),
                                     ("between_anchors" if s["inner_anchor_edges_complete"] else "between_hsp_edges",
                                      s["left_anchor_boundary_0based"], s["right_anchor_boundary_0based"])):
            for label, values in depths.items():
                coverage.append(dict(structure_id=s["structure_id"], contig=c, segment=segment, filter=label,
                                     start_0based=start, end_0based_exclusive=end, length_bp=end-start,
                                     coverage_unit="primary_mate_alignment_read_depth", **_stats(values[c][start:end])))
        for start in range(s["interval_start_0based"], s["interval_end_0based_exclusive"], config["window_bp"]):
            end = min(start+config["window_bp"], s["interval_end_0based_exclusive"])
            windows.append(dict(structure_id=s["structure_id"], contig=c, start_0based=start,
                                end_0based_exclusive=end, **_stats(depths["quality_filtered"][c][start:end])))
    for (sid, side), info in boundary_info.items():
        s = info["structure"]
        complete = s[f"{side}_anchor_inner_edge_complete"]
        boundaries.append(dict(structure_id=sid, contig=s["contig"], side=side,
                               anchor_id=s[f"{side}_anchor_id"], boundary_0based=info["position"],
                               complete_inner_anchor_edge=complete, unique_qnames=info["qnames"],
                               distinct_read_start_orientation=len(info["starts"]),
                               distinct_pair_start_coordinates=len(info["pair_starts"]),
                               forward_primary_mate_alignments=info["forward"], reverse_primary_mate_alignments=info["reverse"],
                               nm_histogram=json.dumps(dict(sorted(info["nm"].items())), sort_keys=True),
                               flank_bp_each_side=config["boundary_flank_bp"], min_mapq=config["min_mapq"], min_baseq=config["min_baseq"],
                               count_unit="unique_QNAME_within_this_candidate_boundary_not_UMI_molecules",
                               boundary_type="complete_inner_anchor_edge" if complete else "partial_HSP_edge_not_gene_boundary"))
    receipt = {"status": "COMPLETE", "counts": dict(counts), "header_validation": header_result,
               "qname_validation": "PASS_contiguous_groups_no_duplicate_primary_mates",
               "boundary_pair_validation": "both_primary_mates_present_with_reciprocal_coordinates_and_strands",
               "sqlite_integrity_check": integrity, "sqlite_validation_seconds": round(sql_seconds, 6),
               "qname_assumption": "Input QNAME must uniquely identify a sequencing fragment across the entire file; indistinguishable reused names cannot establish molecule independence.",
               "coverage_unit": "primary_mate_alignment_read_depth; overlapping mates each contribute",
               "mapq_255": "excluded_from_quality_filtered_evidence_MAPQ_unknown",
               "duplicate_flag": "included_and_counted; no UMI or PCR molecule deduplication",
               "scope": "conditional on complete supplied assembly; no species uniqueness or biological absence claim"}
    return {"coverage": coverage, "boundaries": boundaries, "windows": windows, "receipt": receipt}


def analyze_readback(path, contigs, structures, config, work_dir=None):
    with alignment_records(path) as (header, records, reader):
        result = analyze_records(records, header, contigs, structures, config, work_dir)
    result["receipt"].update(reader)
    return result


def write_readback_tables(out, result):
    from .structure import write_tsv
    write_tsv(Path(out) / "coverage.tsv", result["coverage"], COVERAGE_FIELDS)
    write_tsv(Path(out) / "boundary_support.tsv", result["boundaries"], BOUNDARY_FIELDS)
    write_tsv(Path(out) / "coverage_windows.tsv", result["windows"], WINDOW_FIELDS)


def _fastq_records(path, mate):
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    with opener(path, "rt", encoding="ascii", newline="") as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq, plus, qual = (handle.readline().rstrip("\r\n") for _ in range(3))
            if not header.startswith("@") or not plus.startswith("+") or not seq or len(seq) != len(qual):
                raise ValueError(f"Malformed/truncated FASTQ in {path}")
            if set(seq.upper()) - set("ACGTRYSWKMBDHVN") or any(not 33 <= ord(c) <= 126 for c in qual):
                raise ValueError(f"Invalid FASTQ sequence/quality in {path}")
            tokens = header[1:].split()
            if not tokens:
                raise ValueError("Empty FASTQ QNAME")
            name = tokens[0]
            if name.endswith(("/1", "/2")):
                if not name.endswith(f"/{mate}"):
                    raise ValueError("FASTQ mate suffix disagrees with its input file")
                name = name[:-2]
            if len(tokens) > 1 and tokens[1].startswith(("1:", "2:")):
                if not tokens[1].startswith(f"{mate}:"):
                    raise ValueError("FASTQ CASAVA mate label disagrees with its input file")
            if not name or name == "*":
                raise ValueError("Invalid FASTQ fragment name")
            yield name


def validate_read_pairs(r1, r2, work_dir=None):
    before = time.perf_counter()
    r1, r2 = Path(r1).resolve(), Path(r2).resolve()
    if r1 == r2 or os.path.samefile(r1, r2):
        raise ValueError("R1 and R2 must be distinct files")
    count = 0
    with _unique_names(work_dir) as names:
        for left, right in itertools.zip_longest(_fastq_records(r1, 1), _fastq_records(r2, 2)):
            if left is None or right is None or left != right:
                raise ValueError("Paired FASTQ names/counts differ")
            _register_name(names, left)
            count += 1
            if count % 50000 == 0:
                names.commit()
        names.commit()
        integrity = names.execute("PRAGMA integrity_check").fetchone()[0]
        if integrity != "ok":
            raise RuntimeError("FASTQ QNAME SQLite integrity_check failed")
    if not count:
        raise ValueError("Paired FASTQ input is empty")
    return {"status": "PASS", "paired_fragments": count, "both_streams_complete": True,
            "unique_names": True, "sqlite_integrity_check": integrity,
            "validation_seconds": round(time.perf_counter()-before, 6)}


def map_all_contigs(contigs, r1, r2, out, config, receipt, bin_dir=None, threads=1):
    from .structure import executable, file_record, run_command, utc_now
    out = Path(out)
    r1, r2 = Path(r1).resolve(), Path(r2).resolve()
    receipt["fastq_validation"] = validate_read_pairs(r1, r2, out)
    receipt["reads"] = [file_record(r1), file_record(r2)]
    build, bowtie, samtools = [executable(name, bin_dir) for name in ("bowtie2-build", "bowtie2", "samtools")]
    receipt["mapping_versions"] = {}
    for name, command in (("bowtie2", bowtie), ("samtools", samtools)):
        result = subprocess.run([command, "--version"], capture_output=True, text=True, check=True)
        receipt["mapping_versions"][name] = result.stdout.splitlines()[0]
    index = out / "all_contigs_index"
    run_command([build, "--threads", str(threads), str(contigs), str(index)], out / "index.log", receipt)
    bam = out / "all_contigs_readback.bam"
    mode = config.get("reporting_mode", "k_multi_audit")
    if mode not in ("default_best", "k_multi_audit"):
        raise ValueError("Unknown Bowtie2 reporting_mode")
    reporting_args = ["-k", str(config["bowtie2_k"])] if mode == "k_multi_audit" else []
    argv = [bowtie, "--very-sensitive", "-N", "1", "-L", "20", "--score-min", "L,0,-0.4", *reporting_args,
            "--no-mixed", "--no-discordant", "-I", "0",
            "-X", str(config["max_fragment_bp"]), "--no-unal", "-p", str(threads), "-x", str(index),
            "-1", str(r1), "-2", str(r2)]
    converter = [samtools, "view", "-b", "-o", str(bam), "-"]
    entry = {"argv": argv, "pipe_to": converter, "started_utc": utc_now()}
    receipt["commands"].append(entry)
    receipt["mapping_reporting_mode"] = mode
    receipt["mapping_quality_scope"] = "MAPQ255 unknown and excluded; available MAPQ remains aligner-estimated and assembly-conditional, not experimentally calibrated"
    with (out / "bowtie2.log").open("w") as log, (out / "samtools.log").open("w") as convert_log:
        producer = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=log)
        try:
            consumer = subprocess.Popen(converter, stdin=producer.stdout, stderr=convert_log)
        except BaseException:
            producer.stdout.close()
            producer.terminate()
            producer.wait()
            raise
        producer.stdout.close()
        try:
            consumer_code, producer_code = consumer.wait(), producer.wait()
        except BaseException:
            consumer.terminate()
            producer.terminate()
            consumer.wait()
            producer.wait()
            raise
    entry.update(returncode=producer_code, pipe_returncode=consumer_code, completed_utc=utc_now())
    if consumer_code or producer_code:
        raise RuntimeError(f"Competitive readback failed: bowtie2={producer_code}, samtools={consumer_code}")
    run_command([samtools, "quickcheck", "-v", str(bam)], out / "bam_quickcheck.log", receipt)
    receipt["mapping_reference_scope"] = "every_contig_in_input_FASTA; no selected-target-only index"
    return bam
