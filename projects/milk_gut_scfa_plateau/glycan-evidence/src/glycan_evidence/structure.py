"""Configurable two-anchor geometry and descriptive all-contig readback.

This module never emits an absence or taxonomic-carrier call. Coordinates are
zero-based, half-open; physical left/right are independent of query gene names.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import shutil
import subprocess
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

SCHEMA_VERSION = "glycan-evidence.structure/1"
BLAST_FIELDS = ("qseqid sseqid pident length qlen slen qstart qend sstart send "
                "evalue bitscore gaps qseq sseq").split()
STRUCTURE_FIELDS = [
    "structure_id", "contig", "contig_length_bp", "contig_sequence_sha256",
    "interval_start_0based", "interval_end_0based_exclusive",
    "orientation_to_reference", "left_anchor_id", "right_anchor_id",
    "left_anchor_boundary_0based", "right_anchor_boundary_0based",
    "left_anchor_inner_edge_complete", "right_anchor_inner_edge_complete",
    "inner_anchor_edges_complete", "both_anchors_full_length",
    "between_hsp_edges_bp", "between_anchor_bp", "anchor_screen_pass",
    "left_anchor_query_span_fraction", "right_anchor_query_span_fraction",
    "left_anchor_identity_percent", "right_anchor_identity_percent",
    "interval_sequence_sha256", "candidate_evidence", "structural_classification",
]


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def file_record(path):
    path = Path(path)
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": h.hexdigest()}


def sequence_sha(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def read_fasta(path):
    """Read strict, unique-ID nucleotide FASTA without optional dependencies."""
    records, name, parts = {}, None, []

    def save():
        if name is None:
            return
        sequence = "".join(parts).upper()
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"Empty or invalid nucleotide FASTA record: {name}")
        if name in records:
            raise ValueError(f"Duplicate FASTA identifier: {name}")
        records[name] = sequence

    with Path(path).open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                tokens = line[1:].split()
                if not tokens:
                    raise ValueError("Empty FASTA identifier")
                name, parts = tokens[0], []
            elif name is None:
                raise ValueError("FASTA sequence precedes its header")
            else:
                if any(x.isspace() for x in line):
                    raise ValueError(f"Whitespace inside FASTA sequence: {name}")
                parts.append(line)
    save()
    if not records:
        raise ValueError("FASTA contains no records")
    return records


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def load_config(config=None):
    if config is None:
        config = Path(__file__).parent / "presets" / "h5.json"
    if isinstance(config, (str, Path)):
        config = json.loads(Path(config).read_text(encoding="utf-8-sig"))
    else:
        config = json.loads(json.dumps(config))
    if config.get("schema_version") != "glycan-evidence.structure-config/1":
        raise ValueError("Unsupported structure configuration schema")
    anchors = config.get("anchors", [])
    if len(anchors) != 2 or len({a.get("id") for a in anchors}) != 2:
        raise ValueError("Exactly two distinct anchor IDs are required")
    if {a.get("reference_order") for a in anchors} != {0, 1}:
        raise ValueError("Anchor reference_order must be 0 and 1")
    for anchor in anchors:
        if not isinstance(anchor["id"], str) or not anchor["id"] or any(x.isspace() for x in anchor["id"]):
            raise ValueError("Anchor IDs must be nonempty FASTA identifiers")
        if anchor.get("reference_strand") not in ("+", "-"):
            raise ValueError("Each anchor needs reference_strand '+' or '-'")
    if not isinstance(config.get("max_inter_anchor_bp"), int) or config["max_inter_anchor_bp"] < 0:
        raise ValueError("max_inter_anchor_bp must be a nonnegative integer")
    screen = config["anchor_screen"]
    if not 0 <= screen["min_identity_percent"] <= 100 or not 0 < screen["min_query_span_fraction"] <= 1:
        raise ValueError("Invalid anchor screening thresholds")
    if config.get("decision_rule") != "continuous_evidence_only_no_validated_classification":
        raise ValueError("Only continuous descriptive evidence is implemented")
    rb = config["readback"]
    # Configuration schema v1 before the MAPQ audit used -k reporting. Preserve
    # that explicit effective interpretation instead of silently changing it.
    rb.setdefault("reporting_mode", "k_multi_audit")
    if rb["reporting_mode"] not in ("default_best", "k_multi_audit"):
        raise ValueError("Unknown Bowtie2 reporting_mode")
    for key in ("min_mapq", "min_baseq"):
        if not isinstance(rb.get(key), int) or not 0 <= rb[key] <= 254:
            raise ValueError(f"Invalid {key}")
    for key in ("boundary_flank_bp", "window_bp", "bowtie2_k", "max_fragment_bp"):
        if not isinstance(rb.get(key), int) or rb[key] < 1:
            raise ValueError(f"Invalid {key}")
    if sorted(rb.get("breadth_depths", [])) != [1, 5]:
        raise ValueError("This schema reports breadth at depths 1 and 5")
    blast = config["blast"]
    if blast.get("task") != "blastn" or not math.isfinite(blast["evalue"]) or blast["evalue"] <= 0:
        raise ValueError("Invalid blast task/evalue")
    if not isinstance(blast.get("word_size"), int) or blast["word_size"] < 4:
        raise ValueError("Invalid BLAST word_size")
    return config


def parse_hit(row, contigs, anchors):
    hit = dict(row)
    for key in ("length", "qlen", "slen", "qstart", "qend", "sstart", "send", "gaps"):
        hit[key] = int(hit[key])
    for key in ("pident", "bitscore", "evalue"):
        hit[key] = float(hit[key])
        if not math.isfinite(hit[key]):
            raise ValueError(f"Non-finite BLAST {key}")
    if hit["qseqid"] not in anchors or hit["sseqid"] not in contigs:
        raise ValueError("BLAST identifier not present in input FASTA")
    if hit["qlen"] != len(anchors[hit["qseqid"]]) or hit["slen"] != len(contigs[hit["sseqid"]]):
        raise ValueError("BLAST length disagrees with FASTA")
    if not 1 <= hit["qstart"] <= hit["qend"] <= hit["qlen"]:
        raise ValueError("BLAST query must be searched in the supplied plus orientation")
    if not 1 <= min(hit["sstart"], hit["send"]) <= max(hit["sstart"], hit["send"]) <= hit["slen"]:
        raise ValueError("Invalid BLAST subject coordinates")
    if not 0 <= hit["pident"] <= 100 or hit["length"] < 1:
        raise ValueError("Invalid BLAST alignment metrics")
    if len(hit["qseq"]) != hit["length"] or len(hit["sseq"]) != hit["length"]:
        raise ValueError("BLAST aligned sequences have inconsistent lengths")
    if len(hit["qseq"].replace("-", "")) != hit["qend"] - hit["qstart"] + 1:
        raise ValueError("BLAST query sequence does not match its coordinates")
    if len(hit["sseq"].replace("-", "")) != abs(hit["send"] - hit["sstart"]) + 1:
        raise ValueError("BLAST subject sequence does not match its coordinates")
    expected_query = anchors[hit["qseqid"]][hit["qstart"]-1:hit["qend"]]
    lo, hi = min(hit["sstart"], hit["send"])-1, max(hit["sstart"], hit["send"])
    subject = contigs[hit["sseqid"]][lo:hi]
    if hit["sstart"] > hit["send"]:
        subject = reverse_complement(subject)
    if hit["qseq"].replace("-", "").upper() != expected_query or hit["sseq"].replace("-", "").upper() != subject:
        raise ValueError("BLAST alignment sequence is inconsistent with supplied FASTA")
    hit.update(lo=lo, hi=hi, subject_strand="+" if hit["sstart"] <= hit["send"] else "-",
               query_span_fraction=(hit["qend"]-hit["qstart"]+1)/hit["qlen"])
    return hit


def _inner_complete(hit, anchor):
    # Reference-left anchor points inward at its genomic right edge, regardless
    # of coding strand. Reference-right anchor points at its genomic left edge.
    query_start_is_inner = (anchor["reference_order"] == 0) == (anchor["reference_strand"] == "-")
    if query_start_is_inner:
        return hit["qstart"] == 1 and hit["qseq"][0] != "-" and hit["sseq"][0] != "-"
    return hit["qend"] == hit["qlen"] and hit["qseq"][-1] != "-" and hit["sseq"][-1] != "-"


def find_structures(hits, contigs, config):
    """Return every geometrically coherent HSP pair, plus rejected pairs.

    HSPs are not merged across gaps and the top-scoring pair is not selected.
    Even below-screen or incomplete-inner-edge pairs remain in the output.
    """
    cfg = load_config(config)
    defs = sorted(cfg["anchors"], key=lambda a: a["reference_order"])
    by_id = {a["id"]: [] for a in defs}
    for hit in hits:
        if hit["qseqid"] in by_id:
            by_id[hit["qseqid"]].append(hit)
    structures, sequences, rejected, seen = [], {}, [], set()
    for a in sorted(by_id[defs[0]["id"]], key=lambda h: (h["sseqid"], h["lo"], h["hi"], h["qstart"])):
        for b in sorted(by_id[defs[1]["id"]], key=lambda h: (h["sseqid"], h["lo"], h["hi"], h["qstart"])):
            diagnostic = {"anchor_0_hit_id": a.get("hit_id", ""), "anchor_1_hit_id": b.get("hit_id", ""),
                          "anchor_0_contig": a["sseqid"], "anchor_1_contig": b["sseqid"]}
            orientations = ["+" if h["subject_strand"] == d["reference_strand"] else "-"
                            for h, d in zip((a, b), defs)]
            reason = None
            if a["sseqid"] != b["sseqid"]:
                reason = "anchors_on_different_contigs"
            elif orientations[0] != orientations[1]:
                reason = "inconsistent_anchor_orientation"
            if reason:
                rejected.append(dict(diagnostic, reason=reason))
                continue
            orientation = orientations[0]
            left, right = (a, b) if orientation == "+" else (b, a)
            left_def, right_def = defs if orientation == "+" else list(reversed(defs))
            gap = right["lo"] - left["hi"]
            if gap < 0 or gap > cfg["max_inter_anchor_bp"]:
                rejected.append(dict(diagnostic, reason="overlap_or_order_conflict" if gap < 0 else "distance_exceeds_config"))
                continue
            key = tuple((h["sseqid"], h["lo"], h["hi"], h["qstart"], h["qend"], h["subject_strand"]) for h in (a, b))
            if key in seen:
                continue
            seen.add(key)
            sid = f"structure_{len(structures)+1:06d}"
            sequence = contigs[a["sseqid"]][left["lo"]:right["hi"]]
            if orientation == "-":
                sequence = reverse_complement(sequence)
            sequences[sid] = sequence
            lc, rc = _inner_complete(left, left_def), _inner_complete(right, right_def)
            screen = all(h["pident"] >= cfg["anchor_screen"]["min_identity_percent"] and
                         h["query_span_fraction"] >= cfg["anchor_screen"]["min_query_span_fraction"] for h in (a, b))
            structures.append({
                "structure_id": sid, "contig": a["sseqid"], "contig_length_bp": len(contigs[a["sseqid"]]),
                "contig_sequence_sha256": sequence_sha(contigs[a["sseqid"]]),
                "interval_start_0based": left["lo"], "interval_end_0based_exclusive": right["hi"],
                "orientation_to_reference": orientation,
                "left_anchor_id": left["qseqid"], "right_anchor_id": right["qseqid"],
                "left_anchor_boundary_0based": left["hi"], "right_anchor_boundary_0based": right["lo"],
                "left_anchor_inner_edge_complete": lc, "right_anchor_inner_edge_complete": rc,
                "inner_anchor_edges_complete": lc and rc,
                "both_anchors_full_length": all(h["qstart"] == 1 and h["qend"] == h["qlen"] for h in (a, b)),
                "between_hsp_edges_bp": gap, "between_anchor_bp": gap if lc and rc else "",
                "anchor_screen_pass": screen,
                "left_anchor_query_span_fraction": left["query_span_fraction"],
                "right_anchor_query_span_fraction": right["query_span_fraction"],
                "left_anchor_identity_percent": left["pident"], "right_anchor_identity_percent": right["pident"],
                "interval_sequence_sha256": sequence_sha(sequence),
                "candidate_evidence": "complete_inner_edges" if lc and rc and screen else "unresolved_partial_or_below_screen_anchors",
                "structural_classification": "unresolved_no_validated_classification_rule",
            })
    return structures, sequences, rejected


def executable(name, bin_dir=None):
    candidate = str(Path(bin_dir) / name) if bin_dir else shutil.which(name)
    if not candidate or not Path(candidate).is_file():
        raise FileNotFoundError(f"Required executable not found: {name}")
    return candidate


def run_command(argv, log, receipt):
    item = {"argv": [str(x) for x in argv], "started_utc": utc_now()}
    receipt["commands"].append(item)
    with Path(log).open("w", encoding="utf-8") as handle:
        result = subprocess.run(item["argv"], stdout=handle, stderr=subprocess.STDOUT, check=False)
    item.update(returncode=result.returncode, completed_utc=utc_now())
    if result.returncode:
        raise RuntimeError(f"Command failed ({result.returncode}); see {log}")


def run_structure(contigs, anchors, out, r1=None, r2=None, bam=None, config=None, bin_dir=None, threads=1):
    """Run configurable anchor BLAST and optional descriptive competitive readback.

    ``bam`` also accepts SAM for small portable fixtures. An existing alignment
    must contain ALL assembly contigs in @SQ and contiguous QNAME groups. This
    verifies its declared reference dictionary, not its unknown mapping history.
    """
    from .readback import analyze_readback, map_all_contigs, write_readback_tables
    from . import __version__

    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    receipt = {"schema_version": SCHEMA_VERSION, "status": "RUNNING", "started_utc": utc_now(),
               "package_version": __version__, "python_version": platform.python_version(),
               "commands": [], "structural_classification": "unresolved_no_validated_classification_rule",
               "coverage_unit": "primary_mate_alignment_read_depth_not_fragment_depth",
               "scope": "assembly-conditional sequence evidence; no abundance, species-carrier or functional call"}
    try:
        if not isinstance(threads, int) or threads < 1:
            raise ValueError("threads must be a positive integer")
        if bool(r1) != bool(r2) or (bam is not None and r1 is not None):
            raise ValueError("Supply either both paired FASTQ files, or BAM/SAM, or neither")
        contigs, anchors = Path(contigs).resolve(), Path(anchors).resolve()
        cfg = load_config(config)
        cs, all_queries = read_fasta(contigs), read_fasta(anchors)
        queries = {a["id"]: all_queries[a["id"]] for a in cfg["anchors"]}
        receipt.update(contigs=file_record(contigs), anchor_input=file_record(anchors), config=cfg,
                       config_sha256=hashlib.sha256(json.dumps(cfg, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
                       code=[file_record(__file__), file_record(Path(__file__).with_name("readback.py"))],
                       contig_count=len(cs), anchor_count=len(queries))
        write_json(out / "effective_config.json", cfg)
        query_path = out / "anchor_queries.fasta"
        query_path.write_text("".join(f">{name}\n{seq}\n" for name, seq in queries.items()), encoding="ascii")
        blast = executable("blastn", bin_dir)
        version = subprocess.run([blast, "-version"], capture_output=True, text=True, check=True)
        receipt["blast_version"] = version.stdout.strip()
        raw_hits = out / "anchor_blast.tsv"
        run_command([blast, "-task", cfg["blast"]["task"], "-query", str(query_path), "-subject", str(contigs),
                     "-evalue", str(cfg["blast"]["evalue"]),
                     "-word_size", str(cfg["blast"]["word_size"]), "-num_threads", str(threads),
                     "-max_target_seqs", str(max(5, len(cs))), "-outfmt", "6 " + " ".join(BLAST_FIELDS), "-out", str(raw_hits)],
                    out / "blast.log", receipt)
        with raw_hits.open(encoding="utf-8") as handle:
            hits = [parse_hit(row, cs, queries) for row in csv.DictReader(handle, BLAST_FIELDS, delimiter="\t")]
        for n, hit in enumerate(hits, 1):
            hit["hit_id"] = f"hit_{n:06d}"
        coordinate_names = {"qstart": "query_start_1based", "qend": "query_end_1based_inclusive",
                            "sstart": "subject_alignment_start_1based", "send": "subject_alignment_end_1based_inclusive"}
        hit_rows = [{**{coordinate_names.get(k, k): v for k, v in hit.items() if k not in ("lo", "hi")},
                     "subject_start_0based": hit["lo"], "subject_end_0based_exclusive": hit["hi"]} for hit in hits]
        write_tsv(out / "anchor_hits.tsv", hit_rows, ["hit_id"] + [coordinate_names.get(k, k) for k in BLAST_FIELDS] +
                  ["subject_start_0based", "subject_end_0based_exclusive", "subject_strand", "query_span_fraction"])
        structures, sequences, rejected = find_structures(hits, cs, cfg)
        write_tsv(out / "structures.tsv", structures, STRUCTURE_FIELDS)
        write_tsv(out / "rejected_anchor_pairs.tsv", rejected,
                  ["anchor_0_hit_id", "anchor_1_hit_id", "anchor_0_contig", "anchor_1_contig", "reason"])
        (out / "anchored_sequences.fasta").write_text("".join(f">{name}\n{seq}\n" for name, seq in sequences.items()), encoding="ascii")
        readback = None
        if bam is not None or r1 is not None:
            if bam is None:
                bam = map_all_contigs(contigs, r1, r2, out, cfg["readback"], receipt, bin_dir, threads)
                receipt["alignment_provenance"] = "generated_against_every_supplied_contig"
            else:
                receipt["alignment_provenance"] = "external_alignment_reference_dictionary_verified_mapping_history_unverified"
            readback = analyze_readback(bam, cs, structures, cfg["readback"], work_dir=out)
            write_readback_tables(out, readback)
            receipt["alignment_input"] = file_record(bam)
            receipt["readback"] = readback["receipt"]
        else:
            receipt["readback"] = {"status": "NOT_REQUESTED", "interpretation": "unresolved_without_readback"}
            write_readback_tables(out, {"coverage": [], "boundaries": [], "windows": []})
        report = summarize_evidence(structures, readback)
        write_json(out / "evidence_summary.json", report)
        receipt.update(status="COMPLETE", completed_utc=utc_now(), structure_count=len(structures),
                       rejected_pair_reasons=dict(Counter(r["reason"] for r in rejected)), summary=report)
        receipt["outputs"] = [file_record(p) for p in sorted(out.iterdir()) if p.is_file() and
                              p.suffix in (".json", ".tsv", ".fasta", ".log") and p.name != "receipt.json"]
        write_json(out / "receipt.json", receipt)
        (out / "COMPLETE").write_text(utc_now() + "\n", encoding="ascii")
        return receipt
    except BaseException as exc:
        receipt.update(status="FAILED", failed_utc=utc_now(), error=f"{type(exc).__name__}: {exc}")
        write_json(out / "receipt.json", receipt)
        raise


def summarize_evidence(structures, readback=None):
    summary = {"structural_classification": "unresolved_no_validated_classification_rule", "candidate_count": len(structures),
               "decision_rule": "continuous_evidence_only_no_validated_classification", "candidates": []}
    if not structures:
        summary["observation"] = "unresolved_no_geometrically_paired_anchors"
    elif len(structures) > 1:
        summary["observation"] = "unresolved_multiple_candidate_structures"
    else:
        summary["observation"] = "single_candidate_geometry_observed"
    for structure in structures:
        item = {"structure_id": structure["structure_id"], "candidate_evidence": structure["candidate_evidence"],
                "evidence_level": "anchor_geometry_without_readback", "readback_observation": "unresolved_no_readback"}
        if readback is not None:
            sid = structure["structure_id"]
            coverage = [x for x in readback["coverage"] if x["structure_id"] == sid and
                        x["segment"] == "full_anchored_interval" and x["filter"] == "quality_filtered"]
            bounds = [x for x in readback["boundaries"] if x["structure_id"] == sid]
            if coverage:
                item.update(min_read_depth=coverage[0]["min_read_depth"], breadth_ge_1=coverage[0]["breadth_ge_1"],
                            boundary_unique_qnames=[x["unique_qnames"] for x in bounds])
                observed = coverage[0]["breadth_ge_1"] == 1 and len(bounds) == 2 and all(x["unique_qnames"] > 0 for x in bounds)
                item["readback_observation"] = ("full_coverage_and_each_boundary_observed_descriptive_only" if observed
                                                 else "unresolved_low_or_incomplete_read_evidence")
                item["evidence_level"] = ("read_supported_anchor_interval" if observed and structure["inner_anchor_edges_complete"] and structure.get("anchor_screen_pass", False)
                                          else "partial_read_evidence")
        summary["candidates"].append(item)
    return summary
