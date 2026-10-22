"""Phenotype-blind, source-bound inventory of existing metagenomic assemblies.

No scheduler, networking, FASTQ extraction, phenotype reader, or rule tuning.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import re
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
DICT_FILES = ["FROZEN.json", "TRAITS.tsv", "trait_proteins.faa", "trait_cds.fna", "MODULE_RULES.tsv"]
TB_FIELDS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend", "sstart", "send", "sframe", "evalue", "bitscore", "gaps", "gapopen", "btop"]
ORF_FIELDS = ["gene_id", "sample_id", "contig_id", "start1", "end1", "strand", "phase", "prodigal_id", "partial", "start_type", "cds_length", "protein_length", "cds_sha256", "protein_sha256", "start_codon", "stop_codon", "in_frame", "unambiguous_cds", "valid_start", "start_annotation_matches", "terminal_stop", "internal_stop", "translation_matches", "complete_predicted_orf", "integrity_reason"]
TRAIT_FIELDS = ["sample_id", "trait_id", "status", "complete_gene_ids", "partial_or_unresolved_gene_ids", "ambiguous_gene_ids", "protein_hsp_count", "translated_hsp_count", "translated_high_identity_hsp_count", "assembly_state", "annotation_conflict", "negative_is_absence", "functional_activity_asserted"]
MODULE_FIELDS = ["sample_id", "module_id", "status", "required_traits", "complete_same_contig_ids", "candidate_same_contig_ids", "complete_traits_in_sample", "partial_or_unresolved_traits_in_sample", "missing_candidate_traits", "same_strain_asserted", "complete_utilization_asserted"]
CONTIG_MODULE_FIELDS = ["sample_id", "module_id", "contig_id", "required_traits", "complete_traits", "partial_or_unresolved_traits", "gene_ids", "complete_set_on_contig", "candidate_set_on_contig", "assembly_colocation_only"]
TB_SUMMARY_FIELDS = ["trait_id", "trait_role", "contig_id", "strand", "hsp_count", "frames", "query_union_aa", "query_union_percent", "high_identity_hsp_count", "high_identity_query_union_aa", "high_identity_query_union_percent", "subject_min1", "subject_max1", "subject_span_nt", "largest_uncovered_subject_gap_nt", "hsp_row_numbers", "single_orf_asserted"]
PROTEIN_CONTEXT_FIELDS = ["sample_id", "gene_id", "contig_id", "start1", "end1", "strand", "trait_id", "trait_role", "single_direction_pass"]
ASSIGNMENT_FIELDS = ["gene_id", "assembly_accession", "contig_id", "start0", "end0", "strand", "locus_tag", "old_locus_tag", "protein_id", "trait_id", "trait_role", "mapping_status", "exact_seed_ids", "seed_near_best_ids", "within_genome_near_best_gene_ids", "competitor_ambiguity", "forward_identity", "forward_query_coverage", "forward_subject_coverage", "forward_bitscore", "reverse_identity", "reverse_query_coverage", "reverse_subject_coverage", "reverse_bitscore", "reverse_pass", "reference_instrain_gene_exported", "sequence_sha256", "functional_equivalence_asserted", "claim_ceiling"]
BASES = "TCAG"
CODON_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE11 = dict(zip((a+b+c for a in BASES for b in BASES for c in BASES), CODON_AA))
START11 = {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def seqsha(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def record(path):
    path = Path(path).resolve()
    before = path.stat()
    digest = sha(path)
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"file changed during hashing: {path}")
    return {"path": str(path), "bytes": after.st_size, "sha256": digest}


def verify_record(entry):
    if record(entry["path"]) != {k: entry[k] for k in ("path", "bytes", "sha256")}:
        raise ValueError(f"bound file changed: {entry['path']}")


def rows(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, data, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(data)


def fasta(path):
    sequences, headers = {}, {}
    ident, chunks = None, []
    def save():
        if ident is not None:
            sequence = "".join(chunks).upper()
            if not sequence:
                raise ValueError(f"empty FASTA record: {ident}")
            sequences[ident] = sequence
    with Path(path).open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                header = line[1:]
                ident = header.split()[0]
                if ident in headers:
                    raise ValueError(f"duplicate FASTA id: {ident}")
                headers[ident], chunks = header, []
            elif ident is None:
                raise ValueError("FASTA sequence before header")
            else:
                chunks.append(line)
    save()
    return sequences, headers


def write_fasta(path, sequences):
    with Path(path).open("w", encoding="ascii") as handle:
        for ident, sequence in sequences.items():
            handle.write(f">{ident}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start+80] + "\n")


def attributes(text):
    result = {}
    for item in text.strip().strip(";").split(";"):
        if item:
            key, value = item.split("=", 1)
            if key in result:
                raise ValueError(f"duplicate attribute: {key}")
            result[key] = value
    return result


def revcomp(sequence):
    if not set(sequence) <= set("ACGTRYMKSWBDHVN"):
        raise ValueError("non-IUPAC nucleotide sequence")
    return sequence.translate(str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN"))[::-1]


def translate11(cds, initiator):
    peptide = "".join(CODE11.get(cds[i:i+3], "X") for i in range(0, len(cds)-2, 3))
    if initiator and cds[:3] in START11 and peptide:
        peptide = "M" + peptide[1:]
    return peptide


def inspect_orf(cds, peptide, partial, start_type):
    if partial not in ("00", "01", "10", "11"):
        raise ValueError(f"invalid Prodigal partial flag: {partial}")
    in_frame = len(cds) % 3 == 0
    start, stop = cds[:3], cds[-3:]
    # Edge is a boundary continuation, not an initiator, even for TTG/GTG.
    translated = translate11(cds, start_type != "Edge")
    terminal = in_frame and stop in ("TAA", "TAG", "TGA")
    normalized = translated[:-1] if terminal else translated
    supplied = peptide[:-1] if peptide.endswith("*") else peptide
    internal = "*" in normalized
    checks = {"in_frame": in_frame, "unambiguous_cds": set(cds) <= set("ACGT"),
              "valid_start": start in START11 and start_type != "Edge", "start_annotation_matches": start_type == "Edge" or start_type == start, "terminal_stop": terminal,
              "internal_stop": internal, "translation_matches": normalized == supplied}
    reasons = []
    if partial != "00":
        reasons.append("predicted_partial_" + partial)
    for key in ("in_frame", "unambiguous_cds", "valid_start", "start_annotation_matches", "terminal_stop", "translation_matches"):
        if not checks[key]:
            reasons.append(key + "_false")
    if internal:
        reasons.append("internal_stop")
    return {"start_codon": start, "stop_codon": stop, **checks,
            "complete_predicted_orf": not reasons, "integrity_reason": ";".join(reasons) or "complete_predicted_coding_sequence"}


def parse_predictions(assembly, protein_path, cds_path, gff_path, sample):
    proteins_raw, headers = fasta(protein_path)
    cds, cds_headers = fasta(cds_path)
    if set(proteins_raw) != set(cds):
        raise ValueError("Prodigal protein/CDS identifiers disagree")
    gff = {}
    with Path(gff_path).open(encoding="ascii") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                raise ValueError("invalid Prodigal GFF CDS")
            attr = attributes(fields[8])
            key = attr["ID"]
            if key in gff or fields[6] not in ("+", "-") or fields[7] != "0":
                raise ValueError("duplicate GFF ID or unexpected strand/phase")
            gff[key] = (fields, attr)
    out, genes, proteins, used = [], {}, {}, set()
    for ident, raw_peptide in proteins_raw.items():
        pieces = headers[ident].split(" # ")
        if len(pieces) != 5 or cds_headers[ident] != headers[ident]:
            raise ValueError("Prodigal header format or protein/CDS metadata mismatch")
        start, end, strand = map(int, pieces[1:4])
        attr = attributes(pieces[4])
        fields, gattr = gff[attr["ID"]]
        contig = fields[0]
        if attr["ID"] in used or contig not in assembly or not 1 <= start <= end <= len(assembly[contig]):
            raise ValueError("invalid or duplicated predicted coordinates")
        used.add(attr["ID"])
        if (start, end, strand) != (int(fields[3]), int(fields[4]), 1 if fields[6] == "+" else -1):
            raise ValueError("FASTA/GFF coordinates disagree")
        if any(gattr.get(key) != value for key, value in attr.items()):
            raise ValueError("FASTA/GFF attributes disagree")
        genomic = assembly[contig][start-1:end]
        if strand == -1:
            genomic = revcomp(genomic)
        if genomic != cds[ident]:
            raise ValueError("predicted CDS differs from oriented assembly interval")
        peptide = raw_peptide[:-1] if raw_peptide.endswith("*") else raw_peptide
        if not peptide or not re.fullmatch(r"[A-Z*]+", peptide):
            raise ValueError("invalid predicted peptide")
        proteins[ident] = peptide
        integrity = inspect_orf(genomic, raw_peptide, attr["partial"], attr["start_type"])
        out.append({"gene_id": ident, "sample_id": sample, "contig_id": contig, "start1": start, "end1": end,
                    "strand": strand, "phase": 0, "prodigal_id": attr["ID"], "partial": attr["partial"], "start_type": attr["start_type"],
                    "cds_length": len(genomic), "protein_length": len(peptide), "cds_sha256": seqsha(genomic),
                    "protein_sha256": seqsha(peptide), **integrity})
        genes[ident] = {"assembly_accession": sample, "contig_id": contig, "start0": start-1, "end0": end, "strand": strand,
                        "locus_tag": "", "old_locus_tag": "", "protein_id": ident, "instrain_gene_exported": "False"}
    if used != set(gff):
        raise ValueError("GFF/protein feature sets disagree")
    return out, genes, proteins


def union_length(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(end, merged[-1][1])
        else:
            merged.append([start, end])
    return sum(end-start+1 for start, end in merged)


def parse_tblastn(path, seeds, assembly, traits, mapping_rules):
    hits, grouped = [], defaultdict(list)
    for line_no, line in enumerate(Path(path).read_text(encoding="ascii").splitlines(), 1):
        parts = line.split("\t")
        if len(parts) != len(TB_FIELDS):
            raise ValueError(f"invalid tblastn field count at {line_no}")
        row = dict(zip(TB_FIELDS, parts))
        for field in ("length", "qlen", "slen", "qstart", "qend", "sstart", "send", "sframe", "gaps", "gapopen"):
            row[field] = int(row[field])
        for field in ("pident", "evalue", "bitscore"):
            row[field] = float(row[field])
            if not math.isfinite(row[field]):
                raise ValueError("non-finite tblastn metric")
        if row["qseqid"] not in seeds or row["sseqid"] not in assembly:
            raise ValueError("unbound tblastn sequence")
        if row["qlen"] != len(seeds[row["qseqid"]]) or row["slen"] != len(assembly[row["sseqid"]]):
            raise ValueError("tblastn length mismatch")
        if not (1 <= row["qstart"] <= row["qend"] <= row["qlen"] and 1 <= min(row["sstart"], row["send"]) <= max(row["sstart"], row["send"]) <= row["slen"]):
            raise ValueError("invalid tblastn coordinates")
        if row["sframe"] not in (-3, -2, -1, 1, 2, 3) or (row["send"] >= row["sstart"]) != (row["sframe"] > 0):
            raise ValueError("tblastn frame/strand disagreement")
        if not (0 <= row["pident"] <= 100 and row["evalue"] >= 0 and row["bitscore"] > 0 and row["length"] > 0 and row["gaps"] >= 0 and row["gapopen"] >= 0):
            raise ValueError("invalid tblastn metrics")
        if row["qend"] - row["qstart"] + 1 > row["length"] or not row["btop"]:
            raise ValueError("invalid tblastn span/BTOP")
        row["hsp_row_number"] = line_no
        row["trait_role"] = traits[row["qseqid"]]["role"]
        row["strand"] = 1 if row["sframe"] > 0 else -1
        row["query_coverage_percent"] = 100 * (row["qend"]-row["qstart"]+1) / row["qlen"]
        row["high_identity_hsp"] = row["pident"] >= mapping_rules["minimum_identity_percent"] and row["evalue"] <= mapping_rules["maximum_evalue"]
        hits.append(row)
        grouped[(row["qseqid"], row["sseqid"], row["strand"])].append(row)
    summary = []
    for (seed, contig, strand), group in sorted(grouped.items()):
        strong = [r for r in group if r["high_identity_hsp"]]
        coverage = union_length((r["qstart"], r["qend"]) for r in group)
        strong_cov = union_length((r["qstart"], r["qend"]) for r in strong)
        intervals = sorted((min(r["sstart"], r["send"]), max(r["sstart"], r["send"])) for r in group)
        edge, gaps = intervals[0][1], []
        for start, end in intervals[1:]:
            gaps.append(max(0, start-edge-1))
            edge = max(edge, end)
        first, last = intervals[0][0], max(x[1] for x in intervals)
        summary.append({"trait_id": seed, "trait_role": traits[seed]["role"], "contig_id": contig, "strand": strand,
                        "hsp_count": len(group), "frames": ";".join(map(str, sorted({r['sframe'] for r in group}))),
                        "query_union_aa": coverage, "query_union_percent": 100*coverage/len(seeds[seed]),
                        "high_identity_hsp_count": len(strong), "high_identity_query_union_aa": strong_cov,
                        "high_identity_query_union_percent": 100*strong_cov/len(seeds[seed]),
                        "subject_min1": first, "subject_max1": last, "subject_span_nt": last-first+1,
                        "largest_uncovered_subject_gap_nt": max(gaps, default=0),
                        "hsp_row_numbers": ";".join(str(r["hsp_row_number"]) for r in group), "single_orf_asserted": False})
    return hits, summary


def evidence_tables(sample, assembly_state, traits, modules, orfs, assignments, forward, reverse, translated):
    orf_by_id = {r["gene_id"]: r for r in orfs}
    complete, partial, ambiguous = defaultdict(list), defaultdict(list), defaultdict(list)
    weak = Counter(r["sseqid"] for r in forward) + Counter(r["qseqid"] for r in reverse)
    for hit in assignments:
        if hit["trait_role"] != "positive_enzyme_seed":
            continue
        candidate = hit["mapping_status"] in ("exact_seed_peptide", "high_similarity_candidate", "ambiguous_high_similarity", "ambiguous_identical_seed_peptides")
        if candidate and hit["competitor_ambiguity"]:
            ambiguous[hit["trait_id"]].append(hit["gene_id"])
        elif candidate:
            target = complete if orf_by_id[hit["gene_id"]]["complete_predicted_orf"] else partial
            target[hit["trait_id"]].append(hit["gene_id"])
    tb = Counter(r["qseqid"] for r in translated)
    tbstrong = Counter(r["qseqid"] for r in translated if r["high_identity_hsp"])
    trait_rows = []
    for ident, trait in sorted(traits.items()):
        if trait["role"] != "positive_enzyme_seed":
            continue
        if assembly_state == "empty_assembly":
            status = "assembly_not_evaluable"
        elif complete[ident]:
            status = "complete_ORF_candidate"
        elif partial[ident]:
            status = "partial_or_unresolved_ORF_candidate"
        elif ambiguous[ident]:
            status = "homology_ambiguous_candidate"
        elif tbstrong[ident]:
            status = "translated_candidate_prediction_unresolved"
        elif weak[ident] or tb[ident]:
            status = "weak_hit_only"
        else:
            status = "not_detected_in_existing_assembly"
        trait_rows.append({"sample_id": sample, "trait_id": ident, "status": status,
                           "complete_gene_ids": ";".join(complete[ident]), "partial_or_unresolved_gene_ids": ";".join(partial[ident]),
                           "ambiguous_gene_ids": ";".join(ambiguous[ident]), "protein_hsp_count": weak[ident],
                           "translated_hsp_count": tb[ident], "translated_high_identity_hsp_count": tbstrong[ident],
                           "assembly_state": assembly_state, "annotation_conflict": trait["annotation_conflict"],
                           "negative_is_absence": False, "functional_activity_asserted": False})
    module_rows, contig_rows = [], []
    for module in modules:
        required = set(module["required_all"].split(";"))
        full_traits = required & {t for t, ids in complete.items() if ids}
        partial_traits = required & {t for t, ids in partial.items() if ids}
        by_contig = defaultdict(lambda: {"complete": set(), "partial": set(), "genes": []})
        for kind, source in (("complete", complete), ("partial", partial)):
            for trait in sorted(required):
                for gene in source[trait]:
                    group = by_contig[orf_by_id[gene]["contig_id"]]
                    group[kind].add(trait)
                    group["genes"].append(gene)
        full_contigs, candidate_contigs = [], []
        for contig, group in sorted(by_contig.items()):
            full = required <= group["complete"]
            candidate = required <= group["complete"] | group["partial"]
            if full:
                full_contigs.append(contig)
            if candidate:
                candidate_contigs.append(contig)
            contig_rows.append({"sample_id": sample, "module_id": module["module_id"], "contig_id": contig,
                                "required_traits": ";".join(sorted(required)), "complete_traits": ";".join(sorted(group["complete"])),
                                "partial_or_unresolved_traits": ";".join(sorted(group["partial"])), "gene_ids": ";".join(sorted(group["genes"])),
                                "complete_set_on_contig": full, "candidate_set_on_contig": candidate, "assembly_colocation_only": True})
        if assembly_state == "empty_assembly":
            status = "assembly_not_evaluable"
        elif full_contigs:
            status = "complete_candidates_same_assembly_contig"
        elif required <= full_traits:
            status = "complete_candidates_dispersed_in_sample"
        elif candidate_contigs:
            status = "partial_or_unresolved_candidates_same_contig"
        elif required <= full_traits | partial_traits:
            status = "partial_or_unresolved_candidates_dispersed_in_sample"
        else:
            status = "candidate_set_not_resolved"
        module_rows.append({"sample_id": sample, "module_id": module["module_id"], "status": status,
                            "required_traits": ";".join(sorted(required)), "complete_same_contig_ids": ";".join(full_contigs),
                            "candidate_same_contig_ids": ";".join(candidate_contigs), "complete_traits_in_sample": ";".join(sorted(full_traits)),
                            "partial_or_unresolved_traits_in_sample": ";".join(sorted(partial_traits)),
                            "missing_candidate_traits": ";".join(sorted(required-full_traits-partial_traits)),
                            "same_strain_asserted": False, "complete_utilization_asserted": False})
    return trait_rows, module_rows, contig_rows


def load_mapper(path):
    spec = importlib.util.spec_from_file_location("frozen_trait_mapping", path)
    module = importlib.util.module_from_spec(spec)
    old = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = old
    return module


def contextualize_hits(data, direction, sample, genes, traits, mapping_rules, mapper):
    result = []
    for hit in data:
        gene_id, trait_id = (hit["qseqid"], hit["sseqid"]) if direction == "proteins_to_seeds" else (hit["sseqid"], hit["qseqid"])
        gene = genes[gene_id]
        result.append({**hit, "sample_id": sample, "gene_id": gene_id, "contig_id": gene["contig_id"],
                       "start1": gene["start0"]+1, "end1": gene["end0"], "strand": gene["strand"],
                       "trait_id": trait_id, "trait_role": traits[trait_id]["role"],
                       "single_direction_pass": mapper.qualifies(hit, mapping_rules)})
    return result


def source_inputs(config_path):
    config_path = Path(config_path).resolve()
    config = read_json(config_path)
    keys = {"schema", "inventory", "dictionary", "mapping_source", "mapping_rules", "plan", "tool_receipts", "tools", "threads"}
    if set(config) != keys or config["schema"] != "sequence-integrity.assembly-config.v1":
        raise ValueError("unexpected configuration schema or keys")
    if type(config["threads"]) is not int or config["threads"] < 1 or config["threads"] > 128:
        raise ValueError("threads must be 1..128")
    if set(config["tools"]) != {"prodigal", "diamond", "tblastn", "makeblastdb"}:
        raise ValueError("four explicit executable paths required")
    paths = [config_path, Path(__file__).resolve(), HERE/"ASSEMBLY_RULES.json", Path(config["inventory"]),
             Path(config["mapping_source"]), Path(config["mapping_rules"]), Path(config["plan"])]
    if not isinstance(config["tool_receipts"], list) or not config["tool_receipts"]:
        raise ValueError("explicit tool provenance receipts required")
    paths += [Path(x) for x in config["tool_receipts"]]
    paths += [Path(config["dictionary"])/name for name in DICT_FILES]
    paths += [Path(config["tools"][key]) for key in sorted(config["tools"])]
    if any(not p.is_absolute() for p in paths):
        raise ValueError("all configured paths must be absolute")
    bindings = [record(p) for p in paths]
    rules = read_json(HERE/"ASSEMBLY_RULES.json")
    expected = [(config["mapping_source"], "mapping_source_sha256"), (config["mapping_rules"], "mapping_rules_sha256"),
                (Path(config["dictionary"])/"FROZEN.json", "dictionary_frozen_sha256")]
    for path, key in expected:
        if sha(path) != rules[key]:
            raise ValueError(f"wrong frozen source: {path}")
    if sha(config["tools"]["prodigal"]) != rules["prodigal_sha256"]:
        raise ValueError("Prodigal differs from bound v2.6.3 Linux binary")
    dictionary = Path(config["dictionary"])
    frozen = read_json(dictionary/"FROZEN.json")
    entries = {e["path"]: e for e in frozen["outputs"]}
    for name in DICT_FILES[1:]:
        entry = entries[name]
        if sha(dictionary/name) != entry["sha256"] or (dictionary/name).stat().st_size != entry["bytes"]:
            raise ValueError("frozen dictionary mismatch")
    seeds, _ = fasta(dictionary/"trait_proteins.faa")
    cds, _ = fasta(dictionary/"trait_cds.fna")
    trait_rows = rows(dictionary/"TRAITS.tsv")
    traits = {r["trait_id"]: r for r in trait_rows}
    if len(traits) != len(trait_rows) or set(seeds) != set(traits) or set(cds) != set(traits):
        raise ValueError("seed identifiers disagree")
    if Counter(t["role"] for t in traits.values()) != {"positive_enzyme_seed": 6, "specificity_comparator": 3}:
        raise ValueError("six positive/three comparator seeds required")
    for ident, trait in traits.items():
        if seqsha(seeds[ident]) != trait["protein_sha256"] or seqsha(cds[ident]) != trait["cds_sha256"]:
            raise ValueError("seed sequence mismatch")
    modules = rows(dictionary/"MODULE_RULES.tsv")
    if len(modules) != 4 or len({m["module_id"] for m in modules}) != 4:
        raise ValueError("four frozen modules required")
    for module in modules:
        if module["complete_utilization_call_enabled"] != "False" or any(traits[t]["role"] != "positive_enzyme_seed" for t in module["required_all"].split(";")):
            raise ValueError("invalid module role or claim")
    inventory = read_json(config["inventory"])
    report = inventory["historical_report"]
    if report["status"] != "COMPLETE" or report["samples_n"] != 50:
        raise ValueError("complete original 50-library inventory required")
    checked_rows = inventory["assembly_sha_reverification_and_tools"]["assemblies"]
    checked = {r["sample_id"]: r for r in checked_rows}
    if len(checked_rows) != 50 or len(checked) != len(checked_rows):
        raise ValueError("assembly checks must contain exactly 50 unique libraries")
    samples = []
    for original in report["samples"]:
        item = checked[original["sample"]]
        if not item["matches_historical_report"] or not item["stable_stat"] or item["path"] != original["contig_path"] or item["sha256"] != original["contig_sha256"]:
            raise ValueError("assembly source receipts disagree")
        samples.append({"sample_id": original["sample"], "dyad": original["dyad"], "window": original["window"],
                        "historical_contigs": original["contigs"], "historical_bases": original["bases"],
                        "assembly": {k: item[k] for k in ("path", "bytes", "sha256")}})
    if len(samples) != 50 or len({x["sample_id"] for x in samples}) != 50 or len(checked) != 50:
        raise ValueError("exactly 50 unique libraries required")
    return config, rules, bindings, samples, seeds, traits, modules


def fresh(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=False)
    return path


def complete_run(out, manifest):
    if (out/"COMPLETE").exists() or (out/"FAILED.json").exists():
        raise ValueError("cannot finalize an already completed or failed directory")
    excluded = {"RUN.json", "COMPLETE", "FAILED.json"}
    outputs = []
    for path in sorted(out.rglob("*")):
        if path.is_symlink():
            raise ValueError("output symlinks forbidden")
        if path.is_file() and path.relative_to(out).as_posix() not in excluded:
            rec = record(path)
            rec["path"] = path.relative_to(out).as_posix()
            outputs.append(rec)
    manifest = {**manifest, "status": "COMPLETE", "finished_utc": now(), "outputs": outputs}
    write_json(out/"RUN.json", manifest)
    (out/"COMPLETE").write_text(sha(out/"RUN.json") + "\n", encoding="ascii")
    return manifest


def verify_complete(out):
    out = Path(out)
    manifest = read_json(out/"RUN.json")
    if manifest["status"] != "COMPLETE" or (out/"COMPLETE").read_text().strip() != sha(out/"RUN.json") or (out/"FAILED.json").exists():
        raise ValueError("invalid completion receipt")
    seen = set()
    for rec in manifest["outputs"]:
        relative = Path(rec["path"])
        if relative.is_absolute() or ".." in relative.parts or rec["path"] in seen:
            raise ValueError("unsafe/duplicate output path")
        seen.add(rec["path"])
        path = out/relative
        if path.is_symlink() or path.stat().st_size != rec["bytes"] or sha(path) != rec["sha256"]:
            raise ValueError("completed output changed")
    actual = {p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file() and p.relative_to(out).as_posix() not in {"RUN.json", "COMPLETE"}}
    if actual != seen:
        raise ValueError("unbound completed output")
    return manifest


def prepare(config_path, out):
    out = fresh(out)
    try:
        config, rules, bindings, samples, _, _, _ = source_inputs(config_path)
        for sample in samples:
            verify_record(sample["assembly"])
        prepared = {"schema": "sequence-integrity.assembly-prepared.v1", "created_utc": now(),
                    "config_path": str(Path(config_path).resolve()), "bindings": bindings, "samples": samples,
                    "index_base": 0, "sample_order": "historical_report.samples order; never filtered by evidence", "rules": rules,
                    "analysis_run": False, "phenotype_inputs": []}
        write_json(out/"PREPARED.json", prepared)
        return complete_run(out, {"schema": "sequence-integrity.assembly-prepare-run.v1", "created_utc": now(), "bindings": bindings, "sample_count": 50, "analysis_run": False})
    except Exception as exc:
        write_json(out/"FAILED.json", {"at_utc": now(), "error": str(exc)})
        raise


class Commands:
    def __init__(self, out):
        self.out, self.records = out, []
        (out/"logs").mkdir()

    def run(self, name, argv):
        stdout = self.out/"logs"/(name+".stdout.txt")
        stderr = self.out/"logs"/(name+".stderr.txt")
        entry = {"name": name, "argv": [str(x) for x in argv], "started_utc": now(), "stdout": str(stdout.relative_to(self.out)), "stderr": str(stderr.relative_to(self.out))}
        self.records.append(entry)
        write_json(self.out/"COMMANDS.json", self.records)
        with stdout.open("wb") as out_handle, stderr.open("wb") as err_handle:
            result = subprocess.run(entry["argv"], stdout=out_handle, stderr=err_handle, cwd=self.out, check=False)
        entry.update(returncode=result.returncode, finished_utc=now())
        write_json(self.out/"COMMANDS.json", self.records)
        if result.returncode:
            raise RuntimeError(f"{name} failed ({result.returncode}); see retained logs")


def run_sample(prepared_dir, index, out):
    # Validate preparation before allocating the fresh sample directory.
    preparation = verify_complete(prepared_dir)
    prepared_path = Path(prepared_dir)/"PREPARED.json"
    prepared = read_json(prepared_path)
    if preparation["schema"] != "sequence-integrity.assembly-prepare-run.v1" or prepared["schema"] != "sequence-integrity.assembly-prepared.v1":
        raise ValueError("wrong prepared receipt schema")
    config, rules, bindings, samples, seeds, traits, modules = source_inputs(prepared["config_path"])
    if bindings != prepared["bindings"] or bindings != preparation["bindings"] or samples != prepared["samples"] or rules != prepared["rules"]:
        raise ValueError("prepared sources/configuration changed")
    if type(index) is not int or not 0 <= index < len(samples):
        raise ValueError("sample index must be 0..49")
    sample = samples[index]
    verify_record(sample["assembly"])
    out = fresh(out)
    started = now()
    manifest = {"schema": "sequence-integrity.assembly-sample-run.v1", "status": "RUNNING", "created_utc": started,
                "sample_index": index, "sample": sample, "prepared": record(prepared_path),
                "prepared_run": record(Path(prepared_dir)/"RUN.json"), "bindings": bindings, "rules": rules, "phenotype_inputs": []}
    write_json(out/"RUN.json", manifest)
    try:
        commands = Commands(out)
        tools = config["tools"]
        for name, args in (("prodigal", ["-v"]), ("diamond", ["version"]), ("tblastn", ["-version"]), ("makeblastdb", ["-version"])):
            commands.run(name+"_version", [tools[name], *args])
        assembly, _ = fasta(sample["assembly"]["path"])
        if len(assembly) != sample["historical_contigs"] or sum(map(len, assembly.values())) != sample["historical_bases"]:
            raise ValueError("assembly lengths/count disagree with source report")
        if any(not set(seq) <= set("ACGTRYMKSWBDHVN") for seq in assembly.values()):
            raise ValueError("assembly contains non-IUPAC nucleotide")
        mapper = load_mapper(config["mapping_source"])
        mapping_rules = read_json(config["mapping_rules"])
        orfs, assignments, forward, reverse, tb_hits, tb_summary = [], [], [], [], [], []
        genes = {}
        state = "empty_assembly" if not assembly else "assembly_scanned"
        if assembly:
            commands.run("prodigal", [tools["prodigal"], "-i", sample["assembly"]["path"], "-p", "meta", "-g", "11", "-f", "gff", "-o", out/"predicted.gff", "-a", out/"predicted.raw.faa", "-d", out/"predicted.cds.fna"])
            orfs, genes, proteins = parse_predictions(assembly, out/"predicted.raw.faa", out/"predicted.cds.fna", out/"predicted.gff", sample["sample_id"])
            write_fasta(out/"predicted.normalized.faa", proteins)
            if proteins:
                seed_path = Path(config["dictionary"])/"trait_proteins.faa"
                commands.run("seed_makedb", [tools["diamond"], "makedb", "--in", seed_path, "--db", out/"seeds", "--threads", str(config["threads"])])
                commands.run("protein_makedb", [tools["diamond"], "makedb", "--in", out/"predicted.normalized.faa", "--db", out/"proteins", "--threads", str(config["threads"])])
                for name, query, database in (("proteins_to_seeds", out/"predicted.normalized.faa", out/"seeds"), ("seeds_to_proteins", seed_path, out/"proteins")):
                    commands.run(name, [tools["diamond"], "blastp", "--query", query, "--db", database, "--out", out/(name+".raw.tsv"), "--threads", str(config["threads"]), "--very-sensitive", "--evalue", "0.001", "--max-target-seqs", "0", "--max-hsps", "0", "--outfmt", "6", *mapper.FIELDS])
                forward = mapper.read_hits(out/"proteins_to_seeds.raw.tsv", proteins, seeds)
                reverse = mapper.read_hits(out/"seeds_to_proteins.raw.tsv", seeds, proteins)
                assignments = mapper.classify(seeds, traits, genes, proteins, forward, reverse, mapping_rules)
            else:
                state = "assembly_scanned_no_predicted_orfs"
            commands.run("assembly_makedb", [tools["makeblastdb"], "-in", sample["assembly"]["path"], "-dbtype", "nucl", "-parse_seqids", "-out", out/"assembly_db"])
            commands.run("tblastn", [tools["tblastn"], "-query", Path(config["dictionary"])/"trait_proteins.faa", "-db", out/"assembly_db", "-out", out/"tblastn.raw.tsv", "-num_threads", str(config["threads"]), "-db_gencode", "11", "-evalue", str(rules["tblastn_search_evalue"]), "-seg", rules["tblastn_seg"], "-comp_based_stats", rules["tblastn_comp_based_stats"], "-word_size", str(rules["tblastn_word_size"]), "-matrix", rules["tblastn_matrix"], "-gapopen", str(rules["tblastn_gapopen"]), "-gapextend", str(rules["tblastn_gapextend"]), "-max_target_seqs", str(max(1, len(assembly))), "-outfmt", "6 " + " ".join(TB_FIELDS)])
            tb_hits, tb_summary = parse_tblastn(out/"tblastn.raw.tsv", seeds, assembly, traits, mapping_rules)
        # Stable empty outputs document not-run/no-prediction states without hiding denominators.
        for name in ("proteins_to_seeds.raw.tsv", "seeds_to_proteins.raw.tsv", "tblastn.raw.tsv"):
            if not (out/name).exists():
                (out/name).write_text("", encoding="ascii")
        write_tsv(out/"PREDICTED_ORFS.tsv", orfs, ORF_FIELDS)
        write_tsv(out/"PROTEIN_ASSIGNMENTS.tsv", assignments, ASSIGNMENT_FIELDS)
        for name, data in (("proteins_to_seeds", forward), ("seeds_to_proteins", reverse)):
            annotated = contextualize_hits(data, name, sample["sample_id"], genes, traits, mapping_rules, mapper)
            write_tsv(out/(name+".parsed.tsv"), annotated, mapper.FIELDS + ["qcov", "scov"] + PROTEIN_CONTEXT_FIELDS)
        write_tsv(out/"TBLASTN_HSPS.tsv", tb_hits, TB_FIELDS + ["hsp_row_number", "trait_role", "strand", "query_coverage_percent", "high_identity_hsp"])
        write_tsv(out/"TBLASTN_QUERY_COVERAGE.tsv", tb_summary, TB_SUMMARY_FIELDS)
        trait_rows, module_rows, contig_rows = evidence_tables(sample["sample_id"], state, traits, modules, orfs, assignments, forward, reverse, tb_hits)
        write_tsv(out/"TRAIT_STATUS.tsv", trait_rows, TRAIT_FIELDS)
        write_tsv(out/"MODULE_SAMPLE_STATUS.tsv", module_rows, MODULE_FIELDS)
        write_tsv(out/"CONTIG_MODULE_EVIDENCE.tsv", contig_rows, CONTIG_MODULE_FIELDS)
        verify_record(sample["assembly"])
        for binding in bindings:
            verify_record(binding)
        if record(prepared_path) != manifest["prepared"] or record(Path(prepared_dir)/"RUN.json") != manifest["prepared_run"]:
            raise ValueError("preparation changed during run")
        return complete_run(out, {**manifest, "assembly_state": state, "counts": {"contigs": len(assembly), "orfs": len(orfs), "positive_trait_rows": len(trait_rows), "module_rows": len(module_rows), "tblastn_hsps": len(tb_hits)}, "claims": {"functional_activity": False, "same_strain": False, "complete_utilization": False, "not_detected_means_absent": False}})
    except Exception as exc:
        write_json(out/"FAILED.json", {"at_utc": now(), "sample_id": sample["sample_id"], "error": str(exc)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--config", required=True, type=Path)
    prep.add_argument("--out", required=True, type=Path)
    run = sub.add_parser("run")
    run.add_argument("--prepared", required=True, type=Path)
    run.add_argument("--index", required=True, type=int)
    run.add_argument("--out", required=True, type=Path)
    verify = sub.add_parser("verify")
    verify.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.command == "prepare":
        result = prepare(args.config, args.out)
    elif args.command == "run":
        result = run_sample(args.prepared, args.index, args.out)
    else:
        result = verify_complete(args.out)
    print(json.dumps({"status": result["status"], "schema": result["schema"]}))


if __name__ == "__main__":
    main()
