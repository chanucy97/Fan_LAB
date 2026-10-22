"""Frozen, phenotype-blind reference annotation. Production alignment runs on HPC."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIELDS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen"]
DICT_FILES = ["FROZEN.json", "TRAITS.tsv", "trait_proteins.faa", "trait_cds.fna", "MODULE_RULES.tsv"]
PANEL_FILES = ["BUILD_SUMMARY.json", "proteins.faa", "cds_coordinates.tsv", "panel.fna", "contigs.tsv", "genomes.tsv"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def seqsha(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def indexed(rows, key):
    result = {}
    for row in rows:
        ident = row[key]
        if not ident or ident in result:
            raise ValueError(f"empty/duplicate {key}: {ident}")
        result[ident] = row
    return result


def fasta(path):
    result, key, chunks = {}, None, []
    def save():
        if key is not None:
            sequence = "".join(chunks).upper()
            if not sequence or key in result:
                raise ValueError(f"empty/duplicate FASTA sequence: {key}")
            result[key] = sequence
    with Path(path).open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                key, chunks = line[1:].split()[0], []
            elif key is None:
                raise ValueError("FASTA data before header")
            else:
                chunks.append(line)
    save()
    return result


def truth(value):
    if value not in ("True", "False"):
        raise ValueError(f"invalid boolean: {value}")
    return value == "True"


def load_inputs(dictionary, panel):
    rules = json.loads((HERE / "RULES.json").read_text())
    if sha(dictionary / "FROZEN.json") != rules["dictionary_frozen_sha256"]:
        raise ValueError("unexpected dictionary freeze")
    frozen = json.loads((dictionary / "FROZEN.json").read_text())
    for entry in frozen["outputs"]:
        if entry["path"] in DICT_FILES and sha(dictionary / entry["path"]) != entry["sha256"]:
            raise ValueError("dictionary output differs from original freeze")
    panel_summary = json.loads((panel / "BUILD_SUMMARY.json").read_text())
    for entry in panel_summary["outputs"]:
        name = Path(entry["path"].replace("\\", "/")).name
        if name in PANEL_FILES and sha(panel / name) != entry["sha256"]:
            raise ValueError("reference panel output differs from build receipt")
    seeds = fasta(dictionary / "trait_proteins.faa")
    cds_seeds = fasta(dictionary / "trait_cds.fna")
    traits = indexed(tsv(dictionary / "TRAITS.tsv"), "trait_id")
    genes = indexed(tsv(panel / "cds_coordinates.tsv"), "gene_id")
    proteins = fasta(panel / "proteins.faa")
    contigs = fasta(panel / "panel.fna")
    contig_rows = indexed(tsv(panel / "contigs.tsv"), "contig_id")
    if set(seeds) != set(traits) or set(cds_seeds) != set(traits):
        raise ValueError("dictionary identifiers disagree")
    if len(seeds) != 9 or Counter(x["role"] for x in traits.values()) != {"positive_enzyme_seed": 6, "specificity_comparator": 3}:
        raise ValueError("expected frozen six positive / three comparator dictionary")
    if set(contigs) != set(contig_rows):
        raise ValueError("contig identifiers disagree")
    for ident, sequence in seeds.items():
        if seqsha(sequence) != traits[ident]["protein_sha256"] or seqsha(cds_seeds[ident]) != traits[ident]["cds_sha256"]:
            raise ValueError(f"dictionary sequence mismatch: {ident}")
    for ident, sequence in contigs.items():
        if seqsha(sequence) != contig_rows[ident]["sequence_sha256"] or len(sequence) != int(contig_rows[ident]["length"]):
            raise ValueError(f"contig sequence mismatch: {ident}")
    exported = {key for key, row in genes.items() if truth(row["source_protein_exported"])}
    if exported != set(proteins):
        raise ValueError("protein FASTA and coordinate export flags disagree")
    for ident, sequence in proteins.items():
        row = genes[ident]
        if seqsha(sequence) != row["protein_sequence_sha256"] or len(sequence) != int(row["protein_length"]):
            raise ValueError(f"protein sequence mismatch: {ident}")
    for ident, row in genes.items():
        if row["contig_id"] not in contigs or row["assembly_accession"] != contig_rows[row["contig_id"]]["assembly_accession"]:
            raise ValueError(f"gene/contig/genome mismatch: {ident}")
        if not (0 <= int(row["start0"]) < int(row["end0"]) <= len(contigs[row["contig_id"]])):
            raise ValueError(f"invalid coordinates: {ident}")
        if int(row["start1"]) != int(row["start0"]) + 1 or int(row["end1"]) != int(row["end0"]):
            raise ValueError(f"coordinate conventions disagree: {ident}")
    return seeds, cds_seeds, traits, genes, proteins, contigs, contig_rows


def bindings(dictionary, panel):
    rows = []
    for role, base, names in (("dictionary", dictionary, DICT_FILES), ("panel", panel, PANEL_FILES), ("mapping", HERE, ["map_traits.py", "RULES.json"])):
        for name in names:
            path = base / name
            rows.append({"root": role, "path": name, "bytes": path.stat().st_size, "sha256": sha(path)})
    return rows


def verify_bindings(manifest, dictionary, panel, out):
    if manifest["inputs"] != bindings(dictionary, panel):
        raise ValueError("frozen input, rule, or source-code bytes changed")
    if manifest["rules"] != json.loads((HERE / "RULES.json").read_text()):
        raise ValueError("manifest rule copy differs from bound rule file")
    expected = ["exact_seed_matches.tsv", "nanh2_reference_regions.bed", "nanh2_reference_regions.tsv"]
    if sorted(entry["path"] for entry in manifest["prepared_artifacts"]) != expected:
        raise ValueError("prepared coordinate and exact-match artifacts must all be bound")
    for entry in manifest["prepared_artifacts"]:
        path = out / entry["path"]
        if sha(path) != entry["sha256"] or path.stat().st_size != entry["bytes"]:
            raise ValueError("prepared reference artifact changed")


def exact_matches(seeds, proteins):
    byseq = defaultdict(list)
    for trait_id, sequence in seeds.items():
        byseq[sequence].append(trait_id)
    return {gene_id: sorted(byseq[sequence]) for gene_id, sequence in proteins.items() if sequence in byseq}


def rescue_nanh2(cds_seeds, contigs, genes, contig_rows):
    trait_id = "NanH2_Binf_ATCC15697"
    sequence = cds_seeds[trait_id]
    complement = str.maketrans("ACGT", "TGCA")
    rows = []
    for contig_id, genome_sequence in contigs.items():
        for strand, query in ((1, sequence), (-1, sequence.translate(complement)[::-1])):
            start = genome_sequence.find(query)
            while start >= 0:
                end = start + len(query)
                overlaps = [r for r in genes.values() if r["contig_id"] == contig_id and int(r["start0"]) < end and int(r["end0"]) > start]
                pseudo = [r["gene_id"] for r in overlaps if truth(r["annotated_pseudogene"])]
                rows.append({"region_id": f"{contig_id}__historical_NanH2_{start+1}_{end}", "trait_id": trait_id, "contig_id": contig_id,
                             "assembly_accession": contig_rows[contig_id]["assembly_accession"], "start0": start, "end0": end, "strand": strand,
                             "cds_sha256": seqsha(sequence), "overlapping_gene_ids": ";".join(r["gene_id"] for r in overlaps),
                             "overlapping_pseudogene_ids": ";".join(pseudo), "status": "exact_historical_cds_reference_region",
                             "annotation_conflict": bool(pseudo), "ordinary_instrain_gene_id": "", "cohort_intact_gene_asserted": False})
                start = genome_sequence.find(query, start + 1)
    return rows


def prepare(dictionary, panel, out):
    if (out / "INPUTS.json").exists():
        raise ValueError("prepared directory already exists; use a new output directory")
    out.mkdir(parents=True, exist_ok=True)
    seeds, cds_seeds, traits, genes, proteins, contigs, contig_rows = load_inputs(dictionary, panel)
    rules = json.loads((HERE / "RULES.json").read_text())
    if rules["complete_utilization_call_enabled"] or rules["phenotype_input_allowed"] or rules["functional_equivalence_asserted"]:
        raise ValueError("forbidden biological claim configuration")
    exact = exact_matches(seeds, proteins)
    exact_rows = [{"gene_id": gene, "trait_ids": ";".join(ids), "assembly_accession": genes[gene]["assembly_accession"], "protein_sha256": seqsha(proteins[gene])} for gene, ids in sorted(exact.items())]
    write_tsv(out / "exact_seed_matches.tsv", exact_rows, ["gene_id", "trait_ids", "assembly_accession", "protein_sha256"])
    rescue = rescue_nanh2(cds_seeds, contigs, genes, contig_rows)
    fields = ["region_id", "trait_id", "contig_id", "assembly_accession", "start0", "end0", "strand", "cds_sha256", "overlapping_gene_ids", "overlapping_pseudogene_ids", "status", "annotation_conflict", "ordinary_instrain_gene_id", "cohort_intact_gene_asserted"]
    write_tsv(out / "nanh2_reference_regions.tsv", rescue, fields)
    with (out / "nanh2_reference_regions.bed").open("w", encoding="ascii") as handle:
        for row in rescue:
            handle.write(f"{row['contig_id']}\t{row['start0']}\t{row['end0']}\t{row['region_id']}\t0\t{'+' if row['strand'] == 1 else '-'}\n")
    artifacts = [{"path": name, "bytes": (out / name).stat().st_size, "sha256": sha(out / name)} for name in ["exact_seed_matches.tsv", "nanh2_reference_regions.bed", "nanh2_reference_regions.tsv"]]
    manifest = {"schema": rules["schema"], "created_utc": datetime.now(timezone.utc).isoformat(), "inputs": bindings(dictionary, panel), "prepared_artifacts": artifacts,
                "counts": {"seeds": len(seeds), "proteins": len(proteins), "exact_peptide_matches": len(exact), "nanh2_exact_reference_regions": len(rescue)},
                "rules": rules, "cohort_or_phenotype_inputs": [], "alignment_has_run": False}
    write_json(out / "INPUTS.json", manifest)
    return manifest


def read_hits(path, queries, subjects):
    rows = []
    with Path(path).open(encoding="ascii") as handle:
        for line_no, line in enumerate(handle, 1):
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(FIELDS):
                raise ValueError(f"malformed alignment row {line_no}")
            row = dict(zip(FIELDS, parts))
            for field in ("length", "qstart", "qend", "sstart", "send", "qlen", "slen"):
                row[field] = int(row[field])
            for field in ("pident", "evalue", "bitscore"):
                row[field] = float(row[field])
                if not math.isfinite(row[field]):
                    raise ValueError("non-finite alignment metric")
            if row["qseqid"] not in queries or row["sseqid"] not in subjects:
                raise ValueError("alignment contains unbound sequence identifier")
            if row["qlen"] != len(queries[row["qseqid"]]) or row["slen"] != len(subjects[row["sseqid"]]):
                raise ValueError("alignment sequence length mismatch")
            if not (0 <= row["pident"] <= 100 and row["evalue"] >= 0 and row["bitscore"] > 0 and row["length"] > 0):
                raise ValueError("invalid alignment metrics")
            for prefix in ("q", "s"):
                if not 1 <= row[prefix+"start"] <= row[prefix+"end"] <= row[prefix+"len"]:
                    raise ValueError("invalid protein alignment coordinates")
                span = row[prefix+"end"] - row[prefix+"start"] + 1
                if span > row["length"]:
                    raise ValueError("alignment span exceeds alignment length")
                row[prefix+"cov"] = 100 * span / row[prefix+"len"]
            rows.append(row)
    return rows


def best_pairs(rows):
    pairs = {}
    for row in rows:
        key = (row["qseqid"], row["sseqid"])
        rank = lambda x: (x["bitscore"], x["qcov"], x["scov"], x["pident"], -x["qstart"], -x["sstart"])
        if key not in pairs or rank(row) > rank(pairs[key]):
            pairs[key] = row
    return pairs


def qualifies(row, rules):
    return row is not None and row["pident"] >= rules["minimum_identity_percent"] and row["qcov"] >= rules["minimum_query_coverage_percent"] and row["scov"] >= rules["minimum_subject_coverage_percent"] and row["evalue"] <= rules["maximum_evalue"]


def classify(seeds, traits, genes, proteins, forward_rows, reverse_rows, rules):
    forward, reverse = best_pairs(forward_rows), best_pairs(reverse_rows)
    exact = exact_matches(seeds, proteins)
    by_gene, reverse_by_genome = defaultdict(list), defaultdict(list)
    for (gene_id, trait_id), row in forward.items():
        by_gene[gene_id].append(row)
    for (trait_id, gene_id), row in reverse.items():
        reverse_by_genome[(trait_id, genes[gene_id]["assembly_accession"])].append(row)
    output = []
    for gene_id in sorted(set(by_gene) | set(exact)):
        hits = sorted(by_gene[gene_id], key=lambda r: (-r["bitscore"], r["sseqid"]))
        exact_ids = exact.get(gene_id, [])
        top = hits[0] if hits else None
        trait_id = exact_ids[0] if len(exact_ids) == 1 else top["sseqid"] if top else ""
        row = forward.get((gene_id, trait_id))
        rev = reverse.get((trait_id, gene_id))
        seed_near = sorted(r["sseqid"] for r in hits if r["bitscore"] >= hits[0]["bitscore"] * rules["near_best_bitscore_fraction"])
        genome_hits = reverse_by_genome.get((trait_id, genes[gene_id]["assembly_accession"]), [])
        best_genome_score = max((r["bitscore"] for r in genome_hits), default=0)
        protein_near = sorted(r["sseqid"] for r in genome_hits if r["bitscore"] >= best_genome_score * rules["near_best_bitscore_fraction"])
        ambiguity = len(seed_near) != 1 or seed_near != [trait_id] or protein_near != [gene_id]
        if len(exact_ids) > 1:
            status = "ambiguous_identical_seed_peptides"
        elif len(exact_ids) == 1:
            status = "exact_seed_peptide"
        elif qualifies(row, rules) and qualifies(rev, rules):
            status = "ambiguous_high_similarity" if ambiguity else "high_similarity_candidate"
        else:
            status = "below_operational_threshold"
        annotation = genes[gene_id]
        output.append({"gene_id": gene_id, "assembly_accession": annotation["assembly_accession"], "contig_id": annotation["contig_id"],
                       "start0": annotation["start0"], "end0": annotation["end0"], "strand": annotation["strand"],
                       "locus_tag": annotation["locus_tag"], "old_locus_tag": annotation["old_locus_tag"], "protein_id": annotation["protein_id"],
                       "trait_id": trait_id, "trait_role": traits[trait_id]["role"] if trait_id else "", "mapping_status": status,
                       "exact_seed_ids": ";".join(exact_ids), "seed_near_best_ids": ";".join(seed_near), "within_genome_near_best_gene_ids": ";".join(protein_near),
                       "competitor_ambiguity": ambiguity, "forward_identity": row["pident"] if row else "", "forward_query_coverage": row["qcov"] if row else "",
                       "forward_subject_coverage": row["scov"] if row else "", "forward_bitscore": row["bitscore"] if row else "",
                       "reverse_identity": rev["pident"] if rev else "", "reverse_query_coverage": rev["qcov"] if rev else "",
                       "reverse_subject_coverage": rev["scov"] if rev else "", "reverse_bitscore": rev["bitscore"] if rev else "",
                       "reverse_pass": qualifies(rev, rules), "reference_instrain_gene_exported": annotation["instrain_gene_exported"],
                       "sequence_sha256": seqsha(proteins[gene_id]), "functional_equivalence_asserted": False,
                       "claim_ceiling": traits[trait_id]["claim_ceiling"] if trait_id else "unresolved"})
    return output


def summarize(dictionary, panel, out):
    if (out / "COMPLETE").exists():
        raise ValueError("completed run exists; refusing to rewrite result")
    if (out / "FAILED.json").exists():
        raise ValueError("failed run cannot be summarized; rerun alignment first")
    manifest = json.loads((out / "INPUTS.json").read_text())
    verify_bindings(manifest, dictionary, panel, out)
    receipt = json.loads((out / "ALIGNMENT_RECEIPT.json").read_text())
    if receipt["status"] != "ALIGNMENTS_COMPLETE" or receipt["input_manifest_sha256"] != sha(out / "INPUTS.json"):
        raise ValueError("missing/mismatched successful alignment receipt")
    if receipt["commands_sha256"] != sha(out / "COMMANDS.json"):
        raise ValueError("alignment commands changed after execution")
    if sorted(entry["path"] for entry in receipt["files"]) != ["proteins_to_seeds.tsv", "seeds_to_proteins.tsv"]:
        raise ValueError("alignment receipt must bind both expected outputs exactly once")
    for entry in receipt["files"]:
        if sha(out / entry["path"]) != entry["sha256"]:
            raise ValueError("alignment file changed after execution")
    seeds, _, traits, genes, proteins, _, _ = load_inputs(dictionary, panel)
    forward = read_hits(out / "proteins_to_seeds.tsv", proteins, seeds)
    reverse = read_hits(out / "seeds_to_proteins.tsv", seeds, proteins)
    forward_pairs, reverse_pairs = best_pairs(forward), best_pairs(reverse)
    for gene_id, trait_ids in exact_matches(seeds, proteins).items():
        for trait_id in trait_ids:
            if (gene_id, trait_id) not in forward_pairs or (trait_id, gene_id) not in reverse_pairs:
                raise ValueError("known exact peptide control missing from a search direction")
    rows = classify(seeds, traits, genes, proteins, forward, reverse, manifest["rules"])
    if not rows:
        raise ValueError("empty mapping; expected known exact controls")
    write_tsv(out / "gene_trait_mapping.tsv", rows, list(rows[0]))
    summary = {"status": "COMPLETE", "input_manifest_sha256": sha(out / "INPUTS.json"), "forward_hsps": len(forward), "reverse_hsps": len(reverse),
               "status_counts": dict(Counter(r["mapping_status"] for r in rows)), "trait_mapping_rows": len(rows),
               "reference_proteins_without_reported_seed_hit": len(proteins) - len(rows), "no_hit_is_not_biological_absence": True,
               "functional_equivalence_asserted": False, "cohort_or_phenotype_inputs": [],
               "files": [{"path": name, "sha256": sha(out / name), "bytes": (out / name).stat().st_size} for name in ["proteins_to_seeds.tsv", "seeds_to_proteins.tsv", "gene_trait_mapping.tsv", "nanh2_reference_regions.tsv", "nanh2_reference_regions.bed", "exact_seed_matches.tsv"]]}
    write_json(out / "MAPPING_SUMMARY.json", summary)
    return summary


def run(dictionary, panel, out, diamond, threads):
    manifest = json.loads((out / "INPUTS.json").read_text())
    verify_bindings(manifest, dictionary, panel, out)
    if (out / "COMPLETE").exists():
        raise ValueError("completed run exists; use a new output directory")
    if threads < 1:
        raise ValueError("threads must be positive")
    (out / "FAILED.json").unlink(missing_ok=True)
    log = out / "diamond.log"
    commands = []
    try:
        version = subprocess.check_output([diamond, "version"], text=True).strip()
        for label, source in (("seeds", dictionary / "trait_proteins.faa"), ("proteins", panel / "proteins.faa")):
            commands.append([diamond, "makedb", "--in", str(source), "--db", str(out / label), "--threads", str(threads)])
        config = manifest["rules"]["alignment"]
        for label, query, db in (("proteins_to_seeds", panel / "proteins.faa", "seeds"), ("seeds_to_proteins", dictionary / "trait_proteins.faa", "proteins")):
            commands.append([diamond, "blastp", "--query", str(query), "--db", str(out / db), "--out", str(out / (label + ".tsv")),
                             "--threads", str(threads), "--" + config["sensitivity"], "--evalue", str(config["search_evalue"]),
                             "--max-target-seqs", str(config["max_target_seqs"]), "--max-hsps", str(config["max_hsps"]), "--outfmt", "6", *FIELDS])
        write_json(out / "COMMANDS.json", {"version": version, "commands": commands, "input_manifest_sha256": sha(out / "INPUTS.json")})
        with log.open("w", encoding="utf-8") as handle:
            for command in commands:
                subprocess.run(command, check=True, stdout=handle, stderr=subprocess.STDOUT)
        write_json(out / "ALIGNMENT_RECEIPT.json", {"status": "ALIGNMENTS_COMPLETE", "input_manifest_sha256": sha(out / "INPUTS.json"),
                   "commands_sha256": sha(out / "COMMANDS.json"), "files": [{"path": name, "sha256": sha(out / name)} for name in ["proteins_to_seeds.tsv", "seeds_to_proteins.tsv"]]})
        result = summarize(dictionary, panel, out)
        result["diamond_version"] = version
        result["commands_sha256"] = sha(out / "COMMANDS.json")
        write_json(out / "MAPPING_SUMMARY.json", result)
        (out / "COMPLETE").write_text(sha(out / "MAPPING_SUMMARY.json") + "\n", encoding="ascii")
        return result
    except Exception as exc:
        (out / "COMPLETE").unlink(missing_ok=True)
        write_json(out / "FAILED.json", {"status": "FAILED", "error": str(exc), "input_manifest_sha256": sha(out / "INPUTS.json")})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "run", "summarize", "verify"])
    parser.add_argument("--dictionary", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--diamond", default="diamond")
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    dictionary, panel, out = args.dictionary.resolve(), args.panel.resolve(), args.out.resolve()
    if args.command == "prepare":
        result = prepare(dictionary, panel, out)
    elif args.command == "run":
        result = run(dictionary, panel, out, args.diamond, args.threads)
    elif args.command == "summarize":
        result = summarize(dictionary, panel, out)
    else:
        verify_bindings(json.loads((out / "INPUTS.json").read_text()), dictionary, panel, out)
        load_inputs(dictionary, panel)
        result = {"status": "PASS", "input_manifest_sha256": sha(out / "INPUTS.json")}
    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()
