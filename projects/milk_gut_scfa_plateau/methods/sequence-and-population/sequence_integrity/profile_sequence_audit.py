"""Audit fixed enzyme-reference loci using existing inStrain calls, not new calling.

The output sequence is a masked, unphased substitution consensus. It does not
assess indels, reconstruct an isolate, or certify an intact sample ORF.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import platform
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
LABEL = "REFERENCE_BACKED_SUBSTITUTION_CONSENSUS_INDELS_UNASSESSED_UNPHASED"
BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
# NCBI table 11 shares ordinary codon translations with the standard code.
CODONS = {a+b+c: aa for (a,b,c),aa in zip(
    ((a,b,c) for a in "TCAG" for b in "TCAG" for c in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}
STARTS_11 = {"TTG", "CTG", "ATT", "ATC", "ATA", "ATG", "GTG"}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def record(path):
    path = Path(path)
    return {"path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha(path)}


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def table(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def keyed(rows, field):
    result = {r[field]: r for r in rows}
    if len(result) != len(rows):
        raise ValueError("duplicate key: " + field)
    return result


def boolean(value):
    if value in (True, "True", "true", "1"):
        return True
    if value in (False, "False", "false", "0"):
        return False
    raise ValueError("invalid boolean: " + str(value))


def eq(actual, expected):
    if expected is None:
        if actual not in (None, ""):
            raise ValueError("expected missing value")
    elif not math.isclose(float(actual), float(expected), rel_tol=1e-7, abs_tol=1e-7):
        raise ValueError("numeric source mismatch: " + str((actual, expected)))


def fasta(path):
    result, key = {}, None
    with Path(path).open(encoding="ascii") as stream:
        for line in stream:
            line = line.strip()
            if line.startswith(">"):
                key = line[1:].split()[0]
                if key in result:
                    raise ValueError("duplicate FASTA key")
                result[key] = []
            elif line:
                if key is None:
                    raise ValueError("FASTA lacks header")
                result[key].append(line.upper())
    return {key: "".join(parts) for key, parts in result.items()}


def revcomp(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def translate(sequence):
    if len(sequence) % 3:
        raise ValueError("reference CDS is not codon aligned")
    aa = [CODONS.get(sequence[i:i+3], "X") for i in range(0, len(sequence), 3)]
    if sequence[:3] in STARTS_11:
        aa[0] = "M"
    return "".join(aa)


def write_tsv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=fields, extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


class Bindings:
    def __init__(self):
        self.files = {}

    def add(self, path, expected=None):
        actual = record(path)
        if expected is not None and any(actual[k] != expected[k] for k in ("bytes", "sha256")):
            raise ValueError("source hash/size differs: " + str(path))
        previous = self.files.setdefault(actual["path"], actual)
        if previous != actual:
            raise ValueError("source changed during audit")
        return Path(path)

    def stable(self):
        for path, expected in self.files.items():
            if record(path) != expected:
                raise ValueError("source changed during audit: " + path)


def load_roster(base, bindings, rules):
    code = base / "gene_evidence"
    plan_dir = code / "reference_plan_v2"
    plan_path = bindings.add(plan_dir / "REFERENCE_PLAN.json")
    if sha(plan_path) != rules["gene_plan_sha256"]:
        raise ValueError("unexpected frozen gene plan")
    plan = read_json(plan_path)
    roots = {"panel": base / "reference_panel", "dictionary": base / "trait_dictionary",
             "mapping": base / "trait_mapping/prepared_final", "code": code}
    for entry in plan["inputs"]:
        bindings.add(roots[entry["root"]] / entry["path"], entry)
    for entry in plan["artifacts"]:
        bindings.add(plan_dir / entry["path"], entry)
    for entry in read_json(code / "INSTALLED_API.json")["files"]:
        bindings.add(Path(entry["remote_source"]), entry)
    coordinates = keyed(table(roots["panel"] / "cds_coordinates.tsv"), "gene_id")
    panel = fasta(roots["panel"] / "panel.fna")
    seeds = fasta(roots["dictionary"] / "trait_proteins.faa")
    roster = []
    for row in table(roots["mapping"] / "gene_trait_mapping.tsv"):
        if row["mapping_status"] not in rules["accepted_mapping_states"] or boolean(row["competitor_ambiguity"]):
            continue
        c = coordinates[row["gene_id"]]
        if any(str(c[k]) != str(row[k]) for k in ("contig_id", "assembly_accession", "start0", "end0", "strand")):
            raise ValueError("mapped coordinate disagreement")
        if not boolean(c["instrain_gene_exported"]) or boolean(c["compound_location"]) or int(c["codon_start"]) != 1 or int(c["transl_table"]) != 11:
            raise ValueError("unsupported accepted reference CDS")
        roster.append({"locus_id": row["gene_id"], "trait_id": row["trait_id"], "role": row["trait_role"],
                       "assembly_accession": row["assembly_accession"], "contig_id": row["contig_id"],
                       "start0": int(row["start0"]), "end0": int(row["end0"]), "strand": int(row["strand"]),
                       "annotation_conflict": False, "evidence_type": "ordinary_reference_CDS",
                       "reference_cds_sha256": c["cds_sequence_sha256"]})
    c = read_json(plan_dir / "CURATED_HISTORICAL_CDS.json")
    if c["status"] != "CURATED_REFERENCE_DNA_VERIFIED" or c["trait_id"] != "NanH2_Binf_ATCC15697":
        raise ValueError("curated NanH2 identity changed")
    roster.append({"locus_id": c["region_id"], "trait_id": c["trait_id"], "role": "positive_enzyme_seed",
                   "assembly_accession": c["assembly_accession"], "contig_id": c["contig_id"],
                   "start0": int(c["start0"]), "end0": int(c["end0"]), "strand": int(c["strand"]),
                   "annotation_conflict": True, "evidence_type": "curated_historical_CDS",
                   "reference_cds_sha256": c["cds_sha256"]})
    roster.sort(key=lambda r: r["locus_id"])
    counts = Counter((r["role"], r["evidence_type"]) for r in roster)
    if len(roster) != 16 or counts[("positive_enzyme_seed", "ordinary_reference_CDS")] != 11 or counts[("positive_enzyme_seed", "curated_historical_CDS")] != 1 or counts[("specificity_comparator", "ordinary_reference_CDS")] != 4:
        raise ValueError("fixed 16-locus denominator changed")
    if len({r["trait_id"] for r in roster if r["role"] == "positive_enzyme_seed"}) != 6:
        raise ValueError("six-positive-trait denominator changed")
    for r in roster:
        ref = panel[r["contig_id"]][r["start0"]:r["end0"]]
        coding = ref if r["strand"] == 1 else revcomp(ref)
        if r["strand"] not in (-1, 1) or len(ref) != r["end0"] - r["start0"] or hashlib.sha256(coding.encode()).hexdigest() != r["reference_cds_sha256"]:
            raise ValueError("reference coding sequence/hash mismatch")
        if len(coding) % 3 or r["trait_id"] not in seeds:
            raise ValueError("invalid reference/seed")
    return roster, panel, seeds, plan


def load_profile(base, sample, bindings, plan, rules):
    directory = base / "profiles_v1" / sample
    run_path = bindings.add(directory / "RUN.json")
    run = read_json(run_path)
    marker = bindings.add(directory / "COMPLETE")
    if run["status"] != "COMPLETE" or run["sample"]["sample_id"] != sample or marker.read_text().strip() != run["completed_utc"]:
        raise ValueError("profile incomplete or wrong sample")
    for k, v in plan["rules"]["expected_profile_config"].items():
        if run["config"][k] != v:
            raise ValueError("profile filter drift: " + k)
    if {Path(r["path"]).name for r in run["references"]} != {"panel.fna", "panel.stb", "genes.instrain.fna"}:
        raise ValueError("unexpected profile reference roster")
    for entry in run["references"]:
        bindings.add(base / "reference_panel" / Path(entry["path"]).name, entry)
    gene_dir = base / "gene_evidence/results_v2" / sample
    summary_path = bindings.add(gene_dir / "SUMMARY.json")
    summary = read_json(summary_path)
    if summary["status"] != "COMPLETE" or summary["sample_id"] != sample or summary["reference_plan_sha256"] != rules["gene_plan_sha256"]:
        raise ValueError("gene evidence identity/freeze mismatch")
    if bindings.add(gene_dir / "COMPLETE").read_text().strip() != sha(summary_path):
        raise ValueError("gene evidence completion mismatch")
    for entry in summary["outputs"]:
        bindings.add(gene_dir / entry["path"], entry)
    for entry in summary["profile_inputs"]:
        bindings.add(Path(entry["path"]), entry)
    genes = keyed(table(gene_dir / "gene_evidence.tsv"), "gene_id")
    regions = keyed(table(gene_dir / "region_evidence.tsv"), "region_id")
    carriers = keyed(table(gene_dir / "carrier_evidence.tsv"), "assembly_accession")
    for row in list(genes.values()) + list(regions.values()) + list(carriers.values()):
        if row["sample_id"] != sample:
            raise ValueError("cross-sample evidence row")
    for row in genes.values():
        if row["profile_run_sha256"] != sha(run_path) or row["reference_plan_sha256"] != rules["gene_plan_sha256"]:
            raise ValueError("gene output source mismatch")
    if run["profile_status"] == "NO_ANALYSABLE_PAIRS_NOT_BIOLOGICAL_ABSENCE":
        return run, genes, regions, carriers, None, None, None, {}, None
    if run["profile_status"] != "PROFILE_GENERATED_COVERAGE_ELIGIBILITY_NOT_YET_SCORED":
        raise ValueError("unexpected profile status")
    success = [x for x in run["commands"] if x["stage"] == "instrain_profile" and x["returncode"] == 0]
    if not success:
        raise ValueError("no successful profile command")
    argv = success[-1]["argv"]
    for flag, value in {"--min_read_ani": .95, "--min_mapq": 9, "--min_cov": 5, "--min_freq": .05, "--fdr": 1e-6, "--min_insert": 50, "--max_insert_relative": 3}.items():
        if argv.count(flag) != 1 or float(argv[argv.index(flag)+1]) != value:
            raise ValueError("profile command differs: " + flag)
    if argv.count("--skip_mm_profiling") != 1 or argv[argv.index("--pairing_filter")+1] != "paired_only":
        raise ValueError("profile pair/mm mode differs")
    if Path(argv[argv.index("-o")+1]).resolve() != (directory / "instrain").resolve():
        raise ValueError("wrong profile command output")
    alias, original = directory / (sample + ".bam"), directory / "profile.bam"
    entry = [x for x in run.get("profile_bam_aliases", []) if Path(x["path"]).resolve() == alias.resolve()]
    if Path(argv[2]).resolve() != alias.resolve() or len(entry) != 1 or alias.is_symlink() or original.is_symlink() or not os.path.samefile(alias, original):
        raise ValueError("wrong sample BAM identity")
    st, entry = alias.stat(), entry[0]
    if Path(entry["source"]).resolve() != original.resolve() or entry["same_file"] is not True or (st.st_size, st.st_ino, st.st_dev) != (entry["bytes"], entry["inode"], entry["device"]):
        raise ValueError("sample BAM binding changed")
    api = read_json(base / "gene_evidence/INSTALLED_API.json")
    source = next(x for x in api["files"] if x["path"].endswith("profile/profile_utilities.py"))
    if run["instrain_base_quality_source"]["sha256"] != source["sha256"]:
        raise ValueError("base quality implementation changed")
    outputs = {str(Path(x["path"]).resolve()): x for x in run["outputs"] if not x.get("removed_after_success", False)}
    def upstream(relative):
        path = directory / relative
        expected = outputs.get(str(path.resolve()))
        if expected is None:
            raise ValueError("source absent from original RUN: " + relative)
        return bindings.add(path, expected)
    attrs = keyed(table(upstream("instrain/raw_data/attributes.tsv")), "name")
    if attrs["version"]["value"] != rules["profile_version"] or attrs["covT"]["type"] != "special" or attrs["Rdic"]["type"] != "pickle":
        raise ValueError("unsupported profile schema")
    cov = upstream("instrain/raw_data/" + Path(attrs["covT"]["value"]).name)
    rdic = upstream("instrain/raw_data/" + Path(attrs["Rdic"]["value"]).name)
    if cov.name != "covT.hd5" or rdic.name != "Rdic.pickle":
        raise ValueError("unexpected raw storage type")
    raw = upstream("instrain/raw_data/raw_snp_table.csv.gz")
    cumulative = upstream("instrain/raw_data/cumulative_snv_table.csv.gz")
    info = keyed(table(upstream("instrain/output/instrain_gene_info.tsv")), "gene")
    return run, genes, regions, carriers, cov, raw, cumulative, info, record(rdic)


def select_variants(path, roster, require_frequencies=False):
    if path is None:
        return {}
    intervals = defaultdict(list)
    for row in roster:
        intervals[row["contig_id"]].append((row["start0"], row["end0"]))
    selected = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"scaffold", "position", "ref_base", "A", "C", "T", "G", "con_base", "var_base", "mm", "allele_count", "class", "cryptic", "position_coverage"}
        if require_frequencies:
            required |= {"con_freq", "var_freq", "ref_freq"}
        if not required.issubset(set(reader.fieldnames or [])):
            # inStrain can write an empty dataframe without usable column names.
            remainder = list(reader)
            if not remainder and len(reader.fieldnames or []) <= 1:
                return {}
            raise ValueError("SNV source columns missing")
        for row in reader:
            if int(row["mm"]) != 0:
                raise ValueError("unexpected nonzero mismatch bucket")
            contig, pos = row["scaffold"], int(row["position"])
            if any(start <= pos < end for start, end in intervals.get(contig, [])):
                key = (contig, pos)
                if key in selected:
                    raise ValueError("duplicate candidate variant coordinate")
                selected[key] = row
    return selected


def crosscheck_variant_tables(raw, cumulative):
    if set(raw) != set(cumulative):
        raise ValueError("raw/cumulative candidate SNV roster differs")
    fields = ["ref_base", "A", "C", "T", "G", "con_base", "var_base", "mm", "allele_count", "class", "cryptic", "position_coverage"]
    for key in raw:
        if any(raw[key][f] != cumulative[key][f] for f in fields):
            raise ValueError("raw/cumulative candidate SNV differs")


def call_existing(row, reference, depth):
    """Interpret an original caller row; never call a new alternate allele."""
    counts = {base: int(row[base]) for base in BASES}
    if any(x < 0 for x in counts.values()) or sum(counts.values()) != int(row["position_coverage"]) or sum(counts.values()) != depth:
        raise ValueError("SNV counts do not match filtered depth")
    if row["ref_base"] != reference:
        raise ValueError("SNV reference differs from panel")
    nalleles = int(row["allele_count"])
    if not 0 <= nalleles <= 4:
        raise ValueError("invalid original allele count")
    if row["class"] not in ("SNS", "SNV", "con_SNV", "pop_SNV", "DivergentSite", "AmbiguousReference"):
        raise ValueError("unrecognized original caller class")
    for prefix, base in (("con", row["con_base"]), ("var", row["var_base"]), ("ref", reference)):
        if prefix + "_freq" in row and base in BASES and depth:
            eq(row[prefix + "_freq"], counts[base] / depth)
    mixed = nalleles > 1
    reason = ""
    if boolean(row["cryptic"]):
        reason = "original_call_cryptic"
    elif nalleles < 1 or row["class"] in ("DivergentSite", "AmbiguousReference"):
        reason = "original_caller_unresolved"
    elif row["con_base"] not in BASES:
        reason = "invalid_original_consensus"
    elif sum(x == max(counts.values()) for x in counts.values()) != 1:
        reason = "tied_major_alleles"
    elif counts[row["con_base"]] != max(counts.values()):
        raise ValueError("original consensus is not the most abundant allele")
    return ("N" if reason else row["con_base"]), reason, mixed, counts


def make_consensus(reference, depth, variants, start0=0, strand=1):
    if strand not in (-1, 1) or depth is not None and len(depth) != len(reference):
        raise ValueError("invalid locus dimensions")
    sequence, reasons, mixed_offsets, used = [], [], set(), []
    for offset, base in enumerate(reference):
        pos = start0 + offset
        d = None if depth is None else int(depth[offset])
        reason = "coverage_unavailable" if d is None else "zero_depth" if d == 0 else "depth_below_5" if d < 5 else "ambiguous_reference" if base not in BASES else ""
        call = "N" if reason else base
        row = variants.get(pos)
        if row is not None:
            if d is None:
                raise ValueError("variant on unavailable scaffold")
            existing, caller_reason, mixed, counts = call_existing(row, base, d)
            if not reason:
                reason, call = caller_reason, existing
            if mixed:
                mixed_offsets.add(offset if strand == 1 else len(reference)-1-offset)
            used.append({"position0": pos, "reference_base": base, "depth": d, "consensus_base_forward": call,
                         "mask_reason": reason, "mixed_original_call": mixed, "counts": counts, "source": row})
        sequence.append(call)
        reasons.append(reason)
    forward = "".join(sequence)
    coding = forward if strand == 1 else revcomp(forward)
    return coding, reasons, mixed_offsets, used


def mask_runs(reasons, start0):
    result, lo = [], 0
    while lo < len(reasons):
        hi = lo + 1
        while hi < len(reasons) and reasons[hi] == reasons[lo]:
            hi += 1
        if reasons[lo]:
            result.append({"start0": start0+lo, "end0": start0+hi, "length": hi-lo, "reason": reasons[lo]})
        lo = hi
    return result


def protein_effects(ref_nt, consensus_nt, mixed_offsets):
    ref_aa, sample_aa = translate(ref_nt), translate(consensus_nt)
    changes, stops = [], []
    for i, (r, s) in enumerate(zip(ref_aa, sample_aa)):
        mixed = any(j in mixed_offsets for j in range(3*i, 3*i+3))
        if s != "X" and s != r:
            changes.append({"aa_position1": i+1, "reference_aa": r, "consensus_aa": s, "codon_mixed": mixed})
        if s == "*" and i < len(sample_aa)-1 and r != "*":
            stops.append({"aa_position1": i+1, "codon_mixed": mixed})
    return ref_aa, sample_aa, changes, stops


def seed_positions(reference_protein, seed, rules):
    """Reference-only alignment, never fitted to a sample consensus."""
    ref = reference_protein[:-1] if reference_protein.endswith("*") else reference_protein
    seed = seed[:-1] if seed.endswith("*") else seed
    if ref == seed:
        return "identical_reference_seed", {i: i for i in range(len(ref))}
    import Bio
    from Bio.Align import PairwiseAligner
    if Bio.__version__ != rules["biopython_version"]:
        raise ValueError("unfrozen Biopython version for seed alignment")
    cfg = rules["seed_alignment"]
    aligner = PairwiseAligner(mode=cfg["mode"], match_score=cfg["match_score"], mismatch_score=cfg["mismatch_score"],
                              open_gap_score=cfg["open_gap_score"], extend_gap_score=cfg["extend_gap_score"])
    aligns = iter(aligner.align(ref, seed))
    first = next(aligns)
    if next(aligns, None) is not None:
        return "ambiguous_optimal_reference_seed_alignment", {}
    mapping = {}
    for (r0, r1), (s0, s1) in zip(first.aligned[0], first.aligned[1]):
        if r1-r0 != s1-s0:
            raise ValueError("unexpected aligned block")
        mapping.update({int(r0+i): int(s0+i) for i in range(r1-r0)})
    return "unique_fixed_reference_seed_alignment", mapping


def dense_locus_cov(handle, contig, start, end, reference_length):
    if handle is None or contig + "::0" not in handle:
        return None
    import numpy as np
    data = handle[contig + "::0"][()]
    if data.ndim != 2 or data.shape[0] != 2 or not np.issubdtype(data.dtype, np.integer):
        raise ValueError("unsupported filtered covT format")
    values, positions = data
    if np.any(values <= 0) or np.any(positions < 0) or np.any(positions >= reference_length) or len(np.unique(positions)) != len(positions):
        raise ValueError("invalid sparse filtered depth")
    mask = (positions >= start) & (positions < end)
    result = np.zeros(end-start, dtype=np.int64)
    result[positions[mask]-start] = values[mask]
    return result


VARIANT_FIELDS = ["sample_id", "locus_id", "trait_id", "role", "assembly_accession", "contig_id", "position0", "coding_position1", "strand", "reference_base_forward", "original_consensus_base", "original_variant_base", "consensus_base_forward", "A", "C", "G", "T", "depth", "original_allele_count", "original_class", "original_cryptic", "original_con_freq", "original_var_freq", "original_ref_freq", "mixed_original_call", "mask_reason"]
HOLE_FIELDS = ["sample_id", "locus_id", "trait_id", "role", "contig_id", "start0", "end0", "length", "reason"]
NT_FIELDS = ["sample_id", "locus_id", "trait_id", "role", "coding_position1", "reference_nt", "consensus_nt", "mixed_original_call"]
AA_FIELDS = ["sample_id", "locus_id", "trait_id", "role", "aa_position1", "reference_aa", "consensus_aa", "differs_from_reference", "seed_position1", "seed_aa", "differs_from_seed", "reference_already_differs_from_seed", "seed_alignment_status", "codon_mixed", "candidate_premature_stop"]
TRAIT_FIELDS = ["sample_id", "trait_id", "role", "fixed_locus_denominator", "original_supported_loci", "fully_depth5_callable_loci", "complete_substitution_consensus_loci", "mixed_loci", "premature_stop_candidate_loci", "annotation_conflict_loci", "intact_sample_orf_asserted"]
LOCUS_FIELDS = ["sample_id", "locus_id", "trait_id", "role", "assembly_accession", "contig_id", "start0", "end0", "strand", "length", "evidence_type", "annotation_conflict", "original_gene_status", "original_carrier_breadth_depth5", "original_carrier_background_median_depth", "original_carrier_background_supported_cds_count", "original_carrier_background_eligible", "original_carrier_eligible", "mean_depth", "breadth_depth1", "depth5_callable_bases", "depth5_callable_fraction", "zero_depth_bases", "low_depth_positive_bases", "consensus_N_bases", "consensus_resolved_fraction", "full_depth5_callable", "complete_substitution_consensus", "original_variant_rows", "original_SNS_rows", "original_mixed_rows", "mixed_site_count", "caller_masked_bases", "reference_nt_substitutions", "reference_aa_substitutions", "seed_aa_differences", "seed_alignment_status", "candidate_premature_stop_count", "candidate_premature_stop_aa_positions1", "candidate_premature_stop_mixed_codon_count", "start_codon_state", "terminal_codon_state", "reference_cds_sha256", "consensus_nt_sha256", "consensus_aa_sha256", "sequence_label", "flags", "indels_assessed", "full_length_phasing_assessed", "intact_sample_orf_asserted", "functional_activity_asserted"]
ROSTER_FIELDS = ["locus_id", "trait_id", "role", "assembly_accession", "contig_id", "start0", "end0", "strand", "annotation_conflict", "evidence_type", "reference_cds_sha256"]
TABLES = {"locus_summary.tsv": LOCUS_FIELDS, "variants.tsv": VARIANT_FIELDS, "coverage_holes.tsv": HOLE_FIELDS,
          "nucleotide_substitutions.tsv": NT_FIELDS, "amino_acid_comparisons.tsv": AA_FIELDS,
          "trait_denominators.tsv": TRAIT_FIELDS}


def locus_result(sample, target, reference_forward, depth, variants, evidence, carrier, upstream_gene, seed, rules):
    coding, reasons, mixed, used = make_consensus(reference_forward, depth, variants, target["start0"], target["strand"])
    ref_coding = reference_forward if target["strand"] == 1 else revcomp(reference_forward)
    ref_aa, consensus_aa, aa_changes, stops = protein_effects(ref_coding, coding, mixed)
    n = len(reference_forward)
    mean = None if depth is None else sum(int(x) for x in depth) / n
    breadth1 = None if depth is None else sum(int(x) >= 1 for x in depth) / n
    callable_bases = None if depth is None else sum(int(x) >= 5 for x in depth)
    breadth5 = None if callable_bases is None else callable_bases / n
    ordinary = target["evidence_type"] == "ordinary_reference_CDS"
    prefix = "gene_" if ordinary else ""
    for key, value in (("mean_depth", mean), ("breadth_depth1", breadth1), ("breadth_depth5", breadth5)):
        eq(evidence[prefix+key], value)
    if ordinary and upstream_gene is not None:
        if upstream_gene["scaffold"] != target["contig_id"] or int(upstream_gene["start"]) != target["start0"] or int(upstream_gene["end"])+1 != target["end0"]:
            raise ValueError("upstream gene coordinate mismatch")
        for key, value in (("coverage", mean), ("breadth", breadth1), ("breadth_minCov", breadth5)):
            actual = upstream_gene[key]
            if key == "breadth_minCov" and actual in ("", "nan", "NaN", "NA") and value == 0:
                actual = 0
            eq(actual, value)
    elif ordinary and mean is not None and mean > 0:
        raise ValueError("positive locus has no original gene_info row")
    gene_status = evidence["gene_status" if ordinary else "region_status"]
    if not ordinary and not boolean(evidence["annotation_conflict"]):
        raise ValueError("NanH2 annotation conflict was lost")
    alignment_status, seed_map = seed_positions(ref_aa, seed, rules)
    identity = {"sample_id": sample, **{k: target[k] for k in ("locus_id", "trait_id", "role")}}
    holes = [{**identity, "contig_id": target["contig_id"], **r} for r in mask_runs(reasons, target["start0"])]
    vrows = []
    for v in used:
        src = v["source"]
        vrows.append({**identity, "assembly_accession": target["assembly_accession"], "contig_id": target["contig_id"],
                      "position0": v["position0"], "coding_position1": v["position0"]-target["start0"]+1 if target["strand"] == 1 else target["end0"]-v["position0"],
                      "strand": target["strand"], "reference_base_forward": v["reference_base"],
                      "original_consensus_base": src["con_base"], "original_variant_base": src["var_base"],
                      "consensus_base_forward": v["consensus_base_forward"], **v["counts"], "depth": v["depth"],
                      "original_allele_count": int(src["allele_count"]), "original_class": src["class"], "original_cryptic": boolean(src["cryptic"]),
                      "original_con_freq": src.get("con_freq", ""), "original_var_freq": src.get("var_freq", ""), "original_ref_freq": src.get("ref_freq", ""),
                      "mixed_original_call": v["mixed_original_call"], "mask_reason": v["mask_reason"]})
    nts = [{**identity, "coding_position1": i+1, "reference_nt": r, "consensus_nt": c, "mixed_original_call": i in mixed}
           for i, (r, c) in enumerate(zip(ref_coding, coding)) if c in BASES and r != c]
    aas, seed_diff_count = [], 0 if seed_map else None
    for i, (r, c) in enumerate(zip(ref_aa, consensus_aa)):
        seed_pos = seed_map.get(i)
        seed_aa = seed[seed_pos] if seed_pos is not None else ""
        seed_diff = c != "X" and seed_pos is not None and c != seed_aa
        if seed_diff_count is not None:
            seed_diff_count += int(seed_diff)
        if c != "X" and (r != c or seed_diff):
            aas.append({**identity, "aa_position1": i+1, "reference_aa": r, "consensus_aa": c, "differs_from_reference": r != c,
                        "seed_position1": seed_pos+1 if seed_pos is not None else None, "seed_aa": seed_aa,
                        "differs_from_seed": seed_diff if seed_pos is not None else None,
                        "reference_already_differs_from_seed": r != seed_aa if seed_pos is not None else None,
                        "seed_alignment_status": alignment_status, "codon_mixed": any(x in mixed for x in range(i*3, i*3+3)),
                        "candidate_premature_stop": c == "*" and i < len(consensus_aa)-1 and r != "*"})
    flags = ["indels_unassessed", "not_full_length_phased"]
    if depth is None:
        flags.append("filtered_coverage_unavailable")
    elif breadth5 < 1:
        flags.append("coverage_incomplete")
    caller_masked = sum(bool(r) and r not in ("zero_depth", "depth_below_5", "coverage_unavailable") for r in reasons)
    if caller_masked:
        flags.append("caller_or_reference_unresolved")
    if mixed:
        flags.append("mixed_original_alleles")
    if stops:
        flags.append("premature_stop_candidate")
    if target["annotation_conflict"]:
        flags.append("historical_annotation_conflict")
    if "ambiguous" in alignment_status:
        flags.append("seed_alignment_ambiguous")
    carrier_ok = carrier["breadth_depth5"] not in (None, "") and float(carrier["breadth_depth5"]) >= .5 and boolean(carrier["background_eligible"])
    start_state = "unresolved" if "N" in coding[:3] else "recognized_table11_start" if coding[:3] in STARTS_11 else "non_start_substitution_candidate"
    terminal_state = "unresolved" if "N" in coding[-3:] else "stop_codon" if consensus_aa[-1] == "*" else "terminal_stop_not_present_in_substitution_consensus"
    if start_state == "non_start_substitution_candidate":
        flags.append("start_codon_substitution_candidate")
    if terminal_state == "terminal_stop_not_present_in_substitution_consensus":
        flags.append("terminal_stop_substitution_candidate")
    summary = {**identity, **{k: target[k] for k in ("assembly_accession", "contig_id", "start0", "end0", "strand", "evidence_type", "annotation_conflict", "reference_cds_sha256")},
               "length": n, "original_gene_status": gene_status,
               "original_carrier_breadth_depth5": carrier["breadth_depth5"], "original_carrier_background_median_depth": carrier["background_median_depth"],
               "original_carrier_background_supported_cds_count": carrier["background_supported_cds_count"], "original_carrier_background_eligible": boolean(carrier["background_eligible"]),
               "original_carrier_eligible": carrier_ok, "mean_depth": mean, "breadth_depth1": breadth1,
               "depth5_callable_bases": callable_bases, "depth5_callable_fraction": breadth5,
               "zero_depth_bases": None if depth is None else sum(int(x) == 0 for x in depth),
               "low_depth_positive_bases": None if depth is None else sum(0 < int(x) < 5 for x in depth),
               "consensus_N_bases": coding.count("N"), "consensus_resolved_fraction": (n-coding.count("N"))/n,
               "full_depth5_callable": breadth5 == 1, "complete_substitution_consensus": "N" not in coding,
               "original_variant_rows": len(used), "original_SNS_rows": sum(v["source"]["class"] == "SNS" for v in used),
               "original_mixed_rows": sum(v["mixed_original_call"] for v in used), "mixed_site_count": len(mixed),
               "caller_masked_bases": caller_masked, "reference_nt_substitutions": len(nts), "reference_aa_substitutions": len(aa_changes),
               "seed_aa_differences": seed_diff_count, "seed_alignment_status": alignment_status,
               "candidate_premature_stop_count": len(stops), "candidate_premature_stop_aa_positions1": ";".join(str(r["aa_position1"]) for r in stops),
               "candidate_premature_stop_mixed_codon_count": sum(r["codon_mixed"] for r in stops), "start_codon_state": start_state, "terminal_codon_state": terminal_state,
               "consensus_nt_sha256": hashlib.sha256(coding.encode()).hexdigest(), "consensus_aa_sha256": hashlib.sha256(consensus_aa.encode()).hexdigest(),
               "sequence_label": LABEL, "flags": ";".join(flags), "indels_assessed": False, "full_length_phasing_assessed": False,
               "intact_sample_orf_asserted": False, "functional_activity_asserted": False}
    return summary, vrows, holes, nts, aas, coding, consensus_aa


def denominators(sample, loci):
    grouped = defaultdict(list)
    for row in loci:
        grouped[(row["trait_id"], row["role"])].append(row)
    result = []
    for (trait, role), values in sorted(grouped.items()):
        result.append({"sample_id": sample, "trait_id": trait, "role": role, "fixed_locus_denominator": len(values),
                       "original_supported_loci": sum(r["original_gene_status"] in ("supported", "historical_region_supported_annotation_conflict") for r in values),
                       "fully_depth5_callable_loci": sum(r["full_depth5_callable"] for r in values),
                       "complete_substitution_consensus_loci": sum(r["complete_substitution_consensus"] for r in values),
                       "mixed_loci": sum(r["mixed_site_count"] > 0 for r in values),
                       "premature_stop_candidate_loci": sum(r["candidate_premature_stop_count"] > 0 for r in values),
                       "annotation_conflict_loci": sum(r["annotation_conflict"] for r in values), "intact_sample_orf_asserted": False})
    return result


def run_audit(base, sample, out):
    if out.exists():
        raise FileExistsError("use a fresh output directory")
    out.mkdir(parents=True)
    try:
        import Bio
        import h5py
        import numpy
        bindings = Bindings()
        rules = read_json(bindings.add(HERE / "PROFILE_RULES.json"))
        bindings.add(Path(__file__))
        if Bio.__version__ != rules["biopython_version"]:
            raise ValueError("unfrozen Bio version")
        roster, panel, seeds, plan = load_roster(base, bindings, rules)
        run, genes, regions, carriers, cov_path, raw_path, cumulative_path, gene_info, rdic = load_profile(base, sample, bindings, plan, rules)
        raw = select_variants(raw_path, roster)
        variants = select_variants(cumulative_path, roster, True)
        crosscheck_variant_tables(raw, variants)
        result_tables = {name: [] for name in TABLES}
        nt_fasta, aa_fasta = [], []
        handle = h5py.File(cov_path, "r") if cov_path else None
        try:
            if handle:
                for key in handle:
                    contig, mm = key.rsplit("::", 1)
                    if contig not in panel or int(mm) != 0:
                        raise ValueError("unrecognized covT scaffold/mm bucket")
            for t in roster:
                reference = panel[t["contig_id"]][t["start0"]:t["end0"]]
                depth = dense_locus_cov(handle, t["contig_id"], t["start0"], t["end0"], len(panel[t["contig_id"]]))
                vr = {pos: row for (contig, pos), row in variants.items() if contig == t["contig_id"] and t["start0"] <= pos < t["end0"]}
                ev = genes[t["locus_id"]] if t["evidence_type"] == "ordinary_reference_CDS" else regions[t["locus_id"]]
                if ev["trait_id"] != t["trait_id"] or ev["assembly_accession"] != t["assembly_accession"] or ev["contig_id"] != t["contig_id"]:
                    raise ValueError("target evidence identity mismatch")
                result = locus_result(sample, t, reference, depth, vr, ev, carriers[t["assembly_accession"]], gene_info.get(t["locus_id"]), seeds[t["trait_id"]], rules)
                summary, vs, holes, nts, aas, coding, peptide = result
                result_tables["locus_summary.tsv"].append(summary)
                for name, rows in (("variants.tsv", vs), ("coverage_holes.tsv", holes), ("nucleotide_substitutions.tsv", nts), ("amino_acid_comparisons.tsv", aas)):
                    result_tables[name].extend(rows)
                header = ">" + sample + "|" + t["locus_id"] + " " + LABEL + " role=" + t["role"] + " annotation_conflict=" + str(t["annotation_conflict"])
                nt_fasta.append(header + "\n" + "\n".join(coding[i:i+80] for i in range(0, len(coding), 80)) + "\n")
                aa_fasta.append(header + "\n" + "\n".join(peptide[i:i+80] for i in range(0, len(peptide), 80)) + "\n")
        finally:
            if handle:
                handle.close()
        result_tables["trait_denominators.tsv"] = denominators(sample, result_tables["locus_summary.tsv"])
        if len(result_tables["locus_summary.tsv"]) != 16 or len(result_tables["trait_denominators.tsv"]) != 9:
            raise ValueError("incomplete output denominator")
        for name, rows in result_tables.items():
            write_tsv(out / name, rows, TABLES[name])
        write_tsv(out / "ROSTER.tsv", roster, ROSTER_FIELDS)
        (out / "substitution_consensus.fna").write_text("".join(nt_fasta), encoding="ascii")
        (out / "substitution_consensus.faa").write_text("".join(aa_fasta), encoding="ascii")
        bindings.stable()
        files = list(TABLES) + ["ROSTER.tsv", "substitution_consensus.fna", "substitution_consensus.faa"]
        receipt = {"schema": rules["schema"], "status": "COMPLETE", "sample_id": sample, "completed_utc": datetime.now(timezone.utc).isoformat(),
                   "code_sha256": sha(Path(__file__)), "rules_sha256": sha(HERE / "PROFILE_RULES.json"), "roster_sha256": sha(out / "ROSTER.tsv"),
                   "profile_status": run["profile_status"], "source_files": list(bindings.files.values()), "rdic_bound_not_loaded": rdic,
                   "runtime": {"python": platform.python_version(), "biopython": Bio.__version__, "numpy": numpy.__version__, "h5py": h5py.__version__},
                   "denominators": {"loci": 16, "positive_loci_including_historical": 12, "comparator_loci": 4, "positive_traits": 6, "comparator_traits": 3},
                   "row_counts": {name: len(rows) for name, rows in result_tables.items()},
                   "outputs": [{"path": name, "bytes": (out/name).stat().st_size, "sha256": sha(out/name)} for name in files],
                   "source_files_modified": False, "bam_read": False, "rdic_loaded": False, "phenotype_data_read": False,
                   "new_variant_calling": False, "indels_assessed": False, "full_length_phasing_assessed": False, "intact_sample_orf_asserted": False,
                   "sequence_label": LABEL, "limitations": ["Reference mapping and accepted-read filters constrain detection", "Absent SNV rows are reference-filled only at callable original-caller-resolved positions", "Mixed non-tied consensus follows original caller but does not phase a full-length strain", "No indel evidence or complete sample ORF is established"]}
        (out / "SUMMARY.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
        (out / "COMPLETE").write_text(sha(out / "SUMMARY.json") + "\n", encoding="ascii")
        return {"status": "COMPLETE", "sample_id": sample, "loci": 16, "output": str(out)}
    except Exception as exc:
        (out / "FAILED.json").write_text(json.dumps({"status": "FAILED", "sample_id": sample, "error": type(exc).__name__ + ": " + str(exc), "utc": datetime.now(timezone.utc).isoformat()}, indent=2) + "\n")
        raise


def collect(sample_manifest, audit_root, out):
    if out.exists():
        raise FileExistsError("use a fresh collection directory")
    out.mkdir(parents=True)
    try:
        bindings = Bindings()
        bindings.add(sample_manifest)
        bindings.add(Path(__file__))
        bindings.add(HERE/"PROFILE_RULES.json")
        samples = [r["sample_id"] for r in table(sample_manifest)]
        if len(samples) != 50 or len(set(samples)) != 50:
            raise ValueError("collection requires exactly 50 unique sample IDs")
        combined = {name: [] for name in TABLES}
        sources, signatures, reference_roster = [], set(), None
        nt, aa = [], []
        for sample in samples:
            directory = audit_root / sample
            bindings.add(directory/"SUMMARY.json")
            bindings.add(directory/"COMPLETE")
            summary = read_json(directory / "SUMMARY.json")
            if summary["status"] != "COMPLETE" or summary["sample_id"] != sample or (directory / "COMPLETE").read_text().strip() != sha(directory / "SUMMARY.json"):
                raise ValueError("incomplete/mismatched sample audit")
            if summary["code_sha256"] != sha(Path(__file__)) or summary["rules_sha256"] != sha(HERE / "PROFILE_RULES.json"):
                raise ValueError("collection code/rules drift")
            expected_outputs = set(TABLES) | {"ROSTER.tsv", "substitution_consensus.fna", "substitution_consensus.faa"}
            if len(summary["outputs"]) != len(expected_outputs) or {r["path"] for r in summary["outputs"]} != expected_outputs:
                raise ValueError("collection output roster differs")
            signatures.add((summary["code_sha256"], summary["rules_sha256"], summary["roster_sha256"]))
            if len(signatures) != 1:
                raise ValueError("mixed audit freezes")
            roster = table(directory / "ROSTER.tsv")
            reference_roster = roster if reference_roster is None else reference_roster
            if roster != reference_roster:
                raise ValueError("mixed candidate rosters")
            expected_ids = {r["locus_id"] for r in roster}
            for entry in summary["outputs"]:
                bindings.add(directory / entry["path"], entry)
            for name in TABLES:
                rows = table(directory / name)
                if len(rows) != summary["row_counts"][name] or any(r["sample_id"] != sample for r in rows):
                    raise ValueError("audit table identity/count mismatch")
                if name == "locus_summary.tsv" and (len(rows) != 16 or {r["locus_id"] for r in rows} != expected_ids):
                    raise ValueError("incomplete 16-locus grid")
                combined[name].extend(rows)
            nt.append((directory / "substitution_consensus.fna").read_text())
            aa.append((directory / "substitution_consensus.faa").read_text())
            sources.append(record(directory / "SUMMARY.json"))
        if len(combined["locus_summary.tsv"]) != 800 or len(combined["trait_denominators.tsv"]) != 450:
            raise ValueError("incomplete final denominators")
        for name, rows in combined.items():
            write_tsv(out/name, rows, TABLES[name])
        write_tsv(out / "ROSTER.tsv", reference_roster, ROSTER_FIELDS)
        (out / "substitution_consensus.fna").write_text("".join(nt), encoding="ascii")
        (out / "substitution_consensus.faa").write_text("".join(aa), encoding="ascii")
        bindings.stable()
        files = list(TABLES) + ["ROSTER.tsv", "substitution_consensus.fna", "substitution_consensus.faa"]
        receipt = {"schema": "reference-backed-substitution-audit.collection.v1", "status": "COMPLETE", "samples": samples,
                   "sample_manifest": record(sample_manifest), "input_summaries": sources, "input_files_verified_before_after": list(bindings.files.values()), "completed_utc": datetime.now(timezone.utc).isoformat(),
                   "code_sha256": sha(Path(__file__)), "rules_sha256": sha(HERE/"PROFILE_RULES.json"), "roster_sha256": sha(out/"ROSTER.tsv"),
                   "denominators": {"samples": 50, "loci_per_sample": 16, "locus_rows": 800, "positive_locus_rows": 600, "comparator_locus_rows": 200, "trait_rows": 450},
                   "row_counts": {name: len(rows) for name, rows in combined.items()},
                   "outputs": [{"path": name, "bytes": (out/name).stat().st_size, "sha256": sha(out/name)} for name in files],
                   "phenotype_data_read": False, "indels_assessed": False, "intact_sample_orf_asserted": False, "sequence_label": LABEL}
        (out/"SUMMARY.json").write_text(json.dumps(receipt, indent=2)+"\n")
        (out/"COMPLETE").write_text(sha(out/"SUMMARY.json")+"\n")
        return {"status": "COMPLETE", "samples": 50, "locus_rows": 800}
    except Exception as exc:
        (out/"FAILED.json").write_text(json.dumps({"status": "FAILED", "error": type(exc).__name__+": "+str(exc)}, indent=2)+"\n")
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser("audit")
    one.add_argument("--base", type=Path, required=True)
    one.add_argument("--sample", required=True)
    one.add_argument("--out", type=Path, required=True)
    merge = sub.add_parser("collect")
    merge.add_argument("--sample-manifest", type=Path, required=True)
    merge.add_argument("--audit-root", type=Path, required=True)
    merge.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "audit":
        if not args.sample or any(c not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-" for c in args.sample):
            raise ValueError("invalid sample ID")
        result = run_audit(args.base.resolve(), args.sample, args.out.resolve())
    else:
        result = collect(args.sample_manifest.resolve(), args.audit_root.resolve(), args.out.resolve())
    print(json.dumps(result))


if __name__ == "__main__":
    main()
