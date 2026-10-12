"""Quantify marker/locus evidence from query-name-grouped Bowtie2 paired SAM.

Input contract: Bowtie2 --no-mixed --no-discordant -a; mapped records have
AS and YS tags. Each fragment's records must be contiguous. No biological calls
are made here. Reference annotations are nucleotide BLAST alignments, with
1-based inclusive query/subject coordinates and explicit gapped strings.

Each mate is tested separately against every near-best position. A feature is
resolved only when every such position projects sequenced bases onto that same
feature (>=50 clear bases for either channel). Identity also
requires a single near-best taxon group across the entire fragment. A mate can
support multiple adjacent features. A fragment contributes at most once to a
feature; support from one unambiguous mate suffices. Ambiguous fragments have
candidate evidence but no supporting mate, making the two counts disjoint.

Only one deterministic highest-scoring position per supporting mate contributes
coverage. Overlapping mate projections are unioned per feature/fragment; mean
depth therefore measures fragment depth, not twice-counted overlapping reads.
SAM reverse flags do not reverse reference coordinates; BLAST query and subject
directions determine feature projection. No insertion, deletion, skipped region,
soft clip or unsequenced mate gap contributes coverage.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

MARGINS = (0, 6, 12, 18)
IDENTITY_MIN_BASES = 50
LOCUS_MIN_BASES = 50
BIN_SIZE = 2048
EPSILON = 1e-9
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


@dataclass(frozen=True)
class Feature:
    feature_id: str
    channel: str
    taxon_group: str
    length: int
    eligible_identity: int
    pseudo: int


@dataclass
class Annotation:
    feature_id: str
    projection: dict[int, int]
    priority: tuple


@dataclass(frozen=True)
class Alignment:
    mate: int
    contig: str
    pos: int
    cigar: str
    reverse: bool
    score: float

    def position_key(self):
        return self.mate, self.contig, self.pos, self.cigar, self.reverse

    def best_key(self):
        return -self.score, self.contig, self.pos, self.cigar, self.reverse


def read_tsv(path, required):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        yield from reader


def load_panel(contig_path, feature_path, annotation_path):
    contigs = {}
    for r in read_tsv(contig_path, ("contig", "taxon_group", "accession")):
        if not r["contig"] or not r["taxon_group"] or r["contig"] in contigs:
            raise ValueError("Contigs require unique names and nonempty taxon groups")
        contigs[r["contig"]] = r["taxon_group"]
    features = {}
    for r in read_tsv(feature_path, ("feature_id", "channel", "taxon_group", "length", "eligible_identity", "pseudo")):
        f = Feature(r["feature_id"], r["channel"], r["taxon_group"], int(r["length"]), int(r["eligible_identity"]), int(r["pseudo"]))
        if not f.feature_id or f.feature_id in features or f.channel not in {"identity", "locus"}:
            raise ValueError("Features require unique IDs and identity/locus channels")
        if f.length <= 0 or f.eligible_identity not in {0, 1} or f.pseudo not in {0, 1}:
            raise ValueError(f"Invalid feature metadata: {f.feature_id}")
        if f.channel == "identity" and not f.taxon_group:
            raise ValueError(f"Identity marker lacks a taxon group: {f.feature_id}")
        features[f.feature_id] = f
    bins = defaultdict(list)
    annotations_n = 0
    required = ("contig", "feature_id", "subject_start1", "subject_end1", "query_start1", "query_end1", "qseq", "sseq", "pident")
    for r in read_tsv(annotation_path, required):
        if r["contig"] not in contigs or r["feature_id"] not in features:
            raise ValueError("Annotation refers to an unknown contig or feature")
        qstart, qend = int(r["query_start1"]), int(r["query_end1"])
        sstart, send = int(r["subject_start1"]), int(r["subject_end1"])
        qs, ss = r["qseq"], r["sseq"]
        pident = float(r["pident"])
        if not math.isfinite(pident) or not 0 <= pident <= 100:
            raise ValueError("Invalid annotation pident")
        if len(qs) != len(ss) or not qs or any(q == s == "-" for q, s in zip(qs, ss)):
            raise ValueError("Invalid nucleotide BLAST gapped strings")
        if sum(c != "-" for c in qs) != abs(qend - qstart) + 1 or sum(c != "-" for c in ss) != abs(send - sstart) + 1:
            raise ValueError("BLAST coordinates disagree with gapped-string lengths")
        if min(sstart, send, qstart, qend) < 1 or max(qstart, qend) > features[r["feature_id"]].length:
            raise ValueError("Annotation coordinates exceed feature bounds")
        qstep, sstep = (1 if qend >= qstart else -1), (1 if send >= sstart else -1)
        qpos, spos = qstart, sstart
        projection = {}
        for qbase, sbase in zip(qs, ss):
            if qbase != "-" and sbase != "-":
                projection[spos] = qpos
            if qbase != "-":
                qpos += qstep
            if sbase != "-":
                spos += sstep
        # Resolve overlapping HSP coordinate conflicts deterministically. A single
        # sequenced genomic base cannot inflate one feature via multiple HSPs.
        priority = (-pident, -len(projection), sstart, send, qstart, qend, qs, ss)
        annotation = Annotation(r["feature_id"], projection, priority)
        if projection:
            for bucket in range((min(projection) - 1) // BIN_SIZE, (max(projection) - 1) // BIN_SIZE + 1):
                bins[(r["contig"], bucket)].append(annotation)
        annotations_n += 1
    return contigs, features, bins, annotations_n


def sequenced_reference_positions(pos: int, cigar: str) -> list[int]:
    tokens = CIGAR_RE.findall(cigar)
    if not tokens or "".join(n + op for n, op in tokens) != cigar:
        raise ValueError(f"Invalid mapped CIGAR: {cigar}")
    positions = []
    current = pos
    for n, op in tokens:
        length = int(n)
        if length <= 0:
            raise ValueError(f"Zero-length CIGAR operation: {cigar}")
        if op in {"M", "=", "X"}:
            positions.extend(range(current, current + length))
            current += length
        elif op in {"D", "N"}:
            current += length
    return positions


def project_alignment_maps(alignment, bins, features=None):
    positions = sequenced_reference_positions(alignment.pos, alignment.cigar)
    annotations = {}
    for bucket in {(p - 1) // BIN_SIZE for p in positions}:
        for annotation in bins.get((alignment.contig, bucket), ()):
            annotations[id(annotation)] = annotation
    chosen = defaultdict(dict)
    for annotation in annotations.values():
        target = chosen[annotation.feature_id]
        for p in positions:
            q = annotation.projection.get(p)
            if q is not None and (p not in target or annotation.priority < target[p][0]):
                target[p] = (annotation.priority, q)
    full = {f: {item[1] for item in coordinates.values()} for f, coordinates in chosen.items() if coordinates}
    if features is None:
        return full, full, {}
    locus_ids_by_genomic_base = defaultdict(set)
    for feature_id, coordinates in chosen.items():
        if features[feature_id].channel == "locus":
            for genomic_base in coordinates:
                locus_ids_by_genomic_base[genomic_base].add(feature_id)
    clear, shared = {}, {}
    for feature_id, coordinates in chosen.items():
        if features[feature_id].channel != "locus":
            clear[feature_id] = full.get(feature_id, set())
            continue
        clear_bases, shared_bases = set(), set()
        for genomic_base, (_, feature_base) in coordinates.items():
            target = clear_bases if len(locus_ids_by_genomic_base[genomic_base]) == 1 else shared_bases
            target.add(feature_base)
        if clear_bases:
            clear[feature_id] = clear_bases
        if shared_bases:
            shared[feature_id] = shared_bases
    return full, clear, shared


def project_alignment(alignment, bins):
    """Raw coordinate projection; evidence resolution uses the clear/shared maps."""
    return project_alignment_maps(alignment, bins)[0]


def sam_groups(handle) -> Iterator[tuple[str, list[list[str]]]]:
    current_name, group = None, []
    for line_number, line in enumerate(handle, 1):
        if not line.strip() or line.startswith("@"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise ValueError(f"SAM line {line_number}: fewer than 11 columns")
        name = fields[0]
        if current_name is not None and name != current_name:
            yield current_name, group
            group = []
        current_name = name
        group.append(fields)
    if current_name is not None:
        yield current_name, group


def parse_fragment(name, records, contigs, counters):
    unique = {}
    for fields in records:
        counters["sam_alignment_records"] += 1
        flag = int(fields[1])
        if flag & 0x4:
            counters["unmapped_records"] += 1
            continue
        if flag & 0x800:
            counters["ignored_supplementary_records"] += 1
            continue
        if not flag & 0x1 or not flag & 0x2 or flag & 0x8 or bool(flag & 0x40) == bool(flag & 0x80):
            raise ValueError(f"{name}: mapped SAM record is not a complete concordant paired alignment")
        contig, pos = fields[2], int(fields[3])
        if contig not in contigs or pos < 1:
            raise ValueError(f"{name}: unknown contig or invalid POS")
        tags = {}
        for token in fields[11:]:
            pieces = token.split(":", 2)
            if len(pieces) == 3:
                tags[pieces[0]] = pieces[2]
        if "AS" not in tags or "YS" not in tags:
            raise ValueError(f"{name}: AS and YS tags are required for paired scores")
        score = float(tags["AS"]) + float(tags["YS"])
        if not math.isfinite(score):
            raise ValueError(f"{name}: non-finite paired alignment score")
        a = Alignment(1 if flag & 0x40 else 2, contig, pos, fields[5], bool(flag & 0x10), score)
        key = a.position_key()
        if key in unique:
            counters["duplicate_position_records_removed"] += 1
            if a.score > unique[key].score:
                unique[key] = a
        else:
            unique[key] = a
    alignments = list(unique.values())
    if alignments and {a.mate for a in alignments} != {1, 2}:
        raise ValueError(f"{name}: one mate is missing from the mapped fragment")
    counters["unique_alignment_positions"] += len(alignments)
    return alignments


def assess_fragment(alignments, contigs, features, bins, metrics, resolved_by_margin):
    projections = {a: project_alignment_maps(a, bins, features) for a in alignments}
    best_score = max(a.score for a in alignments)
    for margin in MARGINS:
        near = [a for a in alignments if best_score - a.score <= margin + EPSILON]
        by_mate = {mate: [a for a in near if a.mate == mate] for mate in (1, 2)}
        if any(not rows for rows in by_mate.values()):
            raise ValueError("AS+YS near-best paired scores disagree between mates")
        groups = {contigs[a.contig] for a in near}
        supported = defaultdict(set)
        shared_candidates = defaultdict(set)
        candidates = set()
        for mate_alignments in by_mate.values():
            best = min(mate_alignments, key=Alignment.best_key)
            full_maps = [projections[a][0] for a in mate_alignments]
            maps = [projections[a][1] for a in mate_alignments]
            union = set().union(*(set(m) for m in full_maps))
            for feature_id in union:
                f = features[feature_id]
                if f.channel == "identity" and (not f.eligible_identity or f.pseudo):
                    continue
                minimum = IDENTITY_MIN_BASES if f.channel == "identity" else LOCUS_MIN_BASES
                candidate_minimum = minimum if f.channel == "identity" else 1
                if not any(len(m.get(feature_id, ())) >= candidate_minimum for m in full_maps):
                    continue
                candidates.add(feature_id)
                if f.channel == "locus":
                    # Keep shared evidence in a separate ambiguity audit even if
                    # another >=50-base exclusive segment supports the feature.
                    for a in mate_alignments:
                        shared_candidates[feature_id].update(projections[a][2].get(feature_id, ()))
                if not all(len(m.get(feature_id, ())) >= minimum for m in maps):
                    continue
                if f.channel == "identity" and groups != {f.taxon_group}:
                    continue
                supported[feature_id].update(projections[best][1][feature_id])
        for feature_id in candidates:
            metric = metrics[(feature_id, margin)]
            if shared_candidates.get(feature_id):
                metric["shared_candidate_fragments"] += 1
                metric["shared_candidate_bases"] += len(shared_candidates[feature_id])
            bases = supported.get(feature_id)
            if bases:
                metric["fragments"] += 1
                metric["aligned_bases"] += len(bases)
                metric["covered"].update(bases)
                resolved_by_margin[margin][features[feature_id].channel].add(feature_id)
            else:
                metric["ambiguous_fragments"] += 1


def run(sam, contig_path, feature_path, annotation_path, sample, out):
    contigs, features, bins, annotations_n = load_panel(contig_path, feature_path, annotation_path)
    metrics = {(f, m): {"fragments": 0, "aligned_bases": 0, "ambiguous_fragments": 0, "shared_candidate_fragments": 0, "shared_candidate_bases": 0, "covered": set()}
               for f in features for m in MARGINS}
    counters = Counter()
    resolved = {m: {"identity": set(), "locus": set()} for m in MARGINS}
    handle = sys.stdin if str(sam) == "-" else (gzip.open(sam, "rt", encoding="utf-8") if str(sam).endswith(".gz") else open(sam, encoding="utf-8"))
    try:
        for name, records in sam_groups(handle):
            counters["total_fragments"] += 1
            alignments = parse_fragment(name, records, contigs, counters)
            if not alignments:
                continue
            counters["total_mapped_fragments"] += 1
            assess_fragment(alignments, contigs, features, bins, metrics, resolved)
    finally:
        if handle is not sys.stdin:
            handle.close()
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fields = ["sample", "channel", "feature_id", "taxon_group", "length", "eligible_identity", "pseudo", "margin", "fragments", "covered_bases", "breadth", "aligned_bases", "mean_depth", "ambiguous_fragments"]
    with (out / "gene_metrics.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for feature_id in sorted(features):
            f = features[feature_id]
            for margin in MARGINS:
                metric = metrics[(feature_id, margin)]
                covered_bases = len(metric["covered"])
                writer.writerow({"sample": sample, **vars(f), "margin": margin,
                                 "fragments": metric["fragments"], "covered_bases": covered_bases,
                                 "breadth": covered_bases / f.length, "aligned_bases": metric["aligned_bases"],
                                 "mean_depth": metric["aligned_bases"] / f.length,
                                 "ambiguous_fragments": metric["ambiguous_fragments"]})
    summary = {"sample": sample, "margins": list(MARGINS), "identity_min_sequenced_marker_bases_per_mate": IDENTITY_MIN_BASES,
               "locus_min_exclusive_genomic_bases_per_mate": LOCUS_MIN_BASES,
               **{k: counters[k] for k in ("total_fragments", "total_mapped_fragments", "sam_alignment_records", "unmapped_records", "ignored_supplementary_records", "duplicate_position_records_removed", "unique_alignment_positions")},
               "panel_contigs": len(contigs), "panel_features": len(features), "annotation_HSPs": annotations_n,
               "eligible_identity_markers": sum(f.channel == "identity" and f.eligible_identity and not f.pseudo for f in features.values()),
               "locus_features": sum(f.channel == "locus" for f in features.values()),
               "observed_features_by_margin": {str(m): {channel: len(ids) for channel, ids in channels.items()} for m, channels in resolved.items()},
               "paired_score": "AS + YS; ties included with numerical tolerance 1e-9",
               "coverage_unit": "unique projected feature positions per fragment; overlapping mate bases counted once",
               "ambiguous_fragments": "feature candidates with no unambiguous supporting mate; disjoint from fragments",
               "locus_taxon_group": "feature metadata only; no locus-to-subspecies assignment is inferred",
               "locus_shared_ambiguity": {str(m): {f: {"candidate_fragments": metric["shared_candidate_fragments"], "candidate_projected_bases": metric["shared_candidate_bases"]}
                                          for (f, mm), metric in metrics.items() if mm == m and metric["shared_candidate_fragments"]}
                                          for m in MARGINS},
               "locus_shared_ambiguity_definition": "Unresolved co-annotated genomic bases, excluded from positive coverage/depth. Candidate bases union alternative near-best positions per fragment, so are not a measured depth. This audit may overlap fragments with >=50 independent exclusive bases.",
               "unseen_features": "all feature x margin rows retained as zero evidence; zero is not an absence call",
               "biology_calls_made": False, "qname_contract": "all records for each QNAME must be contiguous"}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sam", required=True, help="SAM, SAM.gz or - for stdin; grouped by QNAME")
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--features", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    try:
        run(args.sam, args.contigs, args.features, args.annotations, args.sample, args.out)
    except (ValueError, KeyError) as error:
        parser.exit(2, f"map_evidence: {error}\n")


if __name__ == "__main__":
    main()
