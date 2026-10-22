"""Readback primitives for one pre-specified candidate; no file-opening CLI.

Inputs are in-memory records/handles supplied by the separately configured
driver. Tests use synthetic records only. Sequence support uses ORIGINAL read
qualities; covT parity uses the original pileup counter after overlap handling.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, asdict
import hashlib
import re


SAMPLE_D = "SAMPLE_FLANK_DELETION"
SAMPLE_R = "SAMPLE_FLANK_RETAINED"
REF_PREFIX = "REFERENCE_RETAINED::"
REF_CONSUMES = frozenset("MDN=X")
QUERY_CONSUMES = frozenset("MIS=X")


def sha_sequence(sequence):
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()


def cigar_parts(cigar):
    parts = [(int(n), op) for n, op in re.findall(r"([0-9]+)([MIDNSHP=X])", cigar)]
    if not parts or "".join(str(n) + op for n, op in parts) != cigar or any(n <= 0 for n, _ in parts):
        raise ValueError("Invalid CIGAR")
    return parts


@dataclass(frozen=True)
class Read:
    sample_id: str
    qname: str
    reference: str
    start0: int
    cigar: str
    sequence: str
    qualities: tuple[int, ...] | None
    flag: int
    mapq: int

    def validate(self):
        if not self.sample_id or not self.qname or not self.reference or self.start0 < 0:
            raise ValueError("Missing read identity or invalid start")
        parts = cigar_parts(self.cigar)
        if sum(n for n, op in parts if op in QUERY_CONSUMES) != len(self.sequence):
            raise ValueError("CIGAR/query length mismatch")
        if set(self.sequence) - set("ACGTN"):
            raise ValueError("Read sequence must be alignment-oriented ACGTN")
        if self.qualities is not None:
            if len(self.qualities) != len(self.sequence) or any(type(q) is not int or not 0 <= q <= 93 for q in self.qualities):
                raise ValueError("Invalid original qualities")
        if not 0 <= self.mapq <= 255 or self.flag < 0:
            raise ValueError("Invalid mapping metadata")
        return parts

    @property
    def end0(self):
        return self.start0 + sum(n for n, op in cigar_parts(self.cigar) if op in REF_CONSUMES)

    @property
    def mate(self):
        return 1 if self.flag & 64 and not self.flag & 128 else 2 if self.flag & 128 and not self.flag & 64 else 0


@dataclass(frozen=True)
class Assay:
    candidate_id: str
    sample_id: str
    reference: str
    start0: int
    end0: int
    gap_intervals0: tuple[tuple[int, int], ...]
    # Four frozen comparison strings are reduced to the same anchor/event window.
    templates: tuple[tuple[str, str], ...]
    anchor_length: int
    min_base_quality: int
    min_mapq: int
    reference_to_candidate_orientation: int = 1

    def validate(self):
        if not self.candidate_id or self.start0 < 0 or self.end0 <= self.start0:
            raise ValueError("Invalid assay identity/interval")
        if self.anchor_length < 1 or self.min_base_quality < 0 or self.min_mapq < 0:
            raise ValueError("Invalid explicit rules")
        if self.reference_to_candidate_orientation not in (-1, 1):
            raise ValueError("Invalid orientation")
        ts = dict(self.templates)
        if len(ts) != len(self.templates) or SAMPLE_D not in ts or SAMPLE_R not in ts:
            raise ValueError("Missing or duplicate comparison contexts")
        if len(ts) != 4 or sum(k.startswith(REF_PREFIX) for k in ts) != 2 or len(set(ts.values())) != 4 or any(not s for s in ts.values()):
            raise ValueError("The four contexts must be nonempty and sequence-distinguishable")
        if any(set(s) - set("ACGT") for s in ts.values()):
            raise ValueError("Unresolved comparison sequence")
        width = self.end0 - self.start0
        if len(ts[SAMPLE_R]) != width or not self.gap_intervals0:
            raise ValueError("Retained template length or gap missing")
        gap_lengths = {b-a for a, b in self.gap_intervals0}
        if len(gap_lengths) != 1 or next(iter(gap_lengths)) <= 0:
            raise ValueError("Inconsistent gap lengths")
        gap_length = next(iter(gap_lengths))
        if len(ts[SAMPLE_D]) != width-gap_length:
            raise ValueError("Deleted template length mismatch")
        for a, b in self.gap_intervals0:
            if not self.start0+self.anchor_length <= a < b <= self.end0-self.anchor_length:
                raise ValueError("Gap overlaps an exact anchor")
            x, y = a-self.start0, b-self.start0
            if ts[SAMPLE_R][:x] + ts[SAMPLE_R][y:] != ts[SAMPLE_D]:
                raise ValueError("Gap is not equivalent in sample-flank comparator")
        if any(len(s) != width for k, s in ts.items() if k != SAMPLE_D):
            raise ValueError("Retained context length mismatch")
        return self


def original_filter_reasons(read, min_mapq):
    reasons = []
    if not read.flag & 1 or read.mate == 0:
        reasons.append("invalid_pair_identity")
    if not read.flag & 2:
        reasons.append("not_proper_pair")
    if read.flag & (4 | 8 | 0x100 | 0x200 | 0x800):
        reasons.append("unmapped_failed_nonprimary_or_mate_unmapped")
    if read.mapq == 255:
        reasons.append("mapq255_unknown")
    elif read.mapq < min_mapq:
        reasons.append("mapq_below_original_minimum")
    return reasons


def _project(read, assay):
    """Projection uses original CIGAR coordinates; never aligns a new sequence."""
    r, q = read.start0, 0
    mapping, deletions, insertions, skips = {}, [], [], []
    for n, op in read.validate():
        if op in "M=X":
            for pos in range(max(r, assay.start0), min(r+n, assay.end0)):
                mapping[pos] = q + pos-r
        elif op == "D" and r < assay.end0 and r+n > assay.start0:
            interval = (max(r, assay.start0), min(r+n, assay.end0))
            if deletions and deletions[-1][1] == interval[0]:
                deletions[-1] = (deletions[-1][0], interval[1])
            else:
                deletions.append(interval)
        elif op == "N" and r < assay.end0 and r+n > assay.start0:
            skips.append((r, r+n))
        elif op == "I" and assay.start0 < r < assay.end0:
            insertions.append((r, q, q+n))
        if op in REF_CONSUMES:
            r += n
        if op in QUERY_CONSUMES:
            q += n
    return mapping, deletions, insertions, skips


def classify_read(read, assay, accepted_qnames):
    """Single-read original-BQ evidence; accepted_qnames is this reference's Rdic."""
    assay.validate()
    if read.sample_id != assay.sample_id or read.reference != assay.reference:
        raise ValueError("Sample/reference mismatch")
    if not isinstance(accepted_qnames, (set, frozenset)):
        raise TypeError("skip-mm Rdic member must be a QNAME set")
    reasons = original_filter_reasons(read, assay.min_mapq)
    accepted = read.qname in accepted_qnames
    direction = "reverse" if bool(read.flag & 16) ^ (assay.reference_to_candidate_orientation == -1) else "forward"
    result = {
        "candidate_id": assay.candidate_id, "sample_id": read.sample_id,
        "qname": read.qname, "reference": read.reference, "start0": read.start0,
        "end0": read.end0, "cigar": read.cigar, "flag": read.flag, "mate": read.mate,
        "mapq": read.mapq, "direction": direction,
        "selected_start_offset_from_gap0": read.start0-assay.gap_intervals0[0][0],
        "accepted_original_rdic": accepted, "original_filter_reasons": reasons,
        "main_record_eligible": accepted and not reasons,
        "sequence_sha256": sha_sequence(read.sequence),
        "qualities_sha256": None if read.qualities is None else hashlib.sha256(bytes(read.qualities)).hexdigest(),
        "status": "NOT_SPANNING_BOTH_ANCHORS", "allele": None,
        "structural_class": None, "matching_contexts": [], "min_original_span_bq": None,
    }
    mapping, deletions, insertions, skips = _project(read, assay)
    a = assay.anchor_length
    left = list(range(assay.start0, assay.start0+a))
    right = list(range(assay.end0-a, assay.end0))
    if any(p not in mapping for p in left+right):
        return result
    lq, rq = mapping[left[0]], mapping[right[0]]
    if [mapping[p] for p in left] != list(range(lq, lq+a)) or [mapping[p] for p in right] != list(range(rq, rq+a)):
        result["status"] = "ANCHOR_CIGAR_INTERRUPTED"
        return result
    if skips:
        result["status"] = "REFERENCE_SKIP_UNRESOLVED"
        return result
    query_end = mapping[right[-1]]+1
    sequence = read.sequence[lq:query_end]
    result["spanning_query_interval0"] = [lq, query_end]
    result["spanning_sequence_sha256"] = sha_sequence(sequence)
    if read.qualities is None:
        result["status"] = "ORIGINAL_QUALITIES_UNAVAILABLE"
        return result
    result["min_original_span_bq"] = min(read.qualities[lq:query_end])
    if result["min_original_span_bq"] < assay.min_base_quality:
        result["status"] = "ORIGINAL_BQ_BELOW_THRESHOLD"
        return result
    if not deletions and not insertions:
        structural = "RETAINED"
    elif not insertions and len(deletions) == 1 and deletions[0] in assay.gap_intervals0:
        structural = "DELETION"
    else:
        structural = "OTHER_CIGAR"
    result["structural_class"] = structural
    eligible_templates = [(k, s) for k, s in assay.templates
                          if (k == SAMPLE_D) == (structural == "DELETION")]
    if structural == "OTHER_CIGAR":
        result["status"] = "OTHER_CIGAR_UNRESOLVED"
        return result
    matches = [k for k, s in eligible_templates if sequence == s]
    result["matching_contexts"] = sorted(matches)
    if len(matches) > 1:
        result["status"] = "AMBIGUOUS_CONTEXT"
    elif len(matches) == 1:
        result["status"], result["allele"] = "ALLELE", matches[0]
    else:
        anchor_match = any(sequence[:a] == s[:a] and sequence[-a:] == s[-a:] for _, s in assay.templates)
        result["status"] = "OTHER_SEQUENCE_WITH_EXACT_ANCHORS" if anchor_match else "ANCHOR_SEQUENCE_MISMATCH"
    return result


def collapse_qnames(calls, channel):
    """Main and diagnostic channels remain separate; never add their denominators."""
    if channel not in ("main", "diagnostic"):
        raise ValueError("Unknown channel")
    if len({(c["candidate_id"], c["sample_id"]) for c in calls}) > 1:
        raise ValueError("Only one frozen sample/candidate per aggregation")
    selected = [c for c in calls if c["main_record_eligible"] == (channel == "main")]
    unique_records = {}
    for c in selected:
        key = (c["qname"], c["reference"], c["start0"], c["cigar"], c["flag"])
        if key in unique_records and unique_records[key] != c:
            raise ValueError("Inconsistent repeated alignment record")
        unique_records[key] = c
    by_name = defaultdict(list)
    for c in unique_records.values():
        by_name[c["qname"]].append(c)
    fragments = []
    for qname, group in sorted(by_name.items()):
        qualified = [c for c in group if c["status"] == "ALLELE"]
        alleles = sorted({c["allele"] for c in qualified})
        # A high-quality fully spanning unmatched context cannot silently support
        # the other mate's exact context. This includes retained-vs-deleted conflict.
        ambiguity = any(c["status"] in ("AMBIGUOUS_CONTEXT", "OTHER_SEQUENCE_WITH_EXACT_ANCHORS", "OTHER_CIGAR_UNRESOLVED", "ANCHOR_SEQUENCE_MISMATCH") for c in group)
        state = "UNRESOLVED"
        representative = None
        allele = None
        if len(alleles) > 1 or (qualified and ambiguity):
            state = "CONFLICT_EXCLUDED"
        elif len(alleles) == 1:
            state, allele = "ADJUDICABLE", alleles[0]
            # Root-approved: a qualifying mate1 is preferred; otherwise mate2.
            representative = min(qualified, key=lambda c: (c["mate"] != 1, c["reference"], c["start0"], c["cigar"], c["flag"]))
        fragments.append({
            "qname": qname, "channel": channel, "state": state, "allele": allele,
            "observed_alleles": alleles, "record_statuses": dict(Counter(c["status"] for c in group)),
            "references": sorted({c["reference"] for c in group}),
            "multiple_reference_mappings": len({c["reference"] for c in group}) > 1,
            "qualifying_mates": sorted({c["mate"] for c in qualified}),
            "representative": representative,
            "record_count_after_duplicate_fetch_removal": len(group),
        })
    return fragments


@dataclass(frozen=True)
class TechnicalRules:
    min_qnames: int
    min_selected_read_starts: int
    min_each_direction: int
    min_sample_flank_deletion_fraction: float
    coexistence_reference_qnames: int

    def validate(self):
        if any(type(x) is not int or x < 1 for x in (self.min_qnames, self.min_selected_read_starts, self.min_each_direction, self.coexistence_reference_qnames)):
            raise ValueError("Invalid technical count rule")
        if not 0 < self.min_sample_flank_deletion_fraction <= 1:
            raise ValueError("Invalid technical fraction rule")


def summarize_fragments(fragments, rules, expected_contexts):
    rules.validate()
    if any(f["channel"] != "main" for f in fragments):
        raise ValueError("Technical gate accepts the main Rdic channel only")
    if len({f["qname"] for f in fragments}) != len(fragments):
        raise ValueError("Fragments must already be QNAME-deduplicated")
    adjudicable = [f for f in fragments if f["state"] == "ADJUDICABLE"]
    by_allele = Counter(f["allele"] for f in adjudicable)
    if len(expected_contexts) != 4 or len(set(expected_contexts)) != 4 or SAMPLE_D not in expected_contexts or SAMPLE_R not in expected_contexts or any(k not in expected_contexts for k in by_allele):
        raise ValueError("Expected exactly the four frozen contexts")
    for context in expected_contexts:
        by_allele[context] += 0
    deletion = [f for f in adjudicable if f["allele"] == SAMPLE_D]
    d, r = len(deletion), by_allele[SAMPLE_R]
    starts = {f["representative"]["selected_start_offset_from_gap0"] for f in deletion}
    reference_starts = {(f["representative"]["reference"], f["representative"]["start0"]) for f in deletion}
    directions = Counter(f["representative"]["direction"] for f in deletion)
    fraction = d/(d+r) if d+r else None
    native_total = sum(v for k, v in by_allele.items() if k.startswith(REF_PREFIX))
    passed = d >= rules.min_qnames and len(starts) >= rules.min_selected_read_starts and all(directions[x] >= rules.min_each_direction for x in ("forward", "reverse")) and fraction is not None and fraction >= rules.min_sample_flank_deletion_fraction
    return {
        "channel": "original_Rdic_main", "qnames_total": len(fragments),
        "qname_states": dict(Counter(f["state"] for f in fragments)),
        "qnames_by_allele_context": dict(sorted(by_allele.items())),
        "sample_flank_deletion_qnames": d, "sample_flank_retained_qnames": r,
        "sample_flank_adjudicable_qnames": d+r, "sample_flank_deletion_fraction": fraction,
        "all_context_adjudicable_qnames": len(adjudicable),
        "all_context_deletion_fraction_descriptive": d/len(adjudicable) if adjudicable else None,
        "reference_like_retained_qnames": native_total,
        "deletion_distinct_selected_read_starts": len(starts),
        "selected_read_start_definition": "alignment_start0_minus_homologous_gap_start0_across_both_references",
        "deletion_selected_start_offsets_from_gap0": sorted(starts),
        "deletion_distinct_reference_start_coordinates_descriptive": len(reference_starts),
        "deletion_selected_read_directions": {x: directions[x] for x in ("forward", "reverse")},
        "technical_label": "READ_SUPPORTED_ASSEMBLY_ALLELE_IN_SAMPLE_FLANK_CONTEXT" if passed else "TECHNICAL_SUPPORT_CRITERIA_NOT_MET",
        "coexistence_label": "COEXISTING_REFERENCE_LIKE_CONTEXT" if native_total >= rules.coexistence_reference_qnames else None,
        "rules": asdict(rules), "function_loss_asserted": False,
        "fixed_deletion_or_single_strain_frequency_asserted": False,
    }


def original_column_counts(column, accepted_qnames):
    """Equivalent to the skip-mm branch of original get_base_counts_mm.

    column.pileups must come from the original pysam call, which already applied
    overlap handling and BQ30. Do not apply raw-read BQ filtering here.
    """
    if not isinstance(accepted_qnames, (set, frozenset)):
        raise TypeError("Expected original skip-mm QNAME set")
    counts = {b: 0 for b in "ACTG"}
    for p in column.pileups:
        if not p.is_del and not p.is_refskip and p.alignment.query_name in accepted_qnames:
            base = p.alignment.query_sequence[p.query_position]
            if base in counts:
                counts[base] += 1
    return counts


def original_split_coverage(bam_handle, reference, start0, end0_inclusive, accepted_qnames, counter=original_column_counts):
    """Lazy adapter for an ALREADY OPEN authorized handle; this module opens no files."""
    if start0 < 0 or end0_inclusive < start0:
        raise ValueError("Invalid original split")
    columns = bam_handle.pileup(reference, truncate=True, max_depth=100000,
                               stepper="nofilter", compute_baq=True,
                               ignore_orphans=True, ignore_overlaps=True,
                               min_base_quality=30, start=start0, stop=end0_inclusive+1)
    for column in columns:
        if not start0 <= column.pos <= end0_inclusive:
            raise ValueError("Pileup escaped requested split")
        counts = counter(column, accepted_qnames)
        yield column.pos, sum(counts.values())


def compare_covt_depth(observed, expected_sparse, start0, end0):
    """Compare an original split, then crop externally. None is unavailable, not 0."""
    if expected_sparse is None:
        return {"status": "EXPECTED_COVT_UNAVAILABLE", "parity": None, "differences": []}
    actual = {}
    for pos, depth in observed:
        if pos in actual or not start0 <= pos < end0 or type(depth) is not int or depth < 0:
            raise ValueError("Invalid observed depth")
        actual[pos] = depth
    if any(not start0 <= p < end0 or type(d) is not int or d <= 0 for p, d in expected_sparse.items()):
        raise ValueError("Invalid expected sparse covT slice")
    diffs = [{"position0":p, "observed":actual.get(p,0), "expected":expected_sparse.get(p,0)}
             for p in range(start0,end0) if actual.get(p,0) != expected_sparse.get(p,0)]
    return {"status":"MATCH" if not diffs else "MISMATCH_BLOCK_EVENT_INTERPRETATION", "parity":not diffs, "differences":diffs}


def original_splits(length, window_length):
    """Exact inStrain 1.11 iterate_splits semantics; both ends inclusive."""
    if type(length) is not int or type(window_length) is not int or length < 1 or window_length < 1:
        raise ValueError("Invalid original split inputs")
    chunks = length//window_length+1
    chunk_length = int(length/chunks)
    start = 0
    for i in range(chunks):
        end = length-1 if i+1 == chunks else start+chunk_length-1
        yield start,end
        start += chunk_length


def assays_from_candidate(candidate, anchor_length, min_base_quality, min_mapq):
    """In-memory adapter for the frozen sequence receipt, not a count executor."""
    c = candidate
    if c["schema"] != "existing-reference-resolution.candidate.v1" or c["candidate_count"] != 1:
        raise ValueError("Unexpected candidate schema")
    observed, comparator = c["observed_allele"], c["reference_length_comparator"]
    if observed["kind"] != "OBSERVED_ASSEMBLY_SEQUENCE" or comparator["kind"] != "IN_SILICO_SAME_FLANK_INSERTION_COMPARATOR_NOT_OBSERVED" or comparator["not_a_repaired_observation"] is not True:
        raise ValueError("Observed/comparator provenance changed")
    for item in [observed, comparator]:
        if sha_sequence(item["sequence"]) != item["sha256"]:
            raise ValueError("Comparison sequence hash mismatch")
    j = observed["junction_in_window0"]
    x,y = comparator["insertion_interval_in_window0"]
    if x != j or comparator["sequence"][:x]+comparator["sequence"][y:] != observed["sequence"]:
        raise ValueError("Same-flank comparator is inconsistent")
    templates = [(SAMPLE_D, observed["sequence"][j-anchor_length:j+anchor_length]),
                 (SAMPLE_R, comparator["sequence"][x-anchor_length:y+anchor_length])]
    for ref in c["reference_comparisons"]:
        seq = ref["reference_window_sequence"]
        if sha_sequence(seq) != ref["reference_window_sha256"]:
            raise ValueError("Native reference sequence hash mismatch")
        a,b = ref["gap_in_window0"]
        templates.append((REF_PREFIX+ref["accession"],seq[a-anchor_length:b+anchor_length]))
    result = []
    for ref in c["reference_comparisons"]:
        if ref["reference_window_orientation"] != "reference_plus" or ref["gap"]["subject_step"] != 1:
            raise ValueError("Candidate adapter currently expects verified plus orientation")
        a,b = ref["deletion_interval0"]
        gaps = tuple((g["start0"],g["end0"]) for g in ref["equivalent_gap_intervals0"])
        assay = Assay(c["candidate_id"],c["sample_id"],ref["subject_id"],a-anchor_length,b+anchor_length,
                      gaps,tuple(templates),anchor_length,min_base_quality,min_mapq)
        assay.validate()
        result.append(assay)
    if len(result) != 2 or len({a.reference for a in result}) != 2:
        raise ValueError("Expected exactly two distinct original reference loci")
    return result
