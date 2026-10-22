"""Condition fixed-panel evidence on an internal block or a boundary placement.

This changes the observation scope, not the retained candidate ambiguity or
the v1 alignment rules. Query and simulated profiles require identical selection.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from itertools import groupby
import json
from pathlib import Path

from .structure_inference import UNEXPLAINED, _config, _metadata, _record, _eligible_pairs

SCHEMA = "glycan-evidence.structure-conditioning.research.v1"
PROTOCOL = "any_near_best_internal_contiguous_block_or_inner_boundary"


def _qualifies(pair, metadata, min_anchor):
    left, right = metadata[pair["candidate"]]["junctions"]
    internal = any(min(hi, right)-max(lo, left) >= min_anchor
                   for record in pair["records"] for lo, hi in record["blocks"])
    boundary = bool(pair["direct_junctions"] or pair["bridged_junctions"])
    return internal, boundary


def _distribution(histogram):
    n = sum(histogram.values())
    return dict(count=n, minimum=min(histogram) if n else None,
                maximum=max(histogram) if n else None,
                mean=sum(k*v for k, v in histogram.items())/n if n else None,
                histogram={str(k): v for k, v in sorted(histogram.items())})


def extract_locus_conditioned(sam_path, candidate_metadata, config=None, *, expected_full_counts=None):
    """Select QNAMEs using any near-best placement and retain their full features.

    An internal block overlaps [left_junction, right_junction) by at least
    min_anchor aligned reference bases without a CIGAR gap. Alternatively any
    near-best placement may have original direct/bridged junction evidence.
    A qualifying QNAME keeps every original near-best candidate and its original
    across-placement junction intersection, including nonqualifying candidates.

    counts and anchor_only_counts partition accepted full_counts. Unexplained
    local QNAMEs remain separate and are absent from the conditioned denominator.
    Coverage for selected fragments includes all their original eligible or
    near-best placements; inner breadth clips aligned blocks to the interval.
    """
    cfg, metadata = _config(config), _metadata(candidate_metadata)
    if cfg["mode"] != "paired_junction":
        raise ValueError("Conditioning requires paired_junction features")
    if any(len(item["junctions"]) != 2 for item in metadata.values()):
        raise ValueError("Conditioning requires exactly two ordered inner boundaries")
    full_counts, counts, anchor_only = Counter(), Counter(), Counter()
    original_audit, reasons, ledger, seen = Counter(), Counter(), [], set()
    arrays = {scope: {c: bytearray(m["length"]) for c, m in metadata.items()}
              for scope in ("original_eligible", "original_near", "selected_eligible", "selected_near")}
    fragments = {scope: Counter() for scope in arrays}
    query_lengths, template_lengths = Counter(), Counter()
    omitted = 0

    def mark(scope, placements):
        fragments[scope].update({p["candidate"] for p in placements})
        for pair in placements:
            for record in pair["records"]:
                for lo, hi in record["blocks"]:
                    arrays[scope][pair["candidate"]][lo:hi] = b"\1"*(hi-lo)

    header = {}
    with Path(sam_path).open(encoding="utf-8") as handle:
        first = None
        for line in handle:
            if line.startswith("@SQ\t"):
                tags = dict(field.split(":", 1) for field in line.rstrip().split("\t")[1:])
                if tags.get("SN") in header or "SN" not in tags or "LN" not in tags:
                    raise ValueError("Invalid/duplicate SQ header")
                header[tags["SN"]] = int(tags["LN"])
            elif line.startswith("@"):
                continue
            elif line.strip():
                first = line
                break
        if header != {c: m["length"] for c, m in metadata.items()}:
            raise ValueError("SAM SQ IDs/lengths differ from candidate metadata")

        def records():
            if first is not None:
                yield _record(first, metadata)
            for line in handle:
                if line.startswith("@"):
                    raise ValueError("SAM header after records")
                if line.strip():
                    yield _record(line, metadata)

        for qname, group in groupby(records(), key=lambda r: r["qname"]):
            if qname in seen:
                raise ValueError("QNAME repeated in noncontiguous groups")
            seen.add(qname)
            rs = list(group)
            original_audit["qname_groups"] += 1
            original_audit["alignment_records"] += len(rs)
            original_audit["secondary_records"] += sum(bool(r["flag"] & 256) for r in rs)
            original_audit["mapq_zero_records_retained"] += sum(r["mapq"] == 0 and not r["flag"] & 4 for r in rs)
            pairs = _eligible_pairs(rs, metadata, cfg, original_audit)
            mark("original_eligible", pairs)
            evidence, internal, boundary, selected = {}, False, False, False
            best = None
            qualifying = []
            if not pairs:
                feature = UNEXPLAINED
                original_audit["unexplained_fragments"] += 1
            else:
                original_audit["accepted_fragments"] += 1
                best = max(p["score"] for p in pairs)
                near = [p for p in pairs if p["score"] >= best-cfg["near_best_delta"]]
                mark("original_near", near)
                for candidate in sorted({p["candidate"] for p in near}):
                    subset = [p for p in near if p["candidate"] == candidate]
                    evidence[candidate] = dict(score=max(p["score"] for p in subset),
                        direct_junctions=sorted(set.intersection(*(set(p["direct_junctions"]) for p in subset))),
                        bridged_junctions=sorted(set.intersection(*(set(p["bridged_junctions"]) for p in subset))),
                        placements=len(subset))
                feature = json.dumps([[c, e["direct_junctions"], e["bridged_junctions"]]
                    for c, e in evidence.items()], separators=(",", ":"), ensure_ascii=True)
                for pair in near:
                    has_internal, has_boundary = _qualifies(pair, metadata, cfg["min_anchor"])
                    internal |= has_internal
                    boundary |= has_boundary
                    if has_internal or has_boundary:
                        qualifying.append(pair["candidate"])
                selected = internal or boundary
                if selected:
                    counts[feature] += 1
                    reasons["internal_and_boundary" if internal and boundary else "internal_only" if internal else "boundary_only"] += 1
                    mark("selected_eligible", pairs)
                    mark("selected_near", near)
                    representative = min((p for p in near if p["score"] == best), key=lambda p:
                        (p["candidate"], tuple((r["start"], r["end"], tuple(r["blocks"]), r["query_length"])
                                              for r in p["records"])))
                    query_lengths.update(r["query_length"] for r in representative["records"])
                    template_lengths[abs(representative["records"][0]["tlen"])] += 1
                else:
                    anchor_only[feature] += 1
                    if any(min(hi, metadata[p["candidate"]]["junctions"][1]) >
                           max(lo, metadata[p["candidate"]]["junctions"][0])
                           for p in near for r in p["records"] for lo, hi in r["blocks"]):
                        reasons["excluded_subthreshold_internal_overlap"] += 1
            full_counts[feature] += 1
            keep = cfg["ledger_mode"] == "all" or (cfg["ledger_mode"] == "accepted" and selected)
            if keep:
                if cfg["max_ledger"] is None or len(ledger) < cfg["max_ledger"]:
                    ledger.append(dict(qname=qname, feature=feature, selected=selected,
                        internal_evidence=internal, boundary_evidence=boundary,
                        qualifying_candidates=sorted(set(qualifying)), candidates=list(evidence),
                        candidate_evidence=evidence, best_pair_score=best))
                else:
                    omitted += 1

    accepted_full = Counter({f: n for f, n in full_counts.items() if f != UNEXPLAINED})
    if counts+anchor_only != accepted_full:
        raise RuntimeError("Conditioned and anchor-only counts do not conserve accepted QNAMEs")
    if sum(full_counts.values()) != original_audit["qname_groups"]:
        raise RuntimeError("Full counts do not conserve observed QNAME groups")
    if expected_full_counts is not None and dict(full_counts) != expected_full_counts:
        raise ValueError("Reconstructed full_counts differ from historical extraction")
    for key in ("qname_groups", "alignment_records", "accepted_fragments", "unexplained_fragments"):
        original_audit.setdefault(key, 0)
    original_coverage, selected_coverage = {}, {}
    for candidate, item in sorted(metadata.items()):
        left, right = item["junctions"]
        for prefix, output in (("original", original_coverage), ("selected", selected_coverage)):
            eligible = arrays[prefix+"_eligible"][candidate]
            near = arrays[prefix+"_near"][candidate]
            row = dict(length=item["length"], eligible_fragments=fragments[prefix+"_eligible"][candidate],
                near_best_fragments=fragments[prefix+"_near"][candidate],
                eligible_base_breadth=sum(eligible)/item["length"], near_best_base_breadth=sum(near)/item["length"])
            if prefix == "selected":
                row.update(inner_length=right-left, eligible_inner_breadth=sum(eligible[left:right])/(right-left),
                           near_best_inner_breadth=sum(near[left:right])/(right-left))
            output[candidate] = row
    technical = dict(selection="one_deterministic_max_score_pair_placement_per_selected_QNAME; ties_by_candidate_and_coordinates",
        scope="conditioned_local_pairs_only; not_whole_library_or_independent_molecules",
        mate_query_lengths=_distribution(query_lengths), absolute_template_lengths=_distribution(template_lengths))
    audit = dict(original_audit, accepted_fragments=sum(counts.values()), unexplained_fragments=0,
                 original_unexplained_fragments=full_counts[UNEXPLAINED], anchor_only_fragments=sum(anchor_only.values()),
                 ledger_rows_omitted_by_cap=omitted, technical_summary=technical)
    return dict(schema=SCHEMA, protocol=PROTOCOL, config=cfg,
        count_unit="one_QNAME_sequencing_fragment", scope="conditional_on_any_near_best_internal_or_boundary_evidence",
        counts=dict(sorted(counts.items())), anchor_only_counts=dict(sorted(anchor_only.items())),
        full_counts=dict(sorted(full_counts.items())), unexplained_fragments=full_counts[UNEXPLAINED],
        audit=audit, original_audit=dict(original_audit), selection_audit=dict(original_accepted_fragments=sum(accepted_full.values()),
            selected_fragments=sum(counts.values()), anchor_only_fragments=sum(anchor_only.values()),
            unexplained_fragments=full_counts[UNEXPLAINED], qualification_reasons=dict(reasons),
            full_counts_match_expected=None if expected_full_counts is None else True),
        candidate_coverage=selected_coverage, original_candidate_coverage=original_coverage,
        inner_breadth_definition="union_of_selected_eligible_or_near_best_aligned_reference_blocks_clipped_to_[left,right)_divided_by_inner_length",
        anchor_only_definition="no_near_best_placement_meets_internal_block_or_boundary_criterion; subthreshold_internal_overlap_reported_separately",
        fragment_ledger=ledger)
