"""Truth-free reconciliation of reference calls and same-library assemblies.

The primary evidence rule is the existing structure.summarize_evidence rule.
Assembly-only and combined reports retain identical assembly candidate sets;
reference calls never create complete sequences or increase assembly support.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import math
from pathlib import Path
import re

SCHEMA = "glycan-evidence.structure-recovery/1"
DEFAULT_CONFIG = {
    "primary": {"full_interval_breadth_ge_1": 1.0, "min_boundary_qnames": 1},
    "strict": {"min_interval_read_depth": 3, "min_boundary_qnames": 3,
               "min_boundary_distinct_pair_starts": 2},
    "scope": "development_evidence_tiers_not_calibrated_classification",
}


def _require(ok, message):
    if not ok:
        raise ValueError(message)


def _integer(value, label, minimum=0):
    _require(type(value) is int and value >= minimum, "Invalid integer: " + label)
    return value


def _number(value, label, low=0, high=math.inf):
    _require(type(value) in (int, float) and math.isfinite(value) and low <= value <= high,
             "Invalid finite number: " + label)
    return value


def _hash(value, label):
    _require(isinstance(value, str) and re.fullmatch(r"[a-f0-9]{64}", value),
             "Invalid SHA256: " + label)
    return value


def validate_file_binding(record):
    """Verify a small local source file against its retained byte/SHA binding."""
    _require(isinstance(record, dict) and isinstance(record.get("path"), str),
             "Missing file binding")
    _integer(record.get("bytes"), "file bytes")
    _hash(record.get("sha256"), "file SHA")
    path = Path(record["path"])
    _require(path.is_file(), "Source file missing: " + str(path))
    raw = path.read_bytes()
    _require(len(raw) == record["bytes"] and hashlib.sha256(raw).hexdigest() == record["sha256"],
             "Source file binding mismatch: " + str(path))


def _length_class(length, by_length):
    if 4897 <= length <= 4899 and "alternative_4897_4899" in by_length.values():
        return "alternative_4897_4899"
    return by_length.get(length, "inner_length_" + str(length))


def recover_structure(case_id, reference_calls, panel, candidates, config=None):
    """Reconcile already-bound evidence; caller verifies source file bindings.

    Candidates originate from this case's ordinary assembly. Positions are
    physical 0-based half-open coordinates; sequence hashes are reference-
    oriented. Unknown sequence hashes and missing base audits remain unknown.
    """
    _require(isinstance(case_id, str) and case_id, "Missing case ID")
    cfg = deepcopy(DEFAULT_CONFIG if config is None else config)
    _require(cfg == DEFAULT_CONFIG, "This version accepts only its frozen evidence tiers")
    _require(isinstance(panel, list) and panel, "Empty reference panel")
    by_length, full_hashes, inner_hashes, classes, panel_ids = {}, set(), set(), set(), set()
    for item in panel:
        pid = item.get("candidate_id")
        _require(isinstance(pid, str) and pid and pid not in panel_ids, "Duplicate/invalid panel ID")
        panel_ids.add(pid)
        length = _integer(item.get("inner_length"), "panel inner length", 1)
        group = item.get("class_id")
        _require(isinstance(group, str) and group, "Missing panel class")
        _require(length not in by_length or by_length[length] == group, "Ambiguous panel length class")
        by_length[length] = group
        classes.add(group)
        full_hashes.add(_hash(item.get("full_sequence_sha256"), "panel full sequence"))
        inner_hashes.add(_hash(item.get("inner_sequence_sha256"), "panel inner sequence"))
    _require(isinstance(reference_calls, list) and all(isinstance(x, str) for x in reference_calls),
             "Invalid reference calls")
    _require(len(reference_calls) == len(set(reference_calls)) and set(reference_calls) <= classes,
             "Duplicate/out-of-panel reference call")
    _require(isinstance(candidates, list), "Candidates must be a list")
    seen, results = set(), []
    for raw in candidates:
        _require(isinstance(raw, dict), "Candidate must be an object")
        row = deepcopy(raw)
        sid = row.get("structure_id")
        _require(isinstance(sid, str) and sid and sid not in seen, "Duplicate/invalid structure ID")
        seen.add(sid)
        _require(row.get("case_id") == case_id, "Cross-library candidate refused")
        _require(row.get("source_kind") == "own_library_ordinary_assembly", "Invalid candidate source")
        bindings = row.get("source_bindings")
        _require(isinstance(bindings, dict) and {"receipt", "structures", "coverage", "boundaries"} <= set(bindings),
                 "Missing candidate source bindings")
        for name, binding in bindings.items():
            _require(isinstance(binding, dict) and isinstance(binding.get("path"), str), "Invalid source binding")
            _integer(binding.get("bytes"), "binding bytes")
            _hash(binding.get("sha256"), name)
        _require(isinstance(row.get("contig"), str) and row["contig"], "Missing contig")
        _hash(row.get("contig_sequence_sha256"), "contig")
        length = _integer(row.get("contig_length"), "contig length", 1)
        start = _integer(row.get("interval_start"), "interval start")
        end = _integer(row.get("interval_end"), "interval end", 1)
        left = _integer(row.get("left_boundary"), "left boundary")
        right = _integer(row.get("right_boundary"), "right boundary")
        _require(0 <= start < left < right < end <= length, "Invalid interval/boundary coordinates")
        _require(row.get("orientation") in ("+", "-"), "Invalid orientation")
        for field in ("anchor_screen_pass", "inner_edges_complete", "outer_anchors_complete"):
            _require(type(row.get(field)) is bool, "Missing Boolean: " + field)
        _require(not row["outer_anchors_complete"] or row["inner_edges_complete"], "Full anchors require complete inner edges")
        inner_length = row.get("inner_length")
        if row["inner_edges_complete"]:
            _require(_integer(inner_length, "inner length", 1) == right-left, "Inner length/geometry mismatch")
        else:
            _require(inner_length is None, "Partial inner edges cannot define exact inner length")
        full_hash = _hash(row.get("interval_sequence_sha256"), "assembled interval")
        inner_hash = row.get("inner_sequence_sha256")
        if row["inner_edges_complete"]:
            _hash(inner_hash, "assembled inner sequence")
        else:
            _require(inner_hash is None, "Partial geometry cannot define inner sequence hash")
        coverage = row.get("coverage")
        _require(isinstance(coverage, dict), "Missing coverage")
        _require(coverage.get("unit") == "primary_mate_alignment_read_depth" and
                 coverage.get("filter") == "quality_filtered" and coverage.get("segment") == "full_anchored_interval",
                 "Coverage scope mismatch")
        _require(coverage.get("start") == start and coverage.get("end") == end, "Coverage coordinate mismatch")
        breadth = _number(coverage.get("breadth_ge_1"), "breadth", 0, 1)
        depth = _integer(coverage.get("min_read_depth"), "min read depth")
        _require((breadth == 1) == (depth >= 1), "Breadth/minimum-depth inconsistency")
        boundaries = row.get("boundaries")
        _require(isinstance(boundaries, dict) and set(boundaries) == {"left", "right"}, "Missing/extra boundary")
        for side, position in (("left", left), ("right", right)):
            bound = boundaries[side]
            _require(bound.get("position") == position and bound.get("evidence") == "direct_read_spanning",
                     "Boundary position/evidence mismatch")
            n = _integer(bound.get("qnames"), "boundary QNAME count")
            starts = _integer(bound.get("distinct_pair_starts"), "boundary pair starts")
            _require(starts <= n, "Pair starts exceed QNAMEs")
        failures = []
        if not row["anchor_screen_pass"]: failures.append("anchor_screen_failed")
        if not row["inner_edges_complete"]: failures.append("partial_inner_anchor_edges")
        if breadth != 1: failures.append("incomplete_quality_filtered_coverage")
        for side in ("left", "right"):
            if boundaries[side]["qnames"] < 1: failures.append(side + "_boundary_unobserved")
        supported = not failures
        strict_failures = list(failures)
        if depth < cfg["strict"]["min_interval_read_depth"]: strict_failures.append("low_interval_depth")
        for side in ("left", "right"):
            if boundaries[side]["qnames"] < 3: strict_failures.append(side + "_low_boundary_qnames")
            if boundaries[side]["distinct_pair_starts"] < 2: strict_failures.append(side + "_low_pair_start_diversity")
        group = _length_class(inner_length, by_length) if inner_length is not None else None
        full_known = full_hash in full_hashes if row["outer_anchors_complete"] else None
        inner_known = inner_hash in inner_hashes if inner_hash is not None else None
        row.update(evidence_level="read_supported_assembly_candidate" if supported else "partial_or_unresolved_assembly_candidate",
                   assembly_supported=supported, strict_supported=not strict_failures,
                   primary_failures=failures, strict_failures=strict_failures,
                   length_class=group, outside_panel_length_class=(group not in classes) if group else None,
                   exact_full_sequence_in_panel=full_known, exact_inner_sequence_in_panel=inner_known,
                   new_sequence_candidate=(not inner_known) if supported and inner_known is not None else None,
                   sequence_comparison_scope="fixed_15_reference_panel_not_global_novelty")
        audit = row.get("base_concordance")
        if audit is not None:
            _require(isinstance(audit, dict), "Invalid base concordance")
            for key in ("positions", "zero_acgt_support_positions", "assembly_base_no_qname_support_positions", "alternative_exceeds_assembly_positions"):
                _integer(audit.get(key), key)
            _require(audit["positions"] == end-start, "Base audit interval mismatch")
            for key in audit:
                if key.endswith("_positions"):
                    _integer(audit[key], key)
            _require(all(audit[x] <= audit["positions"] for x in audit if x.endswith("_positions")), "Invalid base audit denominator")
        results.append(row)
    assembly_classes = sorted({x["length_class"] for x in results if x["assembly_supported"]})
    strict_classes = sorted({x["length_class"] for x in results if x["strict_supported"]})
    ref = set(reference_calls)
    observed = set(assembly_classes)
    return {"schema": SCHEMA, "case_id": case_id, "config": cfg,
            "reference_supported_classes": sorted(ref), "assembly_supported_classes": assembly_classes,
            "strict_assembly_supported_classes": strict_classes,
            "assembly_only_supported_classes": assembly_classes.copy(),
            "reference_only_classes": sorted(ref-observed), "assembly_only_classes": sorted(observed-ref),
            "concordant_classes": sorted(ref & observed),
            "discordance": bool(ref != observed),
            "discordance_scope": "method_evidence_disagreement_not_biological_conflict",
            "candidate_count": len(results), "candidates": results,
            "observation": "assembly_candidates_observed" if results else "unresolved_no_paired_anchor_assembly_candidate",
            "scope": "same_library_assembly_readback; reused_reads; no_independent_validation; no_abundance_or_carrier_call"}
