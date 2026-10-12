"""Source-held-out marker-shape diagnostics; no identity call or calibrated P value.

Counts are per-marker supporting fragments. A fragment can occur in more than
one marker, so their sum is a composite count, not independent sample size.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
import hashlib
import json
import math
import re

SCHEMA = "glycan-evidence.profile.descriptive.v1"
COUNT_UNIT = "per_marker_supporting_fragments; cross_marker_duplicates_possible"
COVERAGE_UNIT = "reference_bases; aligned_bases=paired_fragment_sequenced_base_union"
SCOPE = "descriptive_target_profile_shape; not_identity_or_absence_or_abundance"
ROW_FIELDS = {"feature_id", "length", "fragments", "covered_bases", "aligned_bases"}


def _name(value, field):
    if not isinstance(value, str) or not value.strip() or any(c in value for c in "\r\n\t"):
        raise ValueError(f"Invalid {field}")
    return value


def _integer(value, field, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, str)):
        raise ValueError(f"{field} must be an integer")
    if isinstance(value, str) and not re.fullmatch(r"[0-9]+", value):
        raise ValueError(f"{field} must be an integer")
    value = int(value)
    if value < minimum:
        raise ValueError(f"{field} must be >= {minimum}")
    return value


def _sha(value):
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError("Expected lowercase SHA256")
    return value


@dataclass(frozen=True)
class MarkerPanel:
    """Explicit complete marker universe and mapping context, including margin."""

    target: str
    markers: tuple[tuple[str, int], ...]
    mapping_contract_sha256: str
    margin: int = 0

    def __post_init__(self):
        _name(self.target, "target")
        _sha(self.mapping_contract_sha256)
        object.__setattr__(self, "margin", _integer(self.margin, "margin"))
        markers = tuple((_name(k, "feature_id"), _integer(n, "length", 1))
                        for k, n in self.markers)
        if not markers or len({k for k, _ in markers}) != len(markers):
            raise ValueError("Empty or duplicate marker universe")
        object.__setattr__(self, "markers", tuple(sorted(markers)))

    @property
    def sha256(self):
        body = json.dumps(asdict(self), sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(body.encode()).hexdigest()


@dataclass(frozen=True)
class Source:
    """Caller-supplied curated source/cluster identity; never inferred from reads.

    Aliases and near clones must share cluster_id. This module enforces that
    declaration; it does not itself discover aliases or calculate ANI.
    """

    source_id: str
    cluster_id: str
    positive_targets: tuple[str, ...]
    exposure: str = "development_exposed"

    def __post_init__(self):
        _name(self.source_id, "source_id")
        _name(self.cluster_id, "cluster_id")
        targets = tuple(self.positive_targets)
        if isinstance(self.positive_targets, str) or len(set(targets)) != len(targets):
            raise ValueError("positive_targets must be distinct target names")
        for target in targets:
            _name(target, "positive target")
        object.__setattr__(self, "positive_targets", targets)
        if self.exposure not in {"development_exposed", "prospective_unseen", "unlabelled_application"}:
            raise ValueError("Unknown source exposure")
        if self.exposure == "unlabelled_application" and targets:
            raise ValueError("Unlabelled application cannot supply positive truth")


@dataclass(frozen=True)
class Observation:
    panel: MarkerPanel
    source: Source
    sample: str
    counts: tuple[int, ...]
    covered_bases: tuple[int, ...]
    aligned_bases: tuple[int, ...]
    input_sha256: str

    def __post_init__(self):
        _name(self.sample, "sample")
        _sha(self.input_sha256)
        for field in ("counts", "covered_bases", "aligned_bases"):
            values = tuple(_integer(n, field) for n in getattr(self, field))
            if len(values) != len(self.panel.markers):
                raise ValueError(f"Incomplete {field}; explicit zero rows required")
            object.__setattr__(self, field, values)
        for (_, length), n, covered, aligned in zip(
                self.panel.markers, self.counts, self.covered_bases, self.aligned_bases):
            if covered > length or aligned < covered or aligned < n:
                raise ValueError("Invalid covered/aligned base relationship")
            if bool(n) != bool(covered) or bool(n) != bool(aligned):
                raise ValueError("Zero count and base support disagree")
            if aligned > n * length:
                raise ValueError("Aligned bases exceed per-marker fragment union bound")


def observation_from_rows(panel, source, sample, rows, *, input_sha256):
    """Validate a complete compact grid; missing markers are errors, not zeros.

    Each row has exactly ROW_FIELDS. The caller explicitly selects one sample,
    target and margin from a larger table, and passes a frozen panel definition.
    """
    by_id = {}
    expected = dict(panel.markers)
    for row in rows:
        if set(row) != ROW_FIELDS:
            raise ValueError("Unexpected compact marker row schema")
        key = _name(row["feature_id"], "feature_id")
        if key in by_id:
            raise ValueError(f"Duplicate marker: {key}")
        if key not in expected or _integer(row["length"], "length", 1) != expected[key]:
            raise ValueError(f"Unknown marker or changed length: {key}")
        by_id[key] = row
    if set(by_id) != set(expected):
        raise ValueError("Incomplete marker grid; explicit zero rows required")
    values = [by_id[k] for k, _ in panel.markers]
    return Observation(panel, source, sample,
                       tuple(r["fragments"] for r in values),
                       tuple(r["covered_bases"] for r in values),
                       tuple(r["aligned_bases"] for r in values), input_sha256)


@dataclass(frozen=True)
class ExpectedProfile:
    panel: MarkerPanel
    probabilities: tuple[float, ...]
    training_sources: tuple[Source, ...]
    training_samples: tuple[str, ...]
    training_input_sha256: tuple[str, ...]
    source_probabilities: tuple[tuple[float, ...], ...]

    def __post_init__(self):
        for field in ("probabilities", "training_sources", "training_samples", "training_input_sha256"):
            object.__setattr__(self, field, tuple(getattr(self, field)))
        object.__setattr__(self, "source_probabilities", tuple(tuple(p) for p in self.source_probabilities))
        if not self.training_sources or not self.training_samples:
            raise ValueError("Empty training provenance")
        if len(self.training_samples) != len(set(self.training_samples)):
            raise ValueError("Duplicate training samples")
        if len(self.training_samples) != len(self.training_input_sha256):
            raise ValueError("Incomplete training hashes")
        if len(set(self.training_input_sha256)) != len(self.training_input_sha256):
            raise ValueError("Duplicate training input hashes")
        for value in self.training_input_sha256:
            _sha(value)
        if len(self.source_probabilities) != len(self.training_sources):
            raise ValueError("Incomplete source profiles")
        for values in (self.probabilities, *self.source_probabilities):
            if (len(values) != len(self.panel.markers) or
                    any(not math.isfinite(p) or p < 0 or p > 1 for p in values) or
                    not math.isclose(sum(values), 1, abs_tol=1e-10)):
                raise ValueError("Invalid expected probabilities")
        mean = [sum(p[i] for p in self.source_probabilities) / len(self.source_probabilities)
                for i in range(len(self.panel.markers))]
        if any(not math.isclose(p, q, rel_tol=1e-12, abs_tol=1e-12)
               for p, q in zip(self.probabilities, mean)):
            raise ValueError("Expected probabilities differ from equal-source mean")
        if len({s.source_id for s in self.training_sources}) != len(self.training_sources):
            raise ValueError("Duplicate training sources")
        if len({s.cluster_id for s in self.training_sources}) != len(self.training_sources):
            raise ValueError("Training sources share a cluster; choose one representative")
        for source in self.training_sources:
            if source.exposure != "development_exposed" or self.panel.target not in source.positive_targets:
                raise ValueError("Training requires exposed, curated target-positive sources")

    def to_dict(self):
        return {"schema": SCHEMA, **asdict(self), "panel_sha256": self.panel.sha256,
                "count_unit": COUNT_UNIT, "scope": SCOPE,
                "estimator": "equal_source_mean_of_within_source_pooled_count_fractions",
                "smoothing": "none; unsupported reference probability remains zero"}


def fit_expected(observations):
    """Fit empirical opportunity fractions using target-positive development sources.

    Technical observations pool within a source; sources then receive equal
    weight. Mapping/length/copy effects are absorbed empirically, not separately
    identified. Source replication is not supplied by technical observations.
    """
    observations = tuple(observations)
    if not observations:
        raise ValueError("No target-positive development observations")
    panel = observations[0].panel
    if len({o.sample for o in observations}) != len(observations):
        raise ValueError("Duplicate training sample")
    groups = defaultdict(list)
    for observation in observations:
        if observation.panel != panel:
            raise ValueError("Training panel/mapping/margin mismatch")
        groups[observation.source.source_id].append(observation)
    sources, probabilities = [], []
    for source_id in sorted(groups):
        group = groups[source_id]
        source = group[0].source
        if any(o.source != source for o in group):
            raise ValueError("Conflicting source provenance")
        counts = [sum(o.counts[i] for o in group) for i in range(len(panel.markers))]
        total = sum(counts)
        if total == 0:
            raise ValueError("Training source has no supporting counts")
        sources.append(source)
        probabilities.append(tuple(n / total for n in counts))
    mean = tuple(sum(p[i] for p in probabilities) / len(probabilities)
                 for i in range(len(panel.markers)))
    return ExpectedProfile(panel, mean, tuple(sources),
                           tuple(o.sample for o in observations),
                           tuple(o.input_sha256 for o in observations), tuple(probabilities))


def diagnose(observation, expected):
    """Describe shape with source/cluster overlap forbidden, even in development.

    IID occupancy expectations are mathematical reference values only. Existing
    feature counts are dependent; no null calibration, P value, power, detection
    threshold, identity call or absence inference is produced.
    """
    if observation.panel != expected.panel:
        raise ValueError("Observation and expected panel/mapping/margin mismatch")
    if (observation.source.source_id in {s.source_id for s in expected.training_sources}
            or observation.source.cluster_id in {s.cluster_id for s in expected.training_sources}
            or observation.sample in expected.training_samples
            or observation.input_sha256 in expected.training_input_sha256):
        raise ValueError("Training/observation source, cluster, sample or input hash overlap")
    counts, p = observation.counts, expected.probabilities
    total = sum(counts)
    frequencies = [n / total for n in counts] if total else [0.0] * len(counts)
    missing_mass = sum(q for n, q in zip(counts, p) if n == 0)
    outside = sum(n for n, q in zip(counts, p) if q == 0)
    # Conditional multinomial idealization, explicitly not a sampling model for
    # correlated per-feature counts. log1p protects small opportunity fractions.
    zero_probabilities = [1.0 if total == 0 else
                          (0.0 if q == 1 else math.exp(total * math.log1p(-q))) for q in p]
    iid_missing_mass = sum(q * z for q, z in zip(p, zero_probabilities))
    js = None
    deviance = None
    if total:
        js = sum(0.5 * value * math.log(2 * value / (r + q))
                 for r, q in zip(frequencies, p) for value in (r, q) if value > 0)
        if not outside:
            deviance = max(0.0, 2 * sum(n * math.log(n / (total * q))
                                       for n, q in zip(counts, p) if n))
    length = sum(n for _, n in observation.panel.markers)
    return {
        "schema": SCHEMA, "sample": observation.sample,
        "source_id": observation.source.source_id, "cluster_id": observation.source.cluster_id,
        "source_exposure": observation.source.exposure, "target": observation.panel.target,
        "panel_sha256": observation.panel.sha256, "input_sha256": observation.input_sha256,
        "training_sources": [s.source_id for s in expected.training_sources],
        "training_clusters": [s.cluster_id for s in expected.training_sources],
        "training_source_count": len(expected.training_sources), "source_overlap": False,
        "status": "descriptive_only" if total else "no_support",
        "calibration": "none; iid_occupancy_is_only_a_mathematical_reference",
        "p_value": None, "identity_call": None, "scope": SCOPE,
        "count_unit": COUNT_UNIT, "coverage_unit": COVERAGE_UNIT,
        "composite_count_total": total, "marker_count": len(counts),
        "zero_markers": counts.count(0),
        "weighted_breadth": sum(observation.covered_bases) / length,
        "weighted_aligned_base_depth": sum(observation.aligned_bases) / length,
        "js_divergence_nats": js,
        "total_variation": sum(abs(r - q) for r, q in zip(frequencies, p)) / 2 if total else None,
        "composite_deviance": deviance,
        "deviance_status": "no_support" if not total else
                           ("reference_support_gap" if outside else "descriptive_composite"),
        "reference_zero_probability_count": outside,
        "reference_zero_probability_fraction": outside / total if total else None,
        "zero_opportunity_mass": missing_mass,
        "zero_opportunity_composite_count": total * missing_mass,
        "iid_expected_zero_markers": sum(zero_probabilities),
        "iid_expected_zero_opportunity_mass": iid_missing_mass,
        "zero_opportunity_excess_over_iid": missing_mass - iid_missing_mass,
        "observed_effective_marker_count": 1 / sum(r * r for r in frequencies) if total else None,
        "expected_effective_marker_count": 1 / sum(q * q for q in p),
        "markers": [{"feature_id": key, "length": size, "count": n,
                     "covered_bases": covered, "aligned_bases": aligned,
                     "expected_fraction": q, "expected_composite_count": total * q,
                     "iid_zero_probability": z}
                    for (key, size), n, covered, aligned, q, z in zip(
                        observation.panel.markers, counts, observation.covered_bases,
                        observation.aligned_bases, p, zero_probabilities)],
    }
