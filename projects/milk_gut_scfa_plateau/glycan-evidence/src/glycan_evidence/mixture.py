"""Two fixed empirical marker profiles; descriptive composite fit, not identity.

Marker support can overlap across markers. Conditioning on its composite sum
does not turn it into independent fragments or a calibrated likelihood test.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
import math

from .profile import COUNT_UNIT, COVERAGE_UNIT, MarkerPanel, Observation, Source, _sha

SCHEMA = "glycan-evidence.mixture.descriptive.v1"
SCOPE = "fixed_background_target_marker_support_mixture; not_abundance_identity_or_probability"


@dataclass(frozen=True)
class ComponentProfile:
    panel: MarkerPanel
    role: str
    probabilities: tuple[float, ...]
    training_sources: tuple[Source, ...]
    training_samples: tuple[str, ...]
    training_input_sha256: tuple[str, ...]
    source_probabilities: tuple[tuple[float, ...], ...]

    def __post_init__(self):
        if self.role not in {"background", "target"}:
            raise ValueError("Unknown component role")
        for name in ("probabilities", "training_sources", "training_samples", "training_input_sha256"):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        object.__setattr__(self, "source_probabilities", tuple(tuple(p) for p in self.source_probabilities))
        if not self.training_sources or not self.training_samples:
            raise ValueError("Missing training provenance")
        if len(self.training_samples) != len(self.training_input_sha256):
            raise ValueError("Incomplete training hashes")
        for values in (self.training_samples, self.training_input_sha256):
            if len(values) != len(set(values)):
                raise ValueError("Duplicate training sample/hash")
        for value in self.training_input_sha256:
            _sha(value)
        for key in ("source_id", "cluster_id"):
            if len({getattr(s, key) for s in self.training_sources}) != len(self.training_sources):
                raise ValueError("Duplicate training source/cluster")
        for source in self.training_sources:
            if source.exposure != "development_exposed":
                raise ValueError("Training sources must be declared development exposed")
            positive = self.panel.target in source.positive_targets
            if positive != (self.role == "target"):
                raise ValueError("Component role disagrees with curated source target label")
        if len(self.source_probabilities) != len(self.training_sources):
            raise ValueError("Missing source probabilities")
        for values in (self.probabilities, *self.source_probabilities):
            if (len(values) != len(self.panel.markers) or
                    any(isinstance(p, bool) or not isinstance(p, (int, float)) or
                        not math.isfinite(p) or p < 0 or p > 1 for p in values) or
                    not math.isclose(math.fsum(values), 1, rel_tol=0, abs_tol=1e-12)):
                raise ValueError("Invalid probability vector")
        mean = [math.fsum(p[i] for p in self.source_probabilities) / len(self.source_probabilities)
                for i in range(len(self.panel.markers))]
        if any(not math.isclose(p, q, rel_tol=1e-12, abs_tol=1e-12)
               for p, q in zip(mean, self.probabilities)):
            raise ValueError("Probabilities differ from equal-source mean")

    def to_dict(self):
        return {"schema": SCHEMA, **asdict(self), "panel_sha256": self.panel.sha256,
                "estimator": "equal_source_mean_of_within_source_pooled_marker_counts",
                "smoothing": "none", "scope": SCOPE, "count_unit": COUNT_UNIT}


def component_from_observations(observations, *, role):
    """Pool technical counts within source, then average source fractions equally."""
    observations = tuple(observations)
    if not observations:
        raise ValueError("No component observations")
    panel = observations[0].panel
    groups = defaultdict(list)
    for observation in observations:
        if not isinstance(observation, Observation) or observation.panel != panel:
            raise ValueError("Component panel/mapping/margin mismatch")
        groups[observation.source.source_id].append(observation)
    sources, probabilities = [], []
    for key in sorted(groups):
        group = groups[key]
        if any(o.source != group[0].source for o in group):
            raise ValueError("Conflicting source declarations")
        counts = [sum(o.counts[i] for o in group) for i in range(len(panel.markers))]
        if not sum(counts):
            raise ValueError("Zero-support component cannot define probabilities")
        sources.append(group[0].source)
        probabilities.append(tuple(n / sum(counts) for n in counts))
    mean = tuple(math.fsum(p[i] for p in probabilities) / len(probabilities)
                 for i in range(len(panel.markers)))
    return ComponentProfile(panel, role, mean, tuple(sources),
                            tuple(o.sample for o in observations),
                            tuple(o.input_sha256 for o in observations), tuple(probabilities))


@dataclass(frozen=True)
class QueryMembership:
    """All read-source identities and upstream audit, not the synthetic mixture ID.

    Read reuse must be audited upstream; unlike source overlap, it cannot be
    established from the remapped table hash. Empty reused_training_samples is
    an explicit upstream assertion, not an inference made by this module.
    """
    sources: tuple[Source, ...]
    reused_training_samples: tuple[str, ...]
    audit_sha256: str

    def __post_init__(self):
        object.__setattr__(self, "sources", tuple(self.sources))
        object.__setattr__(self, "reused_training_samples", tuple(self.reused_training_samples))
        _sha(self.audit_sha256)
        if not self.sources or any(not isinstance(s, Source) for s in self.sources):
            raise ValueError("Explicit all-component source membership required")
        if len({s.source_id for s in self.sources}) != len(self.sources):
            raise ValueError("Duplicate query source")
        if len(self.reused_training_samples) != len(set(self.reused_training_samples)):
            raise ValueError("Duplicate reused training sample")


def _log_likelihood(counts, probabilities):
    if any(n and p == 0 for n, p in zip(counts, probabilities)):
        return None
    return math.fsum(n * math.log(p) for n, p in zip(counts, probabilities) if n)


def _js(frequencies, probabilities):
    return max(0.0, math.fsum(0.5 * value * math.log(2 * value / (r + q))
                             for r, q in zip(frequencies, probabilities)
                             for value in (r, q) if value))


def fit_mixture(observation, background, target, membership, *, near_collinear_tv=1e-8):
    """Maximize sum n_i log((1-w)*p_Bi+w*p_Ti), 0<=w<=1.

    Fixed profiles give a concave objective; its derivative is nonincreasing.
    Boundary derivatives or bisection therefore find its global maximum. No
    pseudocount, marker-specific free background, or independent-read test is
    introduced. near_collinear_tv is a recorded numerical identification guard,
    never a biological decision threshold or calibrated uncertainty interval.
    """
    if not isinstance(membership, QueryMembership):
        raise ValueError("Explicit membership audit required")
    if background.role != "background" or target.role != "target":
        raise ValueError("Reversed component roles")
    if observation.panel != background.panel or observation.panel != target.panel:
        raise ValueError("Panel/mapping/margin mismatch")
    if (isinstance(near_collinear_tv, bool) or not isinstance(near_collinear_tv, (int, float)) or
            not math.isfinite(near_collinear_tv) or not 0 <= near_collinear_tv <= 1e-3):
        raise ValueError("Invalid numerical near-collinearity guard")
    observed_ids = {s.source_id for s in membership.sources}
    observed_clusters = {s.cluster_id for s in membership.sources}
    target_ids = {s.source_id for s in target.training_sources}
    target_clusters = {s.cluster_id for s in target.training_sources}
    background_ids = {s.source_id for s in background.training_sources}
    background_clusters = {s.cluster_id for s in background.training_sources}
    if target_ids & background_ids or target_clusters & background_clusters:
        raise ValueError("Background and target source/cluster overlap")
    if (target_ids & observed_ids or target_clusters & observed_clusters or
            observation.sample in target.training_samples or
            observation.input_sha256 in target.training_input_sha256):
        raise ValueError("Target profile is not source/cluster/sample/hash held out")
    reused = set(membership.reused_training_samples)
    if not reused <= set(background.training_samples) or reused & set(target.training_samples):
        raise ValueError("Undeclared training reuse or target training reads reused")
    directly_matched_samples = {
        sample for sample, digest in zip(background.training_samples, background.training_input_sha256)
        if observation.sample == sample or observation.input_sha256 == digest
    }
    if not directly_matched_samples <= reused:
        raise ValueError("Every directly matched background training sample/hash reuse must be declared")
    if reused and not (background_ids & observed_ids or background_clusters & observed_clusters):
        raise ValueError("Declared read reuse without background source overlap")
    counts, pb, pt = observation.counts, background.probabilities, target.probabilities
    total = sum(counts)
    frequencies = tuple(n / total for n in counts) if total else None
    difference = tuple(t - b for b, t in zip(pb, pt))
    tv = math.fsum(abs(d) for d in difference) / 2
    union_gap = sum(n for n, b, t in zip(counts, pb, pt) if b == t == 0)
    h0 = _log_likelihood(counts, pb) if total else None
    result = {
        "schema": SCHEMA, "sample": observation.sample, "target": observation.panel.target,
        "panel_sha256": observation.panel.sha256, "input_sha256": observation.input_sha256,
        "scope": SCOPE, "count_unit": COUNT_UNIT, "coverage_unit": COVERAGE_UNIT,
        "weight_unit": "fraction_of_composite_marker_support_attributed_to_fixed_target_profile",
        "weight_is_not": "cell_or_genome_abundance; posterior_probability; identity_call",
        "p_value": None, "identity_call": None, "calibration": "none; descriptive_composite_objective",
        "membership": asdict(membership), "target_source_cluster_held_out": True,
        "background_source_seen": bool(background_ids & observed_ids or background_clusters & observed_clusters),
        "background_training_reads_reused": bool(reused),
        "reference_sources": {"background": background.to_dict(), "target": target.to_dict()},
        "composite_count_total": total, "marker_count": len(counts),
        "profile_total_variation": tv, "near_collinear_tv_guard": near_collinear_tv,
        "union_reference_gap_count": union_gap,
        "union_reference_gap_fraction": union_gap / total if total else None,
        "background_reference_gap_count": sum(n for n, p in zip(counts, pb) if p == 0),
        "h0_log_composite_likelihood": h0,
        "h0_status": "no_support" if not total else ("finite" if h0 is not None else "zero_reference_probability"),
        "h1_log_composite_likelihood": None, "h1_minus_h0": None,
        "two_times_log_fit_improvement": None,
        "fit_improvement_status": "undefined",
        "target_marker_support_weight": None, "target_composite_support_component": None,
        "background_composite_support_component": None, "weight_identification_range": None,
        "h1_js_divergence_nats": None, "h1_composite_deviance": None,
        "pure_target_js_divergence_nats": _js(frequencies, pt) if total else None,
        "background_js_divergence_nats": _js(frequencies, pb) if total else None,
        "optimization": {"method": "concave_derivative_bisection", "weight_tolerance": 1e-12,
                         "iterations": 0, "smoothing": "none"},
    }
    weight, fitted = None, None
    if not total:
        status = "no_support"
    elif union_gap:
        status = "reference_support_gap"
    elif tv == 0:
        status = "unidentifiable_identical_profiles"
        fitted = pb
        result["weight_identification_range"] = [0.0, 1.0]
    elif tv <= near_collinear_tv:
        status = "unidentifiable_near_collinear_profiles"
        result["weight_identification_range"] = [0.0, 1.0]
    elif all(n == 0 or d == 0 for n, d in zip(counts, difference)):
        status = "unidentifiable_observed_support"
        result["weight_identification_range"] = [0.0, 1.0]
        result["h1_log_composite_likelihood"] = h0
        result["h1_minus_h0"] = result["two_times_log_fit_improvement"] = 0.0
        result["fit_improvement_status"] = "finite_descriptive"
    else:
        def derivative(w):
            values = []
            for n, b, t, d in zip(counts, pb, pt, difference):
                if not n or not d:
                    continue
                q = (1 - w) * b + w * t
                if q == 0:
                    return math.copysign(math.inf, d)
                values.append(n * d / q)
            return math.fsum(values)

        if derivative(0.0) <= 0:
            weight = 0.0
        elif derivative(1.0) >= 0:
            weight = 1.0
        else:
            lo, hi = 0.0, 1.0
            for iteration in range(80):
                mid = (lo + hi) / 2
                if derivative(mid) > 0:
                    lo = mid
                else:
                    hi = mid
                if hi - lo <= 1e-12:
                    break
            weight = (lo + hi) / 2
            result["optimization"]["iterations"] = iteration + 1
        fitted = tuple((1 - weight) * b + weight * t for b, t in zip(pb, pt))
        status = "fitted_descriptive"
        result["target_marker_support_weight"] = weight
        result["target_composite_support_component"] = total * weight
        result["background_composite_support_component"] = total * (1 - weight)
    if fitted is not None:
        h1 = _log_likelihood(counts, fitted)
        result["h1_log_composite_likelihood"] = h1
        result["h1_js_divergence_nats"] = _js(frequencies, fitted)
        saturated = math.fsum(n * math.log(n / total) for n in counts if n)
        result["h1_composite_deviance"] = max(0.0, 2 * (saturated - h1))
        if h0 is not None:
            improvement = max(0.0, h1 - h0)
            result["h1_minus_h0"] = improvement
            result["two_times_log_fit_improvement"] = 2 * improvement
            result["fit_improvement_status"] = "finite_descriptive"
        else:
            result["fit_improvement_status"] = "h0_impossible_due_to_empirical_zero; no_finite_improvement"
    result["status"] = status
    result["markers"] = [
        {"feature_id": key, "length": length, "count": n,
         "observed_fraction": n / total if total else None,
         "background_probability": b, "target_probability": t,
         "union_reference_gap": b == t == 0,
         "fitted_probability": fitted[i] if fitted is not None else None,
         "expected_composite_count": total * fitted[i] if fitted is not None else None,
         "count_residual": n - total * fitted[i] if fitted is not None else None,
         "target_expected_composite_count": total * weight * t if weight is not None else None,
         "background_expected_composite_count": total * (1 - weight) * b if weight is not None else None}
        for i, ((key, length), n, b, t) in enumerate(zip(observation.panel.markers, counts, pb, pt))]
    return result
