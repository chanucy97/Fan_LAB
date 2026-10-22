# Two-component marker-support model (research module)

`glycan_evidence.mixture` fits a fixed empirical target marker profile together
with a fixed, explicitly identified near-relative background profile. The module
is the dev3 development Python API; it is not part of the frozen dev2 artifact
and does not alter frozen identity calls.

## API

```python
from glycan_evidence.mixture import (
    QueryMembership, component_from_observations, fit_mixture,
)

background = component_from_observations(background_observations, role="background")
target = component_from_observations(other_source_target_observations, role="target")
membership = QueryMembership(
    sources=all_actual_read_sources,
    reused_training_samples=audited_reused_training_sample_ids,
    audit_sha256=source_membership_audit_sha256,
)
result = fit_mixture(observation, background, target, membership)
```

Observations, sources and panels use the existing `profile` API. All marker rows,
including observed zeros, must be present. The mapping contract, marker lengths,
target and margin must match. Component training pools counts within a source,
then averages normalized profiles equally across sources. No per-marker free
background or pseudocount is introduced. The background retains its curated
non-target label; it is never relabelled target-positive to fit the target API.

`QueryMembership` must contain every actual source contributing reads to the
query, rather than only a synthetic composition identifier. Query source IDs,
near-clone clusters, sample IDs and input hashes cannot overlap target training.
Known background-source overlap is allowed and explicitly reported. Upstream
read hashes must establish training-read reuse; a changed remapped table hash
does not establish independent reads. The numeric API validates declarations,
but cannot discover undeclared source identities or authenticate raw reads.

`ComponentProfile.to_dict()` retains source-level profiles, training source and
cluster IDs, samples, table hashes, panel, estimator, units and scope. Like the
profile API, it is a validated in-memory object; no unauthenticated JSON loader
or new CLI command is supplied here.

## Objective and solver

For complete marker counts `n_i`, composite sum `N`, fixed background fractions
`p_Bi`, and fixed target fractions `p_Ti`, fit:

```text
q_i(w) = (1-w) p_Bi + w p_Ti,  0 <= w <= 1
ell(w) = sum_i n_i log(q_i(w))
ell'(w) = sum_i n_i (p_Ti-p_Bi) / q_i(w)
ell''(w) = -sum_i n_i (p_Ti-p_Bi)^2 / q_i(w)^2 <= 0
```

The objective is concave. Endpoint derivative checks select boundary maxima;
otherwise monotone derivative bisection finds the global maximum to an absolute
weight tolerance of `1e-12`. The omitted conditional multinomial constant is
identical under H0 and H1. Counts may share fragments across markers: this is a
**descriptive conditional composite objective**, not an independent-fragment
sampling likelihood or calibrated likelihood-ratio test.

`w` is the fitted target fraction of **composite marker support**. `N*w` is its
fitted support component. Neither is a cell/genome relative abundance, a read
fraction, a posterior probability or a new identity call. Relative mapping
retention, marker length/copy differences, strain differences and finite training
sampling affect the two empirical profiles.

H0 fixes `w=0`; H1 permits `[0,1]`. Finite `H1-H0` and twice that value are
descriptive improvement measures without P values, thresholds or confidence
intervals. If H0 puts zero probability on an observed marker, its log objective
is negative infinity: the output uses `None` plus an explicit impossibility
status, never a finite invented score. A training zero may result from finite
sampling; it is not proof of target identity.

## Retained outcomes

| Status | Meaning |
|---|---|
| `fitted_descriptive` | Unique numerical support-weight optimum within the fixed two-profile family |
| `no_support` | Entire observed vector is zero; weight and shape fit undefined |
| `reference_support_gap` | At least one observed count has zero probability in both components; the full objective is impossible for every weight |
| `unidentifiable_identical_profiles` | Components are equal; any weight gives the same fitted spectrum |
| `unidentifiable_near_collinear_profiles` | Component total variation is within the recorded numerical guard, default `1e-8`; no unstable weight reported |
| `unidentifiable_observed_support` | Profiles differ only in markers with no observed counts; the objective is flat in this observation |

For exact nonidentification, `[0,1]` records the allowed weight range. For the
near-collinear guard it records refusal to narrow that range numerically, not a
mathematical confidence interval. The guard is fixed before this analysis and
must not be tuned as a biological classification threshold.

A union-reference gap stops the full-spectrum fit. Unsupported markers are not
dropped, and no supported-subset weight is passed off as a complete fit. Every
marker is returned with observed count, component probabilities and gap flag.
Identifiable fits additionally return fitted probabilities, expected counts,
signed residuals, component expected counts, Jensen–Shannon divergence in nats
and descriptive composite deviance. Zero-count markers remain in these arrays.

## Executed development challenge

The separate reproducible driver is:

```text
python outputs/20261008_tool_development/mixture_model/run_mixture_models.py NEW_OUTPUT_DIRECTORY
```

It authenticates the prior 35-condition collection and complete metrics, the
frozen component/read memberships and the original training tables. The target
profile uses the **other** longum genome at 30× / substitution error 0.005 for
each positive source. Pure M257 controls use the equal-source mean of the two
target profiles. The background uses M257B_A6 at 30× / error 0.005. This source
was already seen, and its 30× reads are reused in nine conditions. These are
known-background development challenges, not independent validation or evidence
that an external false-positive rate has been corrected.

The final recorded run, `mixture_model/run_v1_final`, contains 15 fitted conditions and
20 reference-gap conditions. It retains all 35 conditions and all marker rows.
The research report records fit improvements, reference failures and discrepancies
between fitted weights and known-origin marker-support fractions. The latter are
simulation provenance quantities, not organism abundance estimates.

## Verification

```text
python -m unittest discover -s software/glycan-evidence/tests -p test_mixture.py -v
```

Fifteen mathematical/provenance tests cover a closed-form binomial solution,
boundary optima, global agreement with a dense objective grid, count scaling,
zero inputs, union support gaps, exact/near/sample-specific nonidentification,
source-balanced profiles, mapping contracts, target-source/cluster exclusion and
explicit background read reuse. These fixtures validate implementation behavior;
they are not biological benchmark replicates.
