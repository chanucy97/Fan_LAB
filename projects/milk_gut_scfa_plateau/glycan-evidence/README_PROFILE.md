# Marker-profile diagnostics (development version)

`glycan_evidence.profile` describes whether all marker counts, including zeros,
resemble an empirical target profile from other source genomes. It reports
distribution differences and sampling-opportunity summaries. It does not change
the frozen identity engine or return a new identity/absence call, P value,
abundance estimate or validated detection limit.

## Python API

The module has no third-party dependencies. Its public entry points are:

```python
from glycan_evidence.profile import (
    MarkerPanel, Source, observation_from_rows, fit_expected, diagnose,
)

panel = MarkerPanel(
    target="longum",
    markers=(("marker_a", 1000), ("marker_b", 500)),
    mapping_contract_sha256=verified_mapping_contract_sha256,
    margin=0,
)
training_source = Source("accession_A", "curated_cluster_A", ("longum",))
query_source = Source("accession_B", "curated_cluster_B", ("longum",))

# Both row lists must contain both markers, with explicit zero-valued rows.
train = observation_from_rows(panel, training_source, "training_sample",
                              training_rows, input_sha256=training_table_sha256)
query = observation_from_rows(panel, query_source, "query_sample",
                              query_rows, input_sha256=query_table_sha256)
model = fit_expected([train])
diagnostic = diagnose(query, model)
model_record = model.to_dict()
```

The example uses caller-supplied, verified SHA values and real row lists; it does
not invent sequencing data. `ExpectedProfile` is also a public frozen dataclass
for validated reconstruction of its recorded fields. `to_dict()` adds schema,
panel hash, estimator, unit and scope metadata to those fields. I/O adapters must
verify the declared schema/hash and reconstruct `MarkerPanel` and `Source`
objects before constructing `ExpectedProfile`; passing arbitrary JSON directly
to `diagnose` is unsupported.

### Complete observation schema

Every compact row has exactly these five fields:

| Field | Definition |
|---|---|
| `feature_id` | Identifier in the fixed complete marker universe |
| `length` | Positive reference-marker length |
| `fragments` | Nonnegative integer number of supporting fragments for this marker |
| `covered_bases` | Reference bases covered by at least one supporting fragment |
| `aligned_bases` | Sum of sequenced-base unions of paired fragments within this marker |

Unsequenced insert gaps do not count as covered sequence. A fragment can support
multiple markers: summing `fragments` across markers gives a **composite count**,
not the number of independent molecules. Base coverage, depth and fragment
counts are related views of the same evidence, not independent observations.

`observation_from_rows` rejects missing/duplicate/unknown markers, changed
lengths, invalid integer counts and impossible count/base relationships. Marker
order is canonicalized. Missing rows never become observed zeros. An upstream
adapter must select one sample, target and margin, verify the complete original
table and mapping provenance, then pass the compact rows and actual table hash.
The numeric API itself cannot authenticate raw files from a caller-supplied hash.

### Source and training contracts

`Source(source_id, cluster_id, positive_targets, exposure)` requires explicit
curation. Aliases and near clones must share the same `cluster_id`; the module
does not discover aliases or compute ANI. Allowed exposures are:

- `development_exposed`: already inspected; can train only for a curated positive target.
- `prospective_unseen`: cannot train. This flag alone does not certify an independent experiment.
- `unlabelled_application`: no positive truth labels; can receive descriptive output.

Training and query must differ in source ID, cluster ID, sample ID and input
table hash. This exclusion applies across all read depths and error rates. A
renamed copy of a training table is rejected. An upstream source audit remains
necessary because different reads from the same unrecognized biological source
can have different hashes. Training cannot count two sources from one cluster
as independent representatives.

`fit_expected` pools technical observations within each source, normalizes that
source's full count vector, and averages the resulting probability vectors with
equal source weight. More reads from one source do not give it more weight than
another source. Length, copy number, reference differences and mapping retention
are absorbed in this empirical profile; they are not separately estimated.

The caller must record the training read design and select compatible conditions.
`MarkerPanel` binds the target, marker lengths, mapping contract and margin. It
does not contain an error model or certify matching sequencing platforms.

## Quantities and interpretation

Let `n_i` be each marker's count, `N=sum(n_i)` its composite total, `r_i=n_i/N`,
and `p_i` the fitted target fractions from other sources.

- **Jensen–Shannon divergence**: `0.5*KL(r||m)+0.5*KL(p||m)`, `m=(r+p)/2`, in
  natural-log units. Range 0 to log(2); it is finite even when a reference has
  zero opportunity. Low counts alone can produce a large divergence.
- **Composite deviance**: `2*sum(n_i*log(n_i/(N*p_i)))`. It is a descriptive
  count-weighted discrepancy, not a chi-square statistic with validated degrees
  of freedom. If any positive count has `p_i=0`, return `None` and
  `deviance_status=reference_support_gap` instead of smoothing it away.
- **Zero-opportunity mass**: `Z=sum(p_i for n_i==0)`. This is the fraction of the
  fitted target profile represented by markers with no observed support.
- **IID occupancy reference**: `Z_iid=sum(p_i*(1-p_i)**N)` and
  `sum((1-p_i)**N)` for the expected number of zero markers. These are
  mathematical reference values under independent multinomial draws. The
  observed counts violate that independence and are drawn from variable
  genomes, so these are not calibrated expectations for real sequencing.
- **Zero-opportunity excess**: `Z-Z_iid`, a signed descriptive difference. It
  helps display the contribution of limited sampling, but has no validated
  cutoff, uncertainty interval or P value in this release.

Zero training probabilities remain zero; no pseudocount is selected. The model
retains each source's probability vector and checks that the stored mean equals
their equal-weight average. A reference support gap is evidence that the small
training library does not cover every observed pattern, not proof of a wrong
taxonomy. A single-target profile can also be distorted by true target/near-relative
coexistence; this implementation does not resolve that mixture.

`status` is `no_support` when `N=0`, otherwise `descriptive_only`. `p_value` and
`identity_call` are always `None`. For no support, shape divergence and deviance
are undefined. Low information is displayed through the count/occupancy outputs;
no arbitrary low-depth decision threshold is introduced.

## Executed development analysis

Run the separate workspace analysis from the project root:

```text
python outputs/20261008_tool_development/profile_development/run_profiles.py --out NEW_DIRECTORY --include-cohort
```

The script requires a new output directory, verifies archived input hashes and
frozen mapping provenance, writes the analysis plan before fitting, and leaves
a failed receipt on errors. Completion means that descriptive analysis finished.

The final recorded run is
`outputs/20261008_tool_development/profile_development/run_v1_final/`:

- 90 conditions from 5 **already exposed development genomes**; 180 target rows.
- Each longum-positive source is evaluated against the other positive source,
  using its 30× observations in the same substitution-error stratum. No source
  is trained on another depth of itself.
- Negative sources use the two longum-positive source profiles equally, with
  separate one-training-source sensitivity tables.
- Infantis has only one positive source, so its 18 positive conditions explicitly
  lack a source-held-out model. The other 72 infantis conditions can be described
  against that one source; they do not validate positive sensitivity.
- All 50 retained cohort samples are described for both targets under both
  training error strata (200 rows). Their actual error rates, organism mixture
  and target truth are unknown; the two strata are sensitivity views.

At substitution error 0.005 and 30×, M257B_A6's longum zero-opportunity excess is 55.68 percentage points;
the two true longum sources give 3.46 and 1.46 points in reciprocal source
holdout. At 0.1×, true longum sources give 0.49 and −8.45 points despite sparse
coverage. Yet low-depth iVE-15 JS divergence (0.329) is as large as M257's
high-depth value (0.322). These observations motivate conditioning a future
diagnostic on sampling opportunity; they do not establish a classifier cutoff.
At 30× the two positive-source JS values are equal by symmetry: they compare
the same pair in opposite directions, not two independent replications.

The method report in `profile_development/PROFILE_DEVELOPMENT_REPORT.md` records
the cohort observations, reference gaps and next validation steps. All five
genomes have influenced development and cannot be reused as unseen tests.

## Tests and scope of validation

```text
python -m unittest discover -s software/glycan-evidence/tests -p test_profile.py -v
```

With the package installed (or `software/glycan-evidence/src` on `PYTHONPATH`),
21 tests cover complete grids, zero-versus-missing evidence, count/base
consistency, source/cluster/hash isolation, equal-source weighting, immutable
models, reference-zero support, zero-input behavior and low/high-count occupancy.

Next work requires more target and near-relative sources, independent-fragment
or block-level null simulation, library/GC/reference-deletion challenges,
explicit mixture models and a newly frozen evaluation. Depth–breadth and
probabilistic read-assignment concepts have prior art; this implementation is a
transparent development diagnostic, not an originality or maturity claim.
