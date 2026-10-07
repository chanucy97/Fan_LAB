# Implementation notes — version 1.0.0

Version 1.0.0 preserves the supplied source snapshot without changing model logic. These notes document behavior relevant to running the implementation and interpreting its outputs. File-level SHA-256 hashes are recorded in [source-manifest.json](../source-manifest.json).

## IPCW convention for tied event and censoring times

In [`evaluate_survival.py`](../code/evaluate_survival.py), the censoring survival estimate at a recorded time `t` is updated as `G_after = G_before × (1 − c(t)/n(t))`, where `c(t)` counts censoring and `n(t)` includes all patients still at risk immediately before that time. Cases use `1/G(T−)`; patients observed beyond 36 months use `1/G(36)`.

When an event and a censoring occur at the same recorded time, this implementation includes both in the censoring risk-set denominator. Results can differ from a convention that removes the tied events before calculating the censoring step. This matters when follow-up times are rounded or grouped. Document the time resolution and tie convention, and check the intended convention before comparing IPCW AUC, Brier scores, or decision curves across implementations. Version 1.0.0 retains the supplied convention.

## Shared AUC and Brier exception handling

The `metrics()` function in [`evaluate_survival.py`](../code/evaluate_survival.py) calculates IPCW AUC and Brier score within one `try` block. A `ValueError` from either the censoring-weight calculation or the AUC calculation assigns `NaN` to both metrics.

In particular, a development sample can have estimable censoring weights but no eligible case–control pairs within the same outer fold. AUC is then undefined, while a Brier score could still be calculated. The current implementation returns both as missing. This can also reduce the number of usable bootstrap replicates. Check the metric-specific reason for missing values rather than interpreting a missing Brier score as proof that its probability error is inherently unestimable.

## Categorical-code naming precision

The `Screen.encode()` method in [`analysis_core.py`](../code/analysis_core.py) creates indicator-column names using `f"{level:g}"`, whose default representation uses six significant digits. Different allowed numeric category codes can therefore produce the same column name. For example, `1000001`, `1000002`, and `1000003` all become `1e+06`; later indicators overwrite earlier ones in the output dictionary.

Use distinct small integer codes for categorical variables and maintain an explicit mapping in the clinical schema and data dictionary. Ordinary codes such as `0`, `1`, `2`, and `3` remain distinct. Review existing category codes before running the source snapshot; a change to the encoder should be recorded as a new implementation version.

## Extreme numerical values in absolute-risk conversion

The `breslow_probability()` function in [`attention_core.py`](../code/attention_core.py) subtracts the largest fitting log-risk before exponentiation and clips application exponents. Very large differences between fitting log-risks can still underflow the fitting risk-set sums to zero, making a Breslow baseline increment infinite or undefined. Clipping the application exponent does not repair an invalid fitting baseline.

The Cox probability conversion in [`analysis_core.py`](../code/analysis_core.py) also exponentiates application log-risks directly. Extreme values can exceed floating-point range. The runner checks final prediction finiteness and probability bounds, but these checks do not substitute for inspecting the fitted baseline and the range of risk scores. Investigate numerical saturation, nonfinite baselines, and warnings before using the corresponding absolute-risk estimates. Version 1.0.0 preserves these numerical paths.

## Small partitions and preflight coverage

Preflight validates the manifest, input interface, outer-fold labels, and basic outer-fold event and censoring counts. Full training creates smaller holdout and inner-fold partitions and applies additional model-level checks, including a minimum of five fitting events for `SparseCox`.

A successful preflight therefore does not guarantee that every nested partition is fit for model training. Inspect the chosen fold counts and event distribution before a full run, and use the emitted fitting and held-out identifiers to review a failed partition.

## Scope of interpretation

Missing or unseen categorical values are replaced using the training-partition mode, with counts recorded in the audit. Review these counts when a category is rare or differs across cohorts. Development C-index and AUC use only within-outer-fold pairs; their bootstrap intervals are conditional on the fitted predictions. These are implementation choices that should accompany reported results.

The supplied `requirements-tested-windows.txt` records the submitted environment. Dependency installation and input preflight are distinct from completing the full nested training workflow. Record the actual environment and successful run outputs when using this release for an analysis.
