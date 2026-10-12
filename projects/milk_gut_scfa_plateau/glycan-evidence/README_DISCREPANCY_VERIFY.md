# Independent discrepancy verification

`discrepancy_verify.run_discrepancy_verify(root, sample, expected_freeze_sha256, expected_results_sha256, out)` replays a sealed complete-QC marker-discrepancy result. It reads saved data and numerical witnesses, calls no optimizer, imports no input Python and creates no new classification, q value, probability or abundance estimate. This module is separate from the existing compatibility and local-support entry points and does not alter their scientific rules or any released package.

The two expected SHA-256 values must be supplied from a trusted handoff: one for `root/STAGE2_FREEZE.json`, one for `root/results/<sample>/STAGE2_RESULTS.json`. Computing expected hashes from untrusted inputs in the same command is not authentication. `out` must be a new directory outside the protected `results` and `frozen` input trees. The command-line interface additionally seals its completed small outputs:

```bash
glycan-evidence discrepancy-verify --root /path/to/sealed/discrepancy \
  --sample sample_name --expected-freeze-sha256 TRUSTED_FREEZE_SHA256 \
  --expected-results-sha256 TRUSTED_RESULT_SHA256 --out /path/to/new/verification
glycan-evidence verify --directory /path/to/new/verification
```

The SHA placeholders must be replaced with the 64-character values from the trusted execution record. Install `.[discrepancy]` for the NumPy dependency. The equivalent Python interface is:

```python
from glycan_evidence.discrepancy_verify import run_discrepancy_verify

receipt = run_discrepancy_verify(
    root="/path/to/sealed/discrepancy",
    sample="sample_name",
    expected_freeze_sha256=trusted_freeze_sha256,
    expected_results_sha256=trusted_result_sha256,
    out="/path/to/new/verification",
)
```

The supported input is the existing stage-2 layout: `STAGE2_FREEZE.json`, its complete file manifest, `FREEZE.json` and its manifest, `PREFLIGHT.json`, frozen model/auditor/mapper members, and the sample's `STAGE2_RESULTS.json`, `RECEIPT.json`, `QC_BRIDGE_PASS`, query, recount, three fits, marker residuals and descriptive result. Python members are hash-checked as bytes only. Model manifests are checked transitively. Every consumed file is bounded to 32 MiB, with a 128 MiB total budget; the verifier never follows recorded absolute paths into FASTQs, candidates, SAMs, executables or external directories.

Verification checks include:

- Complete QC receipt and descriptor hashes; raw, post-QC, screened-candidate and mapped-QNAME denominators; recorded raw-read audit; full QC hashes and all-fragment bucket counts; unchanged QC/screen arguments and recorded tool versions; candidate-content audit and frozen panel/mapper binding. The final recorded pair of QC gzip sizes must fit the frozen cap. Historical `--no-unal` omissions remain explicit.
- All 220 ordered markers and all 16 buckets, exact integer addback, count bounds and zero-support/empty-bucket retention. A fragment may support multiple markers, so their summed support is not a fragment total.
- Full, background and local saved simplex witnesses, source/coordinate order, query/training read separation, predicted scaled rates, Hellinger losses, unsupported support, convex gaps and feasible-weight repair. Feasibility and arithmetic tolerance remain `1e-12`, convex-gap tolerance `1e-8`, following the frozen independent auditor. Raw solver weights may be slightly negative within its existing `1e-12` repair tolerance; accepted witness weights, counts and nonnegative quantities may not be negative. Signed residuals are signed quantities.
- Every marker's full/background/local loss contribution and weighted rate-residual decomposition. This retains the existing technical partition interpretation; the partitions are not biological replicates or sampling confidence intervals.
- Descriptive summary denominators and null inference fields. The full empirical support range comes from the sealed model parameters. The current local diagnostic contract uses the existing `[19, 2545]` annotation; no envelope is fitted or applied. A changed scientific schema needs a separately reviewed adapter.

Successful output consists of `receipt.json`, `replayed_fits.json`, `marker_residuals.json` and `COMPLETE`, whose contents are the receipt SHA-256. The receipt binds input hashes, verifier hash, runtime NumPy version and output hashes, and records zero new optimizer/H/q calls. The root CLI can apply the usual small-output artifact seal afterward. Incomplete or changed input raises an exception; output is never overwritten and failure never produces `COMPLETE`.

The receipt verifies **recorded complete-QC and mapping audit evidence plus saved floating arithmetic**. It does not rescan raw or QC FASTQs, rerun mapping, re-prove each QNAME signature, independently establish that historical source data were truthful, or validate biological identity/function. Those historical measurements remain the responsibility of the sealed producer. Real-library replay is a separate acceptance gate from the small synthetic tests.

Runtime requires NumPy (available in the package's compatibility environment); SciPy is not used. The tests use a hand-solvable 220-marker, 16-bucket fixture with raw `N=192`, post-QC `N=160`, 40 candidates and 32 mapped groups. They exercise successful replay, zero support, empty buckets, hash/path/type tampering, denominator substitution, marker/bucket errors, altered witnesses/residuals, forbidden inference fields, output overwrite and mutation before completion.
