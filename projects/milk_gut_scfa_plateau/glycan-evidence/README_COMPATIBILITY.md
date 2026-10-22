# Reference-family compatibility from completed metrics

The `compatibility` command bridges a completed frozen identity-mapping run to
the unchanged, separately distributed reference-family research v0.1 package.
It constructs the full joint220 query from measured counts and execution
receipts. It does not accept a manually supplied read count or FASTQ checksum.

```bash
glycan-evidence compatibility \
  --metrics /run/metrics/sample_A/gene_metrics.tsv \
  --bundle /releases/reference_family_research_v01/package \
  --expected-manifest-sha256 80747c1887301503f7f3647a721d45aeae25857a0c2caeecf6b80a0fc1566a75 \
  --out /new/path/compatibility_A
```

Python entry point:

```python
from glycan_evidence.compatibility_io import run_compatibility
receipt = run_compatibility(metrics, bundle, expected_manifest_sha256, out)
```

The wrapper requires only the standard library. Executing the frozen package
requires its numerical dependencies in the same Python environment:
`numpy==2.5.3` and `scipy==1.18.1`. The expected manifest digest is a release
trust anchor: obtain it from the release record, not a replacement bundle.
The bundle is not downloaded or silently updated by this command.

## Accepted inputs

`metrics` must be the `gene_metrics.tsv` in the sealed
`run/metrics/sample/` layout produced by `glycan-evidence run`. Required sibling
files include the full feature-by-margin table, identity adapter receipt,
completed mapping manifest, screen counts, summary and COMPLETE marker. The
top-level artifact manifest must bind these files. Historical compact tables
without a complete adapter input audit are rejected.

The bridge verifies the unchanged engine, panel/index and recorded mapper
version contract, all feature-by-margin rows and the exact eligible 220-marker
set. It obtains N and the two decompressed-FASTQ SHA-256 values from the
adapter's complete input audit. Both streams must have reached EOF with unique
paired names. The adapter input records must match the mapper's input records;
N must equal the independently recorded no-QC screen total, and candidate and
competitive-QNAME counts must agree with their summaries. Missing marker rows
are rejected; they are never filled with zero.

The initial bridge **requires recorded `skip_qc=true`** and rejects a contradictory
fastp stage. QC state is checked separately because the historical mapping
digest does not include it. This restriction matches the release's evaluated
simulation-to-mapping route. Applying the default QC run to raw clinical reads
does not qualify it for this bridge merely because its mapper hash matches.
Pre-QC reads can be mapped as explicit inputs, but their earlier processing and
real-cohort calibration are not established by this interface. Existing read
hashes are verified as recorded evidence; this bridge does not reread the large
original FASTQs or infer their biological origin.

## Execution and outputs

Before any bundle code is executed, its manifest is checked against the supplied
digest and every expected payload byte is verified. Only verified files are
copied into a fresh `verified_bundle/`; source-directory extras and cached
bytecode are excluded. A separate Python process executes that immutable copy
with environment Python-path overrides and user-site imports disabled. The
process explicitly enables UTF-8 mode (`-X utf8`) on every platform, including
Windows environments whose default text encoding is not UTF-8. The
original package remains unchanged and rechecks its own numerical dependencies,
reference binding, training-read separation and model contracts.

- `query.json`: exact ordered220 counts, actual N and input/mapping bindings.
- `result.json`: unchanged research API result, including Hellinger compatibility
  certificates and the separate KL descriptor.
- `api.stdout.log`, `api.stderr.log`: original subprocess output.
- `receipt.json`, `COMPLETE`: input, copied package, query, result and execution
  bindings; COMPLETE contains the receipt SHA-256.

Existing output directories are refused. A subprocess or contract failure
retains a FAILED receipt without COMPLETE. Zero support, an out-of-domain
query, reference incompatibility and a numerical unresolved result remain
explicit research outcomes rather than becoming identity calls.

## Interpretation

The target family is **longum**; infantis reference profiles belong to the
background family in this release. `identity_call`, P values, confidence
intervals and organism abundance remain null. S is the sum of per-marker
fragment counts and may count one fragment at multiple markers; N is the
separate actual input-pair denominator. A compatible background explanation
does not establish target absence, and passing this software interface does
not create a validated identity classifier or calibrate real clinical data.
