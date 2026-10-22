# Frozen identity engine adapter

The public adapter is `glycan_evidence.identity`. Its scope is the registered
2026-10-07 v2 competitive reference panel and its pilot identity rules for
`infantis` and `longum`. It does not select thresholds, infer biological absence,
estimate absolute abundance, establish glycan use, or link a locus to a host.

## API

```python
from glycan_evidence.identity import DEFAULT_RULES, apply_rules, run_sample

# Linux production execution. The destination must not already exist.
run_sample(r1, r2, sample, panel, out, threads=8, skip_qc=False, bin_dir=None)

# Standard-library-only compact-table verification and calls on any OS.
receipt = apply_rules(metrics_root, ["sample1", "sample2"], DEFAULT_RULES, calls_out)
```

`run_sample` returns the original completed mapping manifest plus an
`adapter_receipt`. The panel must match the registered v2 reference files and
Bowtie2 index hashes. New reference panels require a separately validated preset.
The registered rules require Bowtie2 2.5.4. Production mapping also requires fastp
unless `skip_qc=True` is explicitly selected. Simulated input and real input must
remain distinguishable by their recorded QC policy.

`verify_completed_metrics(metrics_root, sample, rules=DEFAULT_RULES)` is a public
read-only check returning `data`, `provenance`, and `rules`. `load_rules()` and
`verify_engine()` reject modified registered rule or engine files.

## Preservation and provenance

`legacy_v2/run_sample.py`, `map_evidence.py`, and `calibrate_calls.py` are byte-for-byte
copies of `scripts/rebuild20261007_method_v2/`. Original machine-path provenance
is retained in the study archive; the public release manifest and `identity.py`
record the implementation and rule hashes. The shipped `frozen_rules.json` is also
an exact copy. No legacy scientific logic has been edited. The new adapter
adds input guards and publication protection around this fixed implementation.

The registered rule SHA-256 is
`2c2f9650d932b7439873067caf3aad58d8066e6ea86f6228953e7573360dcd25`.
The primary rule remains margin 0, marker breadth >= 0.1, >= 1 supporting paired
fragment per marker, and >= 20 supported eligible markers. Feature metadata,
all four margin grids, panel indices, code and Bowtie2 version must match the
frozen mapping contract. There is no automatic calibration path in the adapter.
Historical calibration code is preserved solely for reproducibility and its
regression tests; calling it directly is not the public registered preset.

## Input and completion checks

- Explicit sample IDs must be nonempty, unique, and safe single path components.
- Completed tables require `COMPLETE`, a successful `run_manifest.json`, a full
  feature-by-margin grid, valid finite values, metric SHA-256, exact feature
  schema, and matching panel/code/version contract. Duplicate JSON keys,
  column names, feature-margin rows and output provenance IDs are rejected.
- The adapter refuses missing samples instead of representing unread inputs as
  valid biological results. All supplied sample IDs must pass before it writes
  a calls result.
- New raw runs fully scan paired four-line FASTQs before mapping. Pair names,
  `/1`/`/2` or Illumina mate labels, sequence/quality lengths and characters,
  matched EOF, and unique normalized QNAMEs are checked. Exact uniqueness uses
  disk-backed SQLite in a transaction, not a growing in-memory name set.
- After mapping, a second disk-backed check rejects a SAM QNAME reappearing
  after its group ended. Arbitrary group order is accepted; non-contiguous
  repeated groups are not. Same-QNAME alternative alignments remain together.
- Work runs in a hidden staging directory. The visible destination is renamed
  into place only after every required adapter check passes. Errors preserve
  diagnostics with `ADAPTER_FAILED.json` and remove the staging completion
  marker. Existing destination paths are never overwritten, even if empty.

The SQLite checks require temporary disk proportional to unique QNAME count and
add full sequential scans. They belong on the production compute host. The
checks record elapsed time; their full-scale runtime has not yet been measured
for this packaged version. The raw SHA fields in the adapter receipt hash the
**entire decompressed FASTQ byte streams**, not compressed files. Legacy
`raw_inputs` retain their original file-stat-only meaning. The original index
checksum cache is reused only with its exact complete index set and unchanged
size/mtime, as in frozen v2; it is not a fresh cryptographic rehash of cached
index contents. Without that cache the engine computes full index hashes.

Original manifests retain their actual staging execution paths. The adapter
receipt records both `execution_directory` and `published_directory`, giving a
single directory-prefix relocation map. Relative artifact basenames and their
hashes remain unchanged. Historical compact results lack the new full FASTQ/SAM
scans: applying rules to those results does **not** retroactively claim QNAME
uniqueness or a raw-file rehash.

## Outputs

`apply_rules` writes:

- `identity_calls.tsv`: sample, target, identity_state, reason,
  supported_markers, eligible_markers, metrics_available, flag_scope_panel,
  margin, marker_breadth, minimum_fragments_per_marker,
  supported_identity_markers, coverage_unit, scope.
- `identity_calls.json`: identical calls with schema and scope.
- `receipt.json`: rule/engine/adapter/input/output hashes and input checks.
- `COMPLETE`: the receipt's SHA-256.

States are exclusively `detected` and `insufficient_evidence`. Coverage is in
reference-feature bases, and depth counts each paired fragment's union per
base, preserving the original treatment of overlapping mates. These are not
unique original molecules, absolute microbial abundance, or gene activity.

## Verification performed

The three frozen engine test files were migrated with package-scoped imports
only; their original assertions and synthetic fixtures were retained. All
44 frozen tests pass. The adapter tests cover threshold boundaries, provenance
drift, corrupt rules/hashes, duplicate identities, incomplete runs, staged
publication, failed output writes, complete FASTQ hashes, paired EOF/name
errors, and non-contiguous SAM QNAMEs. These synthetic software tests are not
estimates of classifier performance.

`tests/legacy_v2/COHORT_REGRESSION.json` records reapplication to all 50 completed
infant samples: 100 identity calls, 10 compatible fields per call, **zero
differences** from the frozen cohort calls. The detected counts remain 3 for
infantis and 24 for longum. This is backward-compatibility evidence, not new
independent biological validation. The fixture feature JSON contains reference
metadata only; its generated evidence is artificial and is not patient data.

The full packaged raw-data path has now also run with real executables on Wuxi:
Slurm `10902332`, immutable package snapshot `snapshot_394e1afde703`. Two previously
inspected simulated inputs (1,740 infantis read pairs and 1,640 breve read pairs)
were rerun under the original skip-QC policy, and the positive input additionally
exercised the real fastp branch. All three CLI runs, completion/provenance checks,
QNAME scans and `verify` commands passed. Every resulting metrics table was
byte-identical to its original frozen counterpart, including the fastp run.
This is two input sets and three execution conditions, not three independent
biological samples or an external test set.

The environment was Python 3.12.14, Bowtie2 2.5.4 and fastp 1.3.7. Each task
requested 4 CPUs and 24 GiB on `intel`; observed wall times were approximately
9.73, 9.69 and 3.17 seconds. Child-process peak RSS was 3.37–3.42 GiB; this is
not an aggregate job memory measurement. Slurm's MaxRSS field was unavailable
for these short jobs. The fetched compact package contained 74 files, all
verified by SHA-256 and size. Evidence and pre-frozen expectations are under
`outputs/20261008_tool_development/integration/identity/`. Full-size clinical
FASTQ runtime, resource scaling and externally independent performance remain
separate release gates.
