# Local reference support

`local-support` combines a whole-spectrum Hellinger fit with a second fit on 63 coordinates that had zero counts in every frozen training-background profile. It uses the unchanged candidate developed in 219 recorded conditions. The target reference family is longum; infantis is included in the background reference set.

The command derives all 220 marker counts, the complete paired-input denominator and read hashes from a sealed metrics run. It has no manual denominator, truth-label, marker-selection or threshold options.

```sh
glycan-evidence local-support \
  --metrics /completed/rf01_m03/identity/metrics/rf01_m03/gene_metrics.tsv \
  --bundle /models/local_candidate/bundle \
  --expected-manifest-sha256 70839d9f6f3879092c82c6779762a313d13210a6eeedb8ee672682302695e989 \
  --out /new/results/rf01_m03_local
```

Use the separately supplied candidate bundle and the frozen NumPy/SciPy dependencies from `.[compatibility]`. The expected hash above identifies the revision 3 interface snapshot. It requires isolated Python execution and rejects unexpected files, directories and links; scientific data and parameters are unchanged from the first snapshot. A different interface revision requires its own verified hash.

Input metrics must record `skip_qc=true` with complete paired-stream validation, unique QNAMEs, two distinct decompressed FASTQ hashes, the frozen mapper contract and a complete ordered marker universe. Explicitly pre-QC reads can be measured with `run --skip-qc`; the denominator must then be those actual input reads. Historical post-QC counts paired with only the original raw-input denominator are rejected. The separate BF39 input-reconstruction experiment does not change this command's input contract.

The output contains the measured full and local support, two numerical fit certificates, frozen working budgets, reference-column order and one observation:

- `conditional_target_local_support`: both fits meet their fixed working budgets and empirical support domains.
- `insufficient_evidence`: at least one condition is unmet; `reasons` identifies low/zero support, an out-of-domain value, a failed fit or numerical uncertainty.

The full empirical support domain is 19–5,944 events and the local domain is 19–2,545 events. Support counts are feature events: one fragment can contribute to more than one marker. Coordinates that are zero in the training background are not universally species-specific. The output does not estimate taxonomic identity or organism abundance; these fields remain null.

The bridge verifies candidate bytes before execution, copies only manifest-listed files, executes the verified copy in isolated Python, checks output/query bindings and rechecks input reports and helper-code hashes before creating `COMPLETE`. Failures retain logs and omit `COMPLETE`; existing output directories are preserved. For direct API use, invoke `python -I -B bundle/local_candidate_api.py ...`; the lower-level API validates the measured query schema but does not replace the sealed-metrics provenance bridge.

The portable candidate's first Linux regression used the already-exposed m03 mixture. It reproduced both prior Hellinger losses exactly (full 0.10004119611158535; local 0.08762122679481574) with two fits and no recalibration. This verifies the interface; it does not add an independent validation sample. The previous dev4 compatibility command and all frozen v2 behavior remain separate.
