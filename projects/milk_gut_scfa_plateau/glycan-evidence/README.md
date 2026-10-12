# GlycanEvidence — dev7 source review

GlycanEvidence is research software for recording separate subspecies-marker and glycan-locus sequence evidence. This code-only export preserves the accepted `0.1.0.dev7` Python implementation. It includes deterministic CDS replay and saved discrepancy-certificate verification. It is not a clinical classifier.

## Included scope

- `cds-replay` reproduces coding-sequence checks from an explicitly bound input package and saved alignments; it does not run a new alignment, assembly or optimization.
- `discrepancy-verify` replays saved complete-QC numerical certificates without new optimizer calls. Expected hashes must come from an independent trusted record.
- The source includes the existing identity, structural, profile and local-support modules, with their original input contracts.
- Unit tests contain small artificial inputs. Historical cohort regression receipts, real replay fixtures and study data are excluded.

## Installation and inspection

The base package and CDS replay require Python 3.11 or later. The pinned NumPy 2.5.3 extras require Python 3.12 or later. Install from this directory:

```bash
python -m pip install .
glycan-evidence --version
glycan-evidence --help
glycan-evidence cds-replay --help
```

For discrepancy verification, install `.[discrepancy]`. BAM-reading functions require the `structure` extra; compatibility fitting has separate dependencies. Details of the two verification commands are in `README_CDS_REPLAY.md` and `README_DISCREPANCY_VERIFY.md`.

## Registered rules and separate inputs

The original `src/glycan_evidence/legacy_v2/frozen_rules.json` is included without modification and is bound to an exact SHA-256 by the unchanged implementation:

`2c2f9650d932b7439873067caf3aad58d8066e6ea86f6228953e7573360dcd25`

Its pilot calibration used two public reference genomes. The saved calibration counts do not estimate population-wide diagnostic accuracy. The non-runtime historical eligibility note mentions internal study IDs to document analysis inclusion decisions; it contains no participant attributes, patient records, measurements, sequences or clinical mapping. This note is retained to preserve the registered rule and reproducibility.

Default rule loading is supported. Identity CLI commands still require an explicit `--rules` argument. Reference panels and complete replay inputs are supplied separately; changing the registered rules or supplying a modified replacement is rejected.

The full historical acceptance suite is not claimed for this code-only export. A small synthetic CDS check can be run without study data:

```bash
python -m unittest discover -s tests -p test_cds.py -v
```

## Data and reuse terms

No participant-level records, study sequences, sample mappings, FASTQ/BAM files, real cohort fixtures or internal infrastructure logs are distributed here. Study data are restricted by institutional requirements; requests are handled by the corresponding author through the manuscript's stated process. Code availability does not imply approval to release the input data.

No code reuse licence has yet been selected for this candidate. This README does not grant an open-source licence. The version label identifies the preserved implementation, not a claim that all manuscript inputs or workflows are included.
