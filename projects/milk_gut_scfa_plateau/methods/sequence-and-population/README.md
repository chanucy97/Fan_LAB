# Sequence and population methods

Release adapter: **public-adapter.v1**.

This package exposes the analysis code used for assembly ORF assessment,
reference-backed substitution consensus, reference-background scoring,
population similarity and read support for one pre-specified structural
candidate. Tests generate synthetic records in temporary directories. Participant
sequences, sample metadata, results and the real candidate configuration are not
included.

## Modules

| Entry point | Operation | Inputs needed for an actual run |
|---|---|---|
| `comparison/compare_strains.py` | Eligibility-gated inStrain comparison and descriptive pair tables | Completed inStrain profiles and bound receipts; reference FASTA/STB/genes; sample/visit metadata; configuration |
| `sequence_integrity/assembly_scan.py` | Prodigal ORF integrity, reciprocal protein matching and translated-DNA HSP inventory | Existing assemblies and inventory receipts; frozen seed dictionary; mapping rules; tool provenance and executable paths |
| `sequence_integrity/profile_sequence_audit.py` | Depth-aware substitution consensus, mixed sites, coverage holes and translated substitutions | Original inStrain 1.11.0 outputs; reference/seed mapping plan; upstream gene evidence and source receipts |
| `sequence_integrity/enzyme_context/core.py` | Target-masked background and local/flank scores against a reference panel | Assembly/reference FASTA, all target/comparator mask intervals, hashes and a case manifest |
| `sequence_integrity/readback_stage2_deploy_v1/driver.py` | Readback of exactly one pre-specified, hash-bound candidate | Reviewed sequence templates/coordinates, source and selection manifests, original BAM/BAI/covT/Rdic, matching runtime and explicit run flag |
| `trait_mapping/map_traits.py` | Reciprocal protein matching/classification dependency | Seed dictionary, reference protein panel and DIAMOND |
| `carrier_context/carrier_context.py` | BTOP/coverage/chain scoring dependency | Supplied sequence/alignment records; no automatic downloads |

The readback driver does **not** discover variants. It compares four
sequence-distinguishable local templates at two supplied reference loci, checks
the original filtered-coverage calculation, then collapses support by sample and
QNAME. A failed coverage-parity check blocks event classification.

The original cohort-specific `summarize_phase1.py` administrative join is omitted:
its trust setup and tests depend on the private sample roster and inventory.
Stage outputs retain their original schemas and integrity receipts.

## Synthetic tests

Use Python 3.10 or later in an isolated environment:

```bash
python -m pip install -r requirements-test.txt
python -B run_synthetic_tests.py
```

The suites run in separate subprocesses because two historical directories use
the module name `core`. They test coverage boundaries, missing/unknown states,
partial ORFs, alternative starts, reverse coordinates, competing homologues,
multiple HSPs, background masks/ties, original read qualities, mate conflicts,
QNAME deduplication, and receipt tampering. Alignment tools and real read files
are not invoked. BioPython/h5py-dependent tests report explicit skips when these
packages are absent.

Synthetic examples are embedded in `test_assembly_scan.py`,
`test_profile_sequence_audit.py`, `enzyme_context/test_core.py`, and the readback
`test_core.py`/`test_driver.py`. A test named `actual_negative_strand` describes
actual coordinate arithmetic on **synthetic** DNA, not a participant record.

## Scientific rules retained

- Reciprocal protein candidates require 90% identity, 90% query and subject
  coverage and E-value ≤1e-10; the near-best bitscore fraction is 0.95. Search
  reports all HSPs; classification does not stitch partial HSPs into a protein.
- Complete predicted ORFs require `partial=00`, an unambiguous in-frame CDS,
  table-11 start and terminal stop, no internal stop and matching translation.
  Partial predictions and translated-DNA-only evidence remain separate.
- Profile consensus uses depth ≥5 and the original caller; missing coverage,
  ambiguous bases and unresolved calls remain unknown. It does not assess indels
  or full-length phasing and does not establish an intact sample ORF.
- Population eligibility requires depth ≥5 over at least 50% of the entire
  reference, including at least 50% jointly comparable sequence. Eligible
  popANI ≥0.99999 is high similarity. Low similarity is not itself replacement;
  missing comparison is not absence.
- Context scoring retains the 1,000-nt mask buffer, 10,000-nt flanks, coverage
  0.8, covered identity 0.95 and unique-margin 0.02 settings. Scores express
  reference-panel affinity, not taxonomic truth or strain identification.
- Readback templates and rules are explicit inputs. Synthetic tests demonstrate
  the original 20-nt exact anchors, original BQ ≥30 and MAPQ ≥10. Technical
  support gates remain explicit in `TechnicalRules`; tests preserve their
  original 10 fragments, 3 starts, at least 2 fragments in each direction,
  0.9 sample-flank deletion fraction and 10 retained-reference fragments for
  the coexistence rule. No biological function is inferred from these gates.

## External dependencies and scope

Actual profile/comparison execution expects **inStrain 1.11.0** and its original
data schemas. Assembly execution expects **Prodigal 2.6.3**, DIAMOND, `tblastn`
and `makeblastdb`. Context execution needs `blastn` and `makeblastdb`. Readback
additionally needs pysam, h5py and the exact source/runtime receipts for the
original profile. All tools, reference genomes and licensed/restricted data must
be supplied separately; this package does not fetch them.

The study CLI integrity wrappers retain the fixed 50-library, 16-locus,
six-positive/three-comparator, four-module and eleven-reference checks where
present. These are cohort design assertions, not biological thresholds. The
pure functions and synthetic fixtures can be used independently; running the
study wrappers requires their complete input schemas and declared provenance.
Do not fabricate receipts to make a different dataset pass. This is a public
adapter of a study workflow, not a general end-to-end discovery product.

`INPUT_SCHEMA.md` describes external bindings. `PROVENANCE.json` records original
workspace-relative source paths, original/public SHA-256 values and every
adaptation. Public copies do not modify the frozen study sources. Changes are
limited to path portability, the explicit single-candidate identity parameter,
rebinding adapted code/configuration hashes, and removing private/deployment
fixtures. Numerical algorithms are preserved.

Quick CLI checks (these do not load data):

```bash
python comparison/compare_strains.py --help
python sequence_integrity/assembly_scan.py --help
python sequence_integrity/profile_sequence_audit.py --help
python sequence_integrity/enzyme_context/core.py --help
python sequence_integrity/readback_stage2_deploy_v1/driver.py --help
```
