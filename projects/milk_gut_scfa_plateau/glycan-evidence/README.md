# GlycanEvidence

GlycanEvidence links taxonomic-marker evidence to local glycan-locus reconstruction in shotgun metagenomes. It records marker support, assembled anchor intervals, direct boundary fragments and per-base sequence concordance as separate evidence layers. This makes alternative locus structures inspectable when a fixed reference panel leaves a local interval unresolved.

This public source snapshot is **0.1.0.dev10**, dated **22 October 2026**. It extends the previously published dev7 source with candidate-structure inference (dev8), joint fragment assignment (dev9) and assembly-assisted structure recovery (dev10). The executed scientific implementation and frozen v2 identity rules are preserved byte-for-byte. The [release record](RELEASE_STATUS.md) describes the software scope; [the project methods map](../METHODS_CODE_MAP.md) connects these modules to companion comparison and sequence-analysis code.

## Method layers

| Layer | Interface | Output |
|---|---|---|
| Identity and locus signals | `run`, `call` | Separate marker/locus metrics and frozen identity evidence calls |
| Reference-guided recruitment and assembly | `assemble` | Recruited paired reads, local contigs and input provenance |
| Observed locus architecture | `structure` | Anchor geometry, complete/partial intervals, coverage and boundary support |
| Candidate compatibility | `structure-infer` | Fragment-compatible candidates and identifiable weight ranges |
| Joint fragment allocation | `structure-assign` | Model-dependent architecture allocations and support bounds |
| Assembly-assisted recovery | `structure-recover` | Reference calls beside measured, read-supported assembly intervals |
| Sequence checks | `bases`, `cds-replay` | Per-position support and coding-sequence comparisons |
| Integrity and evaluation | `verify`, `benchmark` | File checks and fixed-truth evaluation summaries |

The conditioning module is a Python research API, documented in [README_STRUCTURE_CONDITIONING.md](README_STRUCTURE_CONDITIONING.md); it has no standalone CLI. Compatibility, local-support and discrepancy-certificate commands use additional explicitly bound research bundles described in their READMEs.

## Install

The base package requires Python 3.11 or later. Use **Python 3.12 or later** for the pinned NumPy/SciPy inference dependencies. From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install './projects/milk_gut_scfa_plateau/glycan-evidence[inference]'
glycan-evidence --help
```

On Windows, activate with `.venv\Scripts\Activate.ps1`. For BAM processing on a supported platform, install the `structure` extra (`pysam==0.24.1`). Raw-read workflows also require fastp, Bowtie2, samtools, BLAST+ and, for local assembly, MEGAHIT. Reference sequences and indexes are explicit inputs.

```bash
glycan-evidence doctor --bin-dir /path/to/bioinformatics/bin --structural --require
```

This command checks the environment without installing tools or downloading references. Run production sequencing workflows on a suitable Linux compute host. Small synthetic tests and SAM/table interfaces can run locally.

## Test with artificial inputs

```bash
cd projects/milk_gut_scfa_plateau/glycan-evidence
python -m unittest discover -s tests -v
```

The tests construct artificial reads, alignments, tables and bound-file packages. No cohort files or server account is required. Numerical tests use the inference extra above; BAM-specific execution additionally requires pysam.

## Input-specific workflows

For paired reads and a verified marker panel:

```bash
glycan-evidence run --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --sample sample_A --panel /path/to/verified_panel --rules /path/to/frozen_rules.json \
  --out /new/path/sample_A --threads 12 --bin-dir /path/to/bin
```

For existing contigs and paired reads:

```bash
glycan-evidence structure --contigs sample_contigs.fasta --anchors H5_anchors.fasta \
  --r1 recruited_R1.fastq.gz --r2 recruited_R2.fastq.gz \
  --out /new/path/structure_A --threads 12 --bin-dir /path/to/bin
```

For a source-bound roster of ordinary-assembly candidates and reference calls:

```bash
glycan-evidence structure-recover --roster ROSTER.json --case-id sample_A \
  --expected-roster-sha256 YOUR_ROSTER_SHA256 --out /new/path/recovery_A
glycan-evidence verify --directory /new/path/recovery_A
```

Prepare the roster according to [the recovery input schema](README_STRUCTURE_RECOVERY.md). The command verifies its bindings; roster preparation is a separate input-adapter step. Each workflow requires a new output directory to keep failed or previous runs separate.

## Detailed method documents

- [Frozen identity engine](README_ENGINE.md), [recruitment and assembly](README_ASSEMBLY.md), [anchor structures](README_STRUCTURE.md).
- [Fragment inference](README_STRUCTURE_INFERENCE.md), [joint assignment](README_STRUCTURE_ASSIGNMENT.md), [assembly-assisted recovery](README_STRUCTURE_RECOVERY.md).
- [Base concordance](README_BASES.md), [CDS replay](README_CDS_REPLAY.md), [internal-region conditioning](README_STRUCTURE_CONDITIONING.md).
- [Marker profiles](README_PROFILE.md), [mixture diagnostics](README_MIXTURE.md), [reference-family compatibility](README_COMPATIBILITY.md).
- [Local-support bundle interface](README_LOCAL_SUPPORT.md), [discrepancy-certificate verification](README_DISCREPANCY_VERIFY.md).

## Evidence and evaluation scope

The contribution is a traceable separation of identity, local structure and sequence support. Assembly-assisted recovery includes structures already observed by ordinary assembly; the existing 35-condition development comparison recovered 23 correct classes with no additional class calls, matching ordinary-assembly geometry. This is development evidence, not an independent demonstration of superiority over assembly. Earlier held-out sources became exposed development material once inspected and used for later design.

An unresolved interval remains unresolved. Read depth counts aligned mates; boundary QNAME counts represent fragments, not UMI-verified molecules. Local structure alone does not identify its taxonomic carrier or establish substrate use. The frozen identity classifier retains its documented longum/infantis classification failures. Command-specific documents state the corresponding units, inputs and gates.

## Data, citation and reuse

See [Data and code availability](../docs/DATA_CODE_AVAILABILITY.md) and [Code availability](../CODE_AVAILABILITY.md). Study-specific inputs follow the stated data-access route. Cite the software version and an immutable Git commit or release tag. No software reuse licence has yet been selected for this project; public source availability is not represented as a permissive licence.
