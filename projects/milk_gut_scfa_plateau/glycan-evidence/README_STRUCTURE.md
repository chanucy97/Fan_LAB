# Configurable H5 anchor and readback module

This is a development module. It reports observed anchor geometry, read depth,
and boundary-spanning fragments. It does not yet classify structural alleles
with independently calibrated sensitivity or specificity.

## Run

Linux production requires BLAST `blastn`; mapping paired reads also requires
`bowtie2`, `bowtie2-build`, `samtools`, and Python `pysam`. Dependencies are
checked at use time. The reference dictionary/geometry and tiny SAM tests use
the Python standard library. NumPy, if present, accelerates depth summaries;
the standard-library fallback has the same units and thresholds.

```bash
glycan-evidence structure \
  --contigs assembly.fasta --anchors coding_anchor_queries.fasta \
  --r1 cleaned_R1.fastq.gz --r2 cleaned_R2.fastq.gz \
  --out structure_run --threads 8 --bin-dir /path/to/environment/bin
```

Alternatively, `--bam complete_name_collated.bam` accepts an existing complete
alignment. Do not also supply FASTQ. `--bam fixture.sam` supports small portable
fixtures. Without reads/alignment, geometry is reported and readback remains
unresolved. Output directories must not exist, including empty directories;
failed runs preserve their failure receipt and cannot be silently overwritten.

Python entry point:

```python
from glycan_evidence.structure import run_structure

receipt = run_structure(
    contigs="assembly.fasta", anchors="coding_anchor_queries.fasta",
    out="structure_run", r1=None, r2=None,
    bam="complete_name_collated.bam", config=None,
    bin_dir="/path/to/environment/bin", threads=8,
)
```

`config` accepts a JSON path or a dictionary. Omitting it loads the packaged
`src/glycan_evidence/presets/h5.json`. This preset applies only to the specified
H5 anchors. The general code has no sample names, contig identifiers, known
sample coordinates, or server paths.

## Input contract

- FASTA IDs must be unique first-token identifiers. Sequences must be nonempty
  IUPAC nucleotides. Whitespace inside sequence lines, duplicate IDs, and empty
  records fail. Case is normalized for sequence comparison and sequence hashes;
  original FASTA bytes are separately hashed.
- The anchor FASTA must contain both configured queries. Queries are supplied
  in their coding orientation. Other FASTA entries are preserved in the input
  file provenance but not searched. BLAST searches both relative orientations;
  the reported subject strand determines the assembly orientation relative to
  the supplied query. Restricting BLAST to `-strand plus` would incorrectly
  discard reverse-complement anchor matches.
- Config `anchors` contains exactly two IDs, `reference_order` 0 and 1, and a
  `reference_strand` (`+` or `-`) for each query. Order means genomic left to
  right on the chosen standard-reference orientation; it is not numeric gene
  order or coding direction.
- H5 default: genomic left anchor `Blon_2174`, right anchor `Blon_2178`; both
  supplied coding queries are minus-strand relative to the standard reference.
  Therefore a reversed assembled interval has physical left `Blon_2178` and
  physical right `Blon_2174`. These labels appear explicitly in every result.
- Paired FASTQ must have identical normalized names/counts, exactly one pair
  per name, valid four-line records, and complete streams. Trailing `/1` and
  `/2` are recognized. Full gzip reading checks stream termination. Input
  validation uses an exact disk-backed name index rather than sampling.
- External BAM/SAM must include **every** contig name and length in `@SQ`, even
  contigs that have no anchors. Duplicate names or missing decoys fail. This
  checks the declared reference dictionary, not the unknown producer's mapping
  history or the identity of sequences behind a same-name/same-length header.
  That provenance limitation is explicit in the receipt.
- QNAME records must be contiguous. Name-collate the **complete** BAM first if
  it was coordinate-sorted; do not extract only the selected target contig.
  Recurrence of a previously closed group is an error. Within a group, there
  may be only one primary record per mate. Secondary/supplementary records are
  permitted and inspected but do not add depth. A boundary-supporting group
  must include both primary mates with reciprocal coordinates and strands.
- QNAME must identify a sequencing fragment across the whole input, including
  merged lanes/read groups. An indistinguishably reused name is not an
  independently measured molecule. PCR duplicate flags are reported and not
  used to infer original-molecule independence.

## What the algorithm records

All geometrically coherent HSP pairs are retained; the highest scoring pair is
not silently selected. Pairing requires a common contig, consistent orientation
relative to the configured reference, expected physical order, and the configured
maximum intervening distance. Cross-contig, orientation, overlap/order, and
distance rejections are written separately. A partial or below-screen HSP pair
is retained and explicitly labeled unresolved.

The H5 exploratory screen is 85% nucleotide identity and 90% query span. Exact
inter-anchor length is reported only when both inward-facing query endpoints
are present without a terminal subject gap. Missing an outer query endpoint
does not itself move an intact inner boundary. Missing an inner endpoint
produces only an HSP-edge distance, not an exact gene-boundary distance. The
software does not stitch separate HSPs into a gene or infer evolutionary
deletion breakpoints.

For generated alignments, all input contigs compete in one Bowtie2 index with
the `--very-sensitive -N 1 -L 20 --score-min L,0,-0.4
--no-mixed --no-discordant -I 0 -X 1000` H5 defaults. The v2 default uses
Bowtie2's ordinary best-alignment reporting **without `-k` or `-a`**; setting
`-k 1` is not an equivalent way to request that mode. Thread count is configurable.
The complete read alignment is consumed, even when there are no paired anchors.

The separate `h5_k5_audit_v1.json` preset retains `-k 5` multialignment reporting.
Actual Bowtie2 2.5.4 fixtures exposed primary MAPQ 255 in this mode when a
second-best alignment was not sufficiently searched. [SAM's MAPQ definition](https://samtools.github.io/hts-specs/SAMv1.pdf)
specifies that 255 is unavailable, and the [versioned Bowtie2 implementation](https://raw.githubusercontent.com/BenLangmead/bowtie2/v2.5.4/unique.h)
contains this return path. An older `MAPQ >= 20` filter incorrectly included
255; the new code explicitly excludes it. This is a documented scientific
correction, so that filtered output is not expected to reproduce the older
mistake. Raw primary-read evidence and geometry can still be compared separately.
The [Bowtie2 manual](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#reporting)
describes the distinction between default-best and `-k` reporting. Available
MAPQ is an aligner estimate conditional on the supplied assembly, not an
experimentally calibrated guarantee of location or species identity.

Readback has two parallel coverage summaries: primary aligned bases with base
quality at least 20, and the same observations additionally requiring MAPQ at
least 20. MAPQ 255 means unknown and is excluded from the quality-filtered
summary. Insertions, deletions, reference skips, clipping and the unsequenced
fragment gap do not contribute base coverage. Both overlapping mates contribute
read depth: it is explicitly **read depth, not fragment depth**.

Boundary evidence requires a proper opposite-strand same-contig paired
alignment, MAPQ at least 20, an uninterrupted CIGAR M/=X block extending at
least 20 bp on each side, and base quality at least 20 at **every** base in that
40-bp window. Indels/skips break a spanning block. QNAME is counted once per
candidate boundary, even if both mates span it; mate-alignment counts and
distinct read-start/pair-start coordinates are also reported. These are anchor
boundaries, not reconstructed mutational breakpoints.

`evidence_level=read_supported_anchor_interval` means a screened candidate has
complete inward anchor edges, at least one quality-filtered observation at
every position of the assembled interval, and at least one observed fragment
across each boundary. This is a descriptive label. All continuous depths,
breadths and fragment counts accompany it; one weakly supported fragment is
not presented as a validated positive classification. Multiple candidates,
partial anchors and incomplete/low coverage remain explicit. The separate
`structural_classification` is `unresolved_no_validated_classification_rule`:
no independent classifier is implemented and no absence call is emitted.

## Outputs and audit

| File | Meaning |
| --- | --- |
| `effective_config.json` | The complete effective preset, including quality filters and scope |
| `anchor_queries.fasta`, `anchor_blast.tsv`, `anchor_hits.tsv` | Actual searched queries, raw HSPs, parsed validated HSPs |
| `structures.tsv` | Every coherent anchor pair; actual contig and reference-oriented interval SHA-256 |
| `rejected_anchor_pairs.tsv` | Rejected HSP-pair geometry with reason |
| `anchored_sequences.fasta` | Observed intervals oriented to the configured reference |
| `coverage.tsv` | Whole interval and inter-anchor/HSP-edge coverage; breadth at depths 1/5, minimum, median and mean read depth |
| `boundary_support.tsv` | Physical left/right anchors, boundary coordinates, QNAME and start-coordinate observations |
| `coverage_windows.tsv` | Quality-filtered depth in configurable windows (H5 default 50 bp) |
| `evidence_summary.json` | Observational evidence level and unresolved classification, without species or phenotype claims |
| `receipt.json`, `COMPLETE` | Versioned inputs/config/code/outputs, executable versions, commands and success/failure |

All coordinate fields explicitly use zero-based starts and exclusive ends.
Named boundary coordinates lie between the preceding and following base. Raw
BLAST fields `qstart`, `qend`, `sstart`, `send` in `anchor_blast.tsv` retain
BLAST's one-based inclusive convention and declared field order. Parsed
`anchor_hits.tsv` explicitly names these one-based fields and adds
`subject_start_0based` / `subject_end_0based_exclusive` genomic coordinates.
Physical left/right labels never imply a fixed gene name.

The QNAME index resides in a temporary SQLite file in the output directory (or
the supplied readback work directory). Exact uniqueness and SQLite integrity
checks are timed and recorded; the private temporary index is removed when
closed. Large alignment and FASTQ files are fully read and hashed. Runtime is
therefore proportional to the whole input, not just the selected locus. Failed
commands or malformed inputs write `status=FAILED`; they never write `COMPLETE`.
Scientific results are small TSV/JSON files; large BAM/index files remain with
the compute run.

## Verification scope

Local standard-library tests cover reference/reverse orientation, mixed query
coding strands, outer/inner partial anchors, different contigs, wrong order,
multiple candidates, distance limits, sequence-faithful BLAST parsing, output
protection, coverage units, overlapping mates, low/empty coverage, MAPQ 0/255,
per-base quality, indels, missing/inconsistent mates, incomplete reference
dictionaries, malformed CIGAR/qualities, duplicate/non-contiguous names and
truncated FASTQ.

On the compute host, the following enables actual BLAST/Bowtie2/samtools/pysam
tiny integration tests. Without these tools they are explicitly skipped.

```bash
export GLYCAN_EVIDENCE_TEST_BIN_DIR=/path/to/environment/bin
python -m unittest discover -s tests -p test_structure.py -v
python -m unittest discover -s tests -p test_readback.py -v
```

These fixtures check default-best forward and reversed assemblies, identical
competing contigs, retained unknown-MAPQ behavior under k5, and a bad BAM header.
They test implementation behavior, not external
accuracy. Existing cohort structures provide regression evidence; they are not
a new independent evaluation set. Exact recall of unseen alleles, carrier
assignment, false-positive/false-negative rates, substrate specificity, and
generalization beyond the H5 preset remain separate validation tasks.
