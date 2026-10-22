# Position-level base concordance

`run_base_concordance(readback_receipt, structures, out)` reports every position
of each anchored candidate using the same primary-alignment, MAPQ20/BQ20 and
same-QNAME consensus rules as the frozen 2026-10-09 research audit. It does not
remap reads, modify the assembly or classify variants.

```python
from glycan_evidence.base_concordance import run_base_concordance

receipt = run_base_concordance(
    "structure/receipt.json",
    "structure/structures.tsv",
    "base_concordance",  # must not already exist
)
```

Use the original Linux compute host for production BAMs. Install the package's
`structure` extra there (`pysam==0.24.1`). No optional dependency is needed for a
zero-candidate run without an alignment. The function accepts strings or `Path`
objects and returns the completed output receipt.

## Input provenance

The source must be a completed `glycan-evidence.structure/1` run with a readable,
nonempty sibling `COMPLETE` marker. Historical stand-alone audit receipts are
not silently converted to this contract. The supplied structure table must match
the bytes and SHA-256 of the unique `structures.tsv` entry in `receipt.outputs`.
A relocated copy of that table is allowed; its content must remain identical.
The original contig FASTA and BAM paths bound by the receipt must be accessible,
and their bytes and SHA-256 are rechecked. Candidate contig and oriented interval
sequence hashes are checked separately. Input file size/mtime changes during the
scan cause failure. The source receipt and completion marker are themselves
hashed into the new receipt.

For an alignment, the source readback must report all of:

- `status=COMPLETE`;
- header dictionary `PASS_all_contig_names_and_lengths`;
- `qname_validation=PASS_contiguous_groups_no_duplicate_primary_mates`;
- `boundary_pair_validation=both_primary_mates_present_with_reciprocal_coordinates_and_strands`;
- `sqlite_integrity_check=ok`.

These checks bind the prior full-file validation to the unchanged BAM. The base
scan independently rejects repeated noncontiguous retained QNAME groups, also
under `python -O`. It is not a second execution of the prior full-file validator.
BAM is required for nonzero candidates; SAM and CRAM are not accepted by this
base audit. Candidate assembly intervals currently support A/C/G/T/N only, matching
the frozen audit's orientation transform; other IUPAC ambiguity is rejected.

A zero-row bound structure table is a valid result. If readback was
`NOT_REQUESTED` and no alignment was supplied, the output is `COMPLETE` with
`evidence_status=not_applicable_no_candidates`, an unresolved interpretation,
`bam_sha_reverified=false` and an explicit `not_revalidated_no_alignment_requested`
record. If a BAM was supplied, its SHA and source validations are still checked.
Neither path invents position rows, and neither means biological absence.

## Outputs and field meanings

The schema is `glycan-evidence.base-concordance/1`. Output contains one
`<structure_id>.positions.tsv.gz` per candidate, lowercase `receipt.json`, and
`COMPLETE` containing the receipt SHA-256. The receipt includes every position
file's byte size/hash, per-candidate flag summaries, filters, input bindings and
counting scope. Failed scans retain a FAILED receipt and no COMPLETE marker.

| Field family | Meaning |
| --- | --- |
| `genomic_position_0based` | Forward assembly-contig coordinate; rows are ascending in this coordinate. |
| `full_position_0based` | Reference-oriented interval coordinate; reverse candidates use `end - 1 - genomic_position`. |
| `assembled_genomic_base` / `assembled_oriented_base` | Assembly base in the original contig / reference-oriented interval. |
| `A/C/G/T/N_read_F/R` | Passing primary alignment observations on each BAM read strand. Overlapping mates each contribute. |
| `A/C/G/T/N_qname_F/R/both` | One agreeing base call per QNAME/position, grouped by contributing read strands. `both` is one QNAME, not two. |
| `conflict_qname_F/R/both` | QNAMEs with incompatible passing base observations; no consensus base is assigned. |
| `DEL_OR_SKIP_read_total`, `LOW_BQ_read_total` | Separate read-level counters, excluded from base consensus. Deletion and reference skip are combined; insertions have no reference-position row. |
| `acgt_qname_total` | Sum of agreeing A/C/G/T QNAME counts, excluding N and conflicts. |
| `assembly_fraction` | Assembly-base QNAME support divided by that A/C/G/T denominator; blank for a zero denominator or ambiguous assembly base. |
| `zero_acgt_support`, `assembly_base_no_qname_support` | Literal zero-support flags. The second includes uncovered positions. |
| `alternative_exceeds_assembly` | An A/C/G/T alternative has more QNAME support than the assembled base; ties are not flagged. |
| `alternate_ge2_qnames_both_orientations` | Alternatives supported by at least two agreeing QNAMEs and contributing reads in both orientations. |

All nucleotide count columns and alternative letters remain in genomic/BAM
orientation, including reverse candidates. F/R refer to read alignment strands,
not the candidate's direction. Only the explicitly oriented base and interval
coordinate are normalized. Consumers displaying alternatives in reference
orientation must complement their letters as well.

MAPQ 255 and MAPQ below 20 are excluded; BQ below 20 is counted separately.
Secondary and supplementary alignments are excluded. Duplicate-flagged primary
alignments are included, preserving the audit rule. Same-QNAME consensus avoids
counting two agreeing overlapping mates as two fragments; it is not PCR/UMI
deduplication or proof of original-molecule independence. N/conflict counts and
inspection flags are descriptive. Fractions are not calibrated variant calls,
taxonomic identity, absolute abundance or functional evidence.

## Small regression tests

```text
python -m unittest discover -s tests -p test_base_concordance.py -v
```

The suite checks source binding, missing markers, malformed inputs, zero-candidate
completion and preservation of existing output directories. Actual tiny BAM
fixtures additionally check overlap consensus, conflicts, reverse orientation,
MAPQ/BQ and primary filters, deletion/skip/insertion handling, alternate flags,
and noncontiguous QNAME rejection with optimized Python. Those BAM tests are
explicitly skipped when pysam is unavailable; skips are not successful BAM
validation. The original research audit and its results remain unchanged.

The base reader walks CIGAR coordinates explicitly: M/=/X consume query and
reference, I/S consume query only, D/N consume reference only, and H/P consume
neither. This corrects the pinned pysam 0.24.1 `get_aligned_pairs` padding
behavior, which advances the query index for P. Fixtures cover both a partial
interval that previously could silently report the wrong bases and a complete
interval that could fail at the final base, including base-quality positions.
Unsupported operators, nonpositive lengths, invalid internal clipping and
query/reference length inconsistencies are rejected. This is a coordinate
correction; the evidence filters, QNAME consensus and output columns are
unchanged. Archived first-build sources and previously frozen audit outputs
remain separately available.
