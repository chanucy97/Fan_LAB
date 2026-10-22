# Assembly-assisted structure evidence

`structure_recovery` reconciles fixed-reference model calls with candidate
intervals reconstructed by the same library's ordinary assembly. It preserves
measured structures outside the reference panel and reports the supporting
coverage, direct boundary fragments and available per-base concordance.

The primary rule reproduces the existing `structure.summarize_evidence`
observational rule. This is a research evidence module, with no calibrated
structural classifier, taxonomic carrier assignment or biological novelty test.

## Verified command-line interface

```powershell
glycan-evidence structure-recover --roster ROSTER.json --case-id BF_80_D90 --expected-roster-sha256 <sha256> --out new_result_directory
glycan-evidence verify --directory new_result_directory
```

The command requires a new output directory and the expected roster SHA-256.
It verifies the bound source files, checks that the reference result belongs to
the selected case and contains the roster's exact reference calls, runs the
recovery core, and seals the output artifacts. A failed run has no COMPLETE
marker. It does not regenerate the adapter's normalized records from raw reads;
use the source-bound study adapter for that transformation.

The installed dev10 wheel passed 105 targeted tests and actual command-line
replays for an assembly-supported case and a case without an assembled interval.
These checks establish executable behavior; measured performance is reported in
the dated study evaluation.

## Python interface

```python
from glycan_evidence.structure_recovery import recover_structure

result = recover_structure(
    case_id=case["case_id"],
    reference_calls=case["reference_calls"],
    panel=roster["panel"],
    candidates=case["assembly_candidates"],
)
```

The core is deterministic and does not modify input records. It validates
identifiers, duplicate candidates, source labels, provenance field shapes,
coordinates, finite values, count denominators and exact sequence hash formats.
`validate_file_binding(record)` separately verifies an actual local file's size
and SHA-256. A source adapter must reconstruct normalized records from their
bound tables; the presence of a JSON hash field alone does not prove that
normalization was performed correctly.

The existing-data adapter is
`outputs/20261016_structure_recovery/assembly_design/build_roster.py` in the
study workspace. Its `ROSTER.json` is the complete, truth-free 82-case example.

## Roster schema

Top-level keys:

- `schema`: `truth-free-existing-assembly-roster.v1`.
- `source_inventory`: local file binding `{path, bytes, sha256}`.
- `panel_sources`: bindings for the frozen reference metadata and FASTA.
- `panel`: `{candidate_id, class_id, inner_length, full_sequence_sha256,
  inner_sequence_sha256}` records.
- `cases`: records containing `case_id`, `kind`, `v2_source`,
  `reference_calls`, `assembly_candidates`, `historical_read_supported_ids`
  and `input_inventory_status`.

`kind` and the opaque case ID are export metadata; they do not set thresholds.
`historical_read_supported_ids` verifies replay of the old rule and is not an
input to the core decision.

Each assembly candidate contains:

- `case_id`, `source_kind="own_library_ordinary_assembly"`, `structure_id`.
- `source_bindings` with `receipt`, `structures`, `coverage`, `boundaries`,
  `sequences`, and optional `base_audit` file bindings.
- `contig`, `contig_length`, `contig_sequence_sha256`.
- `interval_start`, `interval_end`, `left_boundary`, `right_boundary` in
  physical, zero-based, half-open coordinates; `orientation` is `+` or `-`.
- `anchor_screen_pass`, `inner_edges_complete`, `outer_anchors_complete`.
- Exact `inner_length` and `inner_sequence_sha256` when both inward edges
  are complete; otherwise both fields are null. The observed reference-
  oriented interval always has `interval_sequence_sha256`.
- `coverage`: `unit="primary_mate_alignment_read_depth"`,
  `filter="quality_filtered"`, `segment="full_anchored_interval"`,
  matching `start/end`, `breadth_ge_1` and `min_read_depth`.
- `boundaries.left/right`: matching `position`,
  `evidence="direct_read_spanning"`, integer `qnames` and
  `distinct_pair_starts`.
- `base_concordance`: null or a source-bound per-base summary with
  `positions`, `zero_acgt_support_positions`,
  `assembly_base_no_qname_support_positions` and
  `alternative_exceeds_assembly_positions`; additional position counts must
  also be nonnegative integers bounded by the interval length.

## Fixed evidence tiers

The primary tier requires the original anchor screen, complete inner anchor
edges, quality-filtered coverage of every base in the observed interval, and
at least one directly spanning QNAME at each physical boundary. The descriptive
output label is `read_supported_assembly_candidate`.

The strict tier additionally requires minimum interval read depth 3, at least
3 boundary QNAMEs per side, and at least 2 distinct pair-start coordinates per
side. Both tiers are exported. Their thresholds are development evidence rules,
not empirically calibrated diagnostic cutoffs.

Read depth counts aligned mates; QNAMEs identify fragments within the input,
not UMI-verified original molecules. Boundary counts originate from all-contig
readback using MAPQ and baseQ at least 20, excluding unknown MAPQ 255, and a
direct uninterrupted 20-base flank on either side. Unsequenced insert gaps do
not count as direct boundary coverage.

## Interpretation and comparison

`reference_supported_classes` retains the original model calls.
`assembly_supported_classes` contains the supported measured intervals, and
is exactly equal to `assembly_only_supported_classes`. Reference-only calls
remain reference-only and never create an assembled sequence. The strict
assembly tier and all unresolved/partial candidates remain visible.

Length class and exact sequence comparisons are separate. The H5
4897–4899-base alias applies only when that alias is present in the supplied
panel. `new_sequence_candidate` means a supported inner sequence differs from
every sequence in the fixed panel; it is not a global novelty claim. Exact
full-sequence comparison is unknown when an outer anchor is incomplete.

`discordance` reports differing method evidence sets. It does not establish
biological conflict, strain replacement, gene loss or functional change.
Base-concordance counts are preserved without editing the assembly sequence.

For fair evaluation, report ordinary assembly geometry, same-gated assembly
readback and the reference model separately. An assembly-assisted module
inherits structures already found by assembly; those are not gains over
assembly. Reusing the assembly reads for readback establishes consistency
within the observed library, not independent validation.

## Existing-data verification

The adapter verifies the frozen v2 COMPLETE/STATE/RESULT chain, same-case
contig and paired-read bindings, readback table/FASTA hashes, geometry,
sequence orientation and each source audit. It never opens truth sequences,
component labels, diet or phenotype records. All 82 cases retain their original
denominators; missing anchors are unresolved rather than absent.

The core has 17 targeted tests. The 82-case primary tier exactly replays the
historical read-supported candidates. Results remain exposed development
evidence; neither these tests nor the case replay establish superiority over
ordinary assembly.
