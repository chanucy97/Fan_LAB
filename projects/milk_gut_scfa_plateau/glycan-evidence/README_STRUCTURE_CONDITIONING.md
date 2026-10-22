# Internal-locus conditioning research API

`glycan_evidence.structure_conditioning.extract_locus_conditioned` selects local sequencing fragments before the existing assignment model. This is a development API with a research runner; there is no public conditioning CLI command.

The observation rule retains a QNAME if any original near-best paired placement has either a continuous aligned block overlapping the internal interval by at least `min_anchor` bases, or original direct/paired-bridge support for an inner boundary. Paired eligibility, near-best score range, mapping identity and anchor requirements come from `structure_inference` unchanged. A qualifying placement can retain a fragment even if another compatible placement has only anchor evidence.

## API and input contract

```python
from glycan_evidence.structure_conditioning import extract_locus_conditioned
from glycan_evidence.structure_assignment import fit_assignment

config = {
    "mode": "paired_junction", "near_best_delta": 6,
    "min_aligned_fraction": 0.9, "min_map_identity": 0.9,
    "min_anchor": 15, "max_insert": 1000,
    "ledger_mode": "none", "max_ledger": None,
}
# metadata[candidate_id] = {
#     "length": interval_length,
#     "junctions": [left_inner_boundary, right_inner_boundary],
#     "class_id": architecture_class,
# }
query = extract_locus_conditioned(
    "query.sam", metadata, config,
    expected_full_counts=historical_query_counts,
)

# Apply EXACTLY the same selector and configuration to EVERY calibration SAM.
selected_profiles = {}
for candidate, sam_path in calibration_sams.items():
    calibration = extract_locus_conditioned(
        sam_path, metadata, config,
        expected_full_counts=historical_calibration_counts[candidate],
    )
    if not calibration["counts"]:
        raise ValueError(f"Empty selected profile: {candidate}")
    selected_profiles[candidate] = calibration["counts"]
assert set(selected_profiles) == set(metadata)

fit = fit_assignment(
    query["counts"], selected_profiles, metadata,
    {"tv_slack": 0.02, "numerical_tolerance": 1e-9},
)
```

Input SAM must be QNAME-collated, with exact candidate IDs/lengths in its `@SQ` header and exactly two ordered inner boundaries per candidate. Repeated noncontiguous QNAMEs fail validation. A QNAME is a sequencing-fragment unit, not a deduplicated molecule. The existing paired extractor retains eligible secondary and MAPQ-zero alignments.

The complete original feature is preserved for each selected QNAME: all original near-best candidates remain, including nonqualifying compatible candidates, and each candidate retains its original across-placement boundary intersection. Selection cannot turn possible boundary evidence into exclusive evidence.

## Output and denominator

The JSON-safe return contains selected `counts`, excluded `anchor_only_counts`, reconstructed `full_counts`, separate `unexplained_fragments`, selected and original audits, selected/original coverage, and optional fragment ledger. These identities hold:

```text
counts + anchor_only_counts = original accepted feature counts
selected N + excluded N + unexplained N = original QNAME groups
full_counts = expected_full_counts, when provided
```

`anchor_only_counts` is an operational label for fragments failing the selection rule. It includes internal overlaps shorter than `min_anchor`; these are separately counted as `excluded_subthreshold_internal_overlap`. It does not mean every excluded alignment has zero internal bases. Unexplained QNAMEs are audited separately and never silently absorbed by a model background component.

Both profiles and query counts are conditional on this selection event. Assignment weights and their bounds describe selected local fragments; they are not directly comparable with weights fitted on the original denominator or with organism abundance. Empty selected queries remain unresolved. An empty selected calibration profile is an explicit invalid calibration, not grounds for silently dropping a candidate or reusing its original profile.

Coverage includes the original eligible/near-best whole-interval breadth, selected whole-interval breadth, and selected `eligible_inner_breadth` / `near_best_inner_breadth`. Inner breadth is the union of aligned reference blocks clipped to `[left,right)`, divided by the internal length. For a selected QNAME, all its original placements contribute to the corresponding coverage. A pure-flank pair can provide bridge evidence with no sequenced internal bases; bridge support and internal breadth remain separate measures.

## Frozen research replay

The actual runner and protocol are in `outputs/20261016_structure_recovery/locus_conditioning/`. It reuses the 15 original calibration SAMs, all 32 existing benchmark SAMs and all 50 cohort SAMs. It receives no benchmark truth or cohort assembly sequence and performs no remapping. On Wuxi, from the staged run directory:

```bash
sbatch --parsable run.sbatch calibrate
# Only after calibration COMPLETE and Slurm COMPLETED / 0:0:
sbatch --parsable --array=0-81%4 run.sbatch query
```

The checked-in `stage.py` and `submit_queries.py` implement these steps with one-time intent/submission locks; they deliberately refuse automatic reruns into existing directories. The study-bound runner requires Slurm, original source paths, frozen hashes and sealed historical result receipts. It is not a general input-discovery CLI.

For every query, the runner verifies source v1/v2 completion receipts, result hashes and the recorded SAM hash; it reproduces historical full feature counts and original coverage exactly. For calibration, it binds the current SAM bytes to this replay and verifies original profile hashes and full-count reproduction; the original calibration receipt did not contain historical SAM byte hashes. New profile counts, FASTA/metadata hashes, core hashes, config, results and completion markers are bound together. Failed runs revoke their COMPLETE marker.

The frozen replay completed all 82 cases with 57 extractor/assignment/conditioning tests passing on the compute environment. Selected and excluded counts conserve all source QNAMEs; all 656 class-level observed boundary-support dictionaries remain unchanged. The original numeric decision gates were retained for development comparison and have not been recalibrated for the new observation distribution. The independent evaluation tables report the full case set and fixed breadth grid; this retrospective revision is not a new independent accuracy benchmark.

Independent scoring found three additional cohort library calls, but the sole additional benchmark call was incompatible with its exact source architecture: a 4895 bp source was assigned to the 4897–4899 bp class. Correct benchmark class calls remained 19, with incompatible calls increasing from zero to one. The original v2 plus assembly result remains the primary analysis (23 correct, zero extra); conditioning plus assembly yielded 23 correct and one extra. This API is retained as an optional diagnostic, not promoted to the default classifier.

Authoritative artifacts: `PROTOCOL.md`, `FROZEN.json`, `COLLECTION.json`, `COMPARISON_VERIFICATION.json`, `ALL82_COMPARISON.tsv`, `ALL656_CLASS_COMPARISON.tsv`, `STRONG_BASELINES.tsv`, and `REPORT_ZH.md` in the run directory. The original v1/v2 core and result files are unchanged.
