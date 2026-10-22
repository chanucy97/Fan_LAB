# Input contract — public-adapter.v1

All paths are explicit local CLI arguments. Inputs are UTF-8 JSON unless noted. Cases and candidates must have stable unique identifiers. Supply evaluation truth only to the evaluator, never to the inference/calibration process. `examples/make_synthetic_inputs.py` creates a fully artificial complete v1 example.

## Common panel and gates

`panel_metadata.json` is a dictionary keyed by candidate ID. Each value has integer `length`, `inner_start` (or `inner_start0`), `inner_end` (or `inner_end0`) and a two-element `junctions` list equal to the inner boundaries. Coordinates are zero-based, end-exclusive, with `0 < inner_start < inner_end < length`. `class_id` is the runtime group alias. The evaluator derives the fixed comparison class from the inner length and requires aliases/groups to agree with full candidate membership.

`gate_v1.json`/`gate_v2.json` can be direct gate dictionaries or nested under `evaluation_gate`. The checked-in examples preserve the source numerical gates. They are method parameters, not a published private candidate roster. Do not tune them after examining test outcomes and describe the result as a frozen evaluation.

## V1 inventory and local history

Inventory fields:

- `case_count`: number of all cases, including no-call/negative cases.
- `positive_alleles`: number of positive interval truth alleles, not unique classes.
- `cases`: list of objects with `case_id`, `fixture_set`, `components` and `error_rate`; each component supplies `positive` and `coverage`.

The retained `fixture_set` value `represented_stress24` selects the source format containing separate full/inner truth FASTA. It is a compatibility label, not a requirement for 24 cases. Other labels select the single-source truth-structure-table format. This transport adapter does not create or select a study fixture set.

History is `{"cases": [{"case_id": "...", "files": {...}}]}`. Each file entry contains `exists`, `text` and `sha256` of the UTF-8 text bytes. Optional `path` and `bytes` are descriptive provenance. Required keys:

| Key | Content |
| --- | --- |
| `truth_full` | FASTA of source-truth anchored intervals |
| `truth_inner` | FASTA with exactly the same truth IDs, for `represented_stress24` |
| `truth_structure_table` | One-row TSV for the other truth format |
| `structure_table` | TSV of assembly-derived anchored intervals; may be absent only with an explicit no-candidate receipt |
| `anchored_fasta` | Assembly-derived interval FASTA; IDs must match structure table |
| `task_receipt` | JSON text with `candidate_count=0` and `geometry_state` beginning `NOT_APPLICABLE` when geometry is absent |

Each structure row contains `structure_id`, `interval_start_0based`, `interval_end_0based_exclusive`, `left_anchor_boundary_0based`, `right_anchor_boundary_0based`, `orientation_to_reference` (`+`/`-`), `interval_sequence_sha256`, `inner_anchor_edges_complete` and `both_anchors_full_length`. Optional `between_anchor_bp` is checked against the extracted inner length. Interval FASTA is already oriented as recorded by the original output. For `-`, inner slicing uses `end-right : end-left`; full/inner exact comparisons canonicalize reverse complements. Do not pass un-oriented whole-contig sequence as an interval.

## V1 result per case

`results-root/<case_id>/RESULT.json` contains:

- `case_id`.
- `extraction.audit.accepted_fragments`; `ledger_rows_omitted_by_cap` must be zero.
- `extraction.config.ledger_mode` must retain the ledger, not `none`.
- `extraction.fragment_ledger`: complete rows with unique `qname`, `candidates` and `candidate_evidence`. Each candidate evidence has `direct_junctions` and `bridged_junctions`. Its boundary positions must be in that candidate's `junctions`. Candidate evidence must already intersect that candidate's compatible placements. The evaluator then intersects across compatible candidates.
- `extraction.candidate_coverage`: all panel IDs, each with `length`, `eligible_base_breadth`, `near_best_base_breadth`, `eligible_fragments` and `near_best_fragments`. Near-best evidence must not exceed eligible evidence.
- `fit`, `paired_only_fit`: original structure-inference fit dictionaries. Essential fields are `status`, `accepted_fragment_count`, `bestfit_tv`, `bestfit_weights`, `config.candidate_groups`, and jointly optimized `group_ranges`. Fit status is `fitted_descriptive` or `no_accepted_evidence`. A full group's range must be present; summing separate candidate bounds is invalid.
- `decision`, `paired_only_decision`, `reference_baseline`: complete outputs from `analysis_support.summarize`. The independent evaluator recomputes and cross-checks them, including every baseline.

## V2 additions

The v1 extraction additionally provides `counts`: a dictionary from serialized feature JSON to integer QNAME counts. A feature is a JSON list of `[candidate, direct_junctions, bridged_junctions]` entries (or `[candidate]` for the boundary-free form). `__unexplained__` is recorded separately. These counts and the ledger must describe the same accepted fragments and support.

`profiles.json` maps every candidate to feature-count dictionaries from the fixed calibration. Profiles are normalized by the evaluator; unexplained entries are retained separately, and profile-only compatibility is tracked rather than converted into observed boundary support.

`v2-results-root/<case_id>/RESULT_V2.json` contains `case_id`, `assignment_fit`, `decision_v2`, `reference_baseline_v2` and `no_profile_decisions`. `assignment_fit` is the untrimmed result of `glycan_evidence.structure_assignment.fit_assignment`, including configuration, observed support, two no-profile models, group/boundary bounds, total-variation values and solver validation. Do not replace this dictionary with point estimates. Runtime call dictionaries are produced by `assignment_support.summarize_assignment`.

`--v1-evaluation` points to a complete immutable v1 `EVALUATION.json`. The caller must preserve the same case IDs, panel and observations. Missing v1/v2 files fail instead of reducing the case count.

## Outputs

V1 writes `EVALUATION.json`, `HISTORY.json`, `METHOD_CASES.tsv`, `CLASS_EVIDENCE.tsv`, `METHOD_SUMMARY.tsv` and `SOURCE_MANIFEST.json`. V2 writes `EVALUATION_V2.json`, `METHOD_CASES_V2.tsv`, `CLASS_EVIDENCE_V2.tsv`, `METHOD_SUMMARY_V2.tsv` and `SOURCE_MANIFEST.json`.

`recovered_classes`, `missed_classes`, `additional_false_classes`, ambiguous calls and no calls remain distinct. `assembly_exact_full`/`assembly_exact_inner` count truth allele IDs, whereas length-class set recovery can collapse multiple same-class alleles. Source manifests record hashes of caller-supplied local files; users should review their generated manifests before publishing paths or data.
