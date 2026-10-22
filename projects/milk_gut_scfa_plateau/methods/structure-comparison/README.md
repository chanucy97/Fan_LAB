# Structure comparison and independent evaluation

This `public-adapter.v1` package evaluates reference-architecture calls, expanded-reference baselines and ordinary assembly on explicit, separate endpoints. It includes the original support summarizers, two independent evaluators and 21 original synthetic tests. It does not perform read mapping, assembly, model fitting or sample selection.

The retained historical comparison recovered **23 of 35 truth-class occurrences**, the same count as ordinary assembly geometry. That result does not establish superiority over ordinary assembly. It is historical context, not a result produced by the synthetic example or rerun during this public adaptation. Exact full-sequence recovery, exact internal-sequence recovery and architecture-class recovery are different endpoints.

## Files and scope

| File | Purpose |
| --- | --- |
| `analysis_support.py` | Runtime descriptive support summaries and breadth/matched-support baselines |
| `evaluate_benchmark.py` | Independently rebuild distinct-QNAME boundary evidence, cross-check runtime calls, and score reference calls and ordinary assembly geometry/exact sequences |
| `assignment_support.py` | Runtime summaries of coupled-assignment bounds and support-only comparators |
| `evaluate_assignment.py` | Independently rebuild observed masks and no-profile bounds, check allocation units and score v2 against immutable v1 observations |
| `test_evaluate_*.py` | Original 13 benchmark and 8 assignment synthetic checks; assignment tests import the sibling public runtime |
| `examples/make_synthetic_inputs.py` | Generate one entirely artificial positive case for a complete v1 CLI smoke run |
| `INPUT_SCHEMA.md` | Input fields, coordinates, denominators and output meanings |
| `PROVENANCE.json` | Workspace-relative source paths, original/public SHA-256, changes and explicitly omitted historical runners |
| `VALIDATION.json` | Actual public-adapter validation results and untested scope |

No cohort table, cohort FASTA, private case inventory, candidate selection file, archived cohort result, credential or server location is included. The example uses invented identifiers and constructed strings; it contains no participant or reference-genome sequence. The actual historical benchmark artifacts are not included, so this directory alone does not reproduce the historical 23/35 comparison.

## Requirements and quick check

Use Python 3.11 or later. Evaluators and the v1 example use only the Python standard library. The assignment tests additionally call `glycan_evidence.structure_assignment`, which requires the sibling `glycan-evidence` package and its `inference` dependencies. From this directory:

```sh
python -m pip install -e "../../glycan-evidence[inference]"
python -m unittest discover -s . -p "test_evaluate_*.py" -v
python evaluate_benchmark.py --help
python evaluate_assignment.py --help
python examples/make_synthetic_inputs.py --output demo_inputs
python evaluate_benchmark.py --inventory demo_inputs/inventory.json --panel-metadata demo_inputs/panel_metadata.json --gate-config examples/gate_v1.json --results-root demo_inputs/results --history-json demo_inputs/history.json --output demo_evaluation
```

`requirements.txt` explains the dependency split. The original runtime's dependency pins remain in its `pyproject.toml`; this adapter does not silently change them. See the validation record for the versions actually exercised. Both evaluator output directories must be new. Example generation also refuses an existing output directory.

The artificial example is a format/integration smoke test. Its preconstructed fit is not an optimization benchmark, realistic read simulation, independent external validation or performance estimate. Both synthetic truth and synthetic prediction are deliberately supplied to exercise serialization and cross-checks.

## Evaluate caller-supplied results

The first CLI requires a local `--history-json` archive with exact text and SHA-256 bindings. This public adapter removes the original SSH fetch helper; it never fetches truth or remote files. Full accepted fragment ledgers are required. Missing cases, duplicate ledger QNAMEs, truncated ledgers, incomplete groups and inconsistent runtime summaries fail evaluation rather than shrinking the denominator.

After producing a complete v1 evaluation, run the coupled-assignment evaluator with caller-owned files:

```sh
python evaluate_assignment.py --inventory inventory.json --panel-metadata panel_metadata.json --profiles profiles.json --gate-config examples/gate_v2.json --v1-evaluation evaluation_v1/EVALUATION.json --v1-results-root v1_cases --v2-results-root v2_cases --output evaluation_v2
```

Panel membership, calibration profiles, truth archives and runtime `RESULT.json`/`RESULT_V2.json` must be prepared separately. Neither evaluator supplies missing evidence. The two support modules expose library functions rather than separate command-line entry points.

## Endpoints and interpretation

- Distinct QNAMEs are counted once. Every compatible candidate must support the same boundary before evidence is class-exclusive. Read names are opaque; origin labels are never parsed.
- Breadth is the maximum single-member breadth, not coverage stitched across different references. Six breadth-only baselines and six matched-support baselines remain visible at 0.50, 0.80 and 0.95.
- Fixed inner-length classes are retained from the source: 4309, 4897–4899, 824 and otherwise exact other lengths. They are comparison categories, not species, strains, functions or globally novel structures. Unsupported lengths are not reassigned to the nearest category.
- Model weight ranges are descriptive ranges under a fixed total-variation budget, not confidence intervals. V2 allocation lower bounds are fractional model allocations, not observed exclusive-fragment counts.
- Abstention does not count as new recovery. An empty call in a negative control does not establish biological absence.
- Assembly sequence comparisons use the archived oriented intervals and their SHA-256. Partial outer anchors can still leave a complete inner interval. Reference-template calls are never credited as reconstructed sequences.
- `junction_profile_ablation` removes junction categories from fitted profiles while retaining the boundary evidence gate. It is not a completely junction-free analysis.

## Adaptation boundary

Evaluation arithmetic and the original tests were retained. The public changes remove SSH/history transport, replace a workspace-relative test import with the installed runtime import, normalize line endings and reject duplicate local history IDs. The configuration examples expose the evaluation gates without copying the original candidate-membership dictionary. `PROVENANCE.json` records exact hashes.

The fixed 82-case recovery driver, positional readback audit, one-case source-chain check and historical simulation/baseline/verification runners are not included. They depend on archived study-specific rosters, receipts, paths and external tools. Their omission is intentional and is not a claim that the entire historical cohort workflow has been made standalone. Their source hashes and reasons are listed in provenance. The ordinary-assembly geometry and exact-sequence endpoint functions already present in `evaluate_benchmark.py` remain available without replacing them with a new recovery algorithm.
