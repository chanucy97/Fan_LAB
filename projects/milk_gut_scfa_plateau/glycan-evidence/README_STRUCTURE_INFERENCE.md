# Experimental local-structure inference (dev8)

`structure-infer` groups proper paired alignments by fragment, retains near-best candidate alternatives, and records direct versus paired support across configured boundaries. It fits the observed feature distribution to explicit candidate profiles by sparse linear programming, reporting local signal weights and descriptive feasible ranges. The frozen v2 identity engine and existing structure commands retain their behavior.

The input SAM must be generated against exactly the supplied candidate dictionary, contain reciprocal concordant mates, retain secondary alignments, and keep each QNAME contiguous. The runner used for the 2026-10-15 evaluation reports all Bowtie2 alignments; MAPQ zero is retained. A fragment contributes one feature count regardless of the number of candidates it matches. Expected profiles must be generated with the same mapper, filters, candidate set and feature encoding. The core does not silently infer these provenance conditions from an arbitrary SAM.

Install the `inference` optional dependencies. Use:

```text
glycan-evidence structure-infer --sam alignments.sam --candidates candidates.json --profiles profiles.json --config inference_config.json --out new_result_directory
```

`candidates.json` maps candidate IDs to length, zero-based interbase `junctions`, and optional metadata. `profiles.json` maps each identical candidate ID to expected feature counts. `inference_config.json` contains `extraction` and `fit` dictionaries; every output binds the actual inputs. The study runner/configuration in `outputs/20261015_structure_inference` also handles reference-panel freezing, simulation-derived expected profiles, matched baselines and cohort execution.

Outputs describe compatibility with supplied local sequences. Distinct compatible candidates do not by themselves establish a mixture. Weight ranges at optimal total variation plus a specified slack are descriptive, not confidence intervals or organism abundances. Profile mismatch, read support and boundary evidence must be inspected together. New benchmark and cohort results are reported separately from code-test success; dev8 is a research development implementation.
