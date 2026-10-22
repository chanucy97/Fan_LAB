# External input contracts

Version: public-adapter.v1. Paths below are logical input names. Supply local
paths when constructing manifests; no original server path is distributed.
Most record bindings are `{path, bytes, sha256}`. A hash identifies content;
it does not establish consent, sample identity or scientific validity.

## Assembly

`assembly_scan.py prepare --config CONFIG.json --out NEW_DIRECTORY` requires
exactly these configuration keys:

```json
{
  "schema": "sequence-integrity.assembly-config.v1",
  "inventory": "USER_ABSOLUTE_PATH_TO_INVENTORY",
  "dictionary": "USER_ABSOLUTE_PATH_TO_DICTIONARY_DIRECTORY",
  "mapping_source": "USER_ABSOLUTE_PATH_TO_PUBLIC_map_traits.py",
  "mapping_rules": "USER_ABSOLUTE_PATH_TO_RULES.json",
  "plan": "USER_ABSOLUTE_PATH_TO_REVIEWED_PLAN",
  "tool_receipts": ["USER_ABSOLUTE_PATH_TO_TOOL_RECEIPT"],
  "tools": {
    "prodigal": "USER_ABSOLUTE_PATH_TO_PRODIGAL",
    "diamond": "USER_ABSOLUTE_PATH_TO_DIAMOND",
    "tblastn": "USER_ABSOLUTE_PATH_TO_TBLASTN",
    "makeblastdb": "USER_ABSOLUTE_PATH_TO_MAKEBLASTDB"
  },
  "threads": 4
}
```

This is a **schema illustration**, not a runnable case. Inventory requires a
`historical_report` with `status`, `samples_n` and `samples`, plus
`assembly_sha_reverification_and_tools.assemblies`; sample IDs, assembly paths,
sizes, SHA values and contig/base counts must agree. Dictionary files are
`FROZEN.json`, `TRAITS.tsv`, `trait_proteins.faa`, `trait_cds.fna` and
`MODULE_RULES.tsv`. Their frozen hashes and counts are checked by
`source_inputs()`. The original dictionary and Prodigal binary hash constraints
are retained; this package does not pretend a synthetic dictionary is the study
dictionary. Synthetic tests directly exercise the algorithm with toy records.

## Profile consensus and population comparison

The profile audit uses `--base`, `--sample` and `--out`. Under the supplied base
it expects `reference_panel/`, `trait_dictionary/`, `trait_mapping/prepared_final/`,
`gene_evidence/reference_plan_v2/`, `gene_evidence/INSTALLED_API.json`,
`gene_evidence/results_v2/<sample>/`, and `profiles_v1/<sample>/`. Consult
`load_roster()` and `load_profile()` for exact receipt and table fields. It uses
existing gene evidence; it does not rerun that upstream workflow. Its frozen
reference plan hash is retained and is an external dependency.

Comparison configuration requires `prepared_dir`, `sample_manifest`,
`within_pairs`, `genome_metadata`, `reference_fasta`, `reference_stb`,
`reference_genes`, `profiles_root`, `output_root`, `instrain`, `samtools` and
`threads`. Sample metadata uses `sample_id`, `dyad_id`, `window`, `nominal_day`.
The comparison `metadata()` tests demonstrate small synthetic input objects
with `enforce_cohort=False`; the study `prepare` command retains its full-cohort
assertions. Nominal visit spacing is not promoted to actual sampling dates.

## Enzyme context

Input schema: `enzyme-context.inputs.v1`. Required top-level content includes:

- `cohort_denominator`, `candidate_count`, `contig_count`, `mask_orf_count`;
- `dependency` record for the bundled `carrier_context.py`, `config_sha256`,
  `reference_fasta` record, `references`, `subjects`, `source_receipts`, `cases`;
- each case: `case_id`, `sample_id`, `contig`, `contig_length`, `contig_sha256`,
  `assembly` record, positive `genes`, and all accepted `mask_orfs` (including
  specificity comparators);
- each gene/mask: `gene_id`, `trait_id`, `trait_role`, `mapping_status`,
  `start0`, `end0`, `strand`, `cds_sha256`, `protein_sha256`, `partial`,
  `complete_predicted_orf`;
- each subject: `subject_id`, `accession`, `length`, `sequence_sha256`.

Intervals are 0-based half-open. The pure-synthetic `ContextTests.manifest()` and
`ContextTests.case()` constructors show the schema. `test_one_contig_one_command`
uses a mocked aligner, not BLAST or participant assemblies. The bundled
`carrier_context/config.json` retains all numerical rules, uses executable names
resolved on PATH, and is itself SHA-bound. Runtime bindings hash the resolved
executables, and commands use those bound absolute executable paths rather than
resolving PATH again. Its scientific configuration is not silently rewritten.

## One pre-specified readback candidate

The CLI requires `--rules`, `--rules-sha256`, `--out`, and the explicit
`--root-authorized-run` flag. The flag authorizes execution of a locally reviewed
configuration; it does not grant access to anyone else's data.

Rules require schema `single-candidate-readback-rules.public-adapter.v1`, status
`FROZEN_FOR_SINGLE_CANDIDATE_READBACK`, `candidate_count: 1`, a nonempty
`sample_id`, `candidate_id`, `implementation_sha256`, `inputs`, `single_read`,
`technical_support`, `assays`, and `parity_windows`. `inputs` binds the basename
and SHA of five separately supplied metadata files: `candidate`, `bindings`,
`filter_sources`, `frozen_candidate_manifest`, and `frozen_selected_cases`.
The selected case must match the candidate and sample across all records.

`assays_from_candidate()` specifies the sequence contract: observed local DNA,
an explicitly labelled in-silico same-flank insertion comparator, two distinct
reference windows, SHA values, junction/gap coordinates and equivalent gap
intervals. It currently requires the verified reference-plus orientation and
four distinguishable templates. This is not arbitrary indel discovery.

Bindings require original profile BAM/BAI, accepted Rdic, covT, attributes,
scaffold lengths, source/runtime records, reference lengths and original split
definitions. Fresh output, SHA/stat checks, primitive-only Rdic decoding,
unchanged filter/runtime and exact coverage parity are enforced before counts.
The real case rules and candidate sequences are intentionally absent. The
readback tests build toy `Read`, `Assay` and `TechnicalRules` objects in memory;
they demonstrate every numerical support boundary without reading real files.
