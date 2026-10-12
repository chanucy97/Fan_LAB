# Coding-sequence replay

`cds-replay` makes the executed Step6 coding checks available through the installed package. It reads saved alignment evidence and candidate sequences, verifies their provenance, and reproduces CDS completeness and exact nucleotide/protein comparisons. All libraries in the declared input ledger remain represented, including libraries with no candidate. Coding status is reported separately from paired-anchor geometry.

```bash
glycan-evidence cds-replay \
  --input-manifest /path/to/input/input_manifest.json \
  --expected-manifest-sha256 ACTUAL_VERIFIED_MANIFEST_SHA256 \
  --out /new/path/cds_result

glycan-evidence verify --directory /new/path/cds_result
```

Use the actual manifest hash from the independently supplied release/fixture record. Input paths are relative to the manifest directory so that the bound package can move between machines. A changed input or a mismatched manifest is rejected. Recorded historical execution paths are provenance text; the command does not execute them.

The input schema is `glycan-evidence.cds-replay-input/1` with the `h5-cds-step6-v3` protocol. It preserves the fixed reference/query scope and accepted HSP-selection and coding rules. A different gene set or protocol requires an explicit new contract rather than editing a reference file while reusing this version. The fixture preparation and contract are provided with the release acceptance material.

The output includes candidate and all-library summary tables, localized nucleotide sequences (`.fna`) and predicted protein sequences (`.faa`), plus execution/provenance records. The top-level CLI seals the report and sequence files together. `verify` detects changes to either the tables or the sequence files. An existing output directory is not reused for a corrected attempt.

The command uses the Python standard library and does not require BLAST or Biopython at replay time. Upstream alignments must already exist. Replaying a completed analysis establishes reproducibility of its coding checks; it adds no new metagenomic reads, independent truth data, biological replicates or functional measurements. Scientific results and their denominators remain those of the bound input package.
