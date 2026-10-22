# Reference-guided recruitment assembly

`glycan_evidence.assembly.run_assembly(r1, r2, sample, recruit_index,
reference_fasta, recruitment_manifest, out, megahit, *, threads=12,
bin_dir=None, expected_qc_pairs=None, memory_bytes=24000000000,
min_contig_len=500, validation_tmpdir=None)` runs complete FASTQ validation, frozen v2 fastp QC and
permissive paired recruitment, then MEGAHIT on every recruited pair.

This is **reference-guided recruitment assembly**, not complete metagenome
de novo assembly. Results are restricted by the explicitly supplied recruitment
reference universe. The executable path, input paths and reference prefix are
arguments; no sample, server or locus coordinates are embedded in the module.

The recruitment manifest has schema `glycan-evidence.recruitment-index.v1`,
status `VERIFIED`, a `reference_fasta` size/SHA-256 record, unique `source_genomes`,
six `index_files` records (`suffix`: 1, 2, 3, 4, rev.1, rev.2; bytes and SHA-256),
and an `index_sequence_verification` record with status PASS and
`canonical_sequence_sha256`. The manifest-generation audit compares the actual
Bowtie2-inspected index sequences with the source FASTA; the adapter verifies
the explicitly bound source and all six index contents on every invocation.
The manifest is a declared trusted configuration, not an automatically inferred
claim that a reference set is exhaustive. Symlink index inputs are read and
hashed through to their actual content, without copying the reference database.

The cohort preset freezes the existing 21-genome recruitment panel. External
candidate reference genomes must not be added during the frozen cohort run.
Existing frozen v2 and previous exploratory results remain read-only.

Input validation scans the complete paired FASTQ streams with matched EOF and
exact disk-backed QNAME uniqueness, and records SHA-256 of decompressed bytes.
An explicit `validation_tmpdir` can place the temporary name database on verified
node-local scratch. The directory must already exist; its free space and actual
path are recorded. Without it, validation uses the output staging directory.
QC always runs; there is no skip-QC flag in this adapter. When an original QC
pair count is supplied, it must match exactly. Both reads are retained when
either primary mate maps, preserving the original frozen recruitment semantics.

MEGAHIT uses explicit threads, byte memory limit and minimum contig length.
The cohort comparison uses the prior MEGAHIT 1.2.9 parameters: 12 threads,
24,000,000,000 bytes and 500 bp. The declared Slurm request additionally leaves
memory for paired validation, QC and operating-system buffers.

Output contains `all_contigs.fasta`, `contig_metadata.json`, recruited paired
FASTQs, QC/recruitment manifests and logs, MEGAHIT outputs, a versioned receipt,
and `COMPLETE`. The latter binds the receipt SHA. Metadata records every contig,
including sequence hashes. Zero recruitment is a successful processed result
with an empty FASTA and `unresolved_no_recruited_reads`; no-contig assembly is
`unresolved_no_contigs`. Neither means biological absence. Failed subprocesses,
changed inputs or QC mismatches leave a failed staging receipt and no published
completed directory. Existing paths are never overwritten.

Prior exploratory assemblies may be reused only after explicit raw-source,
index, code, parameter, tool-version and output hash comparison. Such rows are
labelled `legacy_reused_preflight_not_run`; their historical input checks are
not silently upgraded to the new complete FASTQ preflight.
