"""Reference-guided paired recruitment and assembly with explicit provenance."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

from . import __version__
from .identity import (_json, _new_stage, _sample_id, validate_paired_fastq, verify_engine)
from .legacy_v2 import run_sample as frozen

SCHEMA = "glycan-evidence.assembly.v1"
SCOPE = "reference_guided_recruitment_assembly; provided_reference_universe; not whole_metagenome_de_novo"


def stamp():
    return datetime.now(timezone.utc).isoformat()


def record(path, checksum=True):
    return frozen.file_record(path, checksum=checksum)


def write(path, value):
    temporary = Path(path).with_name(Path(path).name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf8")
    temporary.replace(path)


def verify_recruitment_manifest(recruit_index, reference_fasta, recruitment_manifest):
    """Verify an explicit source/index binding; no implicit reference discovery."""
    manifest_path = Path(recruitment_manifest).resolve()
    manifest = _json(manifest_path)
    if manifest.get("schema") != "glycan-evidence.recruitment-index.v1" or manifest.get("status") != "VERIFIED":
        raise ValueError("Unsupported or unverified recruitment index manifest")
    source = Path(reference_fasta).resolve()
    source_record = manifest.get("reference_fasta", {})
    if source.stat().st_size != source_record.get("bytes") or record(source)["sha256"] != source_record.get("sha256"):
        raise ValueError("Recruitment reference FASTA differs from frozen manifest")
    if not manifest.get("source_genomes") or len(manifest["source_genomes"]) != len(set(manifest["source_genomes"])):
        raise ValueError("Missing or duplicate recruitment source genome IDs")
    check = manifest.get("index_sequence_verification", {})
    if check.get("status") != "PASS" or not re.fullmatch(r"[a-f0-9]{64}", check.get("canonical_sequence_sha256", "")):
        raise ValueError("Missing index-versus-reference sequence verification")
    prefix = Path(recruit_index).absolute()
    files = sorted(prefix.parent.glob(prefix.name + ".*.bt2*"))
    expected = manifest.get("index_files", [])
    suffixes = {"1", "2", "3", "4", "rev.1", "rev.2"}
    if len(files) != 6 or len(expected) != 6:
        raise ValueError("Recruitment requires exactly six Bowtie2 index members")
    entries = {}
    for item in expected:
        suffix = item.get("suffix", "")
        if suffix in entries or suffix not in suffixes:
            raise ValueError("Duplicate or invalid recruitment index member")
        entries[suffix] = item
    actual, extensions = [], set()
    for path in files:
        match = re.fullmatch(re.escape(prefix.name) + r"\.(1|2|3|4|rev\.1|rev\.2)\.(bt2|bt2l)", path.name)
        if not match:
            raise ValueError("Unexpected recruitment index member")
        suffix, extension = match.groups()
        extensions.add(extension)
        current = record(path)
        if current["bytes"] != entries[suffix].get("bytes") or current["sha256"] != entries[suffix].get("sha256"):
            raise ValueError(f"Recruitment index content hash mismatch: {path.name}")
        actual.append({**current, "suffix": suffix, "extension": extension})
    if len(extensions) != 1:
        raise ValueError("Mixed Bowtie2 index formats")
    return {"manifest": record(manifest_path), "reference_fasta": record(source),
            "source_genomes": manifest["source_genomes"], "index_sequence_verification": check,
            "index_files": actual}


def fasta_metadata(path):
    metadata, seen = [], set()
    current, description, pieces = None, None, []

    def finish():
        if current is None:
            return
        sequence = "".join(pieces).upper()
        if not sequence or not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", sequence):
            raise ValueError(f"Empty/invalid assembled contig sequence: {current}")
        metadata.append({"contig": current, "length": len(sequence), "description": description,
                         "sequence_sha256": hashlib.sha256(sequence.encode("ascii")).hexdigest()})

    with Path(path).open(encoding="ascii") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                description = line[1:]
                fields = description.split()
                if not fields or fields[0] in seen:
                    raise ValueError("Missing or duplicate assembled contig ID")
                current = fields[0]
                seen.add(current)
                pieces = []
            elif current is None:
                raise ValueError("Contig FASTA sequence occurs before its header")
            else:
                pieces.append(line)
        finish()
    return metadata


def run_assembly(r1, r2, sample, recruit_index, reference_fasta, recruitment_manifest, out, megahit,
                 *, threads=12, bin_dir=None, expected_qc_pairs=None,
                 memory_bytes=24000000000, min_contig_len=500, validation_tmpdir=None):
    """Run complete QC, broad paired recruitment and MEGAHIT on recruited reads.

    No sample IDs, contig IDs, genomic coordinates, server roots or automatically
    chosen references are encoded here. Zero recruitment is a valid unresolved
    result. The requested output must not already exist.
    """
    sample = _sample_id(sample)
    if type(threads) is not int or threads < 1 or type(memory_bytes) is not int or memory_bytes < 1:
        raise ValueError("threads and memory_bytes must be positive integers")
    if type(min_contig_len) is not int or min_contig_len < 1:
        raise ValueError("min_contig_len must be a positive integer")
    if expected_qc_pairs is not None and (type(expected_qc_pairs) is not int or expected_qc_pairs < 0):
        raise ValueError("expected_qc_pairs must be a nonnegative integer")
    engine_hashes = verify_engine()
    index = verify_recruitment_manifest(recruit_index, reference_fasta, recruitment_manifest)
    bowtie = frozen.command_prefix("bowtie2", bin_dir)
    fastp = frozen.command_prefix("fastp", bin_dir)
    assembler = Path(megahit).resolve()
    if not assembler.is_file():
        raise FileNotFoundError(f"MEGAHIT launcher not found: {assembler}")
    destination, stage = _new_stage(out)
    receipt = {"schema": SCHEMA, "package_version": __version__, "sample": sample,
               "status": "RUNNING", "scope": SCOPE, "coverage_unit": "recruited_paired_fragments; contig_nucleotide_bases",
               "started_utc": stamp(), "execution_directory": str(stage), "published_directory": str(destination),
               "commands": [], "recruitment_index": index, "legacy_engine_sha256": engine_hashes,
               "adapter_source": record(__file__), "megahit_launcher": record(assembler),
               "parameters": {"threads": threads, "memory_bytes": memory_bytes, "min_contig_len": min_contig_len,
                              "qc": "frozen_v2_complete_fastp_stream", "recruitment": "frozen_v2_permissive_paired_screen",
                              "expected_qc_pairs": expected_qc_pairs}}
    started = time.monotonic()
    path = stage / "receipt.json"
    write(path, receipt)
    try:
        validation_dir = Path(validation_tmpdir).resolve() if validation_tmpdir else stage
        if not validation_dir.is_dir():
            raise ValueError("Explicit validation_tmpdir must be an existing directory")
        receipt["validation_temporary_directory"] = str(validation_dir)
        receipt["validation_disk_free_bytes_at_start"] = shutil.disk_usage(validation_dir).free
        receipt["raw_input_validation"] = validate_paired_fastq(r1, r2, validation_dir)
        write(path, receipt)
        panel = stage / "recruitment_panel"
        panel.mkdir()
        for member in index["index_files"]:
            (panel / ("features." + member["suffix"] + "." + member["extension"])).symlink_to(member["path"])
        state = frozen.RunState(stage / "recruitment", sample, False, True)
        state.data.update(analysis_scope=SCOPE, raw_inputs=receipt["raw_input_validation"]["input_stats"],
                          code=[record(Path(frozen.__file__))], recruitment_index_manifest=index["manifest"])
        frozen.capture_version(state, "bowtie2", bowtie)
        frozen.capture_version(state, "fastp", fastp)
        qc_started = time.monotonic()
        state.data["screen"] = frozen.screen(state, Path(r1).resolve(), Path(r2).resolve(), panel, threads, bowtie, fastp)
        receipt["qc_recruitment_elapsed_seconds"] = time.monotonic() - qc_started
        receipt["screen"] = state.data["screen"]
        if expected_qc_pairs is not None and receipt["screen"]["total_fragments"] != expected_qc_pairs:
            raise ValueError("Repeated complete QC count differs from explicitly expected count")
        state.data.update(status="COMPLETE", completed_utc=stamp())
        state.save()
        receipt["recruitment_manifest"] = record(state.out / "run_manifest.json")
        receipt["tool_versions"] = state.data["versions"]
        reads = [state.out / f"candidate_R{mate}.fastq.gz" for mate in (1, 2)]
        receipt["assembly_input"] = [record(read) for read in reads]
        contigs = stage / "all_contigs.fasta"

        def command(argv, logname):
            item = {"argv": list(map(str, argv)), "started_utc": stamp(), "returncode": None}
            receipt["commands"].append(item)
            write(path, receipt)
            with (stage / logname).open("x", encoding="utf8") as log:
                result = subprocess.run(item["argv"], stdout=log, stderr=subprocess.STDOUT)
            item.update(returncode=result.returncode, completed_utc=stamp())
            write(path, receipt)
            result.check_returncode()

        command([assembler, "--version"], "megahit_version.txt")
        receipt["tool_versions"]["megahit"] = (stage / "megahit_version.txt").read_text().strip()
        if receipt["screen"]["candidate_fragments"]:
            command([assembler, "-1", reads[0], "-2", reads[1], "-o", stage / "megahit",
                     "-t", str(threads), "-m", str(memory_bytes), "--min-contig-len", str(min_contig_len)], "assembly.log")
            source = stage / "megahit/final.contigs.fa"
            if not source.is_file():
                raise ValueError("MEGAHIT completed without final.contigs.fa")
            shutil.copyfile(source, contigs)
        else:
            contigs.write_bytes(b"")
        metadata = fasta_metadata(contigs)
        write(stage / "contig_metadata.json", metadata)
        lengths = sorted((row["length"] for row in metadata), reverse=True)
        total, cumulative, n50 = sum(lengths), 0, 0
        for length in lengths:
            cumulative += length
            if cumulative >= total / 2:
                n50 = length
                break
        receipt["assembly_summary"] = {"contigs": len(metadata), "bases": total,
                                       "max_length": max(lengths, default=0), "N50": n50}
        receipt["assembly_state"] = "assembled" if metadata else "unresolved_no_recruited_reads" if not receipt["screen"]["candidate_fragments"] else "unresolved_no_contigs"
        if receipt["raw_input_validation"]["input_stats"] != [record(p, False) for p in (r1, r2)]:
            raise ValueError("Raw input stat changed during recruitment/assembly")
        receipt["outputs"] = [record(contigs), record(stage / "contig_metadata.json")]
        receipt.update(status="COMPLETE", completed_utc=stamp(), elapsed_seconds=time.monotonic() - started)
        write(path, receipt)
        (stage / "COMPLETE").write_text(record(path)["sha256"] + "\n", encoding="ascii")
        stage.rename(destination)
        return receipt
    except BaseException as error:
        receipt.update(status="FAILED", failed_utc=stamp(), elapsed_seconds=time.monotonic() - started,
                       error_type=type(error).__name__, error=str(error))
        write(path, receipt)
        marker = stage / "COMPLETE"
        if marker.exists():
            marker.unlink()
        raise
