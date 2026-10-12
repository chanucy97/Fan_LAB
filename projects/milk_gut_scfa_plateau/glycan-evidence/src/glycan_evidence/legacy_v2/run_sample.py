"""Auditable paired-read screen and competitive genome mapping for one sample.

Normal runs stream fastp's complete interleaved QC reads into a permissive paired
feature screen. --skip-qc bypasses only fastp (for simulated input). Candidate
selection retains both full primary SAM sequences whenever either mate maps.
Competitive mapping and evidence quantification are identical in both modes.
Candidate FASTQs are retained by default; --no-keep-candidates removes only the
two candidate files created by this invocation after successful quantification.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVNacgtrymkswbdhvn", "TGCAYRKMSWVHDBNtgcayrkmswvhdbn")
INDEX_PARTS = {"1", "2", "3", "4", "rev.1", "rev.2"}
INDEX_NAME = re.compile(r"^(features|genomes)\.(1|2|3|4|rev\.1|rev\.2)\.(bt2|bt2l)$")


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def file_record(path, checksum=False):
    path = Path(path).resolve()
    stat = path.stat()
    record = {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if checksum:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        record["sha256"] = digest.hexdigest()
    return record


def index_records(panel):
    """Reuse build hashes only for an exact, unchanged 12-file index set.

    index_sha256.json schema: {"files": [{"name": "features.1.bt2",
    "bytes": integer, "mtime_ns": integer, "sha256": 64 hex chars}, ...]}.
    The build stage owns this cache. A malformed/stale cache is an error, not a
    reason to hash a different index silently. A missing cache uses direct SHA256.
    """
    panel = Path(panel)
    paths = sorted(path for prefix in ("features", "genomes") for path in panel.glob(prefix + ".*.bt2*"))
    parsed = {}
    for path in paths:
        match = INDEX_NAME.fullmatch(path.name)
        if not match or not path.is_file():
            raise ValueError(f"Unexpected Bowtie2 index member: {path.name}")
        parsed[path.name] = match.groups()
    if len(parsed) != 12:
        raise ValueError("Expected exactly 12 Bowtie2 index files: six features and six genomes members")
    for prefix in ("features", "genomes"):
        members = [value for value in parsed.values() if value[0] == prefix]
        if {value[1] for value in members} != INDEX_PARTS or len({value[2] for value in members}) != 1:
            raise ValueError(f"Incomplete or mixed-format {prefix} index")
    cache_path = panel / "index_sha256.json"
    if not cache_path.exists():
        rows = [dict(file_record(path, checksum=True), name=path.name, checksum_source="direct_sha256_no_build_cache") for path in paths]
        return rows, {"checksum_source": "direct_sha256_no_build_cache"}
    cache = json.loads(cache_path.read_text(encoding="utf-8-sig"))
    entries = cache.get("files") if isinstance(cache, dict) else None
    if not isinstance(entries, list) or len(entries) != 12:
        raise ValueError("Build index checksum cache must contain exactly 12 file entries")
    cached = {}
    for entry in entries:
        if not isinstance(entry, dict) or set(("name", "bytes", "mtime_ns", "sha256")) - entry.keys():
            raise ValueError("Invalid build index checksum cache entry")
        name = entry["name"]
        if not isinstance(name, str) or name in cached or name not in parsed:
            raise ValueError("Duplicate or unexpected name in build index checksum cache")
        if type(entry["bytes"]) is not int or type(entry["mtime_ns"]) is not int or entry["bytes"] < 1:
            raise ValueError("Invalid size/timestamp in build index checksum cache")
        if not isinstance(entry["sha256"], str) or not re.fullmatch(r"[0-9a-fA-F]{64}", entry["sha256"]):
            raise ValueError("Invalid SHA256 in build index checksum cache")
        cached[name] = entry
    if set(cached) != set(parsed):
        raise ValueError("Build checksum cache names do not match the complete current index set")
    rows = []
    for path in paths:
        current, entry = file_record(path), cached[path.name]
        if current["bytes"] != entry["bytes"] or current["mtime_ns"] != entry["mtime_ns"]:
            raise ValueError(f"Build checksum cache is stale for {path.name}; rebuild/verify the cache explicitly")
        rows.append(dict(current, name=path.name, sha256=entry["sha256"].lower(), checksum_source="build_cache_exact_set_stat_verified"))
    return rows, {"checksum_source": "build_cache_exact_set_stat_verified", "cache": file_record(cache_path, checksum=True)}


def command_prefix(tool, bin_dir=None):
    if bin_dir:
        base = Path(bin_dir) / tool
        candidates = [base, base.with_suffix(".exe")]
        path = next((p for p in candidates if p.is_file()), None)
        if path is None:
            raise FileNotFoundError(f"Tool not found in --bin-dir: {tool}")
        return [str(path.resolve())]
    path = shutil.which(tool)
    if not path:
        raise FileNotFoundError(f"Required tool not found: {tool}")
    return [path]


class RunState:
    def __init__(self, out, sample, skip_qc, keep_candidates):
        self.out = Path(out).resolve()
        self.out.mkdir(parents=True, exist_ok=True)
        protected = ("COMPLETE", "run_manifest.json", "candidate_R1.fastq.gz", "candidate_R2.fastq.gz", "competitive.sam.gz", "gene_metrics.tsv")
        if any((self.out / name).exists() for name in protected):
            raise FileExistsError("Output already contains pipeline artifacts; use a new output directory")
        self.data = {"sample": sample, "status": "RUNNING", "started_utc": timestamp(),
                     "skip_qc": skip_qc, "keep_candidates": keep_candidates, "commands": [], "versions": {"python": sys.version},
                     "raw_hash_policy": "stat only; large raw FASTQ files are not hashed", "screen": None}
        self.save()

    def save(self):
        pending = self.out / "run_manifest.json.tmp"
        pending.write_text(json.dumps(self.data, indent=2) + "\n", encoding="utf-8")
        pending.replace(self.out / "run_manifest.json")

    def command(self, stage, command):
        item = {"stage": stage, "argv": list(map(str, command)), "started_utc": timestamp(), "returncode": None}
        self.data["commands"].append(item)
        self.save()
        return item

    def finish_command(self, item, returncode):
        item.update(returncode=returncode, completed_utc=timestamp())
        self.save()


def check_returncode(state, item, process):
    code = process.wait()
    state.finish_command(item, code)
    if code:
        raise RuntimeError(f"{item['stage']} exited with status {code}; inspect its log")


def stop_process(process):
    if process.poll() is None:
        process.terminate()
        try:
            process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            process.kill()
            process.wait(timeout=10)


def capture_version(state, tool, prefix):
    command = prefix + ["--version"]
    item = state.command(tool + "_version", command)
    result = subprocess.run(command, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=30)
    state.finish_command(item, result.returncode)
    text = (result.stdout + result.stderr).strip()
    (state.out / (tool + "_version.log")).write_text(text + "\n", encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"{tool} --version exited with status {result.returncode}")
    state.data["versions"][tool] = text[:4096]
    state.save()


def grouped_primary_sam(lines):
    name, rows = None, []
    for line in lines:
        if not line.strip() or line.startswith("@"):
            continue
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 11:
            raise ValueError("Screen SAM has fewer than 11 fields")
        if name is not None and fields[0] != name:
            yield name, rows
            rows = []
        name = fields[0]
        rows.append(fields)
    if name is not None:
        yield name, rows


def restore_read(fields):
    flag, sequence, quality = int(fields[1]), fields[9], fields[10]
    if sequence == "*" or quality == "*" or not sequence or len(sequence) != len(quality):
        raise ValueError("Candidate mate lacks a complete sequence/quality; cannot reconstruct paired QC reads")
    if "H" in fields[5]:
        raise ValueError("Hard-clipped candidate cannot preserve the full QC read")
    if fields[5] != "*":
        operations = re.findall(r"(\d+)([MIDNSHP=X])", fields[5])
        if "".join(n + op for n, op in operations) != fields[5]:
            raise ValueError("Candidate has an invalid CIGAR")
        query_length = sum(int(n) for n, op in operations if op in {"M", "I", "S", "=", "X"})
        if query_length != len(sequence):
            raise ValueError("Candidate sequence length disagrees with its full query-consuming CIGAR")
    if flag & 0x10:
        sequence = sequence.translate(COMPLEMENT)[::-1]
        quality = quality[::-1]
    return sequence, quality


def select_candidates(sam_lines, r1_path, r2_path):
    counts = {"total_fragments": 0, "candidate_fragments": 0, "both_mates_mapped": 0,
              "one_mate_mapped": 0, "unmapped_fragments": 0, "ignored_secondary_or_supplementary_records": 0,
              "candidate_R1_bases": 0, "candidate_R2_bases": 0}
    with gzip.open(r1_path, "wt", encoding="ascii", compresslevel=1, newline="\n") as r1, gzip.open(r2_path, "wt", encoding="ascii", compresslevel=1, newline="\n") as r2:
        for name, records in grouped_primary_sam(sam_lines):
            mates = {}
            for fields in records:
                flag = int(fields[1])
                if flag & (0x100 | 0x800):
                    counts["ignored_secondary_or_supplementary_records"] += 1
                    continue
                if not flag & 0x1 or bool(flag & 0x40) == bool(flag & 0x80):
                    raise ValueError("Screen requires exactly two primary paired SAM records per QNAME")
                mate = 1 if flag & 0x40 else 2
                if mate in mates:
                    raise ValueError("Duplicate primary mate in screen SAM; QNAME grouping or input names are invalid")
                mates[mate] = fields
            if set(mates) != {1, 2}:
                raise ValueError("Screen SAM is missing a primary mate; paired reads cannot be silently discarded")
            counts["total_fragments"] += 1
            mapped = sum(not int(fields[1]) & 0x4 for fields in mates.values())
            if not mapped:
                counts["unmapped_fragments"] += 1
                continue
            # Restore both before writing either, including an unmapped partner.
            reads = {mate: restore_read(fields) for mate, fields in mates.items()}
            counts["candidate_fragments"] += 1
            counts["both_mates_mapped" if mapped == 2 else "one_mate_mapped"] += 1
            for mate, handle in ((1, r1), (2, r2)):
                sequence, quality = reads[mate]
                handle.write(f"@{name}/{mate}\n{sequence}\n+\n{quality}\n")
                counts[f"candidate_R{mate}_bases"] += len(sequence)
    return counts


def screen(state, r1, r2, panel, threads, bowtie, fastp):
    common = ["--very-sensitive", "-N", "1", "-L", "20", "--score-min", "L,0,-0.4", "--reorder", "-k", "1", "-p", str(threads), "-x", str(panel / "features")]
    screen_command = bowtie + common + (["--interleaved", "-"] if fastp else ["-1", str(r1), "-2", str(r2)])
    processes = []
    fastp_process = None
    with (state.out / "screen_bowtie2.log").open("w", encoding="utf-8") as btlog, (state.out / "fastp.log").open("w", encoding="utf-8") as qclog:
        try:
            if fastp:
                qc_command = fastp + ["-i", str(r1), "-I", str(r2), "--stdout", "--thread", str(threads),
                                      "--detect_adapter_for_pe", "--qualified_quality_phred", "20", "--unqualified_percent_limit", "40", "--length_required", "75",
                                      "--json", str(state.out / "fastp.json"), "--html", str(state.out / "fastp.html")]
                qc_item = state.command("fastp", qc_command)
                fastp_process = subprocess.Popen(qc_command, stdout=subprocess.PIPE, stderr=qclog)
                processes.append((qc_item, fastp_process))
            item = state.command("feature_screen", screen_command)
            process = subprocess.Popen(screen_command, stdin=fastp_process.stdout if fastp_process else subprocess.DEVNULL,
                                       stdout=subprocess.PIPE, stderr=btlog, text=True, encoding="ascii")
            processes.append((item, process))
            if fastp_process:
                fastp_process.stdout.close()
            counts = select_candidates(process.stdout, state.out / "candidate_R1.fastq.gz", state.out / "candidate_R2.fastq.gz")
            process.stdout.close()
            # Both processes must succeed even when the resulting SAM is valid.
            for process_item, child in reversed(processes):
                check_returncode(state, process_item, child)
            return counts
        finally:
            for process_item, child in processes:
                stop_process(child)
                if child.stdout is not None and not child.stdout.closed:
                    child.stdout.close()
                if process_item["returncode"] is None:
                    state.finish_command(process_item, child.returncode)


def competitive_mapping(state, panel, threads, bowtie, candidates):
    destination = state.out / "competitive.sam.gz"
    if candidates == 0:
        with gzip.open(destination, "wt", encoding="ascii", compresslevel=1) as handle:
            handle.write("@HD\tVN:1.6\tSO:queryname\n")
        state.data["competitive_mapping"] = "skipped: no screen candidates; header-only SAM retained"
        state.save()
        return destination
    command = bowtie + ["-a", "--very-sensitive", "-N", "1", "-L", "20", "--score-min", "L,0,-0.4",
                         "--no-mixed", "--no-discordant", "-I", "0", "-X", "1000", "--no-unal", "-p", str(threads),
                         "-x", str(panel / "genomes"), "-1", str(state.out / "candidate_R1.fastq.gz"), "-2", str(state.out / "candidate_R2.fastq.gz")]
    item = state.command("competitive_genome_mapping", command)
    with (state.out / "competitive_bowtie2.log").open("w", encoding="utf-8") as log:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=log)
        try:
            with gzip.open(destination, "wb", compresslevel=1) as compressed:
                shutil.copyfileobj(process.stdout, compressed, length=1024 * 1024)
            process.stdout.close()
            check_returncode(state, item, process)
        finally:
            stop_process(process)
            if process.stdout is not None and not process.stdout.closed:
                process.stdout.close()
            if item["returncode"] is None:
                state.finish_command(item, process.returncode)
    return destination


def quantify(state, panel, sam, sample):
    script = Path(__file__).with_name("map_evidence.py").resolve()
    command = [sys.executable, str(script), "--sam", str(sam), "--contigs", str(panel / "panel_contigs.tsv"),
               "--features", str(panel / "feature_manifest.tsv"), "--annotations", str(panel / "annotations.tsv"), "--sample", sample, "--out", str(state.out)]
    item = state.command("map_evidence", command)
    with (state.out / "map_evidence.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    state.finish_command(item, result.returncode)
    if result.returncode:
        raise RuntimeError(f"map_evidence exited with status {result.returncode}; inspect map_evidence.log")
    for name in ("gene_metrics.tsv", "summary.json"):
        if not (state.out / name).is_file() or not (state.out / name).stat().st_size:
            raise RuntimeError(f"Evidence quantification omitted {name}")
    json.loads((state.out / "summary.json").read_text(encoding="utf-8"))


def run(r1, r2, sample, panel, out, threads=1, skip_qc=False, bin_dir=None, keep_candidates=True):
    state = RunState(out, sample, skip_qc, keep_candidates)
    try:
        if threads < 1:
            raise ValueError("--threads must be positive")
        r1, r2, panel = Path(r1).resolve(), Path(r2).resolve(), Path(panel).resolve()
        if r1 == r2:
            raise ValueError("R1 and R2 must be distinct inputs")
        state.data["raw_inputs"] = [file_record(r1), file_record(r2)]
        metadata = [panel / name for name in ("panel_contigs.tsv", "feature_manifest.tsv", "annotations.tsv")]
        state.data["panel_metadata"] = [file_record(path, checksum=True) for path in metadata]
        state.data["panel_index_files"], state.data["index_checksum_verification"] = index_records(panel)
        state.data["code"] = [file_record(Path(__file__), checksum=True), file_record(Path(__file__).with_name("map_evidence.py"), checksum=True)]
        state.save()
        bowtie = command_prefix("bowtie2", bin_dir)
        fastp = None if skip_qc else command_prefix("fastp", bin_dir)
        capture_version(state, "bowtie2", bowtie)
        if fastp:
            capture_version(state, "fastp", fastp)
        state.data["screen"] = screen(state, r1, r2, panel, threads, bowtie, fastp)
        (state.out / "screen_counts.json").write_text(json.dumps(state.data["screen"], indent=2) + "\n", encoding="utf-8")
        state.save()
        sam = competitive_mapping(state, panel, threads, bowtie, state.data["screen"]["candidate_fragments"])
        quantify(state, panel, sam, sample)
        if not keep_candidates:
            for name in ("candidate_R1.fastq.gz", "candidate_R2.fastq.gz"):
                candidate = state.out / name
                # Only this invocation's exact output files may be removed.
                if candidate.resolve().parent != state.out:
                    raise RuntimeError("Candidate cleanup target escaped the output directory")
                candidate.unlink()
        state.data["status"] = "COMPLETE"
        state.data["completed_utc"] = timestamp()
        state.data["outputs"] = [file_record(state.out / name, checksum=name != "competitive.sam.gz")
                                  for name in ("gene_metrics.tsv", "summary.json", "screen_counts.json", "competitive.sam.gz")]
        state.save()
        (state.out / "COMPLETE").write_text(state.data["completed_utc"] + "\n", encoding="ascii")
        return state.data
    except BaseException as error:
        state.data.update(status="FAILED", failed_utc=timestamp(), error_type=type(error).__name__, error=str(error))
        state.save()
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("r1", "r2", "sample", "panel", "out"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--skip-qc", action="store_true")
    parser.add_argument("--bin-dir")
    parser.add_argument("--keep-candidates", action=argparse.BooleanOptionalAction, default=True)
    args = parser.parse_args()
    try:
        run(**vars(args))
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.exit(2, f"run_sample failed: {error}\n")
    print(f"Completed evidence outputs: {Path(args.out).resolve()}")


if __name__ == "__main__":
    main()
