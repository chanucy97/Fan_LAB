"""Create one wholly artificial format smoke fixture; no genomic source data."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from analysis_support import summarize
from evaluate_benchmark import RECOMMENDED_GATE


def build(output):
    output = Path(output)
    if output.exists():
        raise FileExistsError("Refuse to replace example directory")
    output.mkdir(parents=True)
    inner = "G" * 4309
    full = "AC" * 10 + inner + "TC" * 10
    cid, candidate, group = "synthetic_case", "synthetic_candidate", "canonical_4309"
    metadata = {candidate: {"length": len(full), "inner_start": 20,
                "inner_end": 4329, "junctions": [20, 4329], "class_id": group}}
    extraction = {"config": {"ledger_mode": "accepted"},
        "audit": {"accepted_fragments": 20, "ledger_rows_omitted_by_cap": 0},
        "fragment_ledger": [{"qname": "synthetic_q" + str(i), "candidates": [candidate],
            "candidate_evidence": {candidate: {"direct_junctions": [20, 4329], "bridged_junctions": []}}}
            for i in range(20)],
        "candidate_coverage": {candidate: {"length": len(full), "eligible_base_breadth": 1.0,
            "near_best_base_breadth": 1.0, "eligible_fragments": 20, "near_best_fragments": 20}}}
    fit = {"status": "fitted_descriptive", "accepted_fragment_count": 20, "bestfit_tv": 0.0,
           "config": {"candidate_groups": {group: [candidate]}},
           "group_ranges": {group: [1.0, 1.0]}, "bestfit_weights": {candidate: 1.0}}
    decision, baselines = summarize(extraction, fit, metadata, RECOMMENDED_GATE)
    result = {"case_id": cid, "extraction": extraction, "fit": fit, "paired_only_fit": fit,
              "decision": decision, "paired_only_decision": decision, "reference_baseline": baselines}
    row = {"structure_id": "synthetic_interval", "interval_start_0based": 100,
        "interval_end_0based_exclusive": 100+len(full), "left_anchor_boundary_0based": 120,
        "right_anchor_boundary_0based": 4429, "orientation_to_reference": "+",
        "interval_sequence_sha256": hashlib.sha256(full.encode("ascii")).hexdigest(),
        "between_anchor_bp": len(inner), "inner_anchor_edges_complete": True, "both_anchors_full_length": True}
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=list(row), delimiter="\t", lineterminator="\n")
    writer.writeheader(); writer.writerow(row)
    texts = {"truth_full": ">synthetic_truth\n"+full+"\n",
             "truth_inner": ">synthetic_truth\n"+inner+"\n",
             "structure_table": stream.getvalue(),
             "anchored_fasta": ">synthetic_interval\n"+full+"\n"}
    history = {"schema": "artificial-format-fixture.v1", "cases": [{"case_id": cid, "files": {
        key: {"exists": True, "text": text, "bytes": len(text.encode("utf-8")),
              "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest()} for key, text in texts.items()}}]}
    inventory = {"case_count": 1, "positive_alleles": 1, "cases": [{"case_id": cid,
        "fixture_set": "represented_stress24", "components": [{"positive": True, "coverage": 1}], "error_rate": 0}]}
    for name, value in [("inventory.json", inventory), ("panel_metadata.json", metadata), ("history.json", history)]:
        (output/name).write_text(json.dumps(value, indent=2)+"\n", encoding="utf-8")
    target = output/"results"/cid
    target.mkdir(parents=True)
    (target/"RESULT.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    return {"status": "CREATED", "scope": "artificial_format_smoke_only", "case_count": 1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.output)))


if __name__ == "__main__":
    main()
