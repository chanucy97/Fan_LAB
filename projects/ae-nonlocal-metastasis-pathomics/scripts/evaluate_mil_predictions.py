#!/usr/bin/env python
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import pandas as pd

from metrics import binary_metrics


def pick_score_column(df: pd.DataFrame, requested: str | None) -> str:
    if requested:
        return requested
    for col in ["mil_prob_1", "prob_1", "prob", "score", "y_score"]:
        if col in df.columns:
            return col
    raise ValueError("No probability column found. Pass --score-col explicitly.")


def normalize_columns(df: pd.DataFrame, path: str, score_col: str | None) -> pd.DataFrame:
    col = pick_score_column(df, score_col)
    rename = {col: "score"}
    if "true" in df.columns and "label" not in df.columns:
        rename["true"] = "label"
    if "case_id" in df.columns and "patient_id" not in df.columns:
        rename["case_id"] = "patient_id"
    out = df.rename(columns=rename).copy()
    required = {"patient_id", "label", "score"}
    missing = required.difference(out.columns)
    if missing:
        raise ValueError(f"{path}: missing required columns {sorted(missing)}")
    if "slide_id" not in out.columns:
        out["slide_id"] = out["patient_id"].astype(str)
    if "fold" not in out.columns:
        out["fold"] = "unknown"
    if "group" not in out.columns:
        out["group"] = "test"
    out["source_file"] = path
    return out[["patient_id", "slide_id", "label", "fold", "group", "score", "source_file"]]


def main() -> None:
    parser = argparse.ArgumentParser(description="Evaluate WSI/MIL prediction exports.")
    parser.add_argument("--predictions-glob", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--score-col", default=None)
    parser.add_argument("--patient-aggregation", choices=["mean", "max"], default="mean")
    args = parser.parse_args()

    paths = sorted(glob.glob(args.predictions_glob, recursive=True))
    if not paths:
        raise FileNotFoundError(args.predictions_glob)
    rows = [normalize_columns(pd.read_csv(path), path, args.score_col) for path in paths]
    slide_df = pd.concat(rows, ignore_index=True)
    slide_df["label"] = pd.to_numeric(slide_df["label"], errors="raise").astype(int)
    slide_df["score"] = pd.to_numeric(slide_df["score"], errors="raise")

    agg_func = "mean" if args.patient_aggregation == "mean" else "max"
    patient_df = (
        slide_df.groupby(["patient_id", "fold", "group"], as_index=False)
        .agg(label=("label", "first"), score=("score", agg_func), n_slides=("slide_id", "nunique"))
    )

    fold_rows = []
    for (fold, group), part in patient_df.groupby(["fold", "group"], dropna=False):
        row = {"fold": fold, "group": group}
        row.update(binary_metrics(part["label"], part["score"]))
        fold_rows.append(row)
    summary = binary_metrics(patient_df["label"], patient_df["score"])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    slide_df.to_csv(args.output_dir / "slide_predictions.csv", index=False)
    patient_df.to_csv(args.output_dir / "patient_predictions.csv", index=False)
    pd.DataFrame(fold_rows).to_csv(args.output_dir / "fold_metrics.csv", index=False)
    pd.DataFrame([summary]).to_csv(args.output_dir / "summary_metrics.csv", index=False)


if __name__ == "__main__":
    main()
