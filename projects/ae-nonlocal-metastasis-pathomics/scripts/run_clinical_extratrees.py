#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from metrics import binary_metrics


def parse_cols(text: str) -> list[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def make_pipeline(numeric_cols: list[str], categorical_cols: list[str], seed: int) -> Pipeline:
    pre = ColumnTransformer(
        transformers=[
            ("num", Pipeline([("impute", SimpleImputer(strategy="median")), ("scale", StandardScaler())]), numeric_cols),
            ("cat", Pipeline([("impute", SimpleImputer(strategy="most_frequent")), ("onehot", OneHotEncoder(handle_unknown="ignore"))]), categorical_cols),
        ],
        remainder="drop",
    )
    clf = ExtraTreesClassifier(
        n_estimators=800,
        criterion="gini",
        max_features="sqrt",
        class_weight="balanced",
        random_state=seed,
        n_jobs=-1,
    )
    return Pipeline([("preprocess", pre), ("model", clf)])


def main() -> None:
    parser = argparse.ArgumentParser(description="Run clinical ExtraTrees baseline using predefined patient folds.")
    parser.add_argument("--clinical-table", required=True, type=Path)
    parser.add_argument("--fold-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--id-col", default="patient_id")
    parser.add_argument("--label-col", default="label")
    parser.add_argument("--numeric-cols", required=True)
    parser.add_argument("--categorical-cols", default="")
    parser.add_argument("--seed", type=int, default=2024)
    args = parser.parse_args()

    numeric_cols = parse_cols(args.numeric_cols)
    categorical_cols = parse_cols(args.categorical_cols)
    df = pd.read_csv(args.clinical_table)
    folds = pd.read_csv(args.fold_table)
    df[args.id_col] = df[args.id_col].astype(str)
    folds[args.id_col] = folds[args.id_col].astype(str)
    df = df.merge(folds[[args.id_col, args.label_col, "fold", "group"]], on=[args.id_col, args.label_col], how="inner")
    if df.empty:
        raise ValueError("No overlap between clinical table and fold table")

    pred_rows = []
    metric_rows = []
    for fold in sorted(df["fold"].unique()):
        train = df[(df["fold"] == fold) & (df["group"] == "train")]
        test = df[(df["fold"] == fold) & (df["group"] == "test")]
        model = make_pipeline(numeric_cols, categorical_cols, args.seed + int(fold))
        model.fit(train[numeric_cols + categorical_cols], train[args.label_col].astype(int))
        score = model.predict_proba(test[numeric_cols + categorical_cols])[:, 1]
        fold_pred = test[[args.id_col, args.label_col, "fold", "group"]].copy()
        fold_pred["score"] = score
        fold_pred["model"] = "clinical_extratrees"
        pred_rows.append(fold_pred)
        row = {"fold": fold, "model": "clinical_extratrees"}
        row.update(binary_metrics(fold_pred[args.label_col], fold_pred["score"]))
        metric_rows.append(row)

    pred = pd.concat(pred_rows, ignore_index=True)
    summary = {"model": "clinical_extratrees"}
    summary.update(binary_metrics(pred[args.label_col], pred["score"]))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pred.rename(columns={args.id_col: "patient_id", args.label_col: "label"}).to_csv(args.output_dir / "patient_predictions.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(args.output_dir / "fold_metrics.csv", index=False)
    pd.DataFrame([summary]).to_csv(args.output_dir / "summary_metrics.csv", index=False)


if __name__ == "__main__":
    main()
