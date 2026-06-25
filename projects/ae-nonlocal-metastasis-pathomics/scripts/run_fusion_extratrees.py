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
    clf = ExtraTreesClassifier(n_estimators=800, max_features="sqrt", class_weight="balanced", random_state=seed, n_jobs=-1)
    return Pipeline([("preprocess", pre), ("model", clf)])


def main() -> None:
    parser = argparse.ArgumentParser(description="Run clinical + WSI probability fusion under predefined patient folds.")
    parser.add_argument("--clinical-table", required=True, type=Path)
    parser.add_argument("--wsi-predictions", required=True, type=Path)
    parser.add_argument("--fold-table", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--id-col", default="patient_id")
    parser.add_argument("--label-col", default="label")
    parser.add_argument("--numeric-cols", required=True)
    parser.add_argument("--categorical-cols", default="")
    parser.add_argument("--wsi-score-col", default="score")
    parser.add_argument("--seed", type=int, default=2024)
    args = parser.parse_args()

    numeric_cols = parse_cols(args.numeric_cols)
    categorical_cols = parse_cols(args.categorical_cols)
    clinical = pd.read_csv(args.clinical_table)
    wsi = pd.read_csv(args.wsi_predictions)
    folds = pd.read_csv(args.fold_table)
    for frame in [clinical, wsi, folds]:
        frame[args.id_col] = frame[args.id_col].astype(str)
    wsi_feature = wsi[[args.id_col, "fold", "group", args.wsi_score_col]].rename(columns={args.wsi_score_col: "wsi_score"})
    df = clinical.merge(folds[[args.id_col, args.label_col, "fold", "group"]], on=[args.id_col, args.label_col], how="inner")
    df = df.merge(wsi_feature, on=[args.id_col, "fold", "group"], how="inner")
    if df.empty:
        raise ValueError("No overlap among clinical table, WSI predictions, and folds")

    model_numeric = numeric_cols + ["wsi_score"]
    pred_rows = []
    metric_rows = []
    for fold in sorted(df["fold"].unique()):
        train = df[(df["fold"] == fold) & (df["group"] == "train")]
        test = df[(df["fold"] == fold) & (df["group"] == "test")]
        model = make_pipeline(model_numeric, categorical_cols, args.seed + int(fold))
        model.fit(train[model_numeric + categorical_cols], train[args.label_col].astype(int))
        score = model.predict_proba(test[model_numeric + categorical_cols])[:, 1]
        fold_pred = test[[args.id_col, args.label_col, "fold", "group", "wsi_score"]].copy()
        fold_pred["score"] = score
        fold_pred["model"] = "clinical_wsi_fusion_extratrees"
        pred_rows.append(fold_pred)
        row = {"fold": fold, "model": "clinical_wsi_fusion_extratrees"}
        row.update(binary_metrics(fold_pred[args.label_col], fold_pred["score"]))
        metric_rows.append(row)

    pred = pd.concat(pred_rows, ignore_index=True)
    summary = {"model": "clinical_wsi_fusion_extratrees"}
    summary.update(binary_metrics(pred[args.label_col], pred["score"]))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pred.rename(columns={args.id_col: "patient_id", args.label_col: "label"}).to_csv(args.output_dir / "patient_predictions.csv", index=False)
    pd.DataFrame(metric_rows).to_csv(args.output_dir / "fold_metrics.csv", index=False)
    pd.DataFrame([summary]).to_csv(args.output_dir / "summary_metrics.csv", index=False)


if __name__ == "__main__":
    main()
