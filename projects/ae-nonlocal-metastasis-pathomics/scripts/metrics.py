from __future__ import annotations

import math
from typing import Any

import numpy as np
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)


def binary_metrics(y_true: Any, y_score: Any, threshold: float = 0.5) -> dict[str, float]:
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score).astype(float)
    y_pred = (y_score >= threshold).astype(int)
    out = {
        "n": float(len(y_true)),
        "positive": float(np.sum(y_true == 1)),
        "negative": float(np.sum(y_true == 0)),
        "threshold": float(threshold),
        "accuracy": float(accuracy_score(y_true, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, y_pred)),
        "macro_f1": float(f1_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro_recall": float(recall_score(y_true, y_pred, average="macro", zero_division=0)),
        "macro_precision": float(precision_score(y_true, y_pred, average="macro", zero_division=0)),
        "brier": float(brier_score_loss(y_true, y_score)),
    }
    out["roc_auc"] = float(roc_auc_score(y_true, y_score)) if len(np.unique(y_true)) == 2 else math.nan
    out["pr_auc"] = float(average_precision_score(y_true, y_score)) if len(np.unique(y_true)) == 2 else math.nan
    return out


def summarize_by_group(df, group_cols, label_col="label", score_col="score"):
    rows = []
    for key, part in df.groupby(group_cols, dropna=False):
        if not isinstance(key, tuple):
            key = (key,)
        row = dict(zip(group_cols, key))
        row.update(binary_metrics(part[label_col], part[score_col]))
        rows.append(row)
    return rows
