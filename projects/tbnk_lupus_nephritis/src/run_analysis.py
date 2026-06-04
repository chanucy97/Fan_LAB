from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


FEATURES = [
    "cd3_pct",
    "cd4_pct",
    "cd8_pct",
    "b_pct",
    "nk_pct",
    "lymph_abs",
    "cd3_abs",
    "cd4_abs",
    "cd8_abs",
    "b_abs",
    "nk_abs",
    "cd4_cd8_ratio",
]

GROUPS = ["LN", "NS", "MN", "IgAN", "AAV"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Analyze de-identified routine TBNK lymphocyte subset data for LN versus non-LN immune-mediated nephropathies."
    )
    parser.add_argument("--input", required=True, help="Path to a de-identified CSV file.")
    parser.add_argument("--outdir", required=True, help="Directory for output tables.")
    parser.add_argument("--random-state", type=int, default=42, help="Random seed for cross-validation.")
    return parser.parse_args()


def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"analysis_group", "year", *FEATURES}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")
    df = df[df["analysis_group"].isin(GROUPS)].copy()
    df["target_ln"] = (df["analysis_group"] == "LN").astype(int)
    for col in FEATURES + ["year"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def make_pipeline(model) -> Pipeline:
    return Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", model),
        ]
    )


def fdr_bh(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    q = np.empty_like(p)
    prev = 1.0
    n = len(p)
    for rank in range(n - 1, -1, -1):
        idx = order[rank]
        val = p[idx] * n / (rank + 1)
        prev = min(prev, val)
        q[idx] = min(prev, 1.0)
    return q.tolist()


def mannwhitneyu_approx(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    """Two-sided Mann-Whitney U test using a normal approximation.

    This avoids requiring scipy for the lightweight public example workflow.
    For final manuscript analyses, scipy.stats.mannwhitneyu or equivalent
    validated statistical software can be used.
    """
    x = pd.Series(x).dropna().astype(float)
    y = pd.Series(y).dropna().astype(float)
    n1, n2 = len(x), len(y)
    if n1 == 0 or n2 == 0:
        return float("nan"), float("nan")
    combined = pd.concat([x, y], ignore_index=True)
    ranks = combined.rank(method="average")
    r1 = ranks.iloc[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2
    mean_u = n1 * n2 / 2
    tie_counts = combined.value_counts().to_numpy()
    tie_term = np.sum(tie_counts**3 - tie_counts)
    n = n1 + n2
    var_u = n1 * n2 / 12 * ((n + 1) - tie_term / (n * (n - 1))) if n > 1 else float("nan")
    if not np.isfinite(var_u) or var_u <= 0:
        return float(u1), float("nan")
    z = (u1 - mean_u) / math.sqrt(var_u)
    p_value = math.erfc(abs(z) / math.sqrt(2))
    return float(u1), float(p_value)


def cohort_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for group, sub in df.groupby("analysis_group"):
        row = {"analysis_group": group, "n": len(sub), "ln_target": int(sub["target_ln"].sum())}
        for feature in FEATURES:
            row[f"{feature}_median"] = sub[feature].median()
            row[f"{feature}_iqr"] = sub[feature].quantile(0.75) - sub[feature].quantile(0.25)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("analysis_group")


def ln_vs_nonln_tests(df: pd.DataFrame) -> pd.DataFrame:
    ln = df[df["target_ln"] == 1]
    comp = df[df["target_ln"] == 0]
    rows = []
    for feature in FEATURES:
        ln_values = ln[feature].dropna()
        comp_values = comp[feature].dropna()
        u_stat, p_value = mannwhitneyu_approx(ln_values, comp_values)
        rows.append(
            {
                "feature": feature,
                "ln_n": len(ln_values),
                "nonln_n": len(comp_values),
                "ln_median": ln_values.median(),
                "nonln_median": comp_values.median(),
                "delta_ln_minus_nonln": ln_values.median() - comp_values.median(),
                "mannwhitney_u": u_stat,
                "p_value": p_value,
            }
        )
    q_values = fdr_bh([r["p_value"] for r in rows])
    for row, q in zip(rows, q_values):
        row["fdr_q"] = q
    return pd.DataFrame(rows).sort_values("fdr_q")


def model_comparison(df: pd.DataFrame, random_state: int) -> pd.DataFrame:
    models = {
        "logistic_l2": LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear"),
        "random_forest": RandomForestClassifier(
            n_estimators=300,
            max_depth=4,
            min_samples_leaf=5,
            class_weight="balanced",
            random_state=random_state,
        ),
        "gradient_boosting": GradientBoostingClassifier(random_state=random_state, learning_rate=0.05, n_estimators=120, max_depth=2),
    }
    scoring = {
        "roc_auc": "roc_auc",
        "average_precision": "average_precision",
        "balanced_accuracy": "balanced_accuracy",
    }
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=random_state)
    rows = []
    for name, model in models.items():
        pipe = make_pipeline(model)
        scores = cross_validate(pipe, df[FEATURES], df["target_ln"], cv=cv, scoring=scoring, n_jobs=1)
        row = {"model": name, "n": len(df), "ln_n": int(df["target_ln"].sum()), "nonln_n": int((1 - df["target_ln"]).sum())}
        for metric in scoring:
            values = scores[f"test_{metric}"]
            row[f"{metric}_mean"] = float(np.mean(values))
            row[f"{metric}_sd"] = float(np.std(values, ddof=1))
        rows.append(row)
    return pd.DataFrame(rows).sort_values("roc_auc_mean", ascending=False)


def permutation_importance_table(df: pd.DataFrame, random_state: int) -> pd.DataFrame:
    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear"))
    pipe.fit(df[FEATURES], df["target_ln"])
    result = permutation_importance(pipe, df[FEATURES], df["target_ln"], scoring="roc_auc", n_repeats=30, random_state=random_state)
    return pd.DataFrame(
        {
            "feature": FEATURES,
            "roc_auc_decrease_mean": result.importances_mean,
            "roc_auc_decrease_sd": result.importances_std,
        }
    ).sort_values("roc_auc_decrease_mean", ascending=False)


def holdout_by_year(df: pd.DataFrame) -> pd.DataFrame:
    if df["year"].nunique() < 2:
        return pd.DataFrame([{"note": "Temporal holdout was skipped because fewer than two years are present."}])
    test_year = int(df["year"].max())
    train = df[df["year"] < test_year]
    test = df[df["year"] == test_year]
    if train["target_ln"].nunique() < 2 or test["target_ln"].nunique() < 2:
        return pd.DataFrame([{"note": "Temporal holdout was skipped because train or test set has only one class."}])
    pipe = make_pipeline(LogisticRegression(max_iter=5000, class_weight="balanced", solver="liblinear"))
    pipe.fit(train[FEATURES], train["target_ln"])
    prob = pipe.predict_proba(test[FEATURES])[:, 1]
    pred = (prob >= 0.5).astype(int)
    return pd.DataFrame(
        [
            {
                "train_year_max": test_year - 1,
                "test_year": test_year,
                "train_n": len(train),
                "test_n": len(test),
                "roc_auc": roc_auc_score(test["target_ln"], prob),
                "average_precision": average_precision_score(test["target_ln"], prob),
                "balanced_accuracy_at_0_5": balanced_accuracy_score(test["target_ln"], pred),
            }
        ]
    )


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_data(input_path)
    cohort_summary(df).to_csv(outdir / "cohort_summary.csv", index=False)
    ln_vs_nonln_tests(df).to_csv(outdir / "ln_vs_nonln_tests.csv", index=False)
    model_comparison(df, args.random_state).to_csv(outdir / "model_comparison.csv", index=False)
    permutation_importance_table(df, args.random_state).to_csv(outdir / "permutation_importance.csv", index=False)
    holdout_by_year(df).to_csv(outdir / "temporal_holdout.csv", index=False)
    print(f"Analysis complete. Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
