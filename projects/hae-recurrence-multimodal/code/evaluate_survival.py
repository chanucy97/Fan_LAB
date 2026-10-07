"""Generic horizon survival evaluation; never fits prediction models."""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.utils import concordance_index

from analysis_core import HORIZON, SEED, write_json

FUSIONS = ["feature_fusion", "score_fusion", "equal_rank_mean", "attention_fusion"]


def censor_weights(frame, horizon=HORIZON):
    t, e = frame.rfs_time_months.to_numpy(float), frame.rfs_event.to_numpy(int)
    grid = np.unique(t)
    before, after, g = [], [], 1.0
    for at in grid:
        before.append(g)
        g *= 1 - np.sum((t == at) & (e == 0)) / np.sum(t >= at)
        after.append(g)
    before, after = np.asarray(before), np.asarray(after)
    cases, controls = (e == 1) & (t <= horizon), t > horizon
    w = np.zeros(len(frame))
    denominators = before[np.searchsorted(grid, t[cases])]
    index = np.searchsorted(grid, horizon, side="right") - 1
    gh = after[index] if index >= 0 else 1.0
    if not cases.any() or not controls.any() or gh <= 0 or np.any(denominators <= 0):
        raise ValueError("36-month metric is not estimable in this sample")
    w[cases], w[controls] = 1 / denominators, 1 / gh
    return cases, controls, w, float(gh)


def weighted_auc(risk, cases, controls, weights, folds=None):
    """Cumulative/dynamic IPCW AUC of the declared risk score (higher=worse)."""
    r = np.asarray(risk)
    a, b = r[cases, None], r[None, controls]
    pair_weight = weights[cases, None] * weights[None, controls]
    if folds is not None:
        folds = np.asarray(folds)
        pair_weight *= folds[cases, None] == folds[None, controls]
    if pair_weight.sum() <= 0:
        raise ValueError("No comparable case-control pairs")
    return float((((a > b) + .5 * (a == b)) * pair_weight).sum() / pair_weight.sum())


def within_fold_cindex(frame, risk):
    """Count only comparable pairs predicted by the same fitted model."""
    t, e, r = frame.rfs_time_months.to_numpy(), frame.rfs_event.to_numpy(), np.asarray(risk)
    same = frame.cv_fold.to_numpy()[:, None] == frame.cv_fold.to_numpy()[None, :]
    comparable = (e[:, None] == 1) & ((t[:, None] < t[None, :]) |
                  ((t[:, None] == t[None, :]) & (e[None, :] == 0))) & same
    n = comparable.sum()
    if n == 0:
        return float("nan")
    good = (r[:, None] > r[None, :]) + .5 * (r[:, None] == r[None, :])
    return float(good[comparable].sum() / n)


def metrics(frame, risk, probability, development_group="train"):
    try:
        c = within_fold_cindex(frame, risk) if frame.group.iloc[0] == development_group else float(
            concordance_index(frame.rfs_time_months, -np.asarray(risk), frame.rfs_event))
    except ZeroDivisionError:
        c = np.nan
    try:
        cases, controls, w, _ = censor_weights(frame)
        p = np.asarray(probability)
        a = weighted_auc(risk, cases, controls, w,
                         frame.cv_fold.to_numpy() if frame.group.iloc[0] == development_group else None)
        b = float(np.sum(w * (cases.astype(float) - p) ** 2) / len(frame))
    except ValueError:
        a, b = np.nan, np.nan
    return np.asarray([c, a, b])


def resample_index(frame, rng, development_group="train"):
    if frame.group.iloc[0] == development_group:
        chunks = [np.flatnonzero(frame.cv_fold.to_numpy() == f) for f in sorted(frame.cv_fold.unique())]
        return np.concatenate([rng.choice(c, len(c), replace=True) for c in chunks])
    return rng.integers(0, len(frame), len(frame))


def calibration_and_dca(frame, p, model, group):
    cases, controls, w, gh = censor_weights(frame)
    # Tertiles are descriptive display groups; never used to update predictions.
    cuts = np.unique(np.quantile(p, [1 / 3, 2 / 3]))
    bins = np.searchsorted(cuts, p, side="right")
    calibration = []
    for b in np.unique(bins):
        mask = bins == b
        sub = frame.loc[mask]
        km = KaplanMeierFitter().fit(sub.rfs_time_months, sub.rfs_event)
        at = np.searchsorted(km.timeline, HORIZON, side="right") - 1
        interval = km.confidence_interval_.iloc[at]
        supported = bool((sub.rfs_time_months > HORIZON).any())
        calibration.append({"group": group, "model": model, "bin": int(b), "n": len(sub),
                            "mean_predicted_36m": float(np.mean(p[mask])),
                            "km_observed_36m": 1 - float(km.predict(HORIZON)) if supported else np.nan,
                            "ci_low": 1 - float(interval.iloc[1]) if supported else np.nan,
                            "ci_high": 1 - float(interval.iloc[0]) if supported else np.nan,
                            "at_risk_36m": int((sub.rfs_time_months >= HORIZON).sum()),
                            "followup_support": supported})
    dca = []
    for threshold in np.arange(.05, .51, .01):
        positive = p >= threshold
        nb = (np.sum(w * cases * positive) - threshold / (1 - threshold) * np.sum(w * controls * positive)) / len(frame)
        all_nb = (np.sum(w * cases) - threshold / (1 - threshold) * np.sum(w * controls)) / len(frame)
        dca.append({"group": group, "model": model, "threshold": float(threshold),
                    "net_benefit": float(nb), "treat_all": float(all_nb), "treat_none": 0.0,
                    "censor_survival_36m": gh})
    return calibration, dca


def evaluate(predictions, out, replicates=2000, development_group="train", groups=None, master_seed=SEED):
    out = Path(out)
    names = sorted(c[5:] for c in predictions if c.startswith("risk_"))
    rows, folds, cal, dca, comparisons, support = [], [], [], [], [], []
    metric_names = ["Harrell_C", "IPCW_AUC_36m", "IPCW_Brier_36m"]
    pairs = list(combinations(FUSIONS, 2)) + [(f, c) for f in FUSIONS for c in names if c.startswith("stage_")]
    groups = list(groups) if groups is not None else list(predictions.group.drop_duplicates())
    for gi, group in enumerate(groups):
        frame = predictions.loc[predictions.group.eq(group)].reset_index(drop=True)
        cases, controls, w, gh = censor_weights(frame)
        support.append({"group": group, "n": len(frame), "events": int(frame.rfs_event.sum()),
                        "events_by_36m": int(cases.sum()), "observed_beyond_36m": int(controls.sum()),
                        "censored_by_36m": int(np.sum(w == 0)), "G_36m": gh,
                        "max_IPCW": float(w.max())})
        point, boot = {}, {}
        for name in names:
            r, p = frame[f"risk_{name}"].to_numpy(), frame[f"prob36_{name}"].to_numpy()
            point[name] = metrics(frame, r, p, development_group)
            rng = np.random.default_rng(master_seed + gi)
            values = []
            for _ in range(replicates):
                index = resample_index(frame, rng, development_group)
                values.append(metrics(frame.iloc[index], r[index], p[index], development_group))
            boot[name] = np.asarray(values)
            for k, metric in enumerate(metric_names):
                valid = boot[name][:, k][np.isfinite(boot[name][:, k])]
                if len(valid) < .9 * replicates:
                    raise ValueError(f"Insufficient valid bootstrap replicates: {group}/{name}/{metric}")
                rows.append({"group": group, "model": name, "metric": metric, "estimate": point[name][k],
                             "ci_low": float(np.quantile(valid, .025)), "ci_high": float(np.quantile(valid, .975)),
                             "valid_replicates": len(valid), "evaluation": "within_outer_fold_pairs" if group == development_group and k in (0, 1) else "cohort",
                             "prediction_input": "risk_score" if k in (0, 1) else "36m_probability",
                             "CI_scope": "conditional_on_fitted_predictions"})
            a, b = calibration_and_dca(frame, p, name, group)
            cal.extend(a)
            dca.extend(b)
            if group == development_group:
                for fold in sorted(frame.cv_fold.unique()):
                    mask = frame.cv_fold.eq(fold)
                    estimates = metrics(frame.loc[mask], r[mask], p[mask], development_group)
                    folds.append({"model": name, "fold": fold, "n": int(mask.sum()),
                                  **dict(zip(metric_names, estimates))})
        for a, b in pairs:
            for k, metric in enumerate(metric_names[:2]):
                delta = point[a][k] - point[b][k]
                diffs = boot[a][:, k] - boot[b][:, k]
                diffs = diffs[np.isfinite(diffs)]
                # Approximate null-centred bootstrap; nonzero Monte Carlo correction.
                p = (1 + np.sum(np.abs(diffs - delta) >= abs(delta))) / (len(diffs) + 1)
                comparisons.append({"group": group, "model_a": a, "model_b": b, "metric": metric,
                                    "difference_a_minus_b": delta, "ci_low": float(np.quantile(diffs, .025)),
                                    "ci_high": float(np.quantile(diffs, .975)), "approx_bootstrap_p": float(p),
                                    "interpretation": "descriptive_internal" if group == development_group else "validation_comparison"})
    comp = pd.DataFrame(comparisons)
    comp["holm_p"] = np.nan
    # One family across both validation cohorts, all pairwise contrasts and both metrics.
    idx = comp.index[comp.group.ne(development_group)]
    ordered = comp.loc[idx, "approx_bootstrap_p"].sort_values()
    adjusted = np.minimum(1, np.maximum.accumulate(ordered.to_numpy() * np.arange(len(ordered), 0, -1)))
    comp.loc[ordered.index, "holm_p"] = adjusted
    for filename, records in [("performance", rows), ("development_fold_performance", folds),
                              ("calibration_36m", cal), ("decision_curve_36m", dca), ("followup_support", support)]:
        pd.DataFrame(records).to_csv(out / f"{filename}.csv", index=False)
    comp.to_csv(out / "paired_model_comparisons.csv", index=False)
    write_json(out / "evaluation_methods.json", {
        "horizon": HORIZON, "bootstrap_replicates": replicates,
        "development_group": development_group,
        "groups": groups,
        "development_cindex": "Comparable-pair weighted within-fold C-index; no comparisons across different fitted Cox models",
        "AUC": "Cumulative/dynamic IPCW AUC of the SAME declared risk score used by C-index; development case-control pairs restricted to same outer fold, pooled by IPCW pair weight",
        "equal_rank_mean": "Fixed mean score for C-index/AUC; separately trained nonnegative Cox calibration only for probability/Brier/calibration/DCA; zero calibration slope never substitutes a constant score for discrimination",
        "IPCW": "Cohort-specific marginal censoring Kaplan-Meier; event weights G(T-), control weights G(36); independent censoring assumed",
        "development_CIs": "Fixed predictions, resampling within outer folds; conditional descriptive intervals, do not include full model-development variability",
        "validation_CIs": "Patient bootstrap with fixed trained models and re-estimated censoring distribution in every replicate",
        "comparisons": f"Same bootstrap indices for every model; approximate null-centred bootstrap tests; all {len(idx)} validation tests form one Holm family",
        "calibration": "Descriptive tertiles, observed probability from KM with pointwise log-log intervals; no model recalibration on validation cohorts",
        "DCA": "IPCW net benefit at 36 months; thresholds .05-.50 are an exploratory grid, not clinically validated treatment thresholds",
        "report_all_four_fusions": True, "model_selection_by_validation_performance": False})
