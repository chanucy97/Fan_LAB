"""Training-only preprocessing, Cox models, and nested stacking building blocks."""
from __future__ import annotations

import hashlib
import json
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceWarning, ConvergenceError
from scipy.linalg import LinAlgWarning
from sklearn.model_selection import StratifiedKFold
from clinical_schema import validate_clinical, category_specs, fit_category, apply_category

SEED = 20260727
DIMENSIONS = ["clinical", "radiomics", "habitat", "deep_2p5D"]
KEY = ["ID", "group", "cv_fold", "rfs_event", "rfs_time_months"]
OUTCOMES = ["rfs_time_months", "rfs_event"]
HORIZON = 36.0
COX_PRECISION = 1e-9
LASSO_ZERO_TOL = 1e-6  # abs(beta * training SD): numerical zero, not a p-value cutoff


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def fit_cox_stable(frame, penalty, l1):
    """Numerical step-size retries preserve the exact statistical objective."""
    failures = []
    for step in (.25, .1, .025):
        try:
            with warnings.catch_warnings():
                for category in (ConvergenceWarning, RuntimeWarning, LinAlgWarning):
                    warnings.simplefilter("error", category)
                model = CoxPHFitter(penalizer=penalty, l1_ratio=l1).fit(
                    frame, *OUTCOMES, fit_options={"step_size": step, "max_steps": 1000,
                                                 "precision": COX_PRECISION})
            return model, {"step_size": step, "failed_numerical_attempts": failures}
        except (ConvergenceError, ConvergenceWarning, RuntimeWarning, LinAlgWarning, ArithmeticError) as error:
            failures.append({"step_size": step, "reason": str(error)[:300]})
    raise ConvergenceError(f"Cox numerical retries exhausted without changing penalty: {failures}")


def percentile(reference, values):
    """Midrank empirical CDF: identical mapping for fitting and application rows."""
    reference = np.sort(np.asarray(reference, float))
    values = np.asarray(values, float)
    return (np.searchsorted(reference, values, side="left") +
            np.searchsorted(reference, values, side="right")) / (2 * len(reference))


class Screen:
    """All encodings, imputation, correlation filtering and screening fit locally."""
    def fit(self, raw, outcomes, categorical, cap):
        self.columns = list(raw.columns)
        if not isinstance(categorical, dict):
            if categorical:
                raise ValueError("Categorical predictors require explicit schema specifications")
            categorical = {}
        self.category_specs = categorical
        self.categorical = [c for c in self.columns if c in categorical]
        self.category_states = {}
        self.levels, self.fill = {}, {}
        for c in self.columns:
            a = raw[c].dropna()
            if not len(a):
                if c in self.categorical:
                    raise ValueError(f"No observed categorical training values: {c}")
                continue
            if c in self.categorical:
                state = fit_category(raw[c], c, self.category_specs[c])
                self.category_states[c] = state
                self.fill[c] = state["mode"]
                self.levels[c] = state["observed_codes"]
            else:
                self.fill[c] = float(a.median())
        x = self.encode(raw)
        self.nonconstant = x.columns[x.var(ddof=0) > 1e-10].tolist()
        x = x[self.nonconstant]
        if not x.shape[1]:
            raise ValueError("No nonconstant training features")
        corr = x.rank().corr().abs()
        self.keep = []
        for c in x:
            if not self.keep or (corr.loc[c, self.keep] < .90).all():
                self.keep.append(c)
        x = x[self.keep]
        self.mean, self.sd = x.mean(), x.std(ddof=0)
        x = (x - self.mean) / self.sd
        scores, failed = [], []
        for c in x:
            try:
                m, _ = fit_cox_stable(pd.concat([x[[c]], outcomes], axis=1), .01, 0)
                p = float(m.summary.loc[c, "p"])
                if np.isfinite(p):
                    scores.append((p, c))
            except (ValueError, ArithmeticError, ConvergenceWarning) as error:
                failed.append({"feature": c, "reason": str(error)[:180]})
            except Exception as error:
                # lifelines convergence exceptions are version-specific; retain the reason.
                from lifelines.exceptions import ConvergenceError
                if not isinstance(error, ConvergenceError):
                    raise
                failed.append({"feature": c, "reason": str(error)[:180]})
        self.selected = [c for _, c in sorted(scores)[:cap]]
        if not self.selected:
            raise ValueError("No stable univariable training fit")
        self.audit = {"raw_columns": len(raw.columns), "encoded_nonconstant": len(self.nonconstant),
                      "after_correlation": len(self.keep), "selected": self.selected,
                      "screening_p_values_are_ranking_only": True, "failed_univariate": failed,
                      "categorical_fitting_states": self.category_states}
        return self

    def encode(self, raw):
        out = {}
        self.last_application_audit = {}
        for c, fill in self.fill.items():
            a = raw[c].fillna(fill)
            if c in self.categorical:
                a, audit = apply_category(raw[c], c, self.category_specs[c], self.category_states[c])
                self.last_application_audit[c] = audit
                for level in self.levels[c]:
                    out[f"{c}=={level:g}"] = a.eq(level).astype(float)
            else:
                out[c] = a.astype(float)
        return pd.DataFrame(out, index=raw.index)

    def transform(self, raw):
        x = self.encode(raw)[self.keep]
        return ((x - self.mean) / self.sd)[self.selected]


class SparseCox:
    def fit(self, x, outcomes, cap=10, penalty=.05, l1=1.0, select=True):
        if not np.isfinite(x.to_numpy()).all() or outcomes.rfs_event.sum() < 5:
            raise ValueError("Invalid training predictors or fewer than 5 events")
        self.selected = x.columns[x.std(ddof=1) > 1e-12].tolist()
        x = x[self.selected]
        selection_solver = None
        selection_effects = {}
        if select and self.selected:
            first, selection_solver = fit_cox_stable(pd.concat([x, outcomes], axis=1), penalty, l1)
            ranked = (first.params_ * x.std(ddof=1)).abs().sort_values(ascending=False)
            selection_effects = ranked.to_dict()
            limit = min(cap, max(1, int(outcomes.rfs_event.sum() // 5)))
            self.selected = ranked[ranked > LASSO_ZERO_TOL].index[:limit].tolist()
        refits = []
        while self.selected:
            self.model, solver = fit_cox_stable(pd.concat([x[self.selected], outcomes], axis=1), penalty, l1)
            refits.append(solver)
            retained = ((self.model.params_ * x[self.selected].std(ddof=1)).abs() > LASSO_ZERO_TOL)
            if not select or retained.all():
                break
            self.selected = retained.index[retained].tolist()
        if not self.selected:
            self.model = None
            self.reference = np.zeros(len(outcomes))
            time = outcomes.rfs_time_months.to_numpy()
            event = outcomes.rfs_event.to_numpy()
            self.null_h36 = float(sum(np.sum((time == t) & (event == 1)) / np.sum(time >= t)
                                     for t in np.unique(time[(event == 1) & (time <= HORIZON)])))
            self.training_frame = outcomes.copy()
            self.audit = {"selected": [], "coefficients": {}, "baseline_only": True,
                          "baseline_cumulative_hazard_36m": self.null_h36, "events": int(event.sum()),
                          "zero_tolerance_standardized": LASSO_ZERO_TOL, "solver_precision": COX_PRECISION,
                          "selection_standardized_abs_coefficients": selection_effects,
                          "selection_solver": selection_solver, "refits": refits}
            return self
        self.reference = self.model.predict_log_partial_hazard(x[self.selected]).to_numpy()
        self.training_frame = pd.concat([x[self.selected], outcomes], axis=1).copy()
        self.audit = {"selected": self.selected, "coefficients": self.model.params_.to_dict(),
                      "events": int(outcomes.rfs_event.sum()), "penalty": penalty, "l1_ratio": l1,
                      "selection_solver": selection_solver, "final_solver": solver,
                      "refits": refits, "baseline_only": False,
                      "zero_tolerance_standardized": LASSO_ZERO_TOL, "solver_precision": COX_PRECISION,
                      "selection_standardized_abs_coefficients": selection_effects}
        return self

    def predict(self, x):
        if self.model is None:
            return np.zeros(len(x)), np.full(len(x), .5), np.full(len(x), -np.expm1(-self.null_h36))
        raw = self.model.predict_log_partial_hazard(x[self.selected]).to_numpy()
        base = self.model.baseline_cumulative_hazard_
        before = base.loc[base.index <= HORIZON]
        h36 = float(before.iloc[-1, 0]) if len(before) else 0.0
        p36 = -np.expm1(-h36 * np.exp(raw))
        return raw, percentile(self.reference, raw), p36


class BaseModels:
    def fit(self, sources, data, ids, settings=None):
        settings = settings or {}
        self.screens, self.models = {}, {}
        y = data.loc[ids, OUTCOMES]
        all_x, self.audit = [], {}
        for dimension, entries in sources.items():
            chunks, screen_audit = [], []
            for name, frame, categories in entries:
                cap = settings.get("clinical_cap", 5) if dimension == "clinical" else settings.get("per_source_cap", 2)
                screen = Screen().fit(frame.loc[ids], y, categories, int(cap))
                self.screens[name] = screen
                chunks.append(screen.transform(frame.loc[ids]))
                screen_audit.append({"source": name, **screen.audit})
            x = pd.concat(chunks, axis=1)
            model = SparseCox().fit(x, y)
            self.models[dimension] = model
            all_x.append(x)
            self.audit[dimension] = {"screens": screen_audit, "cox": model.audit}
        self.models["feature_fusion"] = SparseCox().fit(
            pd.concat(all_x, axis=1), y, cap=int(settings.get("fusion_cap", 10)))
        self.audit["feature_fusion"] = self.models["feature_fusion"].audit
        self.fitting_ids = data.loc[ids, "ID"].tolist()
        return self

    def predict(self, sources, ids):
        scores, predictions, all_x = {}, {}, []
        self.application_category_audit = {}
        for dimension, entries in sources.items():
            x = pd.concat([self.screens[name].transform(frame.loc[ids]) for name, frame, _ in entries], axis=1)
            for name, _, _ in entries:
                if self.screens[name].last_application_audit:
                    self.application_category_audit[name] = self.screens[name].last_application_audit
            raw, rank, p = self.models[dimension].predict(x)
            scores[dimension] = rank
            predictions[dimension] = (raw, rank, p)
            all_x.append(x)
        predictions["feature_fusion"] = self.models["feature_fusion"].predict(pd.concat(all_x, axis=1))
        return pd.DataFrame(scores, index=ids), predictions


def inner_splits(data, ids, seed, n_splits=4):
    ids = np.asarray(ids)
    counts = data.loc[ids, "rfs_event"].value_counts()
    if len(counts) != 2 or counts.min() < n_splits:
        raise ValueError(f"Inner folds require at least {n_splits} events and {n_splits} censored patients")
    for fit, val in StratifiedKFold(n_splits, shuffle=True, random_state=seed).split(ids, data.loc[ids, "rfs_event"]):
        yield ids[fit], ids[val]


class Builder:
    def __init__(self, sources, data, out, settings=None):
        self.sources, self.data, self.out = sources, data, Path(out)
        self.settings = settings or {}
        self.cache = {}
        (self.out / "model_artifacts").mkdir(parents=True, exist_ok=True)

    def base(self, ids):
        key = tuple(sorted(int(i) for i in ids))
        if key not in self.cache:
            print(f"  fitting base models: n={len(ids)}, events={int(self.data.loc[list(key), 'rfs_event'].sum())}", flush=True)
            model = BaseModels().fit(self.sources, self.data, list(key), self.settings)
            name = hashlib.sha256(np.asarray(key, dtype=np.int64).tobytes()).hexdigest()[:16]
            write_json(self.out / "model_artifacts" / f"base_{name}_audit.json",
                       {"fitting_IDs": model.fitting_ids, "models": model.audit})
            # Save every fitted preprocessing pipeline for exact reapplication/audit.
            joblib.dump(model, self.out / "model_artifacts" / f"base_{name}.joblib", compress=3)
            self.cache[key] = model
        return self.cache[key]

    def stack(self, fitting, application, seed, label, n_splits=4):
        fitting, application = np.asarray(fitting), np.asarray(application)
        if set(fitting) & set(application):
            raise ValueError("Stack fitting/application overlap")
        oof = pd.DataFrame(index=fitting, columns=DIMENSIONS, dtype=float)
        folds = []
        for fit, val in inner_splits(self.data, fitting, seed, n_splits):
            model = self.base(fit)
            score, _ = model.predict(self.sources, val)
            oof.loc[val] = score
            folds.append({"fit_IDs": self.data.loc[fit, "ID"].tolist(), "heldout_IDs": self.data.loc[val, "ID"].tolist(),
                          "heldout_category_handling": model.application_category_audit})
        if not np.isfinite(oof.to_numpy()).all():
            raise ValueError("Incomplete inner OOF predictions")
        full = self.base(fitting)
        apply, predictions = full.predict(self.sources, application)
        write_json(self.out / "model_artifacts" / f"{label}_partition_audit.json",
                   {"fitting_IDs": self.data.loc[fitting, "ID"].tolist(),
                    "application_IDs": self.data.loc[application, "ID"].tolist(), "inner_folds": folds,
                    "application_category_handling": full.application_category_audit})
        return oof, apply, predictions, full


def read_sources(paths, data):
    specs = {"clinical": [("clinical", paths["clinical"])],
             "radiomics": [("CTP_global", paths["ct_p_global"]), ("T2WI_global", paths["t2wi_global"])],
             "habitat": [("CTP_habitat", paths["ct_p_habitat"]), ("T2WI_habitat", paths["t2wi_habitat"])],
             "deep_2p5D": [("CTP_2p5D", paths["ct_p_2p5d"]), ("T2WI_2p5D", paths["t2wi_2p5d"])]}
    sources, audit = {}, {}
    forbidden = {"group", "cv_fold", "rfs_event", "rfs_time_months", "event", "time", "label", "outcome", "status", "recurrence", "dfs", "rfs", "os"}
    for dim, entries in specs.items():
        sources[dim] = []
        for name, path in entries:
            frame = pd.read_csv(path, dtype={"ID": str})
            if frame.ID.isna().any() or frame.ID.duplicated().any() or not set(data.ID).issubset(set(frame.ID)):
                raise ValueError(f"ID mismatch: {path}")
            excluded_ids = sorted(set(frame.ID) - set(data.ID))
            schema_audit = validate_clinical(frame) if dim == "clinical" else None
            bad = [c for c in frame if c.lower() in forbidden or c.lower().startswith('risk_')]
            if bad:
                raise ValueError(f"Forbidden predictor columns: {bad}")
            frame = frame.set_index("ID").loc[data.ID].reset_index(drop=True)
            frame.index = data.index
            frame = frame.apply(pd.to_numeric, errors="raise").astype(float)
            if np.isinf(frame.to_numpy()).any():
                raise ValueError(f"Infinite predictors: {path}")
            categories = category_specs(f"{name}::") if dim == "clinical" else {}
            frame.columns = [f"{name}::{c}" for c in frame]
            sources[dim].append((name, frame, categories))
            audit[name] = {"path": str(Path(path).resolve()), "sha256": digest(path), "features": len(frame.columns),
                           "categorical": list(categories), "clinical_schema_audit": schema_audit,
                           "missing_cells": int(frame.isna().sum().sum()),
                           "excluded_nonmanifest_IDs": excluded_ids}
    return sources, audit
