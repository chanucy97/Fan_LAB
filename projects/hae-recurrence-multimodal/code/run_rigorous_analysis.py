"""Generic patient-level nested survival modelling driven by user-supplied inputs."""
from __future__ import annotations

import argparse
import importlib.metadata
import json
import platform
import shutil
import sys
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from sklearn.model_selection import StratifiedShuffleSplit
from lifelines.statistics import proportional_hazard_test

from analysis_core import (Builder, DIMENSIONS, HORIZON, KEY, OUTCOMES, SEED,
                           SparseCox, digest, percentile, read_sources, write_json)
from attention_core import breslow_probability, fit_predict, train, LEARNING_RATE, WEIGHT_DECAY
from evaluate_survival import evaluate
from clinical_schema import SCHEMA, fit_category, apply_category, set_schema

DEFAULT_PNM = {}


def validate_manifest(data, config):
    if set(data.columns) != set(KEY):
        raise ValueError("Manifest must contain exactly: ID, group, cv_fold, rfs_event, rfs_time_months")
    if data.isna().any().any() or data.ID.duplicated().any():
        raise ValueError("Missing manifest values or duplicated IDs")
    development_group = config["development_group"]
    if development_group not in set(data.group):
        raise ValueError(f"Development group is absent from manifest: {development_group}")
    if data.group.nunique() < 2:
        raise ValueError("At least one separate validation group is required")
    if not data.rfs_event.isin([0, 1]).all() or not np.isfinite(data.rfs_time_months).all() or not data.rfs_time_months.gt(0).all():
        raise ValueError("Invalid survival endpoint")
    dev = data.group.eq(development_group)
    n_outer_folds = int(config.get("outer_folds", 5))
    if n_outer_folds < 2 or set(data.loc[dev, "cv_fold"]) != set(range(n_outer_folds)):
        raise ValueError("Development cv_fold values must be consecutive 0..outer_folds-1")
    if not data.loc[~dev, "cv_fold"].eq(-1).all():
        raise ValueError("Validation rows must have cv_fold=-1")
    for f in range(n_outer_folds):
        part = data.loc[dev & data.cv_fold.eq(f)]
        if part.rfs_event.sum() < 2 or part.rfs_event.eq(0).sum() < 2:
            raise ValueError("Outer fold lacks sufficient events/censoring")


def tune_epochs(builder, fitting, seed, label, epochs, patience, inner_folds):
    """Dedicated inner validation; its labels never enter the subtraining bases."""
    ids = np.asarray(fitting)
    splitter = StratifiedShuffleSplit(1, test_size=.25, random_state=seed)
    sub, val = next(splitter.split(ids, builder.data.loc[ids, "rfs_event"]))
    sub, val = ids[sub], ids[val]
    x, v, _, _ = builder.stack(sub, val, seed + 1, label + "_early_stop", inner_folds)
    best, history = train(x, builder.data.loc[sub, OUTCOMES], epochs, seed,
                          tuning=(v, builder.data.loc[val, OUTCOMES]), patience=patience)
    write_json(builder.out / "model_artifacts" / f"{label}_early_stop_history.json",
               {"selected_epochs": best, "inner_fitting_IDs": builder.data.loc[sub, "ID"].tolist(),
                "inner_validation_IDs": builder.data.loc[val, "ID"].tolist(), "history": history})
    return best


def fit_equal_mean(x, y, apply):
    """Equal weights remain fixed; development-only nonnegative Cox calibration."""
    a, b = x.mean(axis=1).to_numpy(), apply.mean(axis=1).to_numpy()
    t, e = y.rfs_time_months.to_numpy(), y.rfs_event.to_numpy()
    def objective(beta):
        r = beta * a
        loss = 0.0
        for at in np.unique(t[e == 1]):
            death = (e == 1) & (t == at)
            loss -= r[death].sum() - death.sum() * logsumexp(r[t >= at])
        return loss / e.sum() + .05 * beta ** 2 / 2
    fit = minimize_scalar(objective, bounds=(0, 50), method="bounded")
    if not fit.success or fit.x > 49.9:
        raise ValueError("Equal-rank mean calibration did not converge within fixed bounds")
    beta = 0.0 if objective(0) <= fit.fun else float(fit.x)
    p, baseline = breslow_probability(beta * a, y, beta * b)
    return (b, percentile(a, b), p), {"slope": beta, "ridge_penalty": .05, "baseline": baseline,
                                                     "weights": [0.25] * 4}


def pnm_predictions(clinical, data, fitting, apply, pnm_spec):
    results, artifacts, audit = {}, {}, {}
    for name, columns in pnm_spec.items():
        fit, test, encoding = {}, {}, {}
        for c in columns:
            a, b = clinical.loc[fitting, c], clinical.loc[apply, c]
            spec = SCHEMA["variables"][c]
            state = fit_category(a, c, spec)
            a, train_audit = apply_category(a, c, spec, state)
            b, apply_audit = apply_category(b, c, spec, state)
            levels = state["observed_codes"]
            encoding[c] = {"levels": levels, "reference": float(levels[0]), "fitting_state": state,
                           "training_handling": train_audit, "application_handling": apply_audit}
            for level in levels[1:]:
                key = f"{c}=={level}"
                fit[key], test[key] = a.eq(level).astype(float), b.eq(level).astype(float)
        x, v = pd.DataFrame(fit, index=fitting), pd.DataFrame(test, index=apply)
        if x.shape[1] == 0:
            p, base = breslow_probability(np.zeros(len(fitting)), data.loc[fitting, OUTCOMES], np.zeros(len(apply)))
            results[name] = (np.zeros(len(apply)), np.full(len(apply), .5), p)
            artifacts[name] = {"baseline_only": base}
            audit[name] = {"encoding": encoding, "baseline_only": True}
        else:
            model = SparseCox().fit(x, data.loc[fitting, OUTCOMES], penalty=.01, l1=0, select=False)
            results[name] = model.predict(v)
            artifacts[name] = {"model": model, "encoding": encoding}
            audit[name] = {"encoding": encoding, "cox": model.audit}
    return results, artifacts, audit


def fit_partition(builder, clinical, fitting, apply, label, seed, epochs, patience, inner_folds, pnm_spec):
    print(f"[{label}] select epochs using development-only inner holdout", flush=True)
    selected_epochs = tune_epochs(builder, fitting, seed, label, epochs, patience, inner_folds)
    print(f"[{label}] regenerate inner OOF bases; refit attention for {selected_epochs} epochs", flush=True)
    x, v, predictions, full_base = builder.stack(fitting, apply, seed + 2, label, inner_folds)
    y = builder.data.loc[fitting, OUTCOMES]
    score = SparseCox().fit(x, y, cap=4)
    predictions["score_fusion"] = score.predict(v)
    predictions["equal_rank_mean"], equal_audit = fit_equal_mean(x, y, v)
    net, predictions["attention_fusion"], weights, reference, baseline = fit_predict(
        x, y, v, selected_epochs, seed + 3)
    pnm, artifacts, pnm_audit = pnm_predictions(clinical, builder.data, fitting, apply, pnm_spec)
    predictions.update(pnm)
    model_dir = builder.out / "model_artifacts"
    joblib.dump({"score_fusion": score, "equal_rank_mean": equal_audit, "pnm": artifacts},
                model_dir / f"{label}_fusion_pnm.joblib", compress=3)
    torch.save(net.state_dict(), model_dir / f"{label}_attention.pt")
    np.save(model_dir / f"{label}_attention_reference.npy", reference)
    write_json(model_dir / f"{label}_fusion_audit.json", {
        "fitting_IDs": builder.data.loc[fitting, "ID"].tolist(), "application_IDs": builder.data.loc[apply, "ID"].tolist(),
        "epochs": selected_epochs, "score_fusion": score.audit, "equal_rank_mean": equal_audit,
        "attention_baseline": baseline, "pnm": pnm_audit,
        "attention_identity": {"order": DIMENSIONS, "position_code": "fixed_4x4_identity_added_before_ReLU",
                               "trainable_parameters": sum(t.numel() for t in net.parameters())},
        "seed_schedule": {"partition": seed, "inner_holdout": seed, "early_stop_base_folds": seed + 1,
                          "stack_base_folds": seed + 2, "attention_refit": seed + 3}})
    if label == "final":
        diagnostic_models = {**full_base.models, "score_fusion": score,
                             **{k: a["model"] for k, a in artifacts.items() if "model" in a}}
        diagnostic_rows = []
        for name, fitted in diagnostic_models.items():
            if fitted.model is None:
                diagnostic_rows.append({"model": name, "diagnostic_status": "baseline_only_no_covariate_PH_test"})
                continue
            try:
                with warnings.catch_warnings(record=True) as captured:
                    test = proportional_hazard_test(fitted.model, fitted.training_frame, time_transform="rank")
                table = test.summary.reset_index().rename(columns={"index": "term"})
                table["model"] = name
                table["warning"] = " | ".join(str(w.message) for w in captured)
                diagnostic_rows.extend(table.to_dict("records"))
            except (ValueError, ArithmeticError) as error:
                diagnostic_rows.append({"model": name, "diagnostic_error": str(error)})
        pd.DataFrame(diagnostic_rows).to_csv(builder.out / "cox_PH_diagnostics_exploratory.csv", index=False)
    x.assign(ID=builder.data.loc[x.index, "ID"]).to_csv(model_dir / f"{label}_inner_oof_scores.csv", index=False)
    return predictions, weights


def run_models(data, sources, clinical, out, config):
    builder = Builder(sources, data, out, config)
    result = data[KEY].copy()
    weights = data[KEY].copy()
    development_group = config["development_group"]
    n_outer_folds = int(config.get("outer_folds", 5))
    master_seed = int(config.get("master_seed", SEED))
    dev = data.index[data.group.eq(development_group)].to_numpy()
    external = data.index[data.group.ne(development_group)].to_numpy()
    partitions = [(f"outer_{f}", data.index[data.group.eq(development_group) & data.cv_fold.ne(f)].to_numpy(),
                   data.index[data.group.eq(development_group) & data.cv_fold.eq(f)].to_numpy(), master_seed + 100 * f) for f in range(n_outer_folds)]
    final_offset = int(config.get("final_seed_offset", 900))
    if final_offset in {100 * f for f in range(n_outer_folds)}:
        raise ValueError("final_seed_offset collides with an outer-fold seed")
    partitions.append(("final", dev, external, master_seed + final_offset))
    for label, fitting, apply, seed in partitions:
        predictions, w = fit_partition(builder, clinical, fitting, apply, label, seed,
                                       config.get("attention_epochs_max", 400),
                                       config.get("attention_patience", 60),
                                       config.get("inner_folds", 4),
                                       config.get("pnm_models", DEFAULT_PNM))
        for name, (raw, rank, p) in predictions.items():
            if not all(np.isfinite(a).all() for a in (raw, rank, p)) or np.any((p < 0) | (p > 1)):
                raise ValueError(f"Invalid predictions: {label}/{name}")
            result.loc[apply, f"risk_{name}"] = raw
            result.loc[apply, f"percentile_{name}"] = rank
            result.loc[apply, f"prob36_{name}"] = p
        weights.loc[apply, DIMENSIONS] = w
        result.to_csv(Path(out) / "predictions_INCOMPLETE.csv", index=False)
    if result.isna().any().any():
        raise ValueError("Missing final predictions")
    result.to_csv(Path(out) / "predictions.csv", index=False)
    weights.to_csv(Path(out) / "attention_weights.csv", index=False)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    package = Path(__file__).resolve().parent
    config_path = Path(args.config).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if not config.get("clinical_schema"):
        raise ValueError("Specify a private clinical_schema in the config")
    set_schema((config_path.parent / config["clinical_schema"]).resolve())
    if config["attention_learning_rate"] != LEARNING_RATE or config["attention_weight_decay"] != WEIGHT_DECAY:
        raise ValueError("Recorded attention optimizer settings differ from the implementation")
    for filename, expected in config.get("code_sha256", {}).items():
        if digest(package / filename) != expected:
            raise ValueError(f"Code differs from locked package: {filename}")
    paths = {}
    for name, spec in config["inputs"].items():
        if spec["location"] not in {"root", "package"}:
            raise ValueError(f"Invalid input location for {name}")
        root = package if spec["location"] == "package" else args.input_root.resolve()
        path = root / spec["relative_path"]
        expected_hash = spec.get("sha256")
        if expected_hash and digest(path) != expected_hash:
            raise ValueError(f"Input checksum mismatch: {name}: {path}")
        paths[name] = path
    required = {"manifest", "clinical", "ct_p_global", "t2wi_global", "ct_p_habitat",
                "t2wi_habitat", "ct_p_2p5d", "t2wi_2p5d"}
    if set(paths) != required:
        raise ValueError(f"Inputs must contain exactly {sorted(required)}")
    manifest_path = paths["manifest"]
    data = pd.read_csv(manifest_path, dtype={"ID": str})
    validate_manifest(data, config)
    sources, input_audit = read_sources(paths, data)
    clinical = pd.read_csv(paths["clinical"], dtype={"ID": str}).set_index("ID").loc[data.ID].reset_index(drop=True)
    versions = {name: importlib.metadata.version(name) for name in ["numpy", "pandas", "scipy", "scikit-learn", "lifelines", "torch", "joblib"]}
    print(json.dumps({"preflight": "PASS", "cohorts": data.group.value_counts().to_dict(), "versions": versions}), flush=True)
    if args.preflight_only:
        return
    out = args.outdir.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output directory must be new/empty; never mix old results")
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy2(manifest_path, out / "manifest_used.csv")
    shutil.copy2(paths["clinical"], out / "clinical_matrix_used.csv")
    shutil.copytree(package, out / "code_snapshot", ignore=shutil.ignore_patterns("__pycache__", "*.zip"))
    shutil.copy2(config_path, out / "config_used.json")
    if config.get("clinical_schema"):
        shutil.copy2((config_path.parent / config["clinical_schema"]).resolve(), out / "clinical_schema_used.json")
    write_json(out / "run_provenance.json", {"inputs": input_audit, "analysis_lock": config,
               "python": sys.version, "platform": platform.platform(), "versions": versions})
    result = run_models(data, sources, clinical, out, config)
    evaluate(result, out, config["bootstrap_replicates"], config["development_group"],
             list(data.group.drop_duplicates()), int(config.get("master_seed", SEED)))
    (out / "RUN_COMPLETE.txt").write_text("Nested training and 36-month evaluation completed. Review diagnostics and clinical assumptions before manuscript use.\n", encoding="utf-8")
    print(f"[complete] {out}", flush=True)


if __name__ == "__main__":
    main()
