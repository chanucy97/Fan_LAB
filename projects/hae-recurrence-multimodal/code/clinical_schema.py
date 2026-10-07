"""Explicit clinical schema and training-only category handling; no name heuristics."""
import json
from pathlib import Path

import numpy as np
import pandas as pd


def load_schema(path=None):
    schema_path = Path(path) if path is not None else Path(__file__).with_name("clinical_schema.example.json")
    schema = json.loads(schema_path.read_text(encoding="utf-8-sig"))
    if not schema.get("variables"):
        raise ValueError("Clinical schema must declare at least one predictor")
    for name, spec in schema["variables"].items():
        if spec["type"] not in {"continuous", "categorical"}:
            raise ValueError(f"Invalid schema type: {name}")
        if spec["type"] == "categorical":
            codes = spec["allowed_codes"]
            if not codes or len(codes) != len(set(codes)):
                raise ValueError(f"Invalid category domain: {name}")
    return schema


SCHEMA = load_schema()


def set_schema(path):
    """Replace the in-memory schema for a user-supplied generic run."""
    SCHEMA.clear()
    SCHEMA.update(load_schema(path))


def numeric_checked(values, name, spec):
    values = pd.to_numeric(values, errors="raise").astype(float)
    if np.isinf(values.to_numpy()).any():
        raise ValueError(f"Infinite clinical values: {name}")
    if spec["type"] == "categorical":
        invalid = values.notna() & ~values.isin(spec["allowed_codes"])
        if invalid.any():
            counts = values.loc[invalid].value_counts().to_dict()
            raise ValueError(f"Invalid clinical codes for {name}: {counts}; allowed={spec['allowed_codes']}")
    return values


def validate_clinical(frame):
    expected = set(SCHEMA["variables"]) | {"ID"}
    missing, extra = expected - set(frame), set(frame) - expected
    if missing or extra or frame.columns.duplicated().any():
        raise ValueError(f"Clinical schema columns mismatch; missing={sorted(missing)}, extra={sorted(extra)}")
    audit = {}
    for c, spec in SCHEMA["variables"].items():
        numeric = numeric_checked(frame[c], c, spec)
        audit[c] = {"type": spec["type"], "missing": int(numeric.isna().sum())}
        if spec["type"] == "categorical":
            audit[c]["allowed_codes"] = spec["allowed_codes"]
    return audit


def category_specs(prefix="clinical::"):
    return {prefix + c: spec for c, spec in SCHEMA["variables"].items() if spec["type"] == "categorical"}


def fit_category(values, name, spec):
    a = numeric_checked(values, name, spec).dropna()
    if a.empty:
        raise ValueError(f"Categorical predictor has no observed fitting values: {name}")
    # pandas mode is sorted: smallest numeric code resolves a tie reproducibly.
    return {"mode": float(a.mode().iloc[0]), "observed_codes": sorted(float(v) for v in a.unique()),
            "allowed_codes": spec["allowed_codes"], "unseen_policy": "training_mode_with_audit"}


def apply_category(values, name, spec, state):
    a = numeric_checked(values, name, spec)
    missing = a.isna()
    unseen = a.notna() & ~a.isin(state["observed_codes"])
    audit = {"missing_count": int(missing.sum()), "unseen_count": int(unseen.sum()),
             "unseen_code_counts": {str(k): int(v) for k, v in a.loc[unseen].value_counts().items()},
             "replacement_training_mode": state["mode"], "rule": state["unseen_policy"]}
    return a.mask(missing | unseen, state["mode"]), audit
