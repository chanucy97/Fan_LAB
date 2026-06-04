from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create synthetic, non-patient example TBNK data.")
    parser.add_argument("--output", default="data/example/synthetic_tbnk_example.csv")
    parser.add_argument("--n", type=int, default=180)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def clipped_normal(rng: np.random.Generator, mean: float, sd: float, n: int, low: float = 0.0) -> np.ndarray:
    return np.clip(rng.normal(mean, sd, n), low, None)


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    groups = rng.choice(["LN", "NS", "MN", "IgAN", "AAV"], size=args.n, p=[0.24, 0.34, 0.18, 0.10, 0.14])
    years = rng.choice([2024, 2025, 2026], size=args.n, p=[0.36, 0.42, 0.22])
    rows = []
    for idx, group in enumerate(groups, start=1):
        is_ln = group == "LN"
        cd4_pct = clipped_normal(rng, 28 if is_ln else 39, 8, 1)[0]
        cd8_pct = clipped_normal(rng, 41 if is_ln else 30, 7, 1)[0]
        nk_pct = clipped_normal(rng, 7 if is_ln else 13, 5, 1)[0]
        b_pct = clipped_normal(rng, 13 if is_ln else 11, 5, 1)[0]
        cd3_pct = np.clip(cd4_pct + cd8_pct + rng.normal(3, 4), 35, 90)
        lymph_abs = clipped_normal(rng, 1350 if is_ln else 1950, 520, 1, 200)[0]
        cd4_abs = clipped_normal(rng, 430 if is_ln else 760, 240, 1, 20)[0]
        cd8_abs = clipped_normal(rng, 610 if is_ln else 570, 230, 1, 20)[0]
        nk_abs = clipped_normal(rng, 95 if is_ln else 220, 120, 1, 5)[0]
        b_abs = clipped_normal(rng, 190 if is_ln else 210, 110, 1, 5)[0]
        cd3_abs = cd4_abs + cd8_abs + clipped_normal(rng, 60, 50, 1, 0)[0]
        rows.append(
            {
                "record_id": f"synthetic_{idx:04d}",
                "analysis_group": group,
                "year": years[idx - 1],
                "cd3_pct": round(cd3_pct, 3),
                "cd4_pct": round(cd4_pct, 3),
                "cd8_pct": round(cd8_pct, 3),
                "b_pct": round(b_pct, 3),
                "nk_pct": round(nk_pct, 3),
                "lymph_abs": round(lymph_abs, 3),
                "cd3_abs": round(cd3_abs, 3),
                "cd4_abs": round(cd4_abs, 3),
                "cd8_abs": round(cd8_abs, 3),
                "b_abs": round(b_abs, 3),
                "nk_abs": round(nk_abs, 3),
                "cd4_cd8_ratio": round(cd4_abs / cd8_abs, 4),
            }
        )
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Wrote synthetic example data to: {out}")


if __name__ == "__main__":
    main()
