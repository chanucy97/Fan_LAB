# Running the analysis

## 1. Install dependencies

From this project directory:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

## 2. Generate synthetic example data

```bash
python src/make_example_data.py --output data/example/synthetic_tbnk_example.csv --n 180 --seed 2026
```

The generated file is synthetic and does not contain patient records.

## 3. Run the example analysis

```bash
python src/run_analysis.py --input data/example/synthetic_tbnk_example.csv --outdir outputs/example_run
```

## 4. Run on a private de-identified dataset

Keep private data outside the Git repository. For example:

```bash
python src/run_analysis.py --input <private_data_dir>/deidentified_tbnk.csv --outdir outputs/private_run
```

Before committing, remove all runtime outputs and check that no private path, patient identifier, clinical table, or mapping file has been staged.

## 5. Reproducibility notes

- The example analysis is intended to verify code execution and file structure.
- Published manuscript results should be reproduced only from the institution-approved analysis dataset.
- Patient-level data are not part of this GitHub project.
