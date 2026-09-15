# Dependencies and checks

The helpers use R or Python and the following module-level dependencies:

| Module | Dependencies |
| --- | --- |
| Single-cell preprocessing | Seurat; harmony when enabled |
| Composition/permutation and repertoire metrics | Base R |
| Paired bulk and enrichment | edgeR, fgsea |
| Protein transformations | Base R |
| Spatial scoring and extent | NumPy, SciPy |

Install only the dependencies needed for your selected helper, following their upstream documentation. No package installation runs automatically. The source workflows span different environments; this release does not supply a complete lockfile or claim that all stages share a single validated environment.

All R files were parsed with R 4.5.1. Selected base-R helpers and the Python spatial helper were checked on artificial small inputs; the latter used Python 3.10.20, NumPy 2.2.6 and SciPy 1.15.3. These checks validate program behaviour only, not study findings. Dependency-heavy Seurat/Harmony and edgeR/fgsea workflows were not numerically rerun for this release. The source analyses and article figures have not been fully rerun.
