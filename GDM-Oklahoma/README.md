# GDM-Oklahoma: Adjusted AI/AN Gestational Diabetes Prevalence in Oklahoma

## Overview

This project develops a Bayesian spatial model (BEACON) that jointly integrates:
- **OSDH surveillance data** — statewide AI/AN GDM counts, subject to underdetection
- **Cherokee Nation clinical data** — gold-standard GDM diagnoses within CN service counties

The model estimates county-level *adjusted* AI/AN GDM prevalence across all 77 Oklahoma counties, correcting for differential sensitivity in OSDH reporting and incorporating Cherokee Nation detection/coverage probability.

## Data

Original data cannot be shared publicly due to Tribal data governance requirements. A data use agreement with Cherokee Nation Public Health is required to access the CN clinical data.

- **Cherokee Nation Public Health**: [cherokee.org/services/health](https://www.cherokee.org/services/health)

## Demo Script

`demo_aian_gdm_oklahoma.R` replicates the full analysis pipeline using **synthetic data** generated with parameters fixed at the posterior means from Peterson et al. (2025). It requires no external data files.

### Requirements

```r
install.packages(c("nimble", "coda", "tidyverse", "tigris", "sf", "spdep", "MASS"))
```

### Usage

```r
source("demo_aian_gdm_oklahoma.R")
```

Runtime: ~20–30 minutes (3 chains × 10,000 iterations). Adjust `n_chains` / `n_iter` in Section 7 to shorten.

### Outputs (saved to `output/`)

| File | Description |
|------|-------------|
| `demo_results_table.csv` | Posterior summaries and convergence diagnostics |
| `demo_map_rho_true.pdf` | Adjusted AI/AN GDM prevalence map |
| `demo_map_phi.pdf` | OSDH sensitivity surface map |
| `demo_map_kappa.pdf` | CN detection/coverage probability map |

## Citation

Peterson et al. (2025). *Journal of the Royal Statistical Society Series C.*
