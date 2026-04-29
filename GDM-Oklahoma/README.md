# GDM-Oklahoma: Adjusted AI/AN Gestational Diabetes Prevalence in Oklahoma

## Overview

This project develops a Bayesian spatial model (BEACON) that jointly integrates:
- **OSDH surveillance data** — statewide AI/AN GDM counts, subject to underdetection
- **Cherokee Nation clinical data** — gold-standard GDM diagnoses within CN service counties

The model estimates county-level *adjusted* AI/AN GDM prevalence across all 77 Oklahoma counties, correcting for differential sensitivity in OSDH reporting and incorporating Cherokee Nation detection/coverage probability. See [`vignette_beacon_gdm_oklahoma.qmd`](vignette_beacon_gdm_oklahoma.qmd) for a full walkthrough of the model and its outputs.

## Data

Original data cannot be shared publicly due to Tribal data governance requirements. A data use agreement with Cherokee Nation Public Health is required to access the CN clinical data.

- **Cherokee Nation Public Health**: [cherokee.org/services/health](https://www.cherokee.org/services/health)

## Files

| File | Description |
|------|-------------|
| `demo_aian_gdm_oklahoma.R` | Full demo: synthetic data generation + NIMBLE model fit (~20–30 min) |
| `make_example_outputs.R` | Quick demo: 1 chain × 1,500 iterations, generates `output/` in ~3–5 min |
| `vignette_beacon_gdm_oklahoma.qmd` | Pedagogical walkthrough of the model with rendered maps and tables |

## Quick Start

**Option A — generate outputs quickly (~3–5 min), then render the vignette:**
```r
source("make_example_outputs.R")
quarto::quarto_render("vignette_beacon_gdm_oklahoma.qmd")
```

**Option B — run the full demo (~20–30 min):**
```r
source("demo_aian_gdm_oklahoma.R")
```

### Requirements

```r
install.packages(c("nimble", "coda", "tidyverse", "tigris", "sf",
                   "spdep", "MASS", "knitr", "quarto"))
```

## Outputs (saved to `output/`)

| File | Description |
|------|-------------|
| `demo_results_table.csv` | Posterior summaries and convergence diagnostics |
| `demo_county_posteriors.csv` | County-level posterior means and 95% CrIs (used by vignette) |
| `demo_map_rho_true.pdf` | Adjusted AI/AN GDM prevalence map |
| `demo_map_phi.pdf` | OSDH sensitivity surface map |
| `demo_map_kappa.pdf` | CN detection/coverage probability map |

## Citation

Peterson et al. (2025). *Journal of the Royal Statistical Society Series C.*
