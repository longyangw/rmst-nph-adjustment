# rmst-nph-adjustment

R code accompanying the manuscript:

> **Wu, L. (2026).** Restricted Mean Survival Time Under Non-Proportional Hazards: How Much Power Can Covariate Adjustment Recover? *Submitted to Statistics in Medicine.*

## Overview

This repository contains:

- **`rmst_diagnostic.R`** — a design-stage diagnostic tool that computes
  the benchmark quantity $C(\tau)$, the signal centroid $\bar{t}_\delta$,
  the signal-retention ratio (SRR), and Proposition 3's upper bound, and
  recommends RMST versus log-rank based on the proportional-hazards (PH)
  reference range.
- Simulation scripts used to produce Tables 2--4 in the paper.

## Quick start

```r
source("rmst_diagnostic.R")

# Delayed benefit (DGM1 in the paper)
rmst_diagnostic(type = "delayed",
                median0 = 12, HR = 0.55,
                tau = 24, delay = 4)

# Crossing hazards (DGM2)
rmst_diagnostic(type = "crossing",
                median0 = 12,
                HR_early = 1.30, HR_late = 0.50,
                tau = 24, changepoint = 4)

# Gamma frailty (DGM4)
rmst_diagnostic(type = "frailty",
                median0 = 12, HR_cond = 0.70,
                tau = 24, theta = 0.5)

# Exponential PH reference
rmst_diagnostic(type = "exponential",
                median0 = 12, HR = 0.70, tau = 24)
```

The function returns a list with `C_tau`, `Delta_tau`, `t_bar_delta`,
`tbar_ratio`, `SRR`, `t_star`, `SRR_upper_bound`, `regime`, and a
`recommendation` string.

## Output interpretation

The diagnostic compares $\bar{t}_\delta/\tau$ to the PH reference range
$[0.18, 0.28]$ (corresponding to clinically relevant HR values 0.50--0.90
under exponential survival; see Appendix A.6 of the paper):

| $\bar{t}_\delta/\tau$ | Regime | RMST vs log-rank |
|------------------------|--------|------------------|
| Within or below PH range (SRR $\ge 0.72$) | (i)  | RMST competitive; interpretability advantage at no power cost |
| Above PH range (SRR $< 0.72$)             | (ii) | RMST faces structural disadvantage; log-rank or WLR preferred |
| $C(\tau) < 0$                              | (iii) | RMST structurally favored over flat benchmark |

## Repository layout

```
rmst-nph-adjustment/
├── README.md                # this file
├── rmst_diagnostic.R        # design-stage diagnostic function
├── simulation/              # full simulation scripts
│   └── ...
└── LICENSE
```

## Reproducibility

Simulation results in the paper were generated with R 4.x using packages
`survival`, `pseudo`, and base R. To reproduce Tables 2--4:

```bash
Rscript simulation/run_main.R
```

The simulation generates synthetic data only; no patient data were used.

## Citation

If you use this code, please cite:

```bibtex
@article{Wu2025RMST,
  author  = {Wu, Longyang},
  title   = {Restricted Mean Survival Time Under Non-Proportional
             Hazards: How Much Power Can Covariate Adjustment Recover?},
  journal = {Statistics in Medicine},
  year    = {2025},
  note    = {Submitted}
}
```

## Contact

Longyang Wu — `lwu@uwaterloo.ca`

## License

MIT
