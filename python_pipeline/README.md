# Python pipeline — port of the R forecasting code

This folder is a self-contained Python port of the R scripts in
`r_pipeline/midas_codes/`. It reads the **same data files**, follows the
**same methodology** (identical estimators, identical forecasting
conventions) and writes **CSV outputs with the same layout**, so results are
directly comparable with the R runs and the downstream summary scripts work
on either set of outputs.

## Requirements

Python 3.9+ with `numpy`, `pandas`, `scipy` (`pip install numpy pandas scipy`).

## Files

| file | port of (R) | what it does |
|---|---|---|
| `common_utils.py` | `common_utils.R` + `fit_r_midas_fast.R` | data loading and transforms (asinh / log / dlog / z-score), publication-aware LF alignment (surveys LOCF), calendar dummies and country holidays, rolling origins, DM-HLN test, exp-Almon weights, FWL-accelerated concentrated NLS, iterated multi-step forecasting (option c.1) |
| `ar_benchmark.py` | `ar_benchmark.R` | AR(1,2,3,7) + dummies benchmark family (`ar_dum`, `ar_dum_gas`, `ar_dum_brent`, `ar_dum_gas_brent`) |
| `r_midas.py` | `r_midas.R` | standard R-MIDAS: Almon on the 14 AR lags of y only, gas/brent linear at lag 1, 256-spec grid |
| `r_midas_extended.py` | `r_midas_extended.R` | extended R-MIDAS: Almon on y AND on each included fuel block |
| `r_midas_surveys_h.py` | `r_midas_surveys_h.R` | older pooled variant (month dummies only, single warm-started optimisation) |
| `ru_midas.py` | `ru_midas.R` | unrestricted RU-MIDAS with 28 within-month position dummies (Cimadomo / Foroni-Rossini style) |
| `summary_table.py` | `summary_table.R` | unified metric table (RMSE / MAE / asinh / bias-corrected blocks, DM-HLN, spec dictionaries) |
| `combination_subperiod_gr.py` | `combination_subperiod_gr.R` | forecast combinations (equal / trimmed / real-time inverse-MSE), subperiod ratios (pre / crisis / post), Giacomini-Rossi fluctuation test |
| `run_all.py` | `run_all.R` | orchestrator |

## How to run

From inside this folder:

```bash
python3 run_all.py                                  # everything, both countries
COUNTRY_NAME=Italy python3 r_midas.py               # one model, one country
STAGES=ar,r_midas,summary python3 run_all.py        # subset of stages
```

Data are read from `../r_pipeline/midas_codes/data/` by default
(override with `V2_DATA_DIR`). Outputs go to `python_pipeline/forecasts/` and
`python_pipeline/results/`, so nothing overwrites the R outputs.

Environment knobs (same names and meanings as the R pipeline):

| variable | default | meaning |
|---|---|---|
| `COUNTRY_NAME` | `ALL` | `Germany`, `Italy` or `ALL` |
| `SPEC_LIMIT` | 256 | truncate the spec grid |
| `N_CORES` | cpu−1 | worker processes (parallel over specs) |
| `WINDOW_DAYS` | 3650 | rolling window length |
| `RMIDAS_DISCOUNT` | 1 | discounted WLS lambda (1 = off) |
| `ORIGIN_LIMIT` | off | **python-only**: use only the first N rolling origins (quick tests / validation) |

## Equivalence with R — validation results (June 2026)

- **AR benchmark** (pure OLS): validated against the existing R CSVs on the
  first 250 Italian origins, all four specs, h = 1 and h = 28: max absolute
  difference **2e-11 EUR/MWh** (machine precision). This jointly validates
  the shared infrastructure: data transforms, publication-aware LF
  alignment, calendar dummies, Easter/holiday sets, rolling-origin indexing,
  the iterated multi-step recursion and the sinh back-transform.
- **R-MIDAS** (concentrated NLS): validated on the first 300 Italian origins
  for specs 000-003: forecast correlations 0.997-1.000, aggregate RMSE
  matching to 2-4 decimals; spec_003 matched to 2e-5. The remaining
  differences are optimizer-path effects: the existing R CSVs were produced
  by the OLD R protocol (uniform random multistart, R's L-BFGS-B), whereas
  this port follows the CURRENT R protocol (warm start + deterministic
  canonical shapes, SciPy L-BFGS-B with `ftol = factr * eps`). On windows
  where the concentrated SSE is flat in theta, individual forecasts can
  differ visibly; aggregates do not.
- **RU-MIDAS** (pure OLS): the port is internally exact (the position-based
  coefficient lookup in the iterated forecast reproduces the design-row
  product to 0.0), but it differs from the existing
  `ru_midas_it_h*.csv` files by ~0.4-1 EUR/MWh median. Those CSVs were
  written on May 27, BEFORE the May 29 revision of `common_utils.R` that the
  later AR/R-MIDAS runs (which match to 1e-11) used; the largest single
  deviation falls on Epiphany, pointing to a holiday/calendar change between
  the two R versions. If exact cross-checks on RU-MIDAS are needed,
  re-run `ru_midas.R` first to refresh its CSVs with the current
  `common_utils.R`.
- The 256 `spec_xxx` identifiers encode the same variable subsets in both
  pipelines (the spec number is the integer value of the flag vector
  `gas, brent, lf1..lf6` with `gas` as the most significant bit).

## Runtime notes

The NLS families use the same FWL acceleration as the updated R code: the
cross products are computed once per rolling window and each objective
evaluation is O(K^2 Q^2), independent of the window length. With all 256
specs and ~3,980 origins expect roughly 0.5–2 h per country per family on 8
cores; use `SPEC_LIMIT` / `ORIGIN_LIMIT` for quick experiments. `ru_midas.py`
is the slowest family (large saturated OLS designs), exactly as in R.
