# r_pipeline — Reverse-MIDAS electricity price forecasting (R)

R pipeline for forecasting daily electricity prices in Germany and Italy with
monthly macroeconomic variables, using the Reverse-MIDAS framework of Foroni,
Guerin and Marcellino (2018, IJF) with an exponential-Almon restriction on
the autoregressive block over the lag set J = {1, 2, 3, 7}.

## Layout

```
data/                 input CSVs (daily prices + monthly macro, DE and IT)
midas_codes/          the pipeline
  common_utils.R        shared utilities (data, release-aware LF alignment,
                        calendar/holidays, DM-HLN, exp-Almon, forecasting)
  fit_r_midas_fast.R    FWL-accelerated concentrated NLS estimator
  ar_benchmark.R        AR(1,2,3,7) benchmark family (ar_dum = benchmark)
  r_midas.R             standard R-MIDAS (Almon on price lags only)
  r_midas_extended.R    extended R-MIDAS (Almon also on fuel blocks)
  ru_midas.R            unrestricted RU-MIDAS (28 position dummies)
  summary_table*.R      metric tables (RMSE/MAE, asinh, bias-corrected, DM-HLN)
  combination_subperiod_gr.R  forecast combinations, subperiod evaluation,
                        Giacomini-Rossi (2010) fluctuation test
  optimizer_comparison_r_midas.R  optimizer diagnostic (optional)
  run_all.R             stage-selectable orchestrator
  run_ar_r_midas.R      AR + standard R-MIDAS subset driver
  data/                 input CSVs used by the scripts
  not useful files/     superseded legacy scripts (kept for reference)
```

## How to run

```bash
cd midas_codes
Rscript run_all.R                      # default: ar -> r_midas -> summary -> combination
COUNTRY_NAME=Italy SPEC_LIMIT=2 N_CORES=1 Rscript run_ar_r_midas.R   # smoke test
STAGES="ar,r_midas,ru_midas,summary_all" Rscript run_all.R           # other stages
```

Env knobs: `COUNTRY_NAME` (Germany | Italy | ALL), `N_CORES`, `SPEC_LIMIT`,
`WINDOW_DAYS`, `RMIDAS_DISCOUNT`, `LF_RELEASE_MODE` (see
`RELEASE_DATE_ASSUMPTIONS.md` in `midas_codes/` for the publication-date
conventions of the low-frequency variables).

Outputs are written to `midas_codes/forecasts/` and `midas_codes/results/`;
both are gitignored because they are large and fully regenerable.

A line-by-line Python port of this pipeline lives in the companion
`python_pipeline/` folder at the root of this repository.
