# r-midas — Reverse-MIDAS electricity price forecasting

Twin R and Python pipelines for forecasting daily electricity prices in
Germany and Italy with monthly macroeconomic variables, using the
Reverse-MIDAS framework of Foroni, Guerin and Marcellino (2018, IJF) with an
exponential-Almon restriction on the autoregressive block over the lag set
J = {1, 2, 3, 7}.

The two pipelines implement the same methodology end to end (same data,
same estimators, same forecasting conventions, same CSV output layouts), so
results are directly comparable and either one can be used as the reference
implementation. The Python port has been validated against the R outputs:
the OLS families match to machine precision (~1e-11 EUR/MWh) and the
nonlinear R-MIDAS families agree up to optimizer path effects.

## Layout

```
r_pipeline/           the R pipeline (see r_pipeline/README.md)
  data/                 input CSVs (daily prices + monthly macro, DE and IT)
  midas_codes/          scripts, orchestrators, release-date assumptions
python_pipeline/      line-by-line Python port (see python_pipeline/README.md)
```

## Quick start

R (requires only base R + parallel):

```bash
cd r_pipeline/midas_codes
Rscript run_all.R          # default: ar -> r_midas -> summary -> combination
```

Python (requires numpy, pandas, scipy):

```bash
cd python_pipeline
python3 run_all.py         # same default stages
```

Both honour the same environment knobs: `COUNTRY_NAME` (Germany | Italy |
ALL), `N_CORES`, `SPEC_LIMIT`, `WINDOW_DAYS`, `RMIDAS_DISCOUNT`,
`LF_RELEASE_MODE`, `FAMILIES`, `STAGES` (and `ORIGIN_LIMIT`, Python only,
for quick partial runs).

Models covered: AR(1,2,3,7) benchmark family, standard R-MIDAS (Almon on the
price lags only, fuels linear at lag 1), extended R-MIDAS (Almon also on the
fuel blocks), unrestricted RU-MIDAS with 28 within-month position dummies,
plus summary tables (levels / asinh / MAE / feasible bias-corrected, DM-HLN
tests), real-time forecast combinations, subperiod evaluation and the
Giacomini-Rossi (2010) fluctuation test.

Generated outputs (`forecasts/`, `results/`) are gitignored in both
pipelines because they are large and fully regenerable.
