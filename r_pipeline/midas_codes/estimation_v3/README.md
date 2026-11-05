# estimation_v3 — theta2 free, ARX free-lag control, density evaluation

Run of the full v2 pipeline (`run_all.R`, default stages) with the settings
introduced in October 2026:

* **exp-Almon box** `theta in [-8, 8]^2` (`RMIDAS_THETA2_MAX=8`, now the default):
  the former `theta2 <= 0` restriction is dropped so U-shaped lag profiles
  (large lag 1, dip at lags 2-3, weekly revival at lag 7) are admissible.
  7-start L-BFGS-B per window (warm start + flat, two U-shapes, fast/slow decay, hump).
* **theta diagnostics** saved per rolling window in `results/coefficients/r_midas_<cc>.csv`:
  theta1, theta2, implied phi_1 phi_2 phi_3 phi_7 = delta * b(j; theta), min root modulus,
  at_bound, flag_delta_small, flag_degenerate, optim convergence code.
* **ols_check stage** (`ols_coef_check.R`): plain-OLS coefficients of every regressor
  (full sample / subperiods / rolling, Newey-West HAC) -> `results/ols_coef_check_<cc>.csv`.
* **ARX family** (`arx_free_lags.R`): same 256 specs with the four AR coefficients free (OLS);
  control separating the value of the monthly indicators from the lag restriction.
* **density stage** (`density_eval.R`): log predictive likelihood, CRPS, ACPS
  (Iacopini-Ravazzolo-Rossini, c = 0.05/0.25/0.5/0.75/0.95) and PIT calibration of the
  Gaussian-asinh predictive densities -> `results/density_scores_<cc>.csv`,
  `results/density_top5_per_country_h.csv`. Predictive sds in `forecasts/*_sd.csv`.
* Window 3650 days, horizons 1/7/14/21/28, data as of 2026-05-15 (weekend fix NOT yet applied).

Scripts are symlinks to `../`; outputs (`forecasts/`, `results/`, log) live here.
Baseline (theta2 <= 0) results remain in the parent folder; `../theta2relaxed/` holds the
first relaxed run (5 starts, no ARX / density).
