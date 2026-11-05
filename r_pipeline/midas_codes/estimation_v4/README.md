# estimation_v4 — v4 data (real weekends), normalised asinh, IESI, ARX + R-MIDAS + R-MIDAS-L7, DM/GW tests, density, CSSED

Full pipeline run (`run_all.R`, default stages) implementing README_v4.md sections 1, 2.1–2.6
and 2.8 (2.7 `DROP_IMPUTED` not run). Started 22 Oct 2026.

**Data** (`data -> ../data_v4`, rebuilt 22 Oct 2026 with `build_daily_prices.py` + `build_dataset_v4.py`,
epftoolbox `DE.csv` included):
* Germany: real day-ahead prices every day from **2012-01-09** (epftoolbox EPEX-DE to 2015-01-04,
  ENTSO-E after; MAD on the 2015–2017 overlap 0.002); 2005-01-01 → 2012-01-08 Datastream weekdays with
  weekends carried from Friday (estimation windows only; every evaluated target is real).
  `elec_source`: real_all_days 5106, datastream_weekday 1830, carried_forward 734.
* Italy: GME PUN every day; monthly survey = **iesi_it** (ISTAT economic sentiment; `LF_IT_SURVEY=istat_bci_it`
  for the manufacturing-climate robustness run).
* No interpolation anywhere; monthly PMI gaps carried forward.

**Methods**
* Normalised asinh inside each rolling window (`VST_MODE=normalized`): y = asinh((P − a)/b), a = median, b = MAD
  of the window; forecasts a + b sinh(ŷ) (median) and a + b sinh(ŷ) e^{s²/2} (mean, `__bc`); (a, b) per origin in
  `forecasts/vst_params_<cc>.csv`.
* exp-Almon box [−8, 8]² (θ₂ free), 7 starts incl. two U-shapes; θ̂, implied φ_j, boundary / identification flags per window.
* Families: AR, ARX (free lags {1,2,3,7}), R_MIDAS (Almon on {1,2,3,7}), R_MIDAS_L7 (Almon on 1..7), real-time
  combinations of each family.
* Tests: DM-HLN and Giacomini–White (instruments 1, d_{t−h}) on squared / absolute errors, EUR and transformed scale,
  full sample and pre / crisis / post — vs the AR benchmark (`results/tests_vs_benchmark_<cc>.csv`, Table A) and vs
  the daily-only twin (`results/tests_vs_twin_<cc>.csv`, Table B; `tests_vs_twin_summary.csv`; `tests_l7_*` subsets).
* Density: log score, CRPS, ACPS (c = 0.05 … 0.95), PIT / 90% coverage; DM + GW vs benchmark
  (`density_scores_<cc>.csv`) and vs twin (`density_tests_vs_twin_<cc>.csv`).
* CSSED plots: `results/cssed/cssed_<cc>_h<hh>.png` (+ csv).
* `results/ols_coef_check_<cc>.csv`, `ols_ar7_check_<cc>.csv`: plain-OLS coefficient checks (HAC).

Scripts are symlinks to `../`; outputs live here. `estimation_v3/` (old data, plain asinh) is untouched.
