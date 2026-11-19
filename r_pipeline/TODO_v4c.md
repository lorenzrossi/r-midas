# TODO v4c — code changes after the v3 rewrite of Chapter 2

The current setup (decided in November 2026; the German 2005–2011 data with Friday prices carried over
the weekend are kept as they are for now, so there is no imputed-day correction):
- Italy: log prices, results in `estimation_v4_it_log`.
- Germany: window-normalised asinh, results in `estimation_v4` and `estimation_v4_crps`.
- The three families (ARX, R-MIDAS, R-MIDAS-L7) contain only the **252 specifications with at least one monthly
  indicator**. The daily-only alternatives are the AR models with gas and/or Brent. There are no "twin"
  comparisons anywhere in the paper.
- The paper tables and figures are produced by the Python scripts in
  `Forecasting electricity prices/tables_v3/scripts/` (see section 5).

Do the items in this order. Item 0 changes how every future run is set up, so do it first.

---

## 0. FIRST: run only the MIDAS specifications with monthly indicators
The daily-only ("twin") specifications of the three families, i.e. the 4 per family with no monthly
indicator (`lf_vars` empty: `spec_000`, `spec_064`, `spec_128`, `spec_192`), are no longer used anywhere in
the paper and **can be discarded**. From now on, every run of the ARX, R-MIDAS and R-MIDAS-L7 families
estimates only the **252 specifications that contain at least one low-frequency indicator**. The daily-only
reference models are the AR models with gas and/or Brent (`ar_benchmark.R`), which stay as they are.

- [ ] Add a switch `SPEC_REQUIRE_LF` (default **1**) to the spec builder used by `arx_free_lags.R`,
      `r_midas.R` and `r_midas_l7.R`. After the spec dictionary is built, keep only the rows with a
      non-empty monthly set:
```r
SPEC_REQUIRE_LF <- as.integer(Sys.getenv("SPEC_REQUIRE_LF", "1"))
if (SPEC_REQUIRE_LF == 1L) specs <- specs[lengths(specs_lf) > 0]   # or: dict[dict$lf_vars != "", ]
```
      Keep the original `spec_id`s (do not renumber), so the old and new files stay comparable and
      `tables_v3/scripts` keep working.
- [ ] `combination_subperiod_gr.R`: the combinations then average the 252 specifications automatically
      (item 4a becomes a check rather than a change).
- [ ] `summary_table_ar_r_midas.R`, `tests_dm_gw.R`, `density_eval.R`: drop the "vs twin" outputs (Table B,
      `tests_vs_twin_*`, `tests_l7_vs_twin_*`, `density_tests_vs_twin_*`, `crps_twin_summary_*`). With
      `SPEC_REQUIRE_LF=1` there is no twin to compare with, so the code must skip these blocks instead of
      failing. Gate them behind `REPORT_TWINS=1` if you want to keep the option.
- [ ] Existing runs (`estimation_v4`, `estimation_v4_it_log`) do **not** need to be redone: the table scripts
      already drop the twins (`load()` in `common.py`) and recompute the combinations (`comb252.py`). Apply
      the switch to all **new** runs (item 3 and any future rerun).

---

## 1. Density forecasts: rescale the predictive variance (method 1)

### Idea
The model variance $\sigma^2_{h,t} = \sigma^2\sum_{s<h}\psi_s^2$ ignores two things: the uncertainty of the
parameters, and the movements of the regressors held fixed at the origin (gas, Brent, the monthly
indicators). The realised $h$-step errors of the model contain both. Rescale $\sigma_{h,t}$ with the ratio
of the realised to the model-implied second moment, computed **in real time**. Only errors whose target
date is on or before the origin are used. Do this for both italy and germany

For model $m$, horizon $h$ and origin $t$:

$$u_\tau = g(P_{\tau+h}) - \hat y_{\tau+h|\tau}\quad\text{(error on the transformed scale, with the }a_\tau,b_\tau\text{ of origin }\tau),$$

$$c^2_{h,t} = \frac{\sum_{\tau\in W_t} u_\tau^2}{\sum_{\tau\in W_t}\sigma^2_{h,\tau}},\qquad
W_t = \{\text{the last } R \text{ origins } \tau \text{ with } \tau+h \le t\},\qquad
\tilde\sigma_{h,t} = c_{h,t}\,\sigma_{h,t}.$$

Notes on the formula:
- The density stays Gaussian on the transformed scale. Only its width changes.
- Using the ratio rather than replacing $\sigma_{h,t}$ by the root mean squared error keeps the
  window-to-window variation of $\sigma_{h,t}$ (the variance rises when the estimation window becomes
  more volatile).
- The numerator is the mean squared error, not the variance, so a persistent bias also widens the density.

Defaults:
- Window $R = 365$ origins.
- At least $R_{\min} = 180$ resolved errors, otherwise $\tilde\sigma$ is NA.

The first $h + 180$ origins of each horizon (to about mid-2015) therefore drop out of the density
evaluation. They drop out for **all** models, including the benchmark, so DM/GW compare the same days.

**No re-estimation is needed.** Everything uses the existing forecasts and `_sd` files, inside the density
stage.

### Code (`density_eval.R`)
1. Knobs, near the top:
```r
DENSITY_RECAL <- as.integer(Sys.getenv("DENSITY_RECAL", "0"))   # 1 = rescale sigma_h with past h-step errors
RECAL_WINDOW  <- as.integer(Sys.getenv("RECAL_WINDOW", "365"))  # R
RECAL_MIN     <- as.integer(Sys.getenv("RECAL_MIN", "180"))     # R_min
RESULT_DIR    <- file.path(script_dir, Sys.getenv("RESULT_SUBDIR", "results"))   # replaces the fixed "results"
```
2. New function (it can go in `common_utils.R`):
```r
# Rows must be ordered by origin date; for a fixed h the target dates are then increasing too.
recal_sigma <- function(point, sigma, a, b, y, od, td,
                        R = RECAL_WINDOW, Rmin = RECAL_MIN) {
  mu <- vst_z(point, a, b)
  u  <- vst_z(y, a, b) - mu                       # h-step error on the transformed scale
  ok <- is.finite(u) & is.finite(sigma) & sigma > 0
  e2 <- ifelse(ok, u^2, 0); s2 <- ifelse(ok, sigma^2, 0); n1 <- as.numeric(ok)
  ce <- c(0, cumsum(e2)); cs <- c(0, cumsum(s2)); cn <- c(0, cumsum(n1))
  j  <- findInterval(as.numeric(od), as.numeric(td))   # rows 1..j have target <= origin (resolved)
  i0 <- pmax(0L, j - R)                                # window = rows (i0+1):j
  n  <- cn[j + 1] - cn[i0 + 1]
  c2 <- (ce[j + 1] - ce[i0 + 1]) / (cs[j + 1] - cs[i0 + 1])
  out <- sigma * sqrt(c2)
  out[n < Rmin | !is.finite(out)] <- NA_real_
  out
}
```
3. In `run_country()`, right after reading `ar`:
   - order the rows by origin: `o <- order(ar$origin_date); ar <- ar[o, ]` (then match `ars`, `fc`, `fcs`
     to the same keys, as now);
   - `od <- as.Date(ar$origin_date)`.
4. Make `score_block()` take `od, td` and rescale at the top:
```r
score_block <- function(point, sigma, a, b, y, umin, umax, od = NULL, td = NULL) {
  if (DENSITY_RECAL == 1L) sigma <- recal_sigma(point, sigma, a, b, y, od, td)
  ...   # unchanged
}
```
   Pass `od, td` in the two calls: `S_bench <- score_block(..., od = od, td = td)` and inside
   `parallel_lapply(...)`.
5. Add a counter so the log shows the rescaling at work. Print the median $c_{h,t}$ per family and horizon:
   values above 1 mean the original densities were too narrow.
6. Twin outputs: wrap the `twin_rows` / `density_tests_vs_twin_*` / `crps_twin_summary_*` blocks in
   `if (Sys.getenv("REPORT_TWINS", "0") == "1")`. They are not used by the paper any more.

### Unit test (top of `density_eval.R`, skipped when `DENSITY_RECAL = 0`)
```r
local({ old <- VST_MODE; VST_MODE <<- "level"; set.seed(1)
  n <- 3000; od <- as.Date("2015-01-01") + 0:(n-1); td <- od + 7
  y <- rnorm(n, 0, 2); s <- rep(1, n)            # true sd 2, model says 1
  st <- recal_sigma(rep(0, n), s, rep(0, n), rep(1, n), y, od, td, R = 365, Rmin = 180)
  stopifnot(abs(median(st, na.rm = TRUE) - 2) < 0.1)
  VST_MODE <<- old })
```

### Runs
```bash
# Germany (forecasts of estimation_v4, normalised asinh)
cd midas_codes/estimation_v4_crps
COUNTRY_NAME=Germany VST_MODE=normalized DENSITY_RECAL=1 RESULT_SUBDIR=results_recal \
  STAGES=density Rscript run_all.R 2>&1 | tee recal_de.log
# Italy (log run)
cd ../estimation_v4_it_log
COUNTRY_NAME=Italy VST_MODE=log DENSITY_RECAL=1 RESULT_SUBDIR=results_recal \
  STAGES=density Rscript run_all.R 2>&1 | tee recal_it.log
```
- Sensitivity: `RECAL_WINDOW=730` (into `results_recal730`). Report it in the appendix if the conclusions
  change.

### Checks
- [ ] Unit test passes; the median $c_{h,t}$ is above 1 and grows with $h$, more for models with gas or
      monthly indicators.
- [ ] 90% coverage close to 0.90 and PIT sd close to 0.289 at every horizon.
- [ ] CRPS and log score, before and after, **on the same days** (drop the burn-in from the old scores too).
      The log score should improve clearly at h = 21–28; the CRPS changes less.
- [ ] Rankings of the models by CRPS: do they still match the RMSE rankings?

---

## 2. Descriptive statistics on `data_v4`
`descriptive_stats.py` (in `r_pipeline/`) still reads `r_pipeline/data/`, the **old** May 2026 files with
interpolated weekends.

- [ ] **Data path:** `D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "midas_codes", "data_v4")`.
- [ ] **Monthly columns:** drop `istat_bci_it` (keep `iesi_it`). Keep `Electricity_Gas_Steam` (= IPI energy).
- [ ] **Sub-periods aligned with the results:**
      `per = [("2005","2020"), ("2021","2022"), ("2023","2025"), ("2015","2025")]`
      (Table `tab:desc-sub`; the new values are already in `corrections_v3.tex`, block 4a. Check them
      against the script output.)
- [ ] **Indicator–price correlation by sub-period** (Table `tab:desc-corr`, block 4b). Reuse this code:
```python
PER = [("2005-01-01","2020-12-31"), ("2015-01-01","2020-12-31"), ("2021-01-01","2022-12-31"), ("2023-01-01","2025-12-31")]
pm = d.elec_price.resample("ME").mean()
m["price"] = pm.reindex(m.index)
for v in IPI: m[v] = np.log(m[v]).diff()            # IPI as monthly log growth, surveys in levels
corr = {v: [m.loc[a:b, [v, "price"]].dropna().corr().iloc[0, 1] for a, b in PER] for v in inds}
```

- [ ] **Moments** (`tab:desc-daily`, `tab:desc-monthly`): rerun `summ()` on the new data. The rows are price,
      asinh(price), normalised asinh, **log(price) for Italy (new)**, gas and Brent. Update the "14 negative
      days" mentions: Germany has 38.
- [ ] **Lead–lag table** (`tab:desc-leadlag`): rerun on `data_v4`. Check whether the conclusion still
      holds (surveys lead, IPI coincident).
- [ ] **Weekday and month profile** (old TODO, now possible with real weekends): mean price by weekday and
      by month relative to the overall mean, for 2012–2025 in Germany and 2005–2025 in Italy.

---

## 3. Germany: full plain-asinh run
The 2013–2014 validation already favours the normalised asinh (RMSE 12–15% lower). A full run makes the
appendix comparison symmetric with Italy.

```bash
cd midas_codes
mkdir -p estimation_v4_de_asinh && cd estimation_v4_de_asinh
for f in ../*.R; do ln -sf "$f" .; done; ln -sfn ../data_v4 data
COUNTRY_NAME=Germany VST_MODE=asinh SPEC_REQUIRE_LF=1 caffeinate -i Rscript run_all.R 2>&1 | tee run_de_asinh.log
cd ..
COMPARE_CC=de VST_BASE=normalized \
VST_RUNS="normalized=estimation_v4,asinh=estimation_v4_de_asinh" \
Rscript compare_vst_runs.R 2>&1 | tee compare_vst_de.log
```
- Optional: a level run (`VST_MODE=level`, folder `estimation_v4_de_level`) for a three-way comparison as
  for Italy. `log` is impossible (negative prices).
- [ ] Check the log: Almon diagnostics (no bound hits), CRPS unit test passed, no NaN forecasts at
      negative-price days.
- [ ] `tables_v3/scripts/tab_vst.py`: add a Germany panel reading `results_compare/compare_vst_de.csv`
      (normalised vs asinh [vs level]), and remove the TODO in `appendix_v3.tex` (A.6).
- [ ] Run it with `SPEC_REQUIRE_LF=1` (item 0), i.e. only the specifications with monthly indicators.

---

## 4. Other code fixes

**a. Combinations over the 252 macro specifications, in R.** The paper's combinations are currently
recomputed in Python (`tables_v3/scripts/comb252.py`) as averages of the 252 specifications with monthly
indicators. `combination_subperiod_gr.R` still averages all 256.
- Needed only for runs made **without** `SPEC_REQUIRE_LF=1` (the existing v4 runs). In
  `combination_subperiod_gr.R` (`comb_eq`, `comb_trim`, around line 385), restrict the forecast matrix
  `Fm` (and `Fbc`) to `spec_id`s whose `lf_vars` in the family dictionary are non-empty:
```r
dict <- read.csv(file.path(FORECAST_DIR, sprintf("%s_spec_dictionary_%s.csv", fam, cc)), stringsAsFactors = FALSE)
macro_ids <- dict$spec_id[!is.na(dict$lf_vars) & dict$lf_vars != ""]
Fm <- Fm[, colnames(Fm) %in% macro_ids, drop = FALSE]
```
- Rerun `STAGES=combination,tests,cssed` in `estimation_v4` and `estimation_v4_it_log`. Check that
  `comb_r_midas_l7_*` matches `comb252.csv` (differences < 1e-8).

**b. SPA and MCS tests** (`tables_v3/scripts/spa.py`, `mcs.py`). They now run on the 766 models (benchmark,
3 AR, 3 × 252, 6 combinations).
- Moving-block bootstrap, B = 1000, block length max(h, 16).
- Run them per country **in the foreground** (`python3 spa.py de 1 7 14 21 28`); background jobs were killed
  on the Cowork VM.
- Optional: a port to R (the `MCS` package for the model confidence set; SPA is about 30 lines with the
  functions already in the scripts).

**c. `compare_vst_runs.R`.** The rbind bug in the tests output is already fixed. Rerun it for Italy to get
`compare_vst_tests_it.csv` (tests between transformations), and for Germany after item 3.

**d. Mean forecasts after item 1 (optional).** The `__bc` forecasts use the model σ² in
$\exp(\sigma^2_h/2)$. With the rescaled $\tilde\sigma$ the mean correction changes. Recompute
`Table tab:app-bc` with $\tilde\sigma$ if you keep that robustness check.

**e. Housekeeping in the project folder.** Delete when convenient:
- `results_v3_preview.pdf`, `_preview_results_v3.tex` (superseded by `chapter2_v3_preview.pdf`);
- `img/_to_delete/`;
- `references_results_v3.bib` (superseded by `references_v3.bib`);
- `tables_v3/tab_app_dailyonly.tex` (no longer used).

---

## 5. Regenerate the paper tables (after any change above)
From `Forecasting electricity prices/tables_v3/scripts/`. The paths are in `common.py`; `sel.json`,
`comb252.csv`, `spa_*.csv` and `mcs_*.csv` are written to `~/an/`, so change that path if you run elsewhere.
```bash
python3 comb252.py        # combinations over the 252 macro specs (must run first)
python3 tab_main.py       # tab_point_de/it + sel.json
python3 tab_win.py        # tab_winners
python3 tab_marg.py       # tab_marginal (pairs where both specs contain monthly indicators)
python3 tab_crps.py       # tab_crps  (results_recal/ after item 1)
python3 fig_gr.py; python3 gr_t.py          # Giacomini-Rossi figures
python3 spa.py de 1 7 14 21 28; python3 spa.py it 1 7 14 21 28
python3 mcs.py de; python3 mcs.py it        # then: python3 mcs_summary.py
python3 tab_spa.py
python3 tab_app.py; python3 tab_vst.py      # appendix tables
```
Then recompile `chapter2_v3_preview.pdf` (`_preview_v3_full.tex`) and re-check every number quoted in
`results_v3.tex`, `appendix_v3.tex` and `corrections_v3.tex`.
