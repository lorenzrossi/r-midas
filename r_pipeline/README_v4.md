# v4: data fixes, normalised asinh, Giacomini–White tests, CSSED plots, subperiod tests

Interim plan agreed on 22 Oct 2026: run the R-MIDAS analysis on the full 2005–2025
sample now, with German 2005–2011 weekends handled provisionally (see 1.1). Each step
says whether it is **DONE** (scripts written and run) or **TO DO** (specification and
code to add to the R pipeline).

Folders used below:
- `DATA` = `Forecasting electricity prices/data` (project repo; raw downloads in `data/raw_prices/`)
- `MC`   = `r-midas/r_pipeline/midas_codes` (this repository)

---

## 0. Which data to use (read this first)

**Use only the v4 files.** The R pipeline must read exactly these four files:

| File | Content |
|---|---|
| `midas_codes/data_v4/prices_de_daily.csv` | daily: `date, elec_price, gas_price, brent_price, ciss, country, elec_imputed, elec_source` |
| `midas_codes/data_v4/prices_it_daily.csv` | same columns (Italy: `elec_source` = `real_all_days` everywhere) |
| `midas_codes/data_v4/dataset_de_m.csv` | monthly: prices, IPI (CG, Mfg, Energy), `pmi_mfg_de, pmi_serv_de, ifo_bci_de` |
| `midas_codes/data_v4/dataset_it_m.csv` | monthly: prices, IPI, `pmi_mfg_it, pmi_serv_it, istat_bci_it, iesi_it` |

They are copies of `DATA/data_v4/`. Rebuild and copy them **only** with the two scripts:
```bash
cd "$HOME/Documents/Claude/Projects/Forecasting electricity prices/data"
python3 build_daily_prices.py && python3 build_dataset_v4.py
cp data_v4/*.csv ../r-midas/r_pipeline/midas_codes/data_v4/
```

**How the pipeline picks them up.** The R scripts read `file.path(script_dir, "data")`.
Run v4 from `midas_codes/estimation_v4/`, where `data` is a symlink to `../data_v4`:
```bash
cd "<repo>/r-midas/r_pipeline/midas_codes"
mkdir -p estimation_v4 && cd estimation_v4
for f in ../*.R; do ln -sf "$f" .; done          # scripts = symlinks to the shared code
ln -sfn ../data_v4 data                          # data = the v4 files
ls -l data                                       # must show: data -> ../data_v4
```
Never run v4 from `midas_codes/` itself or from `estimation_v3/`: there `data` is the old
folder.

**Check before every run (R):** it stops if the wrong files are loaded.
```r
d <- read.csv("data/prices_de_daily.csv"); m <- read.csv("data/dataset_it_m.csv")
stopifnot("elec_source" %in% names(d),                     # only v4 files have this column
          "iesi_it" %in% names(m))                         # only v4 Italy monthly has IESI
# real weekends: from 2015 (2012 with DE.csv) Saturday must almost never equal Friday
d$date <- as.Date(d$date); sat <- which(as.POSIXlt(d$date)$wday == 6 & d$date >= as.Date("2015-01-05"))
stopifnot(mean(d$elec_price[sat] == d$elec_price[sat - 1]) < 0.01)
table(d$elec_source)    # DE now: real_all_days 4014; with DE.csv: 5106 (2012-01-09 onwards)
```
The quickest check is the first line: if `elec_source` is missing, you are on the old data.

**Do NOT use** (all built with interpolated or synthetic weekends, or superseded):

| Location | Files | Why |
|---|---|---|
| `DATA/` (top level) | `prices_de_daily.csv`, `prices_it_daily.csv`, `dataset_de_m.csv`, `dataset_it_m.csv` | old `preprocessing.ipynb` output: weekends and monthly gaps linearly interpolated (uses future values) |
| `DATA/` | `de0525*.csv`, `it0525*.csv`, `Germany_DE*.csv`, `Italy_IT*.csv`, `*Baseload*.csv`, `*_Monthly_Baseload_IPI*.csv` | earlier versions: weekends interpolated or repeated from Friday, Datastream errors |
| `DATA/` | `preprocessing.ipynb` | superseded by `build_daily_prices.py` and `build_dataset_v4.py` |
| `DATA/raw_prices/` | `*_MGP_PrezziZonali_NAT.xlsx` | "NAT" is not the PUN (spurious 2015 spikes); the scripts ignore it |
| `r_pipeline/data/`, `midas_codes/data/` | all four CSVs | the old (v1–v3) inputs, kept only to reproduce `estimation_v3` |

`DATA/elec_daily_de.csv`, `elec_daily_it.csv` and `iesi_monthly.csv` are **intermediate**
files read by `build_dataset_v4.py`. Don't feed them to R directly.
`DATA/raw_prices/` holds the raw downloads (ENTSO-E, GME PUN, and `DE.csv` from
epftoolbox once downloaded). Add any new download there and rerun the two scripts.

**Optional housekeeping, recommended:** move the obsolete files out of the way so they
can't be picked up by mistake:
```bash
cd "$HOME/Documents/Claude/Projects/Forecasting electricity prices/data"
mkdir -p old && mv prices_de_daily.csv prices_it_daily.csv dataset_de_m.csv dataset_it_m.csv \
  de0525*.csv it0525*.csv old/
```

---

## 1. Data

### 1.1 Germany: real prices from 2012, interim rule for 2005–2011 (DONE, except one download)

| Period | Source | Weekends |
|---|---|---|
| 2005-01-01 → 2012-01-08 | Datastream `EEX - Phelix Base` (weekdays only) | **carried forward from Friday** (interim, to be acknowledged in the paper) |
| 2012-01-09 → 2015-01-04 | epftoolbox EPEX-DE hourly (Lago et al., 2021) | real |
| 2015-01-05 → 2025-12-31 | ENTSO-E day-ahead, Sequence 1 (BZN DE-AT-LU to 2018-09-30, DE-LU after) | real |

Scripts in `DATA`:
- `download_epftoolbox_de.py`: downloads `https://zenodo.org/records/4624805/files/DE.csv`
  to `DATA/raw_prices/DE.csv`.
- `build_daily_prices.py`: daily baseload (mean of all hours or quarter-hours of the
  delivery day) from `raw_prices/`, with ENTSO-E taking priority over epftoolbox where
  they overlap (2015–2017). It prints a quality report, including the
  epftoolbox–ENTSO-E gap on the overlap, which should be about 0.

**TO DO (you, once):** the Cowork sandbox blocks zenodo.org, so `DE.csv` is not
downloaded yet. On your Mac, in Terminal:
```bash
cd "$HOME/Documents/Claude/Projects/Forecasting electricity prices/data"
python3 download_epftoolbox_de.py        # or download the URL above by hand into raw_prices/DE.csv
python3 build_daily_prices.py            # check: "epftoolbox vs ENTSO-E ... MAD = ~0"
python3 build_dataset_v4.py
cp data_v4/*.csv ../r-midas/r_pipeline/midas_codes/data_v4/
```
Until then, the files in `data_v4/` use the Datastream + Friday rule up to **2015-01-04**.
After the rerun it applies only up to **2012-01-08**. In both cases every *evaluated*
target (from 2015-01-05) is a real price; only estimation windows contain carried-forward
weekends.

### 1.2 Italy: PUN and IESI (DONE)
- Electricity: GME **PUN**, every calendar day 2005–2025 (`raw_prices/*_PUN*.xlsx`); from
  2025 this is the PUN Index GME. It matches Datastream to about €0.003/MWh on weekdays in
  2005–2016. The `*_NAT` files are **not** the PUN and are not used.
- **IESI** (`iesi_monthly.csv`, from ISTAT `DCSC_IESI_21`, base 2021 = 100) enters the
  monthly file as the new column `iesi_it`:
  - Jan and Feb 2005 = Mar 2005 value (the series starts in 2005-03; acknowledge in the paper);
  - Apr 2020 = Mar 2020 value (survey not run).
- `istat_bci_it` is kept: it is the ISTAT **manufacturing** confidence climate, with
  Apr 2020 = Mar 2020.

### 1.3 Conventions now applied by `build_dataset_v4.py` (DONE)
- **No interpolation anywhere.** Electricity (Germany before 2012), gas, Brent and CISS
  take the last available value over weekends, holidays and gaps, so nothing comes from
  the future. The old notebook interpolated straight lines that used Monday's value.
- **Monthly gaps are carried forward only.** The Datastream PMIs have gaps that used to be
  interpolated; the script lists them each time it runs:
  - DE: 2005-02, 2007-06/11/12, 2009-04/06, 2019-07;
  - IT: 2005-02, 2007-04/12, 2009-12, 2010-05, 2012-05, 2015-03, 2016-12.

  Consider re-downloading the PMIs.
- **New daily columns:** `elec_imputed` (1 = carried forward) and `elec_source`
  (`real_all_days` / `datastream_weekday` / `carried_forward`).
- **Naming for the paper:** `ifo_bci_de` is the **ifo Business Expectations** (next 6
  months), not the composite climate. The column name is kept for code compatibility.

Outputs: `DATA/data_v4/` and `MC/data_v4/` (`prices_{de,it}_daily.csv`,
`dataset_{de,it}_m.csv`).

---

## 2. R pipeline (TO DO)

### 2.1 New run folder
Create `MC/estimation_v4/` like `estimation_v3/`: symlink every `.R` script to `../`,
point `data -> ../data_v4`, and keep outputs local. Do not overwrite `data/` or `estimation_v3/`.

### 2.2 IESI as the Italian survey indicator (`common_utils.R`)
1. `LF_CANDIDATES_IT`: replace `"istat_bci_it"` with `"iesi_it"`. The manufacturing climate
   becomes a robustness run; add an env switch `LF_IT_SURVEY = iesi_it | istat_bci_it`.
2. In `load_country_data()`, update `survey_cols` for Italy the same way.
3. In `usual_lf_release_date()`, add a rule for `iesi_it`: the 27th of the reference month,
   rolled forward to a weekday (same as the current `istat_bci_it` rule).
4. The number of specifications stays at 256 (still 6 monthly candidates).

### 2.3 Normalised asinh inside each rolling window (Uniejewski, Weron & Ziel, 2018)
For each origin t0 with window W = rows `(t0-WINDOW_DAYS+1):t0`:
```
a = median(P_W)        b = mad(P_W)          # R's mad() already scales by 1.4826
y = asinh((P - a) / b)                       # applied to the window AND to every lagged price it uses
P_hat    = a + b * sinh(y_hat)               # median forecast, EUR/MWh
P_hat_bc = a + b * sinh(y_hat) * exp(s2_h/2) # mean forecast (bias-corrected)
```
`a` and `b` use data up to t0 only, so there is no look-ahead. The recursion of the
iterated forecasts stays on the transformed scale with the origin's (a, b).

Helpers for `common_utils.R`:
```r
VST_MODE <- tolower(Sys.getenv("VST_MODE", "normalized"))      # "raw" reproduces v3
vst_params <- function(P) if (VST_MODE == "raw") list(a = 0, b = 1) else
  list(a = median(P, na.rm = TRUE), b = mad(P, na.rm = TRUE))
vst     <- function(P, p) asinh((P - p$a) / p$b)
vst_inv <- function(z, p, s2 = 0) p$a + p$b * sinh(z) * exp(s2 / 2)
```

Where to change the code (line numbers as of 20 Oct 2026):

| File | Now | Change |
|---|---|---|
| `common_utils.R` ~118 | `daily$y_asinh <- asinh(elec_price)` | keep the raw `elec_price`; transform per window |
| `ar_benchmark.R` ~46, ~97, ~152–154 | `y`, `ytr <- y[win_lo:t0]`, `sinh(mu_h)` | build the lag matrix on **raw prices once**, then per origin `vst()` the window rows; back-transform with `vst_inv()` |
| `arx_free_lags.R` ~66, ~205–211 | same pattern | same |
| `r_midas.R` ~138, ~255, ~469–477 | same, plus `history_y_aug` | same; `history_y_aug` must be the vst of past prices with the origin's (a, b) |
| `density_eval.R` ~96–99, 179, 196 | `mu = asinh(point)`, `U = sinh(...)`, `z = asinh(y)` | `mu = asinh((point - a)/b)`; `U = a + b*sinh(mu + sigma*z)`; `z_y = (asinh((y - a)/b) - mu)/sigma`; density Jacobian `1/sqrt((u-a)^2 + b^2)` |
| `summary_table_ar_r_midas.R` ~86, `combination_subperiod_gr.R` ~206 | errors on `asinh(y)` | errors on `asinh((y - a)/b)` with the origin's (a, b), the same for every model at that origin |

- Save (a, b) per origin in `forecasts/vst_params_<cc>.csv` (`origin_date, a, b`) and
  read it in the density and summary stages.
- **Speed tip:** compute `build_ar_lags(P, AR_LAGS)` on raw prices once. Per origin,
  `asinh((Plag[win_rows, ] - a)/b)` is a cheap vectorised operation.
- **Check after the run:** the German persistence should rise. The median δ was 0.58 with
  the raw asinh; compare with the free-lag sum of AR coefficients.

### 2.4 Giacomini–White (2006) conditional test (supervisors' request)
Valid with a fixed-length rolling window and with nested models. Instruments: the
constant and d_{t−h}. Use d_{t−h}, not d_{t−1}: only the loss differential of the
forecast made h days earlier has been realised at origin t.

Add to `common_utils.R`:
```r
# e1 = errors of model 1 (test), e2 = errors of model 2 (reference), ordered by consecutive daily origins.
# d = L(e1) - L(e2): negative mean -> model 1 more accurate. Returns chi^2 stat with q = 2 df.
gw_test <- function(e1, e2, h, loss = c("sqerr", "abserr"), keep = NULL) {
  loss <- match.arg(loss)
  d  <- if (loss == "sqerr") e1^2 - e2^2 else abs(e1) - abs(e2)
  n  <- length(d)
  dl <- c(rep(NA_real_, h), d[seq_len(max(n - h, 0L))])     # instrument d_{t-h}
  Z  <- cbind(1, dl)
  ok <- is.finite(d) & is.finite(dl)
  if (!is.null(keep)) ok <- ok & keep                       # subperiod rows (lags taken from full history)
  d <- d[ok]; Z <- Z[ok, , drop = FALSE]
  T_ <- length(d); q <- ncol(Z)
  if (T_ < 50L) return(list(stat = NA_real_, pvalue = NA_real_, T = T_, dbar = NA_real_))
  Zd <- Z * d
  m  <- colMeans(Zd)
  Om <- crossprod(Zd) / T_
  if (h > 1L) for (l in 1:(h - 1L)) {                       # Newey-West (Bartlett), truncation h-1: PSD
    G  <- crossprod(Zd[(l + 1L):T_, , drop = FALSE], Zd[1:(T_ - l), , drop = FALSE]) / T_
    Om <- Om + (1 - l / h) * (G + t(G))
  }
  stat <- T_ * drop(crossprod(m, solve(Om, m)))
  list(stat = stat, pvalue = pchisq(stat, q, lower.tail = FALSE), T = T_, dbar = mean(d))
}
```
- The statistic has no sign. Report the direction with `dbar` (negative means model 1 is better).
- Compute for squared and absolute errors (EUR/MWh) and for the transformed scale.
- Every comparison reports **both** DM-HLN and GW, for squared and absolute errors, in
  EUR/MWh and on the transformed scale, for the full sample and each subperiod (2.6).

#### Table A: every model against the AR benchmark
- **Reference:** `ar_dum`, the AR(1,2,3,7) with calendar dummies.
- **Rows:** every model, with and without macro information:
  - AR family: `ar_dum_gas`, `ar_dum_brent`, `ar_dum_gas_brent`;
  - all 256 ARX specifications (free lags, including the daily-only ones);
  - all 256 R-MIDAS specifications;
  - the real-time combinations (`comb_*`, both families).
- **Output:** `results/tests_vs_benchmark_<cc>.csv`, with columns `country, h, period, family,
  model, hf_vars, lf_vars, n, rmse, rmse_ratio, mae_ratio, dm_stat, dm_pval, gw_stat,
  gw_pval, gw_dbar` (+ `_mae` and `_asinh` versions).
- **For the paper:** a main-text table per country × horizon with the daily-only ARX, a few
  pre-specified macro models (all surveys / all IPI / all monthly) and the combinations.
  The full file goes in the online appendix.

#### Table B: each macro specification against its daily-only twin
- **Definition:** the twin of a specification with daily regressors H and monthly set
  L ≠ ∅ is the specification **of the same family** with the same H and **no monthly
  regressor**.
  - ARX family: same free lags.
  - R-MIDAS family: same Almon lag structure.
  - In both cases the twin is nested, and the only difference is the monthly information,
    which is what this table isolates.
  - Each family has 4 daily-only twins (H = none, gas, brent, gas+brent) and 252 macro
    specifications. The ARX twin with H = none is the AR benchmark itself.
- **Combinations:** a family's combination against that family's daily-only gas+brent specification.
- **Twin mapping** from the spec dictionaries (`forecasts/<family>_spec_dictionary_<cc>.csv`,
  where an empty `lf_vars` means no monthly regressor):
```r
twin_map <- function(dict) {
  lf <- ifelse(is.na(dict$lf_vars), "", dict$lf_vars)
  hf <- ifelse(is.na(dict$hf_vars), "", dict$hf_vars)
  twins <- data.frame(hf_vars = hf[lf == ""], twin_id = dict$spec_id[lf == ""])
  macro <- data.frame(spec_id = dict$spec_id[lf != ""], hf_vars = hf[lf != ""], lf_vars = lf[lf != ""])
  merge(macro, twins, by = "hf_vars", all.x = TRUE)   # one row per macro spec: spec_id -> twin_id
}
# for each family, country, h: e_m <- y - F[[spec_id]]; e_t <- y - F[[twin_id]] (same forecast file,
# same origins);  dm_hln(e_m, e_t, h);  gw_test(e_m, e_t, h, keep = in_period)
```
- **Output:** `results/tests_vs_twin_<cc>.csv`, with columns `country, h, period, family,
  model, twin, hf_vars, lf_vars, rmse_ratio_vs_twin, dm_stat, dm_pval, gw_stat, gw_pval,
  gw_dbar` (+ `_mae`, `_asinh`).
- **For the paper:** a summary table per country × horizon × period with:
  - the median RMSE ratio against the twin;
  - the share of the 252 macro specifications that beat their twin at 5% (GW and DM);
  - the same for the specifications containing each indicator, which ties into the horizon
    story (surveys against IPI).

  The full file goes in the appendix.
- **Multiple testing:** with 252 comparisons per cell, report the shares, not the single
  best specification's p-value. A Model Confidence Set can be added later.

### 2.5 "CUSUM" plots: cumulative difference in squared forecast errors (Welch & Goyal, 2008)
```r
cssed <- function(e_bench, e_model) cumsum(ifelse(is.finite(e_bench) & is.finite(e_model),
                                                  e_bench^2 - e_model^2, 0))
```
New script `cssed_plots.R`:
- One panel per country × horizon, x-axis = target date. A rising line means the model is
  beating the reference.
- Curves:
  1. ARX daily-only (`ar_dum_gas_brent`) against AR;
  2. the real-time macro combination (`comb_invmse`, R-MIDAS, `R_MIDAS_L7` and ARX
     families) against AR;
  3. each macro combination against its own family's daily-only gas+brent twin
     (ARX against ARX, R-MIDAS against R-MIDAS, L7 against L7);
  4. `R_MIDAS_L7` without macro (`spec_000`) against AR (the lag-structure effect alone).
- Shade the subperiods: pre ≤ 2020-12-31; crisis 2021–2022; post ≥ 2023.
- Outputs: `results/cssed/cssed_<cc>_h<hh>.png` and `results/cssed/cssed_<cc>.csv`.
  Use base R graphics, so there is no new dependency.
- Keep the Giacomini–Rossi (2010) fluctuation test next to the plots; it gives the formal
  critical values.

### 2.6 All metrics by subperiod, with DM and GW inside each subperiod
In `combination_subperiod_gr.R` (`PERIODS`: full / pre / crisis / post), extend
`eval_model()` so that each period reports the following, for **every family**: AR, ARX,
R_MIDAS, **R_MIDAS_L7** and the combinations:
- **Point metrics:** RMSE, MAE, bias-corrected RMSE and transformed-scale RMSE/MAE, each as
  a ratio to the benchmark and to the daily-only twin.
- **Density metrics:** mean CRPS, log score and the 90% coverage, from `density_eval.R`
  output restricted to the same target dates.
- **Tests:** DM-HLN (as now) and `gw_test(..., keep = in_period)`, for **both** Table A
  (against the AR benchmark) and Table B (against the daily-only twin). Build d and d_{t−h} on
  the **full** origin sequence first, then subset with `keep`, so the first h days of a
  period use the instrument from the end of the previous period (it is known at the origin).
- Report the number of observations per period (crisis ≈ 730).

### 2.7 Optional robustness: no carried-forward targets in estimation
Env switch `DROP_IMPUTED = 1`: give rows with `elec_imputed == 1` zero weight in the
estimation windows (the weighted-LS machinery already exists for `RMIDAS_DISCOUNT`).
Evaluation is unaffected (all evaluated targets are real).

### 2.8 R-MIDAS with consecutive lags 1–7 (new family `R_MIDAS_L7`)
**Why.** The exponential Almon curve is designed for a run of consecutive lags. On the gappy
set {1,2,3,7} the smoothness restriction has little meaning. `R_MIDAS_L7` puts the same Almon
polynomial on lags 1,2,…,7, so it uses the same information as the benchmark (prices up to
one week back). Everything else is unchanged:
- the 256 specifications (4 daily sets × 64 monthly sets);
- the θ box [−8, 8]², with 7 starts;
- iterated forecasts, the normalised asinh (2.3) and the calendar dummies.

**Script.** Parametrise `r_midas.R` instead of copying it, so the two families can't drift apart:
1. At the top of `r_midas.R`, add:
```r
.parse_lags <- function(x) if (grepl(":", x)) eval(parse(text = x)) else as.integer(strsplit(x, ",")[[1]])
RM_LAGS   <- .parse_lags(Sys.getenv("RMIDAS_LAGS", "1,2,3,7"))   # "1:7" for the new family
RM_PREFIX <- Sys.getenv("RMIDAS_PREFIX", "r_midas")              # file prefix
RM_FAMILY <- Sys.getenv("RMIDAS_FAMILY", "R_MIDAS")              # label in the summaries
```
2. In `r_midas.R`, replace `AR_LAGS` with `RM_LAGS` everywhere it refers to the Almon block.
   Lines as of 20 Oct 2026:
   - `build_ar_lags` (~141), `AR_LAG_SETS` (~143), `pre_lag_max` (~215), `phi_cols` (~225);
   - `exp_almon_weights_j` (~345, ~370), `iterated_error_variance` (~378), `ar_min_root_modulus` (~391).
   
   The implied coefficients are saved as `phi_1 … phi_7`.
3. Replace the hard-coded `"r_midas_..."` file names with `RM_PREFIX`:
   - spec dictionary (~206), forecasts (~491), sds (~508), coefficients (~526);
   - e.g. `sprintf("%s_%s_h%02d.csv", RM_PREFIX, cc, h)`.
4. New wrapper `r_midas_l7.R`:
```r
# R-MIDAS with exponential Almon on the consecutive daily lags 1..7 (family R_MIDAS_L7)
Sys.setenv(RMIDAS_LAGS = "1:7", RMIDAS_PREFIX = "r_midas_l7", RMIDAS_FAMILY = "R_MIDAS_L7")
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)][1])))
source(file.path(here, "r_midas.R"), chdir = TRUE)
```
5. Starting values: the U-shaped starts already in `fit_r_midas_fast.R`, (−3.5, 0.40) and
   (−2.8, 0.32), also work on 1..7. For example, (−3.5, 0.4) gives weights of about 0.75 on
   lag 1, 0.12 on lag 7 and 0.08 on lag 2. No change is needed.
6. Wire the family into the other stages:
   - `run_all.R`: stage `r_midas_l7 = "r_midas_l7.R"`, placed after `r_midas` in `DEFAULT_ORDER`;
   - `combination_subperiod_gr.R`: add `r_midas_l7` to `FAMILIES`, which gives `comb_r_midas_l7_*`;
   - `summary_table_ar_r_midas.R`: add file and dictionary entries `R_MIDAS_L7`
     (`r_midas_l7_%s_h%02d.csv`, `r_midas_l7_spec_dictionary_%s.csv`) and include it in the
     family list (~258);
   - `density_eval.R`: `FAM_LABEL["r_midas_l7"] <- "R_MIDAS_L7"` and add it to `DENSITY_FAMILIES`.
7. Check after the run (in the log, as for `R_MIDAS`): the share of windows at the box bounds,
   the share with U-shaped weights, the median implied φ₁…φ₇ and δ. Compare φ with the
   free-lag AR(7) estimates from `ols_coef_check.R`, after adding lags 4–6 there.

**`R_MIDAS_L7` goes through every stage, exactly like `R_MIDAS`.** None of these is optional:

| Stage | What it must produce for `R_MIDAS_L7` |
|---|---|
| Estimation (`r_midas_l7.R`) | point forecasts `r_midas_l7_<cc>_h<hh>.csv`, bias-corrected forecasts (`__bc`), predictive sds `r_midas_l7_<cc>_h<hh>_sd.csv`, spec dictionary, coefficients and θ diagnostics (`results/coefficients/r_midas_l7_<cc>.csv`) |
| Normalised asinh (2.3) | automatic, because the same code is shared; (a, b) per origin come from the same `vst_params_<cc>.csv` |
| Forecast combinations | `comb_r_midas_l7_<cc>_h<hh>.csv` (equal, trimmed, inverse-MSE, two-stage inverse-MSE), real-time weights |
| Point summary (`summary_table_ar_r_midas.R`) | RMSE / MAE / transformed-scale / bias-corrected metrics and top-5 tables, with family `R_MIDAS_L7` |
| Tests (2.4, 2.8) | DM-HLN and GW against the AR benchmark (Table L7-A) and against the twin (Table L7-B) |
| Subperiods (2.6) | all metrics and DM/GW tests for full / pre / crisis / post, in both L7 tables |
| Giacomini–Rossi fluctuation test | computed for the L7 specifications and combinations, as for R-MIDAS (`gr_fluctuation_<cc>.csv`) |
| CSSED plots (2.5) | L7 curves in every panel (curves 2–4 above) |
| Density forecasts (`density_eval.R`) | log score, CRPS, ACPS (c = 0.05 … 0.95), PIT and 90% coverage, and DM tests on the scores against the AR benchmark **and** against the L7 twin, full sample and subperiods; `density_scores_<cc>.csv` and `density_top5_per_country_h.csv` include `R_MIDAS_L7` |
| Marginal contribution (`marginal_contribution.py`) | run on the L7 rows too, so the horizon story (surveys against IPI) can be checked for both lag structures |
| Optional robustness (2.7, `DROP_IMPUTED=1`) | rerun the L7 family as well |

For density forecasts: extend the twin logic to scores. For each macro specification, compare
its log score / CRPS with those of its twin (`dm_hln_loss`, and `gw_test_loss(l1, l2, h, keep)`, identical to
`gw_test` but taking the two score series directly, with d = l1 − l2). Write the results to `results/density_tests_vs_twin_<cc>.csv` for
**all three families** (ARX, R_MIDAS, R_MIDAS_L7).

**Table L7-A: every `R_MIDAS_L7` model against the AR benchmark with lags {1,2,3,7}**
- **Reference:** `ar_dum` (AR(1,2,3,7) + calendar), the same as in Table A.
- **Rows:**
  - all 256 `R_MIDAS_L7` specifications, with and without monthly indicators;
  - `comb_r_midas_l7_*`.
- **Caveat:** this table mixes the lag-structure effect and the macro effect, as agreed. The
  four rows with no monthly regressor isolate the lag-structure effect: `spec_000` against
  `ar_dum` is Almon 1–7 against free {1,2,3,7}, with nothing else different.
- **Tests:** DM-HLN and GW (2.4), squared and absolute errors, EUR and transformed scale,
  full sample and subperiods (2.6).
- **Output:** `results/tests_l7_vs_benchmark_<cc>.csv`, with the same columns as Table A.

**Table L7-B: each `R_MIDAS_L7` macro specification against its twin**
- **Twin:** the `R_MIDAS_L7` specification with the same daily regressors and **no monthly
  indicator**, so the same Almon 1–7 lag structure. It is nested; the only difference is the
  monthly information.
  - Use `twin_map()` from 2.4 on `r_midas_l7_spec_dictionary_<cc>.csv`.
  - 252 macro specifications against 4 twins.
  - `comb_r_midas_l7_*` is compared with the L7 daily-only gas+brent specification.
- **Tests and periods:** as in L7-A.
- **Output:** `results/tests_l7_vs_twin_<cc>.csv`, with the same columns as Table B.
- **For the paper:** a summary table per country × horizon × period with:
  - the median RMSE ratio against the twin;
  - the share of macro specifications beating their twin at 5% (GW and DM);
  - the same share by indicator.

  Put it next to Table B, to show whether the macro effect depends on the lag structure.

---

## 3. Run order and checks
1. Section 1 TO DO (epftoolbox download, rebuild `data_v4`).
2. Code changes 2.2–2.6 (and 2.7 if wanted), then in `estimation_v4/`:
   `Rscript run_all.R` (stages: ols_check, ar, arx, r_midas, r_midas_l7, combination, summary, density),
   then `Rscript cssed_plots.R`.
3. Checks:
   - **Benchmark error changes against v3.** Expected, because targets and weekends changed:
     in 2025, interpolated German weekends were €21.7/MWh too high on average.
   - **Weekday profile of the fitted calendar dummies.** Sunday should now be clearly
     negative for Germany.
   - **German δ.** It should be higher than v3's median of 0.58.
   - **GW against DM.** The unconditional GW (constant only) should be close to DM; the
     conditional test adds the d_{t−h} dimension.
   - **CSSED against the subperiod tables.** The jumps should line up with the crisis
     rows of the subperiod table.

## 4. What the paper must state (data section)
- German day-ahead prices:
  - real every day from 2012-01-09 (epftoolbox until 2015-01-04, ENTSO-E after);
  - 2005-01-01 → 2012-01-08 are Datastream weekday prices with weekend values set to the
    preceding Friday;
  - this affects estimation windows only, never evaluated targets.
- Italian prices: GME PUN (PUN Index GME from 2025).
- Fuel prices and CISS: last business-day value carried over weekends and holidays.
- IESI: Jan–Feb 2005 set to Mar 2005 (series start); Apr 2020 = Mar 2020 (no survey).
- ifo variable: Business Expectations (6 months).
