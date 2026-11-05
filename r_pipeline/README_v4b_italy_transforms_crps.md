# v4b: Italy with prices in levels and in logs, and an exact CRPS

Two additions to the v4 pipeline (see `README_v4.md`):

1. **Temporary Italian runs** with the electricity price modelled (a) in **levels** and (b) in
   **logs**, next to the v4 run (normalised asinh).
2. **CRPS for the density forecasts** (supervisors' request): the exact, untruncated CRPS for
   every model, with tables and DM/GW tests against the benchmark and against the twin.

**Scope.** Everything here is about **Italy**. Germany is not re-estimated and keeps the v4
setup (normalised asinh). The only German step is re-scoring the existing v4 German forecasts
with the exact CRPS: the density stage only, no estimation (section 3.2).

> **Caveat (keep in the log, not in the paper).** Picking the transformation after looking at
> the test results is data snooping. These runs are for orientation only. For the paper, choose
> the transformation per country on a validation period before the evaluation sample
> (e.g. 2013–2014 origins), or average across transformations (Uniejewski, Weron & Ziel, 2018).
> See section 6.

> **Naming.** In the current code, `VST_MODE=raw` means **asinh of the raw price** (the v3
> transformation), not "prices in levels". The new modes are called `level` and `log` to avoid
> confusion. My quick check (review document, table 2.2) found that **asinh of the raw price**
> works best for Italy, so I recommend also running `VST_MODE=asinh` (one more Italian run,
> about 2 hours with the Mac awake).

---

## 1. Transformation modes (`common_utils.R`)

Every model is linear and Gaussian on a transformed scale z = g(P). The pipeline needs four
things for each g:
- `vst()`: the transform g;
- `vst_inv()`: its inverse, plus the mean correction;
- `vst_logjac()`: the log-Jacobian, for the log score;
- the quantile function, for the CRPS.

| `VST_MODE` | z = g(P) | median forecast | mean forecast (`__bc`) | log-Jacobian log g′(u) | domain |
|---|---|---|---|---|---|
| `normalized` (v4 default) | asinh((P − a)/b), a = median, b = MAD of the window | a + b sinh(μ) | a + b sinh(μ) e^{σ²/2} | −½ log((u − a)² + b²) | all u |
| `asinh` (alias `raw`, v3) | asinh(P) | sinh(μ) | sinh(μ) e^{σ²/2} | −½ log(u² + 1) | all u |
| `log` | log P | e^{μ} | e^{μ + σ²/2} | −log u | u > 0 only |
| `level` | P | μ | μ | 0 | all u |

The mean formulas assume z ~ N(μ, σ²) with σ² the h-step variance from
`iterated_error_variance()`, which is unchanged.
- asinh: E[sinh z] = sinh(μ) e^{σ²/2}.
- log: the lognormal mean.
- level: no retransformation, so no bias.

Replace the block at lines ~113–150 of `common_utils.R` with:
```r
VST_MODE <- tolower(Sys.getenv("VST_MODE", unset = "normalized"))
if (VST_MODE == "raw") VST_MODE <- "asinh"                       # backward-compatible alias
stopifnot(VST_MODE %in% c("normalized", "asinh", "log", "level"))

vst_params <- function(P) {
  if (VST_MODE != "normalized") return(list(a = 0, b = 1))      # (a, b) only used by normalized
  a <- stats::median(P, na.rm = TRUE)
  b <- stats::mad(P, na.rm = TRUE)                               # already scaled by 1.4826
  if (!is.finite(b) || b <= 0) b <- stats::sd(P, na.rm = TRUE)
  if (!is.finite(b) || b <= 0) b <- 1
  list(a = a, b = b)
}
vst <- function(P, p) switch(VST_MODE,
  normalized = asinh((P - p$a) / p$b),
  asinh      = asinh(P),
  log        = { if (any(P <= 0, na.rm = TRUE)) stop("VST_MODE=log needs strictly positive prices"); log(P) },
  level      = P)
# s2 = h-step variance on the z scale; s2 = 0 gives the median, s2 > 0 the mean.
vst_inv <- function(z, p, s2 = 0) switch(VST_MODE,
  normalized = p$a + p$b * sinh(z) * exp(s2 / 2),
  asinh      = sinh(z) * exp(s2 / 2),
  log        = exp(z + s2 / 2),
  level      = z)
# inverse WITHOUT mean correction, vectorised over matrices (quantiles, grids)
vst_inv0 <- function(z, a, b) switch(VST_MODE,
  normalized = a + b * sinh(z), asinh = sinh(z), log = exp(z), level = z)
# log g'(u): for the log score on the price scale
vst_logjac <- function(u, a, b) switch(VST_MODE,
  normalized = -0.5 * log((u - a)^2 + b^2),
  asinh      = -0.5 * log(u^2 + 1),
  log        = ifelse(u > 0, -log(u), -Inf),
  level      = 0 * u)
# z of a price, vectorised with per-origin (a, b)
vst_z <- function(u, a, b) switch(VST_MODE,
  normalized = asinh((u - a) / b), asinh = asinh(u),
  log = ifelse(u > 0, log(pmax(u, 1e-12)), -Inf), level = u)
```
- `vst_params_by_origin()` and `load_vst_params()` are unchanged. Outside `normalized` they
  store a = 0, b = 1.
- The estimation scripts (`ar_benchmark.R`, `arx_free_lags.R`, `r_midas.R`) already call only
  `vst()` / `vst_inv()`, so **no change is needed there**.
- The new modes are for **Italy only**. Germany stays on `normalized`. (`log` would stop on the
  German negative prices anyway.)
- **Fuel prices** stay in logs in every mode. Optional: in `level` mode, a switch
  `GAS_SCALE=level` would make gas enter in levels too. Not needed for this temporary run.

### 1.1 Replace the hard-coded asinh in the evaluation scripts
These lines compute "transformed-scale" errors. Make them use the model's own scale:

| File (line) | Now | Replace with |
|---|---|---|
| `summary_table_ar_r_midas.R` (~95) | `asinh((y_a - a_k) / b_k)` etc. | `vst_z(y_a, a_k, b_k)` etc. |
| `combination_subperiod_gr.R` (~215–216) | `asinh((yhat - a)/b) - asinh((y - a)/b)` | `vst_z(yhat, a, b) - vst_z(y, a, b)` |
| `tests_dm_gw.R` (~81–82) | `asinh((y - a)/b)`, `asinh((f - a)/b)` | `vst_z(y, a, b)`, `vst_z(f, a, b)` |

The column names (`rmse_asinh`, …) can stay. They now mean "on the model's transformed
scale" and are **not comparable across runs** with different `VST_MODE`. Across runs, compare
only on the EUR scale (RMSE, MAE, CRPS, log score).

### 1.2 Make the evaluation stages respect `COUNTRY_NAME`
`density_eval.R` (l. ~315), `tests_dm_gw.R` (l. ~229), `cssed_plots.R` (l. ~110) and
`summary_table_ar_r_midas.R` (l. ~239, ~255) loop over `c("de", "it")` whatever
`COUNTRY_NAME` says. In an Italy-only folder there are no German forecasts. In `density_eval.R`
that gives `min(numeric(0)) = Inf` for the CRPS range, and German rows full of NA. Add to
`common_utils.R`:
```r
eval_countries <- function() switch(Sys.getenv("COUNTRY_NAME", unset = "ALL"),
  Germany = "de", Italy = "it", c("de", "it"))
```
Then replace `c("de", "it")` with `eval_countries()` in those five places. In `density_eval.R`,
also write the top-5 file per country, so a German-only rerun does not overwrite the Italian
rows:
```r
bp <- file.path(RESULT_DIR, sprintf("density_top5_per_country_h_%s.csv",
                                    paste(eval_countries(), collapse = "_")))
```

### 1.3 Density stage (`density_eval.R`)
Three edits make the density code transformation-agnostic:
```r
# integrate_scores(), line ~99: grid of price values
U <- vst_inv0(outer(mu, rep(1, G)) + outer(sigma, ZGRID), a, b)
# line ~102: z of the outcome
z_y <- (vst_z(y, a, b) - mu) / sigma
# score_block(), line ~141 and ~145
mu <- vst_z(point, a, b)                                    # point = median forecast
S[ok, "logS"] <- dnorm(r$z_y, log = TRUE) - log(sigma[ok]) + vst_logjac(y[ok], a[ok], b[ok])
```
- `point` must be the **median** forecast (the plain column, not `__bc`). This is already the case.
- In `log` mode the grid is strictly positive, so a negative `umin` simply leaves an empty
  left tail.

---

## 2. Exact CRPS (supervisors' request)

### 2.1 What the code computes now
`density_eval.R` already reports a column `CRPS`. It is integrated only over a fixed range
[`umin`, `umax`] = [min(y) − 100, max(y) + 200]. That is a **threshold-weighted CRPS** with an
indicator weight (Gneiting & Ranjan, 2011). It is proper, but it is **not the standard CRPS**:
the predictive mass outside the range is ignored. This matters most for wide densities at
h = 21–28 and for comparisons across transformations.

**Change:** rename the current column to `CRPS_trunc` (kept, since ACPS uses the same range),
and add the exact CRPS as `CRPS`:

  CRPS(F, y) = ∫ (F(u) − 1{u ≥ y})² du = 2 ∫₀¹ (1{y < Q(τ)} − τ)(Q(τ) − y) dτ

The second form is the quantile representation (Laio & Tamea, 2007), which works for any
transformation.

### 2.2 Functions (add to `common_utils.R` or the top of `density_eval.R`)
```r
# closed forms (used for level / log, and as unit tests)
crps_norm  <- function(y, mu, s) { z <- (y - mu) / s; s * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi)) }
crps_lnorm <- function(y, mu, s) {                  # y > 0; Baran & Lerch (2015), as in scoringRules::crps_lnorm
  w <- (log(y) - mu) / s
  y * (2 * pnorm(w) - 1) - 2 * exp(mu + s^2 / 2) * (pnorm(w - s) + pnorm(s / sqrt(2)) - 1)
}
# quantile representation, any VST_MODE (mu, s on the z scale; a, b per origin)
CRPS_K   <- 1000L
CRPS_TAU <- (seq_len(CRPS_K) - 0.5) / CRPS_K       # midpoint rule
CRPS_ZQ  <- qnorm(CRPS_TAU)
crps_quantile <- function(y, mu, s, a, b, chunk = 500L) {
  n <- length(y); out <- numeric(n)
  for (st in seq(1L, n, by = chunk)) {               # chunks keep memory at chunk x K
    i <- st:min(n, st + chunk - 1L)
    Q <- vst_inv0(outer(mu[i], rep(1, CRPS_K)) + outer(s[i], CRPS_ZQ), a[i], b[i])
    I <- Q > y[i]                                    # y < Q(tau)
    out[i] <- 2 * rowMeans((I - rep(CRPS_TAU, each = length(i))) * (Q - y[i]))
  }
  out
}
crps_exact <- function(y, mu, s, a, b) switch(VST_MODE,
  level = crps_norm(y, mu, s),
  log   = crps_lnorm(y, mu, s),
  crps_quantile(y, mu, s, a, b))                     # normalized / asinh: no closed form
```
- **Accuracy:** I checked these formulas against direct numerical integration. The closed forms
  match to machine precision. The quantile rule with K = 1000 is within about 1e-6 (relative)
  for Gaussian, lognormal and normalised-asinh densities.
- **Cost:** about 4 million operations per model-horizon, roughly 20–40 minutes for the whole
  density stage. Use K = 500 for a quick test run.

### 2.3 Wire it into `score_block()` and the tests
1. `SCORE_COLS <- c("logS", "CRPS", "CRPS_trunc", paste0("ACPS_", ACPS_C))`.
   In `GFUNS`, rename the entry `CRPS` to `CRPS_trunc`.
2. In `score_block()`, after the integration:
   `S[ok, "CRPS"] <- crps_exact(y[ok], mu[ok], sigma[ok], a[ok], b[ok])`.
3. In `period_rows()`, both are losses:
   `sign_ <- if (sc %in% c("CRPS", "CRPS_trunc")) 1 else -1`.
4. The twin tests (`density_tests_vs_twin_<cc>.csv`) loop over `logS` and `CRPS`: keep `CRPS`
   (now the exact one) and add `CRPS_trunc` only if wanted.
5. **Unit test** at the start of `density_eval.R`. It stops if the CRPS code is wrong:
```r
local({ old <- VST_MODE
  VST_MODE <<- "level"; stopifnot(abs(crps_quantile(71, 50, 12, 0, 1) - crps_norm(71, 50, 12)) < 1e-3)
  VST_MODE <<- "log";   stopifnot(abs(crps_quantile(95, log(60), .35, 0, 1) - crps_lnorm(95, log(60), .35)) < 1e-3)
  VST_MODE <<- old })
```

### 2.4 CRPS outputs for the paper
New file `results/crps_summary_<cc>.csv`, built from `density_scores_<cc>.csv`, with one row
per country × h × period × family × model. Columns:
- `crps_mean` (EUR/MWh) and `crps_ratio` (against `ar_dum`);
- `crps_skill` = 1 − `crps_ratio`;
- `dm_stat, dm_pval, gw_stat, gw_pval, gw_dbar` (against `ar_dum`);
- the same against the twin, for macro specifications.

Tables:
1. **Main text:** per country × horizon, the CRPS ratio and DM/GW significance for:
   - AR + gas, AR + gas + Brent;
   - the combinations;
   - the pre-specified macro models (all surveys / all IPI / all monthly).
2. **Twin summary:** per country × horizon × period, the median CRPS ratio against the twin and
   the share of macro specifications with a significantly lower CRPS (GW and DM), overall and
   by indicator. This is the density counterpart of Table B.
3. **Calibration:** next to the CRPS tables, report PIT sd, 90% coverage and the KS p-value
   (already computed).

---

## 3. Runs

### 3.1 Italy: the transformation variants

One folder per variant, built like `estimation_v4/`: scripts as symlinks, and `data -> ../data_v4`.
```bash
cd "<repo>/r-midas/r_pipeline/midas_codes"
for v in level log asinh; do
  mkdir -p estimation_v4_it_$v && cd estimation_v4_it_$v
  for f in ../*.R; do ln -sf "$f" .; done
  ln -sfn ../data_v4 data
  cd ..
done

# run each one with the Mac kept awake (about 2-3 h per variant for Italy only)
cd estimation_v4_it_level && COUNTRY_NAME=Italy VST_MODE=level LF_IT_SURVEY=iesi_it \
  caffeinate -i Rscript run_all.R 2>&1 | tee run_it_level.log
cd ../estimation_v4_it_log && COUNTRY_NAME=Italy VST_MODE=log LF_IT_SURVEY=iesi_it \
  caffeinate -i Rscript run_all.R 2>&1 | tee run_it_log.log
cd ../estimation_v4_it_asinh && COUNTRY_NAME=Italy VST_MODE=asinh LF_IT_SURVEY=iesi_it \
  caffeinate -i Rscript run_all.R 2>&1 | tee run_it_asinh.log    # recommended third variant
```
- Keep everything else as in v4: `RMIDAS_THETA2_MAX=8`, `WINDOW_DAYS=3650`, all four families,
  all stages.
- **Before each full run, do a quick test:**
  `COUNTRY_NAME=Italy SPEC_LIMIT=8 STAGES=ar,arx,summary,density`.
  Check that the forecasts are in a sensible EUR range and that the CRPS unit test passes.

### 3.2 CRPS only, in the existing `estimation_v4/` folder (no re-estimation)
The forecasts and the `_sd` files are already there, so only the density stage runs (about
20–40 minutes):
```bash
cd "<repo>/r-midas/r_pipeline/midas_codes/estimation_v4"
COUNTRY_NAME=Germany VST_MODE=normalized STAGES=density Rscript run_all.R 2>&1 | tee crps_de.log
COUNTRY_NAME=Italy   VST_MODE=normalized STAGES=density Rscript run_all.R 2>&1 | tee crps_it.log
```
- **Germany:** this gives the CRPS tables of section 2.4, and nothing else changes.
- **Italy:** this gives the exact CRPS for the normalised v4 run, which is the baseline for
  the comparison in section 4.
- `VST_MODE=normalized` must be set explicitly, so the CRPS is computed with the transformation
  the forecasts were produced with.

---

## 4. Comparing the runs (`compare_vst_runs.R`, new)
Reads the Italian results of `estimation_v4` (normalised) and of `estimation_v4_it_{level,log,asinh}`.
For each h and each of `ar_dum`, `ar_dum_gas`, `comb_eq` (each family), the best pre-specified
macro model and the L7 no-macro model, it writes `results_compare/compare_vst_it.csv` with:
- RMSE, MAE (EUR), mean exact CRPS, mean log score, 90% coverage, PIT sd, per transformation;
- for each model, **DM-HLN and GW tests between transformations** (e.g. `ar_dum` level against
  `ar_dum` log). The errors share the same targets, so `dm_hln(e_level, e_log, h)` and
  `gw_test(...)` apply directly; for CRPS and log score use `dm_hln_loss` / `gw_test_loss`;
- the same by subperiod (pre / crisis / post).

Comparisons across transformations are valid only on the EUR scale: RMSE, MAE, CRPS, and the
log score (with the Jacobian, as in 1.3). Do not compare the `_asinh` columns across runs.

---

## 5. Checks
- [ ] `VST_MODE=normalized` reproduces the v4 forecasts exactly (refactoring check).
- [ ] `VST_MODE=asinh` reproduces the v3 transformation (and roughly my quick check: AR RMSE
      13.2 / 26.1 / 35.1 / 42.9 at h = 1/7/14/28).
- [ ] `log` runs for Italy without error (min Italian price 10.66 EUR/MWh).
- [ ] The German CRPS rerun leaves `density_scores_it.csv` and the German point-forecast
      tables untouched (check the file timestamps).
- [ ] CRPS unit tests pass; the exact CRPS is ≥ `CRPS_trunc` for every observation (it adds the
      tails).
- [ ] For `level`, the forecasts and the `__bc` forecasts coincide.
- [ ] Italy with gas, h = 28: does the blow-up under `normalized` (RMSE 58.3) disappear under
      `log` / `asinh`?

## 6. Later (for the paper)
- Choose the Italian transformation (and, for consistency, confirm Germany's) on a validation period before 2015 (e.g. 2013–2014
  origins, by mean CRPS of the AR benchmark), and state the rule. Or average the forecasts (and
  mix the densities) across transformations.
- Report the CRPS in the main text next to the RMSE, and the truncated CRPS / ACPS only in the
  appendix.
