#!/usr/bin/env python3
# =============================================================================
# PYTHON PIPELINE - FORECAST COMBINATION, SUBPERIOD EVALUATION AND
#                   GIACOMINI-ROSSI FLUCTUATION TEST
# Port of r_pipeline/midas_codes/combination_subperiod_gr.R.
# Post-processing only: reads the forecast CSVs (Python or R generated, same
# layout), no re-estimation.  Applied identically to Germany and Italy.
#
# 1. Combinations across the 256 specs: comb_eq (equal weights), comb_trim
#    (20% trimmed mean), comb_invmse (real-time inverse-MSE, Bates-Granger,
#    365-day trailing window of RESOLVED errors, min 60 resolved errors).
# 2. Subperiod evaluation vs ar_dum: full / pre (..2020) / crisis (2021-22) /
#    post (2023..); RMSE-MAE-asinh ratios, DM-HLN on the full sample.
# 3. Giacomini-Rossi (2010) two-sided fluctuation test (mu = 0.3), HAC
#    variance with rectangular kernel truncated at h-1.
#    NOTE: critical values transcribed from GR (2010) Table 1; verify against
#    the published table before final use in the paper.
#
# Outputs: forecasts/comb_<family>_<cc>_h<h>.csv,
#          results/summary_subperiod_<cc>.csv, results/gr_fluctuation_<cc>.csv
# =============================================================================

import os
import numpy as np
import pandas as pd

import common_utils as cu

# NOTE: FAMILIES is read at call time (inside run_country), not at import
# time, so orchestrators can set it after importing this module.
def _families():
    return [f.strip() for f in
            os.environ.get("FAMILIES", "r_midas,r_midas_extended").split(",")
            if f.strip()]


TRAIL_DAYS   = int(os.environ.get("TRAIL_DAYS", "365"))
TOP_N_SPECS  = int(os.environ.get("TOP_N_SPECS", "10"))
MIN_RESOLVED = int(os.environ.get("MIN_RESOLVED", "60"))
TRIM         = float(os.environ.get("TRIM", "0.2"))
GR_MU        = float(os.environ.get("GR_MU", "0.3"))
BENCH_COL    = "ar_dum"


def _comb_include_ar():
    """AR candidates added to the inverse-MSE combination POOL (not to the
    benchmark, which stays ar_dum, and not to comb_eq / comb_trim, which
    remain family-only averages).  With the AR forecasts in the pool the
    real-time inverse-MSE weights migrate towards the pure AR whenever the
    macro block turns harmful and back when it pays off.  No look-ahead.
    Set COMB_INCLUDE_AR="" to disable.  Read at call time so orchestrators
    can set the env var after import."""
    raw = os.environ.get("COMB_INCLUDE_AR", "ar_dum,ar_dum_gas")
    return [v.strip() for v in raw.split(",") if v.strip()]

PERIODS = {
    "full":   ("1900-01-01", "2100-01-01"),
    "pre":    ("1900-01-01", "2020-12-31"),
    "crisis": ("2021-01-01", "2022-12-31"),
    "post":   ("2023-01-01", "2100-01-01"),
}

GR_CV = pd.DataFrame({
    "mu":   np.arange(0.1, 1.0, 0.1),
    "cv05": [3.393, 3.179, 3.012, 2.890, 2.779, 2.634, 2.560, 2.433, 2.331],
    "cv10": [3.170, 2.948, 2.766, 2.626, 2.500, 2.356, 2.252, 2.130, 2.020],
})


def rmse_(e):
    e = e[np.isfinite(e)]
    return np.sqrt(np.mean(e**2)) if len(e) else np.nan


def mae_(e):
    e = e[np.isfinite(e)]
    return np.mean(np.abs(e)) if len(e) else np.nan


def trim_mean_rows(F, trim=TRIM):
    out = np.full(F.shape[0], np.nan)
    for i in range(F.shape[0]):
        x = F[i][np.isfinite(F[i])]
        if x.size == 0:
            continue
        x = np.sort(x)
        k = int(np.floor(trim * x.size))
        out[i] = x[k:x.size - k].mean()
    return out


def combine_invmse(F, y, odates, h, trail=TRAIL_DAYS, min_res=MIN_RESOLVED):
    n, _ = F.shape
    E2 = (F - y[:, None])**2
    ok = np.isfinite(E2)
    E2z = np.where(ok, E2, 0.0)
    CS = np.cumsum(E2z, axis=0)
    CN = np.cumsum(ok, axis=0)
    od = odates.astype("datetime64[D]").astype(int)
    out = np.full(n, np.nan)
    eqm = np.nanmean(F, axis=1)
    for i in range(n):
        b = np.searchsorted(od, od[i] - h, side="right") - 1
        a = np.searchsorted(od, od[i] - h - trail, side="right") - 1
        if b <= a or b < 0:
            out[i] = eqm[i]
            continue
        cnt = CN[b] - (CN[a] if a >= 0 else 0)
        sse = CS[b] - (CS[a] if a >= 0 else 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            mse = np.where(cnt >= min_res, sse / np.maximum(cnt, 1), np.nan)
            w = 1.0 / mse
        fi = F[i]
        use = np.isfinite(fi) & np.isfinite(w) & (w > 0)
        out[i] = (np.sum(w[use] * fi[use]) / np.sum(w[use])
                  if use.any() else eqm[i])
    return out


def gr_fluctuation(e_bench, e_model, tdates, h, mu=GR_MU):
    ok = np.isfinite(e_bench) & np.isfinite(e_model)
    d = e_bench[ok]**2 - e_model[ok]**2        # > 0  ->  model better
    td = np.asarray(tdates)[ok]
    n = len(d)
    if n < 100:
        return None
    m = max(20, int(np.floor(mu * n)))
    dbar = d.mean()
    s = ((d - dbar)**2).mean()
    for k in range(1, h):
        s += 2 * np.sum((d[k:] - dbar) * (d[:-k] - dbar)) / n
    if not np.isfinite(s) or s <= 0:
        return None
    cs = np.concatenate([[0.0], np.cumsum(d)])
    Ft = (cs[m:n + 1] - cs[0:n - m + 1]) / (np.sqrt(s) * np.sqrt(m))
    ends = td[m - 1:n]
    i_max, i_min = int(np.argmax(Ft)), int(np.argmin(Ft))
    mu_eff = m / n
    cv05 = float(np.interp(mu_eff, GR_CV["mu"], GR_CV["cv05"]))
    cv10 = float(np.interp(mu_eff, GR_CV["mu"], GR_CV["cv10"]))
    return {"stat_max": float(Ft[i_max]), "date_max": str(ends[i_max])[:10],
            "stat_min": float(Ft[i_min]), "date_min": str(ends[i_min])[:10],
            "mu": mu_eff, "m": m, "n": n, "cv05": cv05, "cv10": cv10}


def eval_model(label, family, h, yhat, y, ybench, tdates):
    e_m = yhat - y
    e_b = ybench - y
    ea_m = np.arcsinh(yhat) - np.arcsinh(y)
    ea_b = np.arcsinh(ybench) - np.arcsinh(y)
    rows = []
    for pn, (d0, d1) in PERIODS.items():
        sel = (tdates >= np.datetime64(d0)) & (tdates <= np.datetime64(d1))
        if not sel.any():
            continue
        ok_l = sel & np.isfinite(e_m) & np.isfinite(e_b)
        ok_a = sel & np.isfinite(ea_m) & np.isfinite(ea_b)
        if not ok_l.any():
            continue
        if pn == "full":
            dm_lvl = cu.dm_hln(e_m[ok_l], e_b[ok_l], h)
            dm_as = (cu.dm_hln(ea_m[ok_a], ea_b[ok_a], h)
                     if ok_a.any() else {"stat": np.nan, "pvalue": np.nan})
        else:
            dm_lvl = dm_as = {"stat": np.nan, "pvalue": np.nan}
        rows.append({
            "family": family, "model": label, "h": h, "period": pn,
            "n": int(np.sum(ok_l)),
            "rmse": rmse_(e_m[ok_l]), "rmse_bench": rmse_(e_b[ok_l]),
            "rmse_ratio": rmse_(e_m[ok_l]) / rmse_(e_b[ok_l]),
            "mae_ratio": mae_(e_m[ok_l]) / mae_(e_b[ok_l]),
            "rmse_ratio_asinh": rmse_(ea_m[ok_a]) / rmse_(ea_b[ok_a]),
            "dm_stat": dm_lvl["stat"], "dm_pval": dm_lvl["pvalue"],
            "dm_stat_asinh": dm_as["stat"], "dm_pval_asinh": dm_as["pvalue"],
        })
    return rows


def run_country(country):
    cc = "de" if country == "Germany" else "it"
    print(f"\n=== Combination + subperiod + GR fluctuation - {country} ===")
    sum_rows, gr_rows = [], []

    for h in cu.HORIZONS:
        ar_fp = os.path.join(cu.FORECAST_DIR, f"ar_{cc}_h{h:02d}.csv")
        if not os.path.exists(ar_fp):
            print(f"  missing {ar_fp} - skipped")
            continue
        ar = pd.read_csv(ar_fp, parse_dates=["origin_date", "target_date"])

        # AR variants vs ar_dum
        y = ar["y_actual"].to_numpy(float)
        yb = ar[BENCH_COL].to_numpy(float)
        td = ar["target_date"].to_numpy()
        for mcol in [c for c in ar.columns
                     if c not in ("origin_date", "target_date",
                                  "y_actual", BENCH_COL)
                     and not c.endswith("__bc")]:
            ym = ar[mcol].to_numpy(float)
            sum_rows += eval_model(mcol, "AR", h, ym, y, yb, td)
            gr = gr_fluctuation(yb - y, ym - y, td, h)
            if gr:
                gr_rows.append({"family": "AR", "model": mcol, "h": h, **gr})

        for fam in _families():
            fp = os.path.join(cu.FORECAST_DIR, f"{fam}_{cc}_h{h:02d}.csv")
            if not os.path.exists(fp):
                print(f"  missing {fp} - skipped")
                continue
            fc = pd.read_csv(fp, parse_dates=["origin_date", "target_date"])
            ar_pool = [c for c in _comb_include_ar()
                       if c in ar.columns and not c.endswith("__bc")]
            ar_cols = list(dict.fromkeys(
                ["origin_date", "target_date", "y_actual", BENCH_COL]
                + ar_pool))
            mg = ar[ar_cols] \
                .merge(fc.drop(columns=["y_actual"]),
                       on=["origin_date", "target_date"]) \
                .sort_values("origin_date")
            specs = [c for c in mg.columns if c.startswith("spec_")
                     and not c.endswith("__bc")]  # exclude bias-corrected twins
            if not specs:
                continue
            F = mg[specs].to_numpy(float)
            y2 = mg["y_actual"].to_numpy(float)
            yb2 = mg[BENCH_COL].to_numpy(float)
            td2 = mg["target_date"].to_numpy()

            # comb_eq / comb_trim: R-MIDAS family only.
            # comb_invmse_fam: inverse-MSE combination of the family only.
            # comb_invmse: TWO-STAGE.  Stage 1 = comb_invmse_fam.  Stage 2 =
            # inverse-MSE over the SMALL pool {stage-1 combo, AR candidates};
            # with 3 candidates the weights can concentrate on the AR when
            # the macro block hurts (in a 258-strong pool the AR weight is
            # diluted to ~2/258 and the fallback never bites).
            comb_eq = np.nanmean(F, axis=1)
            comb_trim = trim_mean_rows(F)
            comb_invmse_fam = combine_invmse(F, y2,
                                             mg["origin_date"].to_numpy(), h)
            if ar_pool:
                F_pool = np.column_stack(
                    [comb_invmse_fam, mg[ar_pool].to_numpy(float)])
                comb_invmse = combine_invmse(F_pool, y2,
                                             mg["origin_date"].to_numpy(), h)
            else:
                comb_invmse = comb_invmse_fam

            combo = pd.DataFrame({
                "origin_date": mg["origin_date"].dt.strftime("%Y-%m-%d"),
                "target_date": mg["target_date"].dt.strftime("%Y-%m-%d"),
                "y_actual": y2, "comb_eq": comb_eq,
                "comb_trim": comb_trim,
                "comb_invmse_fam": comb_invmse_fam,
                "comb_invmse": comb_invmse})
            cfp = os.path.join(cu.FORECAST_DIR,
                               f"comb_{fam}_{cc}_h{h:02d}.csv")
            combo.to_csv(cfp, index=False)

            full_rmse = {s: rmse_(F[:, j] - y2)
                         for j, s in enumerate(specs)}
            top_specs = sorted(full_rmse, key=full_rmse.get)[:TOP_N_SPECS]
            best_spec = top_specs[0]
            models = {f"{s}_expost": F[:, specs.index(s)]
                      for s in top_specs}
            models.update({
                "comb_eq": comb_eq, "comb_trim": comb_trim,
                "comb_invmse_fam": comb_invmse_fam,
                "comb_invmse": comb_invmse})
            for lbl, ym in models.items():
                sum_rows += eval_model(lbl, fam, h, ym, y2, yb2, td2)
                gr = gr_fluctuation(yb2 - y2, ym - y2, td2, h)
                if gr:
                    gr_rows.append({"family": fam, "model": lbl,
                                    "h": h, **gr})
            print(f"  {fam} h={h:02d}: best ex-post spec {best_spec}; "
                  f"combos written to {os.path.basename(cfp)}")

    sum_df = pd.DataFrame(sum_rows)
    sum_df["country"] = cc
    sp = os.path.join(cu.RESULT_DIR, f"summary_subperiod_{cc}.csv")
    sum_df.to_csv(sp, index=False)
    print(f"  wrote {sp}  (rows = {len(sum_df)})")

    if gr_rows:
        gr_df = pd.DataFrame(gr_rows)
        gr_df["country"] = cc
        gr_df["reject_05"] = np.maximum(
            gr_df["stat_max"], -gr_df["stat_min"]) > gr_df["cv05"]
        gp = os.path.join(cu.RESULT_DIR, f"gr_fluctuation_{cc}.csv")
        gr_df.to_csv(gp, index=False)
        print(f"  wrote {gp}  (rows = {len(gr_df)})")


def main():
    country_arg = os.environ.get("COUNTRY_NAME", "ALL")
    countries = (["Germany", "Italy"] if country_arg.upper() == "ALL"
                 else [country_arg])
    for cn in countries:
        run_country(cn)
    print("\nCombination + subperiod + GR fluctuation done.")


if __name__ == "__main__":
    main()
