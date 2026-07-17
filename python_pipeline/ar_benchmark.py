#!/usr/bin/env python3
# =============================================================================
# PYTHON PIPELINE - AR(1,2,3,7) BENCHMARK FAMILY
# Port of r_pipeline/midas_codes/ar_benchmark.R.
#
# Four specifications per country and horizon h in {1,7,14,21,28}:
#   ar_dum            : AR(1,2,3,7) + calendar dummies              [BENCHMARK]
#   ar_dum_gas        : ar_dum + lag 1 of log(gas)
#   ar_dum_brent      : ar_dum + lag 1 of log(brent)
#   ar_dum_gas_brent  : ar_dum + lag 1 of both fuels
#
# OLS on the rolling window; iterated multi-step (option c.1) with calendar
# refresh at the target row; forecasts back-transformed via sinh().
# Outputs: forecasts/ar_<cc>_h<h>.csv  (same layout as the R script).
# =============================================================================

import os
import time
import numpy as np

import common_utils as cu


def run_ar_country(country):
    cc = "de" if country == "Germany" else "it"
    print(f"\n=== AR benchmark - {country} ===")
    dat = cu.load_country_data(country)
    daily = dat["daily"]
    dates = daily["date"]
    n = len(daily)
    y = daily["y_asinh"].to_numpy()

    ar_mat, ar_cols   = cu.build_ar_lags(y, cu.AR_LAGS)
    cal_mat, cal_cols = cu.build_calendar_dummies(dates, country)
    gas_mat, gas_cols = cu.build_hf_exo_lags(daily["gas_log"].to_numpy(),
                                             "gas", 1)
    brt_mat, brt_cols = cu.build_hf_exo_lags(daily["brent_log"].to_numpy(),
                                             "brent", 1)
    hol_set = cu.country_holidays(dates.dt.year.unique(), country)

    X_full = np.hstack([np.ones((n, 1)), cal_mat, ar_mat, gas_mat, brt_mat])
    all_cols = ["intercept"] + cal_cols + ar_cols + gas_cols + brt_cols
    cidx = {c: i for i, c in enumerate(all_cols)}
    base_cols = ["intercept"] + cal_cols + ar_cols

    specs = {
        "ar_dum":           [],
        "ar_dum_gas":       gas_cols,
        "ar_dum_brent":     brt_cols,
        "ar_dum_gas_brent": gas_cols + brt_cols,
    }

    origins = cu.make_origins(n, cu.WINDOW_DAYS,
                              pre_lag_max=max(cu.AR_LAGS + [cu.Q_HF]),
                              h_max=cu.H_MAX)
    n_orig = len(origins)
    print(f"  rolling origins: {n_orig} "
          f"(from {dates.iloc[origins[0]].date()} "
          f"to {dates.iloc[origins[-1]].date()})")

    results = {}
    t0c = time.time()
    for spec_name, extra in specs.items():
        cols = base_cols + extra
        sel = [cidx[c] for c in cols]
        Xs = X_full[:, sel]
        ar_lag_idx = [cols.index(c) for c in ar_cols]
        yhat_mat = np.full((n_orig, cu.H_MAX), np.nan)
        sigma2_vec = np.full(n_orig, np.nan)
        for oi, t0 in enumerate(origins):
            lo = t0 - cu.WINDOW_DAYS + 1
            Xtr = Xs[lo:t0 + 1]
            ytr = y[lo:t0 + 1]
            ok = np.isfinite(Xtr).all(axis=1) & np.isfinite(ytr)
            if ok.sum() < Xtr.shape[1] + 5:
                continue
            beta, _, rank, _ = np.linalg.lstsq(Xtr[ok], ytr[ok], rcond=None)
            if rank < Xtr.shape[1]:
                continue
            res = ytr[ok] - Xtr[ok] @ beta
            sigma2_vec[oi] = float(np.mean(res**2))
            target_dates = [dates.iloc[t0 + 1 + s] for s in range(cu.H_MAX)]
            frow = Xs[t0 + 1].copy()
            yhat_mat[oi] = cu.iterated_forecast_linear(
                beta, cols, y[:t0 + 1], frow, ar_lag_idx, cu.AR_LAGS,
                target_dates, cu.H_MAX, hol_set)
        results[spec_name] = (yhat_mat, sigma2_vec)
    print(f"  fitting done in {time.time() - t0c:.1f}s")

    y_lvl = daily["elec_price"].to_numpy(float)
    for h in cu.HORIZONS:
        cols = {}
        for nm, (ymat, s2) in results.items():
            mu = ymat[:, h - 1]
            cols[nm] = np.sinh(mu)
            cols[f"{nm}__bc"] = np.sinh(mu) * np.exp(s2 / 2)
        fp = os.path.join(cu.FORECAST_DIR, f"ar_{cc}_h{h:02d}.csv")
        nrows = cu.write_forecast_csv(fp, dates, origins, h, y_lvl, cols)
        print(f"  wrote {fp}  (rows = {nrows})")


def main():
    country_arg = os.environ.get("COUNTRY_NAME", "ALL")
    countries = (["Germany", "Italy"] if country_arg.upper() == "ALL"
                 else [country_arg])
    for c in countries:
        run_ar_country(c)
    print("\nAR benchmark family done.")


if __name__ == "__main__":
    main()
