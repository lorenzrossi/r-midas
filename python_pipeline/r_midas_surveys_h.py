#!/usr/bin/env python3
# =============================================================================
# PYTHON PIPELINE - R-MIDAS SURVEYS-H (OLDER POOLED VARIANT)
# Port of r_pipeline/midas_codes/r_midas_surveys_h.R.
#
# Differences from r_midas.py / r_midas_extended.py (kept on purpose, this is
# the older variant retained for comparability with earlier runs):
#   * Linear base: intercept + 11 month-of-year dummies ONLY (no weekday, no
#     holiday dummies).
#   * Exp-Almon on the AR block AND on each included HF block (gas, brent),
#     like the extended model.
#   * Single warm-started optimisation per origin (n_starts = 1), maxit 800
#     on the first origin of a spec and 200 afterwards.
#
# Outputs: forecasts/r_midas_surveys_<cc>_h<h>.csv and the spec dictionary,
# same layout as the R script.
# =============================================================================

import os
import time
import numpy as np
from concurrent.futures import ProcessPoolExecutor

import common_utils as cu

_CTX = {}


def _build_context(country):
    cc = "de" if country == "Germany" else "it"
    dat = cu.load_country_data(country)
    daily = dat["daily"]
    monthly = dat["monthly"]
    dates = daily["date"]
    n = len(daily)
    y = daily["y_asinh"].to_numpy()

    ar_mat, _    = cu.build_hf_exo_lags(y, "y", cu.Q_HF)
    gas_mat, _   = cu.build_hf_exo_lags(daily["gas_log"].to_numpy(),
                                        "gas", cu.Q_HF)
    brent_mat, _ = cu.build_hf_exo_lags(daily["brent_log"].to_numpy(),
                                        "brent", cu.Q_HF)
    Z_full = {"y_lags": ar_mat, "gas": gas_mat, "brent": brent_mat}

    lf_names = (cu.LF_CANDIDATES_DE if country == "Germany"
                else cu.LF_CANDIDATES_IT)
    lf_mat = np.column_stack([
        cu.make_daily_lf(dates, monthly, v, cu.LF_PUB[v],
                         locf=(cu.LF_PUB[v] == cu.PUB_LAG_SURVEY))
        for v in lf_names])

    mon_mat, mon_cols = cu.build_month_dummies(dates)
    X_base = np.hstack([np.ones((n, 1)), mon_mat])
    base_cols = ["intercept"] + mon_cols

    origins = cu.make_origins(n, cu.WINDOW_DAYS,
                              pre_lag_max=max(cu.AR_LAGS + [cu.Q_HF]),
                              h_max=cu.H_MAX)
    masks = cu.build_spec_grid(cu.HF_CANDIDATES, lf_names, cu.spec_limit())

    return {"cc": cc, "daily": daily, "dates": dates, "n": n, "y": y,
            "Z_full": Z_full, "lf_names": lf_names, "lf_mat": lf_mat,
            "X_base": X_base, "base_cols": base_cols,
            "origins": origins, "masks": masks}


def _init_worker(country):
    global _CTX
    _CTX = _build_context(country)


def fit_one_spec(si):
    c = _CTX
    mk = c["masks"].iloc[si]
    lf_names = c["lf_names"]
    hf_used = [v for v in ("gas", "brent") if mk[v] == 1]
    lf_used = [v for v in lf_names if mk[v] == 1]
    lf_idx = [lf_names.index(v) for v in lf_used]

    origins = c["origins"]
    y = c["y"]
    yhat_mat = np.full((len(origins), cu.H_MAX), np.nan)
    theta_warm = None
    lin_cols = list(c["base_cols"]) + lf_used

    for oi, t0 in enumerate(origins):
        lo = t0 - cu.WINDOW_DAYS + 1
        sl = slice(lo, t0 + 1)
        ywin = y[sl]
        Z_list = {"y_lags": c["Z_full"]["y_lags"][sl]}
        for nm in hf_used:
            Z_list[nm] = c["Z_full"][nm][sl]
        X_lin = (np.hstack([c["X_base"][sl], c["lf_mat"][sl][:, lf_idx]])
                 if lf_used else c["X_base"][sl])

        keep = np.isfinite(X_lin).all(axis=1) & np.isfinite(ywin)
        for nm in Z_list:
            keep &= np.isfinite(Z_list[nm]).all(axis=1)
        if keep.sum() < X_lin.shape[1] + 2 * (1 + len(hf_used)) + 5:
            continue
        Xf = X_lin[keep]
        yf = ywin[keep]
        Zf = {nm: M[keep] for nm, M in Z_list.items()}

        if theta_warm is not None and len(theta_warm) != 2 * len(Zf):
            theta_warm = None
        maxit = 800 if theta_warm is None else 200
        fit = cu.fit_r_midas_fast(yf, Zf, Xf, lin_cols,
                                  theta_init=theta_warm,
                                  optim_maxit=maxit, n_starts=1)
        if fit is None:
            continue
        theta_warm = np.concatenate([fit["theta"][nm] for nm in Zf])

        X_lin_row = (np.concatenate([c["X_base"][t0 + 1],
                                     c["lf_mat"][t0 + 1][lf_idx]])
                     if lf_used else c["X_base"][t0 + 1].copy())
        if not np.isfinite(X_lin_row).all():
            continue
        target_dates = [c["dates"].iloc[t0 + 1 + s] for s in range(cu.H_MAX)]
        yhat_mat[oi] = cu.iterated_forecast_r_midas(
            fit, c["Z_full"], X_lin_row, lin_cols, y[:t0 + 1], t0,
            target_dates, ar_block_name="y_lags", hf_block_names=hf_used,
            h_max=cu.H_MAX, holiday_set=None)
    return si, yhat_mat


def run_country(country):
    global _CTX
    _CTX = _build_context(country)
    c = _CTX
    cc = c["cc"]
    spec_ids = list(c["masks"].index)
    n_cores = cu.get_n_cores()
    print(f"\n=== R-MIDAS pooled (surveys-h, older variant) - {country} "
          f"(n_cores={n_cores}) ===")

    dict_df = cu.spec_dictionary(c["masks"], cu.HF_CANDIDATES, c["lf_names"])
    dp = os.path.join(cu.FORECAST_DIR,
                      f"r_midas_surveys_spec_dictionary_{cc}.csv")
    dict_df.to_csv(dp, index=False)
    origins = c["origins"]
    print(f"  rolling origins: {len(origins)}")

    t0c = time.time()
    results = [None] * len(spec_ids)
    if n_cores > 1 and len(spec_ids) > 1:
        with ProcessPoolExecutor(max_workers=n_cores,
                                 initializer=_init_worker,
                                 initargs=(country,)) as ex:
            for si, ym in ex.map(fit_one_spec, range(len(spec_ids))):
                results[si] = ym
    else:
        for si in range(len(spec_ids)):
            _, ym = fit_one_spec(si)
            results[si] = ym
    print(f"  fitting done in {time.time() - t0c:.1f}s")

    y_lvl = c["daily"]["elec_price"].to_numpy(float)
    for h in cu.HORIZONS:
        cols = {sid: np.sinh(results[si][:, h - 1])
                for si, sid in enumerate(spec_ids)}
        fp = os.path.join(cu.FORECAST_DIR,
                          f"r_midas_surveys_{cc}_h{h:02d}.csv")
        nrows = cu.write_forecast_csv(fp, c["dates"], origins, h, y_lvl, cols)
        print(f"  wrote {fp}  (rows = {nrows})")


def main():
    country_arg = os.environ.get("COUNTRY_NAME", "ALL")
    countries = (["Germany", "Italy"] if country_arg.upper() == "ALL"
                 else [country_arg])
    for cn in countries:
        run_country(cn)
    print("\nR-MIDAS pooled multi-h done.")


if __name__ == "__main__":
    main()
