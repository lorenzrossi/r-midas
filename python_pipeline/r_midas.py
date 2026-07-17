#!/usr/bin/env python3
# =============================================================================
# PYTHON PIPELINE - R-MIDAS (POOLED, EXP-ALMON ON y ONLY)
# Port of r_pipeline/midas_codes/r_midas.R (with fit_r_midas_fast).
#
# Almon polynomial on the AR lag set J = {1, 2, 3, 7} of y (same lags as the
# AR benchmark; the actual j values enter the Almon exponent, model.tex eq.
# (almon)); gas_log_lag1 and brent_log_lag1 enter linearly; LF columns linear
# (publication-aware, surveys LOCF over missing releases); calendar controls
# 6 weekday + 11 month + 1 holiday; intercept.
#
# NLS: concentrated (profiled), L-BFGS-B with theta2 <= 0, warm start from
# the previous origin plus deterministic canonical shape restarts (n_starts
# = 5 like the R script; capped at 1 + 5 shapes).
#
# Iterated multi-step (option c.1): AR-Almon feedback; LF and HF lag-1
# columns frozen at row (origin + 1); calendar refreshed at the target row.
#
# Outputs (same layout as the R script):
#   forecasts/r_midas_<cc>_h<h>.csv
#   forecasts/r_midas_spec_dictionary_<cc>.csv
#   results/coefficients/r_midas_<cc>.csv
# =============================================================================

import os
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

import common_utils as cu

N_STARTS_R_MIDAS = 5

_CTX = {}


def _build_context(country):
    cc = "de" if country == "Germany" else "it"
    dat = cu.load_country_data(country)
    daily = dat["daily"]
    monthly = dat["monthly"]
    dates = daily["date"]
    n = len(daily)
    y = daily["y_asinh"].to_numpy()

    # Almon AR block over the lag set J = {1, 2, 3, 7} (model.tex, same lags
    # as the AR benchmark; the actual j values enter the Almon exponent).
    ar_mat, _ = cu.build_ar_lags(y, cu.AR_LAGS)

    gas_lag1 = np.concatenate([[np.nan], daily["gas_log"].to_numpy()[:-1]])
    brt_lag1 = np.concatenate([[np.nan], daily["brent_log"].to_numpy()[:-1]])

    lf_names = (cu.LF_CANDIDATES_DE if country == "Germany"
                else cu.LF_CANDIDATES_IT)
    lf_mat = np.column_stack([
        cu.make_daily_lf(dates, monthly, v, cu.LF_PUB[v],
                         locf=(cu.LF_PUB[v] == cu.PUB_LAG_SURVEY))
        for v in lf_names])

    cal_mat, cal_cols = cu.build_calendar_dummies(dates, country)
    X_base = np.hstack([np.ones((n, 1)), cal_mat])
    base_cols = ["intercept"] + cal_cols
    hol_set = cu.country_holidays(dates.dt.year.unique(), country)

    origins = cu.make_origins(n, cu.WINDOW_DAYS,
                              pre_lag_max=max(cu.AR_LAGS + [cu.Q_HF]),
                              h_max=cu.H_MAX)
    masks = cu.build_spec_grid(cu.HF_CANDIDATES, lf_names, cu.spec_limit())

    return {"cc": cc, "country": country, "daily": daily, "dates": dates,
            "n": n, "y": y, "ar_mat": ar_mat, "gas_lag1": gas_lag1,
            "brt_lag1": brt_lag1, "lf_names": lf_names, "lf_mat": lf_mat,
            "X_base": X_base, "base_cols": base_cols, "hol_set": hol_set,
            "origins": origins, "masks": masks,
            "discount": cu.discount_lambda()}


def _init_worker(country):
    global _CTX
    _CTX = _build_context(country)


def fit_one_spec(si):
    c = _CTX
    masks = c["masks"]
    sid = masks.index[si]
    mk = masks.iloc[si]
    lf_names = c["lf_names"]
    hf_used = [v for v in ("gas", "brent") if mk[v] == 1]
    lf_used = [v for v in lf_names if mk[v] == 1]
    lf_idx = [lf_names.index(v) for v in lf_used]

    origins = c["origins"]
    y = c["y"]
    n_orig = len(origins)
    coef_cols = ["delta_y"] + lf_names
    yhat_mat = np.full((n_orig, cu.H_MAX), np.nan)
    coef_mat = np.full((n_orig, len(coef_cols)), np.nan)
    sigma2_vec = np.full(n_orig, np.nan)
    theta_warm = None

    lin_cols = list(c["base_cols"]) + lf_used + \
        ([f"{v}_lag1" for v in hf_used])

    for oi, t0 in enumerate(origins):
        lo = t0 - cu.WINDOW_DAYS + 1
        sl = slice(lo, t0 + 1)
        ywin = y[sl]
        Zwin = c["ar_mat"][sl]

        parts = [c["X_base"][sl]]
        if lf_used:
            parts.append(c["lf_mat"][sl][:, lf_idx])
        if "gas" in hf_used:
            parts.append(c["gas_lag1"][sl].reshape(-1, 1))
        if "brent" in hf_used:
            parts.append(c["brt_lag1"][sl].reshape(-1, 1))
        X_lin = np.hstack(parts)

        keep = np.isfinite(X_lin).all(axis=1) & \
            np.isfinite(Zwin).all(axis=1) & np.isfinite(ywin)
        if keep.sum() < X_lin.shape[1] + 3:
            continue
        Xf = X_lin[keep]
        yf = ywin[keep]
        Zf = Zwin[keep]
        if c["discount"] < 1:
            sw = np.sqrt(c["discount"] ** (t0 - np.arange(lo, t0 + 1)))[keep]
            Xf = Xf * sw[:, None]
            yf = yf * sw
            Zf = Zf * sw[:, None]

        if theta_warm is not None and len(theta_warm) != 2:
            theta_warm = None
        fit = cu.fit_r_midas_fast(yf, {"y_lags": Zf}, Xf, lin_cols,
                                  theta_init=theta_warm,
                                  optim_maxit=2000, optim_factr=5e8,
                                  n_starts=N_STARTS_R_MIDAS,
                                  lag_sets={"y_lags": cu.AR_LAGS})
        if fit is None:
            continue
        theta_warm = fit["theta"]["y_lags"]
        sigma2_vec[oi] = fit["sigma2"]

        coef_mat[oi, 0] = fit["delta"]["y_lags"]
        for v in lf_used:
            b = fit["beta_lin"].get(v, np.nan)
            if np.isfinite(b):
                coef_mat[oi, 1 + lf_names.index(v)] = b

        row_parts = [c["X_base"][t0 + 1]]
        if lf_used:
            # Low-frequency information available at the forecast origin.
            row_parts.append(c["lf_mat"][t0][lf_idx])
        if "gas" in hf_used:
            row_parts.append([c["gas_lag1"][t0 + 1]])
        if "brent" in hf_used:
            row_parts.append([c["brt_lag1"][t0 + 1]])
        X_lin_row = np.concatenate([np.atleast_1d(p) for p in row_parts])
        if not np.isfinite(X_lin_row).all():
            continue

        target_dates = [c["dates"].iloc[t0 + 1 + s] for s in range(cu.H_MAX)]
        yhat_mat[oi] = cu.iterated_forecast_r_midas(
            fit, {"y_lags": c["ar_mat"]}, X_lin_row, lin_cols,
            y[:t0 + 1], t0, target_dates,
            ar_block_name="y_lags", hf_block_names=(),
            h_max=cu.H_MAX, holiday_set=c["hol_set"],
            lag_sets={"y_lags": cu.AR_LAGS})
    return si, yhat_mat, coef_mat, sigma2_vec


def run_country(country):
    global _CTX
    _CTX = _build_context(country)
    c = _CTX
    cc = c["cc"]
    masks = c["masks"]
    spec_ids = list(masks.index)
    lf_names = c["lf_names"]
    n_cores = cu.get_n_cores()
    print(f"\n=== R-MIDAS standard (Almon on y, linear HF) - {country} "
          f"(n_cores={n_cores}, n_starts={N_STARTS_R_MIDAS}) ===")

    dict_df = cu.spec_dictionary(masks, cu.HF_CANDIDATES, lf_names)
    dp = os.path.join(cu.FORECAST_DIR, f"r_midas_spec_dictionary_{cc}.csv")
    dict_df.to_csv(dp, index=False)
    print(f"  wrote {dp}  (specs = {len(dict_df)})")
    origins = c["origins"]
    print(f"  rolling origins: {len(origins)} "
          f"(from {c['dates'].iloc[origins[0]].date()} "
          f"to {c['dates'].iloc[origins[-1]].date()})")

    t0c = time.time()
    results = [None] * len(spec_ids)
    if n_cores > 1 and len(spec_ids) > 1:
        with ProcessPoolExecutor(max_workers=n_cores,
                                 initializer=_init_worker,
                                 initargs=(country,)) as ex:
            for si, ym, cm, s2 in ex.map(fit_one_spec, range(len(spec_ids))):
                results[si] = (ym, cm, s2)
    else:
        for si in range(len(spec_ids)):
            _, ym, cm, s2 = fit_one_spec(si)
            results[si] = (ym, cm, s2)
    print(f"  fitting done in {time.time() - t0c:.1f}s")

    y_lvl = c["daily"]["elec_price"].to_numpy(float)
    for h in cu.HORIZONS:
        cols = {}
        for si, sid in enumerate(spec_ids):
            mu = results[si][0][:, h - 1]
            s2 = results[si][2]
            cols[sid] = np.sinh(mu)
            cols[f"{sid}__bc"] = np.sinh(mu) * np.exp(s2 / 2)
        fp = os.path.join(cu.FORECAST_DIR, f"r_midas_{cc}_h{h:02d}.csv")
        nrows = cu.write_forecast_csv(fp, c["dates"], origins, h, y_lvl, cols)
        print(f"  wrote {fp}  (rows = {nrows})")

    coef_cols = ["delta_y"] + lf_names
    blocks = []
    odates = c["dates"].iloc[origins].dt.strftime("%Y-%m-%d").to_numpy()
    for si, sid in enumerate(spec_ids):
        cm = results[si][1]
        blk = pd.DataFrame({"origin_date": odates, "spec_id": sid})
        for j, col in enumerate(coef_cols):
            blk[col] = cm[:, j]
        blocks.append(blk)
    coef_df = pd.concat(blocks, ignore_index=True)
    coef_df = coef_df[~coef_df[coef_cols].isna().all(axis=1)]
    cp = os.path.join(cu.COEF_DIR, f"r_midas_{cc}.csv")
    coef_df.to_csv(cp, index=False)
    print(f"  wrote {cp}  (rows = {len(coef_df)})")


def main():
    country_arg = os.environ.get("COUNTRY_NAME", "ALL")
    countries = (["Germany", "Italy"] if country_arg.upper() == "ALL"
                 else [country_arg])
    for cn in countries:
        run_country(cn)
    print("\nR-MIDAS standard (Almon on y, linear HF) done.")


if __name__ == "__main__":
    main()
