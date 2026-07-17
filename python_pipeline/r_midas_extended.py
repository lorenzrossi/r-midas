#!/usr/bin/env python3
# =============================================================================
# PYTHON PIPELINE - R-MIDAS EXTENDED (EXP-ALMON ON y AND ON EACH HF BLOCK)
# Port of r_pipeline/midas_codes/r_midas_extended.R.
#
# Almon polynomial on the Q=14 daily AR lags of y AND on the Q=14 daily lags
# of each included HF exogenous block (gas_log, brent_log).  LF columns and
# calendar controls enter linearly.
#
# NLS: concentrated, L-BFGS-B with theta2 <= 0 per block, warm start plus
# deterministic canonical shape restarts (n_starts = 5 + 5*(K-1), capped).
#
# Iterated multi-step (option c.1): AR-Almon feedback; HF Almon z values and
# LF columns frozen at row (origin + 1); calendar refreshed at target.
#
# Outputs (same layout as the R script):
#   forecasts/r_midas_extended_<cc>_h<h>.csv
#   forecasts/r_midas_extended_spec_dictionary_<cc>.csv
#   results/coefficients/r_midas_extended_<cc>.csv
# =============================================================================

import os
import time
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

import common_utils as cu


def n_starts_for_K(K):
    return 5 + 5 * (K - 1)


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

    cal_mat, cal_cols = cu.build_calendar_dummies(dates, country)
    X_base = np.hstack([np.ones((n, 1)), cal_mat])
    base_cols = ["intercept"] + cal_cols
    hol_set = cu.country_holidays(dates.dt.year.unique(), country)

    origins = cu.make_origins(n, cu.WINDOW_DAYS,
                              pre_lag_max=max(cu.AR_LAGS + [cu.Q_HF]),
                              h_max=cu.H_MAX)
    masks = cu.build_spec_grid(cu.HF_CANDIDATES, lf_names, cu.spec_limit())

    return {"cc": cc, "country": country, "daily": daily, "dates": dates,
            "n": n, "y": y, "Z_full": Z_full, "lf_names": lf_names,
            "lf_mat": lf_mat, "X_base": X_base, "base_cols": base_cols,
            "hol_set": hol_set, "origins": origins, "masks": masks,
            "discount": cu.discount_lambda()}


def _init_worker(country):
    global _CTX
    _CTX = _build_context(country)


def fit_one_spec(si):
    c = _CTX
    masks = c["masks"]
    mk = masks.iloc[si]
    lf_names = c["lf_names"]
    hf_used = [v for v in ("gas", "brent") if mk[v] == 1]
    lf_used = [v for v in lf_names if mk[v] == 1]
    lf_idx = [lf_names.index(v) for v in lf_used]
    K = 1 + len(hf_used)
    n_st = n_starts_for_K(K)

    origins = c["origins"]
    y = c["y"]
    n_orig = len(origins)
    coef_cols = ["delta_y", "delta_gas", "delta_brent"] + lf_names
    yhat_mat = np.full((n_orig, cu.H_MAX), np.nan)
    coef_mat = np.full((n_orig, len(coef_cols)), np.nan)
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
        if keep.sum() < X_lin.shape[1] + 2 * K + 5:
            continue
        Xf = X_lin[keep]
        yf = ywin[keep]
        Zf = {nm: M[keep] for nm, M in Z_list.items()}
        if c["discount"] < 1:
            sw = np.sqrt(c["discount"] ** (t0 - np.arange(lo, t0 + 1)))[keep]
            Xf = Xf * sw[:, None]
            yf = yf * sw
            Zf = {nm: M * sw[:, None] for nm, M in Zf.items()}

        if theta_warm is not None and len(theta_warm) != 2 * len(Zf):
            theta_warm = None
        fit = cu.fit_r_midas_fast(yf, Zf, Xf, lin_cols,
                                  theta_init=theta_warm,
                                  optim_maxit=2000, optim_factr=5e8,
                                  n_starts=n_st)
        if fit is None:
            continue
        theta_warm = np.concatenate([fit["theta"][nm] for nm in Zf])

        coef_mat[oi, 0] = fit["delta"]["y_lags"]
        if "gas" in hf_used:
            coef_mat[oi, 1] = fit["delta"]["gas"]
        if "brent" in hf_used:
            coef_mat[oi, 2] = fit["delta"]["brent"]
        for v in lf_used:
            b = fit["beta_lin"].get(v, np.nan)
            if np.isfinite(b):
                coef_mat[oi, 3 + lf_names.index(v)] = b

        X_lin_row = (np.concatenate([c["X_base"][t0 + 1],
                                     c["lf_mat"][t0 + 1][lf_idx]])
                     if lf_used else c["X_base"][t0 + 1].copy())
        if not np.isfinite(X_lin_row).all():
            continue

        target_dates = [c["dates"].iloc[t0 + 1 + s] for s in range(cu.H_MAX)]
        yhat_mat[oi] = cu.iterated_forecast_r_midas(
            fit, c["Z_full"], X_lin_row, lin_cols,
            y[:t0 + 1], t0, target_dates,
            ar_block_name="y_lags", hf_block_names=hf_used,
            h_max=cu.H_MAX, holiday_set=c["hol_set"])
    return si, yhat_mat, coef_mat


def run_country(country):
    global _CTX
    _CTX = _build_context(country)
    c = _CTX
    cc = c["cc"]
    masks = c["masks"]
    spec_ids = list(masks.index)
    lf_names = c["lf_names"]
    n_cores = cu.get_n_cores()
    print(f"\n=== R-MIDAS extended, multi-h - {country} (n_cores={n_cores}) ===")

    dict_df = cu.spec_dictionary(masks, cu.HF_CANDIDATES, lf_names)
    dp = os.path.join(cu.FORECAST_DIR,
                      f"r_midas_extended_spec_dictionary_{cc}.csv")
    dict_df.to_csv(dp, index=False)
    print(f"  wrote {dp}  (specs = {len(dict_df)})")
    origins = c["origins"]
    print(f"  rolling origins: {len(origins)}")

    t0c = time.time()
    results = [None] * len(spec_ids)
    if n_cores > 1 and len(spec_ids) > 1:
        with ProcessPoolExecutor(max_workers=n_cores,
                                 initializer=_init_worker,
                                 initargs=(country,)) as ex:
            for si, ym, cm in ex.map(fit_one_spec, range(len(spec_ids))):
                results[si] = (ym, cm)
    else:
        for si in range(len(spec_ids)):
            _, ym, cm = fit_one_spec(si)
            results[si] = (ym, cm)
    print(f"  fitting done in {time.time() - t0c:.1f}s")

    y_lvl = c["daily"]["elec_price"].to_numpy(float)
    for h in cu.HORIZONS:
        cols = {sid: np.sinh(results[si][0][:, h - 1])
                for si, sid in enumerate(spec_ids)}
        fp = os.path.join(cu.FORECAST_DIR,
                          f"r_midas_extended_{cc}_h{h:02d}.csv")
        nrows = cu.write_forecast_csv(fp, c["dates"], origins, h, y_lvl, cols)
        print(f"  wrote {fp}  (rows = {nrows})")

    coef_cols = ["delta_y", "delta_gas", "delta_brent"] + lf_names
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
    cp = os.path.join(cu.COEF_DIR, f"r_midas_extended_{cc}.csv")
    coef_df.to_csv(cp, index=False)
    print(f"  wrote {cp}  (rows = {len(coef_df)})")


def main():
    country_arg = os.environ.get("COUNTRY_NAME", "ALL")
    countries = (["Germany", "Italy"] if country_arg.upper() == "ALL"
                 else [country_arg])
    for cn in countries:
        run_country(cn)
    print("\nR-MIDAS extended multi-h done.")


if __name__ == "__main__":
    main()
