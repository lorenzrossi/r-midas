#!/usr/bin/env python3
# =============================================================================
# PYTHON PIPELINE - SUMMARY TABLE
# Port of r_pipeline/midas_codes/summary_table.R.
#
# Reads every forecast CSV under forecasts/ and writes
#   results/summary_table.csv
#   results/summary_top5_per_country_h.csv
# with the same metric blocks as the R script: rmse / rmse_asinh / mae /
# mae_asinh / rmse_bc (bias-corrected back-transform), each with SE, ratio
# vs the ar_dum benchmark and DM-HLN stat/p-value.
# =============================================================================

import os
import numpy as np
import pandas as pd

import common_utils as cu

FAMILY_PATTERNS = {
    "RU_MIDAS":         "ru_midas_{cc}_h{h:02d}.csv",
    "R_MIDAS":          "r_midas_{cc}_h{h:02d}.csv",
    "R_MIDAS_extended": "r_midas_extended_{cc}_h{h:02d}.csv",
    "R_MIDAS_surveys":  "r_midas_surveys_{cc}_h{h:02d}.csv",
    # Real-time combinations from combination_subperiod_gr (run it first).
    "COMB":             "comb_r_midas_{cc}_h{h:02d}.csv",
}
DICT_FILES = {
    "RU_MIDAS":         "ru_midas_spec_dictionary_{cc}.csv",
    "R_MIDAS":          "r_midas_spec_dictionary_{cc}.csv",
    "R_MIDAS_extended": "r_midas_extended_spec_dictionary_{cc}.csv",
    "R_MIDAS_surveys":  "r_midas_surveys_spec_dictionary_{cc}.csv",
}


def _rmse_block(e_t, e_b, h):
    n = len(e_t)
    if n < 5:
        return dict(rmse=np.nan, rmse_se=np.nan, rmse_ratio=np.nan,
                    dm_stat=np.nan, dm_pval=np.nan)
    rmse = float(np.sqrt(np.mean(e_t**2)))
    rmse_b = float(np.sqrt(np.mean(e_b**2)))
    v_e2 = float(np.var(e_t**2, ddof=1))
    rmse_se = float(np.sqrt(v_e2 / (n * 4 * rmse**2))) if rmse > 0 else np.nan
    dm = cu.dm_hln(e_t, e_b, h, loss="sqerr")
    return dict(rmse=rmse, rmse_se=rmse_se,
                rmse_ratio=rmse / rmse_b if rmse_b > 0 else np.nan,
                dm_stat=dm["stat"], dm_pval=dm["pvalue"])


def _mae_block(e_t, e_b, h):
    n = len(e_t)
    if n < 5:
        return dict(mae=np.nan, mae_se=np.nan, mae_ratio=np.nan,
                    dm_stat=np.nan, dm_pval=np.nan)
    mae = float(np.mean(np.abs(e_t)))
    mae_b = float(np.mean(np.abs(e_b)))
    mae_se = float(np.std(np.abs(e_t), ddof=1) / np.sqrt(n))
    dm = cu.dm_hln(e_t, e_b, h, loss="abserr")
    return dict(mae=mae, mae_se=mae_se,
                mae_ratio=mae / mae_b if mae_b > 0 else np.nan,
                dm_stat=dm["stat"], dm_pval=dm["pvalue"])


def compute_metrics(y_actual, y_hat_test, y_hat_bench, h, y_hat_test_bc=None, y_hat_bench_bc=None):
    y_a = np.asarray(y_actual, float)
    y_t = np.asarray(y_hat_test, float)
    y_b = np.asarray(y_hat_bench, float)
    ok = np.isfinite(y_a) & np.isfinite(y_t) & np.isfinite(y_b)
    y_a, y_t, y_b = y_a[ok], y_t[ok], y_b[ok]

    e_t_p, e_b_p = y_a - y_t, y_a - y_b
    r_lev = _rmse_block(e_t_p, e_b_p, h)
    m_lev = _mae_block(e_t_p, e_b_p, h)

    ya_a, ya_t, ya_b = np.arcsinh(y_a), np.arcsinh(y_t), np.arcsinh(y_b)
    e_t_a, e_b_a = ya_a - ya_t, ya_a - ya_b
    r_a = _rmse_block(e_t_a, e_b_a, h)
    m_a = _mae_block(e_t_a, e_b_a, h)

    # Feasible correction is computed at forecast generation from each
    # origin's rolling-window residual variance. Never estimate it from the
    # full evaluation sample.
    if y_hat_test_bc is None or y_hat_bench_bc is None:
        r_bc = dict(rmse=np.nan, rmse_se=np.nan, rmse_ratio=np.nan,
                    dm_stat=np.nan, dm_pval=np.nan)
    else:
        tbc = np.asarray(y_hat_test_bc, float)[ok]
        bbc = np.asarray(y_hat_bench_bc, float)[ok]
        okbc = np.isfinite(tbc) & np.isfinite(bbc) & np.isfinite(y_a)
        r_bc = _rmse_block(y_a[okbc] - tbc[okbc],
                           y_a[okbc] - bbc[okbc], h)
    sig2_t = np.nan
    sig2_b = np.nan

    return {
        "n": int(len(y_a)),
        "rmse": r_lev["rmse"], "rmse_se": r_lev["rmse_se"],
        "rmse_ratio": r_lev["rmse_ratio"],
        "dm_stat": r_lev["dm_stat"], "dm_pval": r_lev["dm_pval"],
        "rmse_asinh": r_a["rmse"], "rmse_se_asinh": r_a["rmse_se"],
        "rmse_ratio_asinh": r_a["rmse_ratio"],
        "dm_stat_asinh": r_a["dm_stat"], "dm_pval_asinh": r_a["dm_pval"],
        "mae": m_lev["mae"], "mae_se": m_lev["mae_se"],
        "mae_ratio": m_lev["mae_ratio"],
        "dm_mae_stat": m_lev["dm_stat"], "dm_mae_pval": m_lev["dm_pval"],
        "mae_asinh": m_a["mae"], "mae_se_asinh": m_a["mae_se"],
        "mae_ratio_asinh": m_a["mae_ratio"],
        "dm_mae_stat_asinh": m_a["dm_stat"],
        "dm_mae_pval_asinh": m_a["dm_pval"],
        "rmse_bc": r_bc["rmse"], "rmse_bc_se": r_bc["rmse_se"],
        "rmse_bc_ratio": r_bc["rmse_ratio"],
        "dm_bc_stat": r_bc["dm_stat"], "dm_bc_pval": r_bc["dm_pval"],
        "sigma2_asinh_test": sig2_t, "sigma2_asinh_bench": sig2_b,
    }


def process_country_h(cc, h):
    ar_path = os.path.join(cu.FORECAST_DIR, f"ar_{cc}_h{h:02d}.csv")
    if not os.path.exists(ar_path):
        raise FileNotFoundError(ar_path)
    ar = pd.read_csv(ar_path)
    bench = ar["ar_dum"].to_numpy(float)
    rows = []

    def _bc(df, col):
        # Mirror the R behaviour: missing __bc column -> None -> NA metrics.
        return df[col] if col in df.columns else None

    ar_specs = [c for c in ar.columns
                if c not in ("origin_date", "target_date", "y_actual")
                and not c.endswith("__bc")]
    for sp in ar_specs:
        m = compute_metrics(ar["y_actual"], ar[sp], bench, h,
                            _bc(ar, f"{sp}__bc"), _bc(ar, "ar_dum__bc"))
        rows.append({"country": cc, "h": h, "family": "AR", "model": sp, **m})

    for fam, pat in FAMILY_PATTERNS.items():
        fp = os.path.join(cu.FORECAST_DIR, pat.format(cc=cc, h=h))
        if not os.path.exists(fp):
            continue
        mi = pd.read_csv(fp)
        common = ar.merge(mi, on="target_date", suffixes=("_ar", ""))
        if len(common) == 0:
            continue
        bench_a = common["ar_dum"].to_numpy(float)
        prefix = "comb_" if fam == "COMB" else "spec_"
        spec_cols = [c for c in mi.columns
                     if c.startswith(prefix) and not c.endswith("__bc")]
        for sp in spec_cols:
            m = compute_metrics(common["y_actual"], common[sp], bench_a, h,
                                _bc(common, f"{sp}__bc"),
                                _bc(common, "ar_dum__bc"))
            rows.append({"country": cc, "h": h, "family": fam,
                         "model": sp, **m})
    return pd.DataFrame(rows)


def main():
    blocks = []
    for cc in ("de", "it"):
        for h in cu.HORIZONS:
            print(f"processing {cc}  h={h} ...")
            try:
                blocks.append(process_country_h(cc, h))
            except FileNotFoundError as e:
                print(f"  skipping: missing {e}")
    summary = pd.concat(blocks, ignore_index=True)

    summary["hf_vars"] = np.nan
    summary["lf_vars"] = np.nan
    for cc in ("de", "it"):
        for fam, pat in DICT_FILES.items():
            dp = os.path.join(cu.FORECAST_DIR, pat.format(cc=cc))
            if not os.path.exists(dp):
                continue
            d = pd.read_csv(dp).set_index("spec_id")
            sel = (summary["country"] == cc) & (summary["family"] == fam)
            if not sel.any():
                continue
            summary.loc[sel, "hf_vars"] = \
                summary.loc[sel, "model"].map(d["hf_vars"]).to_numpy()
            summary.loc[sel, "lf_vars"] = \
                summary.loc[sel, "model"].map(d["lf_vars"]).to_numpy()

    summary = summary.sort_values(
        ["country", "h", "family", "rmse_ratio"], na_position="last")
    out = os.path.join(cu.RESULT_DIR, "summary_table.csv")
    summary.to_csv(out, index=False)
    print(f"\nwrote {out}  (rows = {len(summary)})")

    best = (summary[np.isfinite(summary["rmse_ratio"])]
            .sort_values("rmse_ratio")
            .groupby(["country", "h"], as_index=False)
            .head(5)
            .sort_values(["country", "h", "rmse_ratio"]))
    bp = os.path.join(cu.RESULT_DIR, "summary_top5_per_country_h.csv")
    best.to_csv(bp, index=False)
    print(f"wrote {bp}  (rows = {len(best)})")


if __name__ == "__main__":
    main()
