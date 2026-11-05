"""Average marginal contribution of each monthly indicator to forecast accuracy,
by country and horizon, across all R-MIDAS specifications.

For every indicator v and every pair of specifications that are identical
except that one of them also contains v (same daily regressors, same other
monthly indicators), compute 100 * [log RMSE(with v) - log RMSE(without v)].
Report the mean over pairs (negative = adding v lowers the RMSE, in %) and
the share of pairs in which adding v helps.

Usage: python3 marginal_contribution.py results/summary_table_ar_r_midas.csv
"""
import sys
import numpy as np, pandas as pd

f = sys.argv[1] if len(sys.argv) > 1 else "results/summary_table_ar_r_midas.csv"
d = pd.read_csv(f)
r = d[d.family == "R_MIDAS"].copy()
r["lfset"] = r.lf_vars.fillna("none").apply(lambda s: frozenset([] if s == "none" else s.split("|")))
r["hf"] = r.hf_vars.astype(str)
rows = []
for (c, h), g in r.groupby(["country", "h"]):
    lr = {(k, ls): np.log(x) for k, ls, x in zip(g.hf, g.lfset, g.rmse)}
    for v in sorted(set().union(*g.lfset)):
        diff = np.array([lr[(k, ls | {v})] - val for (k, ls), val in lr.items()
                         if v not in ls and (k, ls | {v}) in lr]) * 100
        rows.append(dict(country=c, h=h, var=v, pairs=len(diff),
                         mean_pct=diff.mean(), share_helps=(diff < 0).mean()))
o = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print("Mean change in RMSE (%) from adding the indicator (negative = improves)")
print(o.pivot_table(index=["country", "var"], columns="h", values="mean_pct").round(2))
print("\nShare of specification pairs in which adding the indicator lowers the RMSE")
print(o.pivot_table(index=["country", "var"], columns="h", values="share_helps").round(2))
o.to_csv(f.replace(".csv", "_marginal_contribution.csv"), index=False)
