"""Descriptive statistics for the Descriptive Statistics section of Chapter 2 (data_v4).

Run from anywhere:  python3 descriptive_stats.py
Needs: pandas, numpy, scipy, statsmodels.
Reads midas_codes/data_v4/prices_{de,it}_daily.csv and dataset_{de,it}_m.csv (TODO_v4c, item 2):
real weekend prices (Germany from 2012-01-09, Italy throughout), no interpolation.
Every table is printed and also written as CSV to desc_v4/ next to this script, so the
numbers in descriptive_statistics.tex / corrections_v3.tex can be checked line by line.

Tables (names as in the paper):
  tab:desc-daily     moments of the daily series: price, asinh(price), normalised asinh
                     (full-sample median / MAD), log(price) for Italy, gas TTF, Brent
  tab:desc-monthly   moments of the monthly indicators as they enter the models
                     (surveys in levels, IPI in monthly log growth, %)
  tab:desc-sub       price by sub-period: 2005-2020, 2021-2022, 2023-2025, 2015-2025 (out of sample)
  tab:desc-corr      correlation of each indicator with the monthly average price, by sub-period
  tab:desc-leadlag   lead-lag correlations (yoy log change of the monthly mean price vs the
                     indicator at month m-L; surveys in levels, IPI in yoy log growth)
  weekday / month    mean price by weekday and by month relative to the overall mean,
                     Germany 2012-2025 (real weekends), Italy 2005-2025
"""
import os, warnings
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.tsa.stattools import adfuller, kpss, acf
warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "midas_codes", "data_v4")
OUT = os.path.join(HERE, "desc_v4"); os.makedirs(OUT, exist_ok=True)
IPI = ("Consumer_Goods", "Manufacturing", "Electricity_Gas_Steam")
DROP_M = ("elec_price", "gas_price", "brent_price", "ciss", "istat_bci_it")   # istat_bci_it replaced by iesi_it
LAB = {"pmi_mfg_de": "PMI manufacturing", "pmi_serv_de": "PMI services", "ifo_bci_de": "ifo",
       "pmi_mfg_it": "PMI manufacturing", "pmi_serv_it": "PMI services", "iesi_it": "IESI",
       "Consumer_Goods": "IPI consumer goods", "Manufacturing": "IPI manufacturing",
       "Electricity_Gas_Steam": "IPI energy"}

def summ(x):
    x = pd.Series(x).dropna()
    return dict(N=len(x), Mean=x.mean(), SD=x.std(), Min=x.min(), Q05=x.quantile(.05),
                Median=x.median(), Q95=x.quantile(.95), Max=x.max(),
                Skew=stats.skew(x), ExKurt=stats.kurtosis(x),
                ACF1=acf(x, nlags=7)[1], ACF7=acf(x, nlags=7)[7],
                ADF_p=adfuller(x, autolag="AIC")[1],           # H0: unit root
                KPSS_p=kpss(x, regression="c", nlags="auto")[1])  # H0: level stationary (p capped at .01/.10)

def load(c):
    d = pd.read_csv(f"{D}/prices_{c}_daily.csv", parse_dates=["date"]).set_index("date")
    m = pd.read_csv(f"{D}/dataset_{c}_m.csv", parse_dates=["date"]).set_index("date")
    m = m[[v for v in m.columns if v not in DROP_M]]
    return d, m

def show(title, df, name, **kw):
    print(f"\n### {title}"); print(df.round(3).to_string(**kw))
    df.to_csv(os.path.join(OUT, name + ".csv"), **({} if kw.get("index", True) else dict(index=False)))

# ---- tab:desc-daily / tab:desc-monthly -------------------------------------------------------
rows, mrows = [], []
for c in ["de", "it"]:
    d, m = load(c)
    p = d.elec_price
    a, b = p.median(), stats.median_abs_deviation(p, scale="normal")
    series = [("price", p), ("asinh(price)", np.arcsinh(p)), ("asinh((price-med)/MAD)", np.arcsinh((p - a) / b))]
    if (p > 0).all(): series.append(("log(price)", np.log(p)))
    series += [("gas TTF", d.gas_price), ("Brent", d.brent_price)]
    for name, x in series:
        rows.append(dict(country=c, var=name, **summ(x)))
    print(f"{c}: {len(p)} days {p.index.min().date()}..{p.index.max().date()}, negative days = {int((p < 0).sum())}, "
          f"zero days = {int((p == 0).sum())}, real (non-imputed) days = {int((d.elec_imputed == 0).sum())}, "
          f"full-sample median = {a:.2f}, MAD = {b:.2f}")
    for v in m.columns:
        x = 100 * np.log(m[v]).diff() if v in IPI else m[v]
        mrows.append(dict(country=c, var=v, **summ(x)))
show("tab:desc-daily (moments of the daily series)", pd.DataFrame(rows), "desc_daily", index=False)
show("tab:desc-monthly (surveys in levels, IPI in monthly log growth, %)", pd.DataFrame(mrows), "desc_monthly", index=False)

# ---- tab:desc-sub: sub-periods aligned with the results section ---------------------------
per = [("2005", "2020"), ("2021", "2022"), ("2023", "2025"), ("2015", "2025")]
srows = []
for c in ["de", "it"]:
    d, _ = load(c)
    for a, b in per:
        x, g = d.loc[a:b, "elec_price"], d.loc[a:b, "gas_price"]
        srows.append(dict(country=c, period=f"{a}-{b}" + (" (out of sample)" if a == "2015" else ""),
                          N=len(x), mean=x.mean(), sd=x.std(), cv=x.std() / x.mean(), min=x.min(), max=x.max(),
                          neg=int((x < 0).sum()), corr_gas=x.corr(g)))
S = pd.DataFrame(srows); print(); print("### tab:desc-sub (check against corrections_v3.tex, block 4a)")
print(S.round(2).to_string(index=False)); S.to_csv(os.path.join(OUT, "desc_sub.csv"), index=False)

# ---- tab:desc-corr: indicator-price correlation by sub-period (block 4b) ------------------
PER = [("2005-01-01", "2020-12-31"), ("2015-01-01", "2020-12-31"), ("2021-01-01", "2022-12-31"), ("2023-01-01", "2025-12-31")]
crows = []
for c in ["de", "it"]:
    d, m = load(c)
    inds = list(m.columns)
    pm = d.elec_price.resample("ME").mean()
    m = m.copy(); m["price"] = pm.reindex(m.index)
    for v in IPI: m[v] = np.log(m[v]).diff()            # IPI as monthly log growth, surveys in levels
    corr = {v: [m.loc[a:b, [v, "price"]].dropna().corr().iloc[0, 1] for a, b in PER] for v in inds}
    for v in inds:
        crows.append(dict(country=c, indicator=LAB.get(v, v), var=v,
                          **{f"{a[:4]}-{b[:4]}": r for (a, b), r in zip(PER, corr[v])}))
C = pd.DataFrame(crows); print(); print("### tab:desc-corr (check against corrections_v3.tex, block 4b)")
print(C.round(2).to_string(index=False)); C.to_csv(os.path.join(OUT, "desc_corr.csv"), index=False)

# ---- tab:desc-leadlag --------------------------------------------------------------------
# yoy log change of the monthly mean price at month m against each indicator at month m-L
# (surveys in levels, IPI in yoy log growth)
for c in ["de", "it"]:
    d, m = load(c)
    pm = d.elec_price.resample("ME").mean(); pm.index = pm.index.to_period("M")
    m = m.copy(); m.index = m.index.to_period("M")
    lp = np.log(pm.clip(lower=1)); yoy = 100 * (lp - lp.shift(12))
    res = {}
    for v in m.columns:
        z = 100 * (np.log(m[v]) - np.log(m[v].shift(12))) if v in IPI else m[v]
        res[v] = [yoy.corr(z.shift(L)) for L in (0, 1, 2, 3, 6)]
    LL = pd.DataFrame(res, index=["L0", "L1", "L2", "L3", "L6"]).T
    LL.insert(0, "indicator", [LAB.get(v, v) for v in LL.index])
    show(f"tab:desc-leadlag {c} (price yoy vs indicator at m-L)", LL, f"desc_leadlag_{c}")

# ---- weekday and month profile (real weekends: DE from 2012, IT from 2005) ---------------
for c, start in [("de", "2012-01-09"), ("it", "2005-01-01")]:
    d, _ = load(c)
    p = d.loc[start:, "elec_price"]
    base = p.mean()
    wd = p.groupby(p.index.dayofweek).mean(); wd.index = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
    mo = p.groupby(p.index.month).mean()
    prof = pd.DataFrame({"weekday": wd.index.tolist() + [""] * 5, "mean_wd": wd.tolist() + [np.nan] * 5,
                         "rel_wd_pct": (100 * (wd / base - 1)).tolist() + [np.nan] * 5,
                         "month": list(range(1, 13)), "mean_month": mo.tolist(),
                         "rel_month_pct": (100 * (mo / base - 1)).tolist()})
    print(f"\n### weekday / month profile {c} from {start} (overall mean {base:.2f} EUR/MWh; rel = % deviation)")
    print(prof.round(2).to_string(index=False)); prof.to_csv(os.path.join(OUT, f"desc_profile_{c}.csv"), index=False)

# ---- data check: share of days lying exactly on the straight line between their neighbours ----
# (should be ~0 in data_v4: interpolated weekends were the defect of the old files; a day
#  carried forward over a weekend is NOT linear with its neighbours unless both equal)
for c in ["de", "it"]:
    d, _ = load(c)
    p = d.elec_price
    lin = ((p.shift(-1) - 2 * p + p.shift(1)).abs() < 0.015)
    print(c, "share of linearly interpolated days by weekday (Mon=0):",
          lin.groupby(d.index.dayofweek).mean().round(3).tolist(),
          "| imputed (carried-forward) days:", int(d.elec_imputed.sum()))
