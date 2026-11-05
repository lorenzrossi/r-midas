"""Descriptive statistics for the Descriptive Statistics section of Chapter 2.

Run from anywhere:  python3 descriptive_stats.py
Needs: pandas, numpy, scipy, statsmodels.
Reads data/prices_{de,it}_daily.csv and data/dataset_{de,it}_m.csv.
"""
import os, warnings
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.tsa.stattools import adfuller, kpss, acf
warnings.filterwarnings("ignore")
pd.set_option("display.width", 250)
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
IPI = ("Consumer_Goods", "Manufacturing", "Electricity_Gas_Steam")

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
    return d, m

rows, mrows = [], []
for c in ["de", "it"]:
    d, m = load(c)
    p = d.elec_price
    a, b = p.median(), stats.median_abs_deviation(p, scale="normal")
    for name, x in [("price", p), ("asinh(price)", np.arcsinh(p)),
                    ("asinh((price-med)/MAD)", np.arcsinh((p - a) / b)),
                    ("gas TTF", d.gas_price), ("Brent", d.brent_price)]:
        rows.append(dict(country=c, var=name, **summ(x)))
    for v in [v for v in m.columns if v not in ("elec_price", "gas_price", "brent_price", "ciss")]:
        x = 100 * np.log(m[v]).diff() if v in IPI else m[v]
        mrows.append(dict(country=c, var=v, **summ(x)))
print(pd.DataFrame(rows).round(3).to_string())
print(pd.DataFrame(mrows).round(3).to_string())

# Sub-periods
per = [("2005", "2019"), ("2020", "2022"), ("2023", "2025"), ("2015", "2025")]
for c in ["de", "it"]:
    d, _ = load(c)
    out = []
    for a, b in per:
        x, g = d.loc[a:b, "elec_price"], d.loc[a:b, "gas_price"]
        out.append(dict(period=f"{a}-{b}", N=len(x), mean=x.mean(), sd=x.std(), cv=x.std() / x.mean(),
                        min=x.min(), max=x.max(), neg=int((x < 0).sum()), corr_gas=x.corr(g)))
    print(c); print(pd.DataFrame(out).round(2).to_string())

# Lead-lag correlations: yoy log change of the monthly mean price at month m
# against each indicator at month m-L (surveys in levels, IPI in yoy log growth)
for c in ["de", "it"]:
    d, m = load(c)
    pm = d.elec_price.resample("ME").mean(); pm.index = pm.index.to_period("M")
    m.index = m.index.to_period("M")
    lp = np.log(pm.clip(lower=1)); yoy = 100 * (lp - lp.shift(12))
    res = {}
    for v in [v for v in m.columns if v not in ("elec_price", "gas_price", "brent_price", "ciss")]:
        z = 100 * (np.log(m[v]) - np.log(m[v].shift(12))) if v in IPI else m[v]
        res[v] = [yoy.corr(z.shift(L)) for L in (0, 1, 2, 3, 6)]
    print(c); print(pd.DataFrame(res, index=["L0", "L1", "L2", "L3", "L6"]).T.round(2))

# Data check: share of days lying exactly on the straight line between their neighbours
for c in ["de", "it"]:
    d, _ = load(c)
    p = d.elec_price
    lin = ((p.shift(-1) - 2 * p + p.shift(1)).abs() < 0.015)
    print(c, "share of linearly interpolated days by weekday (Mon=0):",
          lin.groupby(d.index.dayofweek).mean().round(2).tolist())
