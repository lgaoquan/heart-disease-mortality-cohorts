# -*- coding: utf-8 -*-
"""
Leave-one-out diagnostics for heterogeneity.

Part of the analysis code for:
  Self-reported heart disease and mortality across four ageing cohorts:
  a pooled analysis of 219,364 participants and 34,571 deaths

Data are not distributed with this repository. See README.md for access instructions.
Set the DATA_ROOT environment variable to the folder containing the harmonised
cohort files and the analysis outputs.
"""
import os

DATA_ROOT = os.environ.get("DATA_ROOT", os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "data"))
# -*- coding: utf-8 -*-
"""对扫描出的高异质性组合做留一法诊断：判断异质性是真梯度还是精度假象。"""
import warnings, itertools
import numpy as np
import pandas as pd
from lifelines import CoxTimeVaryingFitter
warnings.filterwarnings("ignore")

IN = r"DATA_ROOT\分析数据\四库合并长表.parquet"
OUT = r"DATA_ROOT\分析数据"
WINDOW = 5
CANDIDATES = [("hearte", "stroke"), ("hibp", "stroke"), ("smokev", "lunge"),
              ("hearte", "lunge"), ("bmi_cat", "lunge"), ("diab", "stroke")]

long = pd.read_parquet(IN)
long = long[long["iwstat"] == 1].copy()
for c in ["stroke", "diab", "hibp", "hearte", "smokev"]:
    long[c] = long[c].where(long[c].isin([0, 1]))
for c in ["bmi", "agey", "iwy", "adlfive"]:
    long.loc[long[c] < -10, c] = np.nan
long["educl"] = long["educl"].where(long["educl"].isin([1, 2, 3]))
long["female"] = np.where(long["gender"].isin([1, 2]), (long["gender"] == 2).astype(float), np.nan)
long["bmi_cat"] = pd.cut(long["bmi"], [0, 18.5, 25, 30, 100], labels=[0, 1, 2, 3]).astype(float)
long["adlfive"] = long["adlfive"].clip(0, 5)

BASEV = ["smokev", "hibp", "hearte", "bmi_cat", "adlfive"]

def build(g, outcome, expo):
    g = g.sort_values(["pid", "wave"])
    first = g.groupby("pid").first()
    ok = first[first[outcome] == 0]
    b = ok[["iwy", "agey", "female", "educl"] + BASEV].copy()
    b.columns = ["y0", "age0", "female", "educl"] + [f"{c}_b" for c in BASEV]
    recs = []
    for pid, sub in g.groupby("pid"):
        if pid not in b.index:
            continue
        sub = sub[sub["iwy"].notna()]
        if len(sub) < 2:
            continue
        yr, st, ex = sub["iwy"].values, sub[outcome].values, sub[expo].values
        y0 = b.loc[pid, "y0"]
        for i in range(len(sub) - 1):
            t0, t1 = yr[i] - y0, yr[i + 1] - y0
            if t1 <= t0 or t0 < 0 or t0 >= WINDOW or np.isnan(st[i + 1]):
                continue
            recs.append({"pid": pid, "start": t0, "stop": min(t1, WINDOW),
                         "event": int(st[i + 1] == 1), "x": ex[i]})
    if not recs:
        return None
    pp = pd.DataFrame(recs).merge(b, left_on="pid", right_index=True, how="left")
    pp["age"] = pp["age0"] + pp["start"]
    return pp

def meta(d):
    w = 1 / d["se"] ** 2
    muf = (w * d["logHR"]).sum() / w.sum()
    Q = float((w * (d["logHR"] - muf) ** 2).sum()); k = len(d)
    C = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = 1 / (d["se"] ** 2 + tau2)
    m = float((wr * d["logHR"]).sum() / wr.sum()); se = float(np.sqrt(1 / wr.sum()))
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return m, se, I2, tau2, Q

for expo, out in CANDIDATES:
    print("=" * 100)
    print(f"暴露 {expo}  ->  结局 {out}")
    print("=" * 100)
    rows = []
    for coh in ["CHARLS", "ELSA", "HRS", "SHARE"]:
        g = long[long["cohort"] == coh]
        pp = build(g, out, expo)
        if pp is None:
            continue
        cov = ["age", "female", "educl"] + [f"{c}_b" for c in BASEV if c != expo]
        d = pp[["event", "start", "stop", "x"] + cov].dropna()
        d = d[d["stop"] > d["start"]].copy()
        if len(d) < 500 or d["event"].sum() < 30:
            continue
        d["pid2"] = pp.loc[d.index, "pid"].astype("category").cat.codes
        try:
            ctv = CoxTimeVaryingFitter(penalizer=0.01)
            ctv.fit(d[["pid2", "start", "stop", "event", "x"] + cov], id_col="pid2",
                    start_col="start", stop_col="stop", event_col="event")
            hr = float(np.exp(ctv.params_["x"]))
            ci = np.exp(ctv.confidence_intervals_.loc["x"]).values
            rows.append({"队列": coh, "事件": int(d["event"].sum()),
                         "HR": round(hr, 3), "CI": f"{ci[0]:.2f}-{ci[1]:.2f}",
                         "logHR": np.log(hr), "se": (np.log(ci[1]) - np.log(ci[0])) / 3.92})
        except Exception as e:
            print(f"   {coh} 拟合失败 {e}")
    r = pd.DataFrame(rows)
    if len(r) < 3:
        print("  可用队列不足\n"); continue
    m, se, I2, tau2, Q = meta(r)
    print(r.to_string(index=False))
    print(f"\n   合并 HR={np.exp(m):.3f} ({np.exp(m-1.96*se):.3f}-{np.exp(m+1.96*se):.3f})  I²={I2:.1f}%  τ²={tau2:.4f}")
    for drop in r["队列"]:
        dd = r[r["队列"] != drop].reset_index(drop=True)
        m2, se2, I22, t2, q2 = meta(dd)
        print(f"     剔除 {drop:<7} -> HR={np.exp(m2):.3f}  I²={I22:.1f}%")
    print()
