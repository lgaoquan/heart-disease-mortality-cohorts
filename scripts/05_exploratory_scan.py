# -*- coding: utf-8 -*-
"""
Exploratory scan of all exposure-outcome combinations.

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
"""统一扫描：对每个暴露 x 结局组合做队列内时变 Cox + 随机效应 meta + 留一法 I²"""
import os, sys, warnings, itertools
import numpy as np
import pandas as pd
from lifelines import CoxTimeVaryingFitter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_intervals import prepare, intervals, EXPO, BASEV

warnings.filterwarnings("ignore")
IN = r"DATA_ROOT\分析数据\四库合并长表.parquet"
OUT = r"DATA_ROOT\分析数据"
WINDOW = 5
OUTCOMES = ["stroke", "hearte", "cancre", "lunge", "dead"]
LABEL = {"stroke": "自报卒中", "hearte": "自报心脏病", "cancre": "自报癌症",
         "lunge": "自报肺病", "dead": "全因死亡"}
OBJ = {"dead"}                      # 客观结局


def meta(r):
    w = 1 / r["se"] ** 2
    muf = (w * r["logHR"]).sum() / w.sum()
    Q = float((w * (r["logHR"] - muf) ** 2).sum()); k = len(r)
    C = w.sum() - (w ** 2).sum() / w.sum()
    t2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = 1 / (r["se"] ** 2 + t2)
    m = float((wr * r["logHR"]).sum() / wr.sum()); se = float(np.sqrt(1 / wr.sum()))
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return m, se, I2, t2


long = prepare(pd.read_parquet(IN))
print(f"长表 {len(long):,} 行；构造各结局的 person-interval …")
skel = {}
for out in OUTCOMES:
    pp = intervals(long, out, WINDOW)
    skel[out] = pp
    print(f"  {LABEL[out]:<8}{len(pp):>9,} 行   事件 {int(pp['event'].sum()):>7,}")

rows = []
for expo, out in itertools.product(EXPO, OUTCOMES):
    if expo == out or (expo, out) in (("hibp", "hearte"), ("hearte", "hearte")):
        continue
    pp = skel[out]
    xcol = f"x_{expo}"
    if xcol not in pp.columns:
        continue
    cov = ["age", "female_b", "educl_b"] + [f"{c}_b" for c in BASEV if c != expo]
    res = []
    for coh in ["CHARLS", "ELSA", "HRS", "SHARE"]:
        d = pp[(pp["cohort"] == coh)][["event", "start", "stop", xcol, "pid"] + cov].dropna(
            subset=["event", "start", "stop", xcol, "age", "female_b"])
        if len(d) < 500 or d["event"].sum() < 30:
            continue
        d = d.copy()
        d["pid2"] = d["pid"].astype("category").cat.codes
        try:
            ctv = CoxTimeVaryingFitter(penalizer=0.01)
            ctv.fit(d[["pid2", "start", "stop", "event", xcol] + cov], id_col="pid2",
                    start_col="start", stop_col="stop", event_col="event")
            hr = float(np.exp(ctv.params_[xcol]))
            ci = np.exp(ctv.confidence_intervals_.loc[xcol]).values
            res.append({"coh": coh, "hr": hr, "logHR": np.log(hr),
                        "se": (np.log(ci[1]) - np.log(ci[0])) / 3.92,
                        "ev": int(d["event"].sum())})
        except Exception:
            continue
    if len(res) < 3:
        continue
    r = pd.DataFrame(res)
    m, se, I2, t2 = meta(r)
    loo = []
    for drop in r["coh"]:
        dd = r[r["coh"] != drop]
        if len(dd) >= 3:
            _, _, i22, _ = meta(dd); loo.append(i22)
    rows.append({"暴露": expo, "结局": LABEL[out], "结局类型": "客观" if out in OBJ else "自报",
                 "队列数": len(r), "总事件": int(r["ev"].sum()),
                 "合并HR": round(np.exp(m), 3),
                 "CI": f"{np.exp(m-1.96*se):.2f}-{np.exp(m+1.96*se):.2f}",
                 "I2_%": round(I2, 1), "tau2": round(t2, 4),
                 "留一法I2最低": round(min(loo), 1) if loo else np.nan,
                 "HR范围": f"{r['hr'].min():.2f}-{r['hr'].max():.2f}",
                 "各队列HR": " ".join(f"{x.coh}:{x.hr:.2f}({x.ev})" for x in r.itertuples())})

tab = pd.DataFrame(rows).sort_values("I2_%", ascending=False)
pd.set_option("display.width", 300)
print("\n" + "=" * 140)
print("统一扫描结果（按 I² 降序）")
print("=" * 140)
print(tab.to_string(index=False))
tab.to_csv(os.path.join(OUT, "统一扫描结果.csv"), index=False, encoding="utf-8-sig")
print(f"\n已写出：{os.path.join(OUT, '统一扫描结果.csv')}")
print("\n判读标准：")
print("  · 留一法I2最低 > 60%  → 异质性稳健（剔除任一队列都不塌）")
print("  · HR范围跨度 > 1.5倍  → 临床上有意义的差异")
print("  · 两者同时满足且结局为客观 → 可支撑顶刊异质性叙事")
