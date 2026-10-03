# -*- coding: utf-8 -*-
"""
Baseline characteristics and event counts by cohort.

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
"""四库基线特征表 + 事件数可行性评估"""
import os
import pandas as pd
import numpy as np

IN = r"DATA_ROOT\分析数据\四库合并长表.parquet"
OUT = r"DATA_ROOT\分析数据"
long = pd.read_parquet(IN)
print(f"读入长表 {len(long):,} 行")

# ---- 编码清洗 ----
BIN = ["stroke", "diab", "hibp", "hearte", "cancre", "lunge", "smokev", "hosp"]
for c in BIN:
    if c in long.columns:
        long[c] = long[c].where(long[c].isin([0, 1]))          # 负值/特殊码 -> NaN
if "gender" in long.columns:
    long["female"] = (long["gender"] == 2).astype("float")     # Gateway: 1=男 2=女
    long.loc[~long["gender"].isin([1, 2]), "female"] = float("nan")
for c in ["bmi", "adlfive", "orient", "cogtot", "cesd", "cesd_z", "agey", "iwy"]:
    if c in long.columns:
        long.loc[long[c] < -10, c] = float("nan")
long["depressed"] = long["dep_bin"]                            # 统一用二分化口径

# 只保留已访谈（iwstat==1）
d = long[long["iwstat"] == 1].copy()
print(f"已访谈 person-wave：{len(d):,}")

# 基线 = 每人首个已访谈波
d = d.sort_values(["cohort", "pid", "wave"])
REN = ["agey", "stroke", "diab", "hibp", "hearte", "cancre", "lunge",
       "depressed", "cesd", "smokev", "shlt", "bmi", "adlfive", "orient", "cogtot", "iwy", "wave"]
base = d.groupby(["cohort", "pid"], as_index=False).first()
base = base.rename(columns={c: f"{c}_b" for c in REN if c in base.columns})

# 随访时长 = 末次已访谈年份 - 基线年份
last = d.groupby(["cohort", "pid"], as_index=False).last()[["cohort", "pid", "iwy", "wave"]]
last = last.rename(columns={"iwy": "iwy_last", "wave": "wave_last"})
base = base.merge(last, on=["cohort", "pid"], how="left")
base["followup_y"] = base["iwy_last"] - base["iwy_b"]

# 事件：基线无卒中者中，后续任一波出现卒中
stroke_ev = (d[d["stroke"] == 1].groupby(["cohort", "pid"]).size()
             .rename("n_stroke").reset_index())
base = base.merge(stroke_ev, on=["cohort", "pid"], how="left")
base["n_stroke"] = base["n_stroke"].fillna(0)
base["eligible"] = base["stroke_b"] == 0
base["incident"] = np.where(base["eligible"] & (base["n_stroke"] > 0), 1,
                            np.where(base["eligible"], 0, np.nan))

rows = []
for coh, g in base.groupby("cohort"):
    gw = d[d["cohort"] == coh]
    elig = g[g["eligible"]]
    rows.append({
        "队列": coh,
        "基线年份(中位)": int(g["iwy_b"].median()) if g["iwy_b"].notna().any() else None,
        "人数": len(g),
        "person-wave 数": len(gw),
        "基线年龄 均值(SD)": f"{g['agey_b'].mean():.1f} ({g['agey_b'].std():.1f})",
        "女性 %": f"{g['female'].mean()*100:.1f}",
        "随访年数 中位(IQR)": f"{g['followup_y'].median():.0f} "
                              f"({g['followup_y'].quantile(.25):.0f}-{g['followup_y'].quantile(.75):.0f})",
        "卒中 基线患病 %": f"{g['stroke_b'].mean()*100:.1f}",
        "糖尿病 %": f"{g['diab_b'].mean()*100:.1f}",
        "高血压 %": f"{g['hibp_b'].mean()*100:.1f}",
        "心脏病 %": f"{g['hearte_b'].mean()*100:.1f}",
        "抑郁(单条目) %": f"{g['depressed_b'].mean()*100:.1f}",
        "当前吸烟 %": f"{g['smokev_b'].mean()*100:.1f}",
        "BMI 均值": f"{g['bmi_b'].mean():.1f}",
        "ADL 受限 %": f"{(g['adlfive_b'] > 0).mean()*100:.1f}",
        "纳入分析人数": len(elig),
        "卒中事件数": int(elig["incident"].sum()),
        "事件率 %": f"{elig['incident'].mean()*100:.1f}",
        "人均随访年": f"{elig['followup_y'].mean():.1f}",
        "总随访人年": f"{elig['followup_y'].sum():,.0f}",
    })

tab = pd.DataFrame(rows).set_index("队列").T
pd.set_option("display.width", 220)
print("\n" + "=" * 100)
print("四库基线特征与卒中事件（可行性核心表）")
print("=" * 100)
print(tab.to_string())
tab.to_csv(os.path.join(OUT, "基线特征表.csv"), encoding="utf-8-sig")

print("\n" + "=" * 100)
print("已访谈者中各变量缺失率（%）")
print("=" * 100)
cols = ["agey", "female", "bmi", "stroke", "diab", "hibp", "hearte", "depressed",
        "cesd", "adlfive", "smokev", "educl", "edyrs", "orient", "cogtot", "hosp", "iwy"]
cols = [c for c in cols if c in d.columns]
miss = d.groupby("cohort")[cols].apply(lambda x: (x.isna().mean() * 100).round(1))
print(miss.to_string())
miss.to_csv(os.path.join(OUT, "基线缺失率.csv"), encoding="utf-8-sig")

print("\n" + "=" * 100)
print("抑郁口径对比（基线）")
print("=" * 100)
print(base.groupby("cohort").agg(
    单条目阳性率=("depressed_b", lambda s: f"{s.mean()*100:.1f}%"),
    单条目缺失=("depressed_b", lambda s: f"{s.isna().mean()*100:.1f}%"),
    量表分均值=("cesd_b", lambda s: f"{s.mean():.1f}"),
    量表分缺失=("cesd_b", lambda s: f"{s.isna().mean()*100:.1f}%"),
).to_string())
