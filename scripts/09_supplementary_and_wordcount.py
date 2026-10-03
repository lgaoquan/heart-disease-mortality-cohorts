# -*- coding: utf-8 -*-
"""
Supplementary tables and manuscript word counts.

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
"""生成补充材料 + 统计稿件字数"""
import os, re
import pandas as pd

PJ = r"DATA_ROOT"
OUT = os.path.join(PJ, "稿件", "补充材料")
os.makedirs(OUT, exist_ok=True)

# ---------- Table S1：41 组合扫描 ----------
s1 = pd.read_csv(os.path.join(PJ, "分析数据", "统一扫描结果.csv"))
s1 = s1.sort_values("I2_%", ascending=False)
s1.to_csv(os.path.join(OUT, "TableS1_探索性扫描41组合.csv"), index=False, encoding="utf-8-sig")
print(f"Table S1：{len(s1)} 行")

# ---------- Table S2：变量映射 ----------
s2 = pd.read_csv(os.path.join(PJ, "变量映射表.csv"))
s2.to_csv(os.path.join(OUT, "TableS2_四库变量映射.csv"), index=False, encoding="utf-8-sig")
print(f"Table S2：{len(s2)} 行")

# ---------- Table S3：死因覆盖率 ----------
pp = pd.read_parquet(os.path.join(PJ, "分析数据", "person_interval_死亡.parquet"))
rows = []
for c in ["CHARLS", "ELSA", "HRS", "SHARE"]:
    g = pp[pp["cohort"] == c]
    tot = int(g["event_all"].sum())
    cv = int(g["event_cv"].sum())
    non = int(g["event_non"].sum())
    rows.append({"Cohort": c, "All-cause deaths": tot,
                 "Cardiovascular deaths": cv, "Non-cardiovascular deaths": non,
                 "Deaths with cause recorded": cv + non,
                 "Cause coverage, %": round((cv + non) / tot * 100, 1) if tot else None})
s3 = pd.DataFrame(rows)
s3.loc[len(s3)] = {"Cohort": "Total", "All-cause deaths": s3["All-cause deaths"].sum(),
                   "Cardiovascular deaths": s3["Cardiovascular deaths"].sum(),
                   "Non-cardiovascular deaths": s3["Non-cardiovascular deaths"].sum(),
                   "Deaths with cause recorded": s3["Deaths with cause recorded"].sum(),
                   "Cause coverage, %": round(s3["Deaths with cause recorded"].sum() / s3["All-cause deaths"].sum() * 100, 1)}
s3.to_csv(os.path.join(OUT, "TableS3_死因覆盖率.csv"), index=False, encoding="utf-8-sig")
print(f"Table S3：{len(s3)} 行")
print(s3.to_string(index=False))

# ---------- Table S4：缺失率 ----------
long = pd.read_parquet(os.path.join(PJ, "分析数据", "四库合并长表.parquet"))
d = long[long["iwstat"] == 1].copy()
for c in ["stroke", "diab", "hibp", "hearte", "cancre", "lunge", "smokev", "hosp"]:
    if c in d.columns:
        d[c] = d[c].where(d[c].isin([0, 1]))
cols = [c for c in ["agey", "gender", "bmi", "stroke", "diab", "hibp", "hearte", "cancre",
                    "lunge", "smokev", "shlt", "adlfive", "orient", "cogtot", "hosp",
                    "cesd", "educl", "edyrs", "rxheart"] if c in d.columns]
miss = d.groupby("cohort")[cols].apply(lambda x: (x.isna().mean() * 100).round(1)).reset_index()
miss.to_csv(os.path.join(OUT, "TableS4_缺失率.csv"), index=False, encoding="utf-8-sig")
print(f"Table S4：{len(miss)} 行 × {len(cols)} 变量")

# ---------- 字数统计 ----------
ms = open(os.path.join(PJ, "稿件", "Manuscript_EJPC_v1.md"), encoding="utf-8").read()
def wc(txt):
    txt = re.sub(r"\[待填[^\]]*\]", "", txt)
    txt = re.sub(r"^#+ .*$", "", txt, flags=re.M)
    return len(re.findall(r"[A-Za-z][A-Za-z'\-\u2019]*", txt))

ab = ms.split("## Abstract")[1].split("## Lay summary")[0]
lay = ms.split("## Lay summary")[1].split("## 1 Introduction")[0]
body = ms.split("## 1 Introduction")[1].split("## Acknowledgements")[0]
print("\n" + "=" * 60)
print("字数统计（EJPC：正文 ~5,000 词；摘要 ≤250；Lay summary ≤250）")
print("=" * 60)
print(f"  Abstract      : {wc(ab):>5} 词   {'✅' if wc(ab) <= 250 else '❌ 超限'}")
print(f"  Lay summary   : {wc(lay):>5} 词   {'✅' if wc(lay) <= 250 else '❌ 超限'}")
print(f"  Main body     : {wc(body):>5} 词   {'✅' if wc(body) <= 5000 else '❌ 超限'}")
