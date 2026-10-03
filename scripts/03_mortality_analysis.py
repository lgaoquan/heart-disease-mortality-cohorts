# -*- coding: utf-8 -*-
"""
Main analysis: heart disease and all-cause, cardiovascular and non-cardiovascular mortality.

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
"""
S1：心血管死亡结局构造 + 主分析（既往心脏病 -> 心血管死亡 / 非心血管死亡 / 全因死亡）
死因来源：四库 Harmonized EOL 的 ragcod（1=癌症 / 2=心血管 / 3=其他）
输出：DATA_ROOT\\分析数据\\心血管死亡分析.csv
"""
import os, sys, warnings
import numpy as np
import pandas as pd
from lifelines import CoxTimeVaryingFitter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cohort_intervals import prepare, BASEV

warnings.filterwarnings("ignore")
OUT = r"DATA_ROOT\分析数据"
IN = os.path.join(OUT, "四库合并长表.parquet")
WINDOW = 5

EOL = {
    "HRS":    (r"DATA_ROOT\整理后数据库\HRS\数据\Harmonized_HRS_EOL\GH_HRS_EOL_b.dta", "hhidpn"),
    "ELSA":   (r"DATA_ROOT\整理后数据库\ELSA\数据\Wave_0-11\h_elsa_eol_a2.dta", "idauniq"),
    "SHARE":  (r"DATA_ROOT\整理后数据库\SHARE\数据\GH_SHARE_EOL_g.dta", "mergeid"),
    "CHARLS": (r"DATA_ROOT\整理后数据库\CHARLS\数据\Harmonized_CHARLS\H_CHARLS_EOL_a.dta", "ID"),
}

# ---- 读 EOL 死因 ----
causes = []
for coh, (p, idv) in EOL.items():
    import pyreadstat
    _, m = pyreadstat.read_dta(p, metadataonly=True, output_format="dict")
    cols = [c for c in m.column_names if c.lower() in (idv.lower(), "ragcod", "racod_h", "racod_e", "racod_s", "racod_c", "raxyear")]
    df = pyreadstat.read_dta(p, usecols=cols, output_format="pandas")[0]
    ren = {}
    for c in df.columns:
        lc = c.lower()
        if lc == idv.lower():
            ren[c] = "pid"
        elif lc == "ragcod":
            ren[c] = "ragcod"
        elif lc == "raxyear":
            ren[c] = "eol_death_year"
        else:
            ren[c] = "racod_specific"
    df = df.rename(columns=ren)
    df = df[["pid", "ragcod", "eol_death_year"] + (["racod_specific"] if "racod_specific" in df.columns else [])]
    df.insert(0, "cohort", coh)
    df["pid"] = df["pid"].astype(str)
    causes.append(df)
    n_cv = (df["ragcod"] == 2).sum()
    print(f"{coh:<8} EOL {len(df):>6} 条   其中心血管死因 {n_cv:>5} 条")
cause = pd.concat(causes, ignore_index=True)

# ---- 长表 + 区间 ----
long = prepare(pd.read_parquet(IN))
long["pid"] = long["pid"].astype(str)
long = long.merge(cause, on=["cohort", "pid"], how="left")
print(f"\n长表 {len(long):,} 行；已关联死因 {long['ragcod'].notna().sum():,} 行")

# 复刻 intervals 逻辑（额外携带 ragcod）
d = long.sort_values(["cohort", "pid", "wave"]).reset_index(drop=True)
g = d.groupby(["cohort", "pid"], sort=False)
d["iwy_f"] = g["iwy"].ffill()
d["w_at_iyw"] = d["wave"].where(d["iwy"].notna())
d["w_at_iyw"] = d.groupby(["cohort", "pid"], sort=False)["w_at_iyw"].ffill()
d["year_est"] = d["iwy"].fillna(d["iwy_f"] + 2.0 * (d["wave"] - d["w_at_iyw"]))
g = d.groupby(["cohort", "pid"], sort=False)
d["is_base"] = (d["alive_iw"] == 1) & (g["alive_iw"].cumsum() == 1)
b = d[d["is_base"]][["cohort", "pid", "wave", "year_est", "agey", "female", "educl"] + BASEV].copy()
ren = {"wave": "base_wave", "year_est": "base_year", "agey": "age0",
       "female": "female_b", "educl": "educl_b"}
for c in BASEV:
    ren[c] = f"{c}_b"
b = b.rename(columns=ren)
# 清洗治疗类变量（0/1），再前向填充
if "rxheart" in d.columns:
    d["rxheart"] = d["rxheart"].where(d["rxheart"].isin([0, 1]))
g = d.groupby(["cohort", "pid"], sort=False)
for c in ["hearte", "diab", "hibp", "smokev", "cesd_z", "bmi_cat", "adlfive", "rxheart"]:
    if c in d.columns:
        d[f"x_{c}"] = g[c].ffill()
d["n_year"] = g["year_est"].shift(-1)
d["n_dead"] = g["dead"].shift(-1)
d["n_cause"] = g["ragcod"].shift(-1)

d = d.merge(b, on=["cohort", "pid"], how="inner")
d["start"] = d["year_est"] - d["base_year"]
d["stop"] = d["n_year"] - d["base_year"]
d["age"] = d["age0"] + d["start"]
m = (d["wave"] >= d["base_wave"]) & d["start"].notna() & d["stop"].notna()
m &= (d["stop"] > d["start"]) & (d["start"] >= 0) & (d["start"] < WINDOW) & d["n_dead"].notna()
pp = d[m].copy()
pp["stop"] = np.minimum(pp["stop"], WINDOW)
pp = pp[pp["is_base"] == False]
pp["event_all"] = pp["n_dead"].astype(int)
pp["event_cv"] = ((pp["n_dead"] == 1) & (pp["n_cause"] == 2)).astype(int)
pp["event_non"] = ((pp["n_dead"] == 1) & (pp["n_cause"].isin([1, 3]))).astype(int)
for c in BASEV:
    pp[f"{c}_b"] = pp[f"{c}_b"].fillna(-1)
pp["educl_b"] = pp["educl_b"].fillna(0)
pp["female_b"] = pp["female_b"].fillna(-1)
pp["pid"] = pp["pid"].astype(str)

print(f"person-interval {len(pp):,} 行")
print(pp.groupby("cohort").agg(interval=("event_all", "size"), 全因死亡=("event_all", "sum"),
                               心血管死亡=("event_cv", "sum"), 非心血管死亡=("event_non", "sum"),
                               死因缺失=("event_all", lambda s: 0)).to_string())

# ---- 主分析 ----
def meta(r):
    w = 1 / r["se"] ** 2
    muf = (w * r["logHR"]).sum() / w.sum()
    Q = float((w * (r["logHR"] - muf) ** 2).sum()); k = len(r)
    C = w.sum() - (w ** 2).sum() / w.sum()
    t2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = 1 / (r["se"] ** 2 + t2)
    mm = float((wr * r["logHR"]).sum() / wr.sum()); se = float(np.sqrt(1 / wr.sum()))
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return mm, se, I2, t2

rows = []
EXPO = "x_hearte"
for out, label in [("event_cv", "心血管死亡"), ("event_non", "非心血管死亡"), ("event_all", "全因死亡")]:
    cov = ["age", "female_b", "educl_b"] + [f"{c}_b" for c in BASEV if c != "hearte"]
    res = []
    for coh in ["CHARLS", "ELSA", "HRS", "SHARE"]:
        sub = pp[pp["cohort"] == coh]
        dd = sub[["event", "start", "stop", EXPO, "pid"] + cov].copy() if False else \
             sub.assign(event=sub[out])[["event", "start", "stop", EXPO, "pid"] + cov]
        dd = dd.dropna(subset=["event", "start", "stop", EXPO, "age", "female_b"])
        if len(dd) < 500 or dd["event"].sum() < 30:
            continue
        dd = dd.copy()
        dd["pid2"] = dd["pid"].astype("category").cat.codes
        try:
            ctv = CoxTimeVaryingFitter(penalizer=0.01)
            ctv.fit(dd[["pid2", "start", "stop", "event", EXPO] + cov], id_col="pid2",
                    start_col="start", stop_col="stop", event_col="event")
            hr = float(np.exp(ctv.params_[EXPO]))
            ci = np.exp(ctv.confidence_intervals_.loc[EXPO]).values
            res.append({"coh": coh, "hr": hr, "logHR": np.log(hr),
                        "se": (np.log(ci[1]) - np.log(ci[0])) / 3.92,
                        "ev": int(dd["event"].sum())})
        except Exception as e:
            print(f"   {out} {coh} 失败 {e}")
    if len(res) < 3:
        continue
    r = pd.DataFrame(res)
    mm, se, I2, t2 = meta(r)
    loo = []
    for drop in r["coh"]:
        d2 = r[r["coh"] != drop]
        if len(d2) >= 3:
            _, _, i22, _ = meta(d2); loo.append(i22)
    rows.append({"暴露": "既往心脏病", "结局": label, "总事件": int(r["ev"].sum()),
                 "合并HR": round(np.exp(mm), 3), "CI": f"{np.exp(mm-1.96*se):.2f}-{np.exp(mm+1.96*se):.2f}",
                 "I2_%": round(I2, 1), "tau2": round(t2, 4),
                 "留一法I2最低": round(min(loo), 1) if loo else np.nan,
                 "HR范围": f"{r['hr'].min():.2f}-{r['hr'].max():.2f}",
                 "各队列HR": " ".join(f"{x.coh}:{x.hr:.2f}({x.ev})" for x in r.itertuples())})

tab = pd.DataFrame(rows)
pd.set_option("display.width", 260)
print("\n" + "=" * 130)
print("既往心脏病 → 三类死亡结局（时变 Cox，统一 5 年窗口）")
print("=" * 130)
print(tab.to_string(index=False))
tab.to_csv(os.path.join(OUT, "心血管死亡分析.csv"), index=False, encoding="utf-8-sig")
print(f"\n已写出：{os.path.join(OUT, '心血管死亡分析.csv')}")
pp.to_parquet(os.path.join(OUT, "person_interval_死亡.parquet"), index=False)
print(f"person-interval 已存：{os.path.join(OUT, 'person_interval_死亡.parquet')}")
