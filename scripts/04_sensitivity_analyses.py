# -*- coding: utf-8 -*-
"""
Absolute risk (Aalen-Johansen), follow-up window, cause-coverage, treatment and sex analyses.

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
S2：敏感性分析与 EJPC 硬性要求项
  A. 绝对风险（Aalen–Johansen 累积发生率，竞争风险校正）
  B. 随访窗口敏感性（3 / 5 年）
  C. 死因覆盖率敏感性（仅 HRS+SHARE）
  D. 治疗修饰（按 rxheart 分层）
  E. 性别分层（SAGER 要求）
输出：DATA_ROOT\\分析数据\\S2_*.csv
"""
import os, sys, warnings
import numpy as np
import pandas as pd
from lifelines import CoxTimeVaryingFitter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

OUT = r"DATA_ROOT\分析数据"
pp = pd.read_parquet(os.path.join(OUT, "person_interval_死亡.parquet"))
COV = ["age", "female_b", "educl_b", "smokev_b", "hibp_b", "bmi_cat_b", "adlfive_b"]
print(f"读入 person-interval {len(pp):,} 行")


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


def fit_cohort(sub, out, expo="x_hearte", cov=None):
    cov = COV if cov is None else cov
    d = sub.assign(event=sub[out])[["event", "start", "stop", expo, "pid"] + cov]
    d = d.dropna(subset=["event", "start", "stop", expo, "age"])
    if len(d) < 300 or d["event"].sum() < 20:
        return None
    d = d.copy(); d["pid2"] = d["pid"].astype("category").cat.codes
    try:
        ctv = CoxTimeVaryingFitter(penalizer=0.01)
        ctv.fit(d[["pid2", "start", "stop", "event", expo] + cov], id_col="pid2",
                start_col="start", stop_col="stop", event_col="event")
        hr = float(np.exp(ctv.params_[expo]))
        ci = np.exp(ctv.confidence_intervals_.loc[expo]).values
        return {"hr": hr, "logHR": np.log(hr),
                "se": (np.log(ci[1]) - np.log(ci[0])) / 3.92, "ev": int(d["event"].sum())}
    except Exception:
        return None


def run(out, cohorts=("CHARLS", "ELSA", "HRS", "SHARE"), sub_filter=None,
        expo="x_hearte", cov=None):
    res = []
    for coh in cohorts:
        sub = pp[pp["cohort"] == coh]
        if sub_filter is not None:
            sub = sub_filter(sub)
        if len(sub) == 0:
            continue
        r = fit_cohort(sub, out, expo, cov)
        if r:
            r["coh"] = coh
            res.append(r)
    if len(res) < 2:
        return None
    r = pd.DataFrame(res)
    mm, se, I2, t2 = meta(r)
    return {"合并HR": round(np.exp(mm), 3),
            "CI": f"{np.exp(mm-1.96*se):.2f}-{np.exp(mm+1.96*se):.2f}",
            "I2_%": round(I2, 1), "队列数": len(r), "总事件": int(r["ev"].sum()),
            "各队列HR": " ".join(f"{x.coh}:{x.hr:.2f}({x.ev})" for x in r.itertuples())}


# ================= A. 绝对风险（Aalen–Johansen CIF）=================
def aj_cif(sub, event_col, comp_col, tmax=5.0, step=0.5):
    """Aalen–Johansen 累积发生率：返回 (时间点, CIF_event, CIF_competing)"""
    times = np.arange(step, tmax + 1e-9, step)
    s = sub.dropna(subset=["start", "stop", event_col, comp_col])
    cif_e, cif_c = [], []
    cum_e = cum_c = 0.0
    S = 1.0
    for t in times:
        at_risk = (s["stop"] >= t - 1e-9).sum()
        d_e = ((s["stop"] <= t) & (s["stop"] > t - step - 1e-9) & (s[event_col] == 1)).sum()
        d_c = ((s["stop"] <= t) & (s["stop"] > t - step - 1e-9) & (s[comp_col] == 1)).sum()
        if at_risk > 0:
            cum_e += S * d_e / at_risk
            cum_c += S * d_c / at_risk
            S *= (1 - (d_e + d_c) / at_risk)
        cif_e.append(cum_e); cif_c.append(cum_c)
    return times, np.array(cif_e), np.array(cif_c)


rows_abs = []
for coh in ["CHARLS", "ELSA", "HRS", "SHARE"]:
    s = pp[pp["cohort"] == coh]
    for lvl, lab in [(0, "无心脏病"), (1, "有心脏病")]:
        g = s[s["x_hearte"] == lvl]
        if len(g) < 100:
            continue
        t, cif_cv, cif_other = aj_cif(g, "event_cv", "event_non")
        _, cif_all, _ = aj_cif(g, "event_all", "event_all")
        n = g["pid"].nunique()
        rows_abs.append({"队列": coh, "暴露组": lab, "人数": n,
                         "5年心血管死亡累积发生率%": round(cif_cv[-1] * 100, 2),
                         "5年非心血管死亡%": round(cif_other[-1] * 100, 2),
                         "5年全因死亡%": round((cif_cv[-1] + cif_other[-1]) * 100, 2)})

abs_tab = pd.DataFrame(rows_abs)
print("\n" + "=" * 110)
print("A. 绝对风险（Aalen–Johansen 累积发生率，5 年）")
print("=" * 110)
print(abs_tab.to_string(index=False))

# 绝对风险差
ard = []
for coh in abs_tab["队列"].unique():
    g = abs_tab[abs_tab["队列"] == coh]
    if len(g) < 2:
        continue
    hi = g[g["暴露组"] == "有心脏病"].iloc[0]
    lo = g[g["暴露组"] == "无心脏病"].iloc[0]
    ard.append({"队列": coh,
                "心血管死亡 绝对风险差%": round(hi["5年心血管死亡累积发生率%"] - lo["5年心血管死亡累积发生率%"], 2),
                "全因死亡 绝对风险差%": round(hi["5年全因死亡%"] - lo["5年全因死亡%"], 2)})
ard_tab = pd.DataFrame(ard)
print("\n绝对风险差（有心脏病 − 无心脏病，5 年）：")
print(ard_tab.to_string(index=False))
abs_tab.to_csv(os.path.join(OUT, "S2_绝对风险.csv"), index=False, encoding="utf-8-sig")
ard_tab.to_csv(os.path.join(OUT, "S2_绝对风险差.csv"), index=False, encoding="utf-8-sig")

# ================= B/C =================
# 3 年窗口：重新截断 stop（不能按 stop<=3 过滤——波间隔 2 年，过滤会只剩第一批区间）
pp3 = pp[pp["start"] < 3.0].copy()
pp3["stop"] = np.minimum(pp3["stop"], 3.0)
pp3 = pp3[pp3["stop"] > pp3["start"]]
print(f"\n3 年窗口 person-interval {len(pp3):,} 行（原 {len(pp):,}）")
_pp_backup = pp
pp = pp3
rows = []
for out, olab in [("event_all", "全因死亡"), ("event_cv", "心血管死亡"), ("event_non", "非心血管死亡")]:
    r = run(out)
    if r:
        rows.append({"分析": "窗口=3年（重新截断）", "结局": olab, **r})
pp = _pp_backup

for out, olab in [("event_all", "全因死亡"), ("event_cv", "心血管死亡"), ("event_non", "非心血管死亡")]:
    r = run(out)
    if r:
        rows.append({"分析": "主分析（5年窗口）", "结局": olab, **r})
    # 仅 HRS+SHARE（死因覆盖敏感性）：两库即可
    r2 = run(out, cohorts=("HRS", "SHARE"))
    if r2:
        rows.append({"分析": "仅 HRS+SHARE", "结局": olab, **r2})

tab = pd.DataFrame(rows)
print("\n" + "=" * 140)
print("B/C. 敏感性分析")
print("=" * 140)
print(tab.to_string(index=False))

# ================= D. 治疗修饰（三分类暴露）=================
# 有心脏病者中 x_hearte 恒定，无法分层估计；改为三分类暴露：
#   0 = 无心脏病  |  1 = 有心脏病·未服药  |  2 = 有心脏病·服药
pp = _pp_backup.copy()
pp["x_hd3"] = np.where(pp["x_hearte"] == 0, 0,
                       np.where(pp["x_rxheart"] == 1, 2,
                                np.where(pp["x_rxheart"] == 0, 1, np.nan)))
rows_d = []
for out, olab in [("event_all", "全因死亡"), ("event_cv", "心血管死亡")]:
    res = []
    for coh in ["CHARLS", "HRS", "SHARE"]:          # ELSA 的 rxheart 100% 缺失
        sub = pp[(pp["cohort"] == coh) & pp["x_hd3"].notna()]
        d = sub.assign(e0=(sub["x_hd3"] == 1).astype(float),
                       e1=(sub["x_hd3"] == 2).astype(float),
                       event=sub[out])[["event", "start", "stop", "e0", "e1", "pid"] + COV]
        d = d.dropna(subset=["event", "start", "stop", "age"])
        if len(d) < 300 or d["event"].sum() < 20:
            continue
        d = d.copy(); d["pid2"] = d["pid"].astype("category").cat.codes
        try:
            ctv = CoxTimeVaryingFitter(penalizer=0.01)
            ctv.fit(d[["pid2", "start", "stop", "event", "e0", "e1"] + COV], id_col="pid2",
                    start_col="start", stop_col="stop", event_col="event")
            ci0 = np.exp(ctv.confidence_intervals_.loc["e0"]).values
            ci1 = np.exp(ctv.confidence_intervals_.loc["e1"]).values
            res.append({"队列": coh, "事件": int(d["event"].sum()),
                        "有心脏病未服药HR": round(float(np.exp(ctv.params_["e0"])), 3),
                        "CI0": f"{ci0[0]:.2f}-{ci0[1]:.2f}",
                        "有心脏病服药HR": round(float(np.exp(ctv.params_["e1"])), 3),
                        "CI1": f"{ci1[0]:.2f}-{ci1[1]:.2f}"})
        except Exception as e:
            print(f"   D {coh} {out} 失败 {e}")
    if res:
        t = pd.DataFrame(res)
        t.insert(0, "结局", olab)
        print("\n" + "=" * 120)
        print(f"D. 治疗修饰（三分类暴露，参照组=无心脏病）：{olab}")
        print("=" * 120)
        print(t.to_string(index=False))
        t.to_csv(os.path.join(OUT, f"S2_治疗修饰_{olab}.csv"), index=False, encoding="utf-8-sig")

# ================= E. 性别分层 =================
COV_NOFEM = [c for c in COV if c != "female_b"]
rows_e = []
for lvl, lab in [(0, "男性"), (1, "女性")]:
    for out, olab in [("event_all", "全因死亡"), ("event_cv", "心血管死亡")]:
        r = run(out, sub_filter=lambda s, l=lvl: s[s["female_b"] == l], cov=COV_NOFEM)
        if r:
            rows_e.append({"亚组": lab, "结局": olab, **r})

if rows_e:
    t = pd.DataFrame(rows_e)
    print("\n" + "=" * 130)
    print("E. 性别分层（SAGER 要求）")
    print("=" * 130)
    print(t.to_string(index=False))
    t.to_csv(os.path.join(OUT, "S2_性别分层.csv"), index=False, encoding="utf-8-sig")
else:
    print("\nE. 性别分层：无可用结果")

tab.to_csv(os.path.join(OUT, "S2_敏感性分析.csv"), index=False, encoding="utf-8-sig")
print(f"\n已写出：S2_绝对风险.csv / S2_绝对风险差.csv / S2_敏感性分析.csv / S2_性别分层.csv")
