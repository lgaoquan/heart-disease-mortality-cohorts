# -*- coding: utf-8 -*-
"""
Generate Figure 1-3 and Table 1-2.

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
S3：图表成品（图内文字用英文，符合 EJPC 投稿要求）
  Figure 1 流程图  Figure 2 森林图  Figure 3 累积发生率曲线
  Table 1 基线特征表（按暴露分组）  Table 2 主分析+敏感性汇总
输出：DATA_ROOT\\图表\\
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from lifelines import CoxTimeVaryingFitter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")

OUT = r"DATA_ROOT\分析数据"
FIG = r"DATA_ROOT\图表"
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 300, "savefig.bbox": "tight"})

COH = ["CHARLS", "ELSA", "HRS", "SHARE"]
COH_FULL = {"CHARLS": "CHARLS (China)", "ELSA": "ELSA (England)",
            "HRS": "HRS (USA)", "SHARE": "SHARE (Europe)"}
COV = ["age", "female_b", "educl_b", "smokev_b", "hibp_b", "bmi_cat_b", "adlfive_b"]

long = pd.read_parquet(os.path.join(OUT, "四库合并长表.parquet"))
long["female"] = np.where(long["gender"].isin([1, 2]), (long["gender"] == 2).astype(float), np.nan)
# 二分类变量清洗：负值缺失码 -> NaN（否则均值会被污染，ELSA 吸烟率曾算出 -40.7%）
for c in ["stroke", "diab", "hibp", "hearte", "cancre", "lunge", "smokev", "hosp"]:
    if c in long.columns:
        long[c] = long[c].where(long[c].isin([0, 1]))
for c in ["bmi", "agey", "iwy", "adlfive", "cesd_z"]:
    if c in long.columns:
        long.loc[long[c] < -10, c] = np.nan
long["educl"] = long["educl"].where(long["educl"].isin([1, 2, 3]))
pp = pd.read_parquet(os.path.join(OUT, "person_interval_死亡.parquet"))
for c in ["x_hearte", "x_diab", "x_hibp", "x_smokev", "x_cesd_z", "x_bmi_cat", "x_adlfive", "x_rxheart"]:
    if c in pp.columns and c != "x_cesd_z":
        pp[c] = pp[c].where(~pp[c].isin([-1])) if c.endswith("_cat") else pp[c]

# ================= Figure 1：流程图 =================
flow = []
for c in COH:
    g = long[long["cohort"] == c]
    n_all = g["pid"].nunique()
    ever = g[g["iwstat"] == 1]["pid"].nunique()
    first = g[g["iwstat"] == 1].sort_values(["pid", "wave"]).groupby("pid").first()
    n_alive_base = int((first["stroke"].notna()).sum())
    n_base_notdead = int(first["dead"].sum() == 0) if "dead" in first.columns else n_alive_base
    n_pp = pp[pp["cohort"] == c]["pid"].nunique()
    flow.append({"cohort": c, "N_total": n_all, "N_interviewed": ever,
                 "N_baseline": n_alive_base, "N_analysed": n_pp,
                 "N_excluded": n_all - n_pp})
flow = pd.DataFrame(flow)
print("流程图数字：")
print(flow.to_string(index=False))

fig, ax = plt.subplots(figsize=(7.4, 5.0))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(cy, title, sub, fc, ec, bold=False):
    h = 1.02
    ax.add_patch(FancyBboxPatch((2.15, cy - h / 2), 5.7, h, boxstyle="round,pad=0.05",
                                fc=fc, ec=ec, lw=1.1))
    ax.text(5, cy + 0.20, title, ha="center", va="center", fontsize=9,
            weight="bold" if bold else "normal")
    ax.text(5, cy - 0.22, sub, ha="center", va="center", fontsize=8.2,
            color="#333", weight="bold" if bold else "normal")

ax.text(5, 9.55, "Four harmonized ageing cohorts", ha="center", fontsize=10.5, weight="bold")
ax.text(5, 9.05, f"Pooled source population: N = {flow['N_total'].sum():,}",
        ha="center", fontsize=8.6, color="#555")
box(8.0, "Interviewed at ≥ 1 wave", f"N = {flow['N_interviewed'].sum():,}",
    "#eef3fa", "#3b6ea5")
box(6.55, "Alive at baseline interview", f"N = {flow['N_baseline'].sum():,}",
    "#eef3fa", "#3b6ea5")
box(5.1, "Analysed (≥ 2 waves, complete covariates)",
    f"N = {flow['N_analysed'].sum():,}   |   deaths = {int(pp['event_all'].sum()):,}",
    "#e8f5ec", "#2e7d4f", bold=True)
for y0, y1 in [(7.49, 7.06), (6.04, 5.61)]:
    ax.add_patch(FancyArrowPatch((5, y0), (5, y1), arrowstyle="-|>",
                                 mutation_scale=10, color="#666", lw=0.9))

y0 = 2.55
for i, c in enumerate(COH):
    x = 0.28 + i * 2.44
    r = flow[flow["cohort"] == c].iloc[0]
    ax.add_patch(FancyBboxPatch((x, y0), 2.25, 1.5, boxstyle="round,pad=0.05",
                                fc="white", ec="#888", lw=0.8))
    ax.text(x + 1.125, y0 + 1.20, COH_FULL[c], ha="center", fontsize=8.2, weight="bold")
    ax.text(x + 1.125, y0 + 0.83, f"source {r['N_total']:,}", ha="center",
            fontsize=7.4, color="#555")
    ax.text(x + 1.125, y0 + 0.52, f"analysed {r['N_analysed']:,}", ha="center", fontsize=7.8)
    d = int(pp[pp["cohort"] == c]["event_all"].sum())
    ax.text(x + 1.125, y0 + 0.21, f"deaths {d:,}", ha="center", fontsize=7.8, color="#b03030")
ax.add_patch(FancyArrowPatch((5, 4.55), (5, 4.10), arrowstyle="-|>",
                             mutation_scale=10, color="#666", lw=0.9))
fig.savefig(os.path.join(FIG, "Figure1_flow.png"), dpi=300)
fig.savefig(os.path.join(FIG, "Figure1_flow.pdf"))
plt.close(fig)
print("Figure 1 完成")

# ================= Table 1：基线特征（按暴露分组）=================
bl = long[long["iwstat"] == 1].sort_values(["cohort", "pid", "wave"]).groupby(["cohort", "pid"]).first().reset_index()
bl["hd"] = bl["hearte"]
rows = []
for c in COH:
    g = bl[bl["cohort"] == c]
    for lvl, lab in [(0, "No heart disease"), (1, "Heart disease")]:
        s = g[g["hd"] == lvl]
        rows.append({
            "Cohort": c, "Group": lab, "N": len(s),
            "Age, mean (SD)": f"{s['agey'].mean():.1f} ({s['agey'].std():.1f})",
            "Female, %": f"{s['female'].mean()*100:.1f}" if s["female"].notna().any() else "—",
            "Education (ISCED), %": " / ".join(
                f"{s['educl'].eq(i).mean()*100:.0f}" for i in (1, 2, 3)) if s["educl"].notna().any() else "—",
            "Hypertension, %": f"{s['hibp'].mean()*100:.1f}",
            "Diabetes, %": f"{s['diab'].mean()*100:.1f}",
            "Current smoker, %": f"{s['smokev'].mean()*100:.1f}",
            "BMI, mean": f"{s['bmi'].mean():.1f}",
            "ADL limited, %": f"{(s['adlfive'] > 0).mean()*100:.1f}",
            "Depressive symptoms, %": f"{s['dep_bin'].mean()*100:.1f}",
        })
t1 = pd.DataFrame(rows)
t1.to_csv(os.path.join(FIG, "Table1_基线特征.csv"), index=False, encoding="utf-8-sig")
print("\nTable 1 完成（%d 行）" % len(t1))

# ================= Figure 2：森林图 =================
def fit(sub, out, expo="x_hearte", cov=None):
    cov = COV if cov is None else cov
    d = sub.assign(event=sub[out])[["event", "start", "stop", expo, "pid"] + cov].dropna(
        subset=["event", "start", "stop", expo, "age"])
    if len(d) < 300 or d["event"].sum() < 20:
        return None
    d = d.copy(); d["pid2"] = d["pid"].astype("category").cat.codes
    try:
        ctv = CoxTimeVaryingFitter(penalizer=0.01)
        ctv.fit(d[["pid2", "start", "stop", "event", expo] + cov], id_col="pid2",
                start_col="start", stop_col="stop", event_col="event")
        hr = float(np.exp(ctv.params_[expo])); ci = np.exp(ctv.confidence_intervals_.loc[expo]).values
        return hr, ci[0], ci[1], int(d["event"].sum())
    except Exception:
        return None


def meta(r):
    w = 1 / r["se"] ** 2; muf = (w * r["logHR"]).sum() / w.sum()
    Q = float((w * (r["logHR"] - muf) ** 2).sum()); k = len(r)
    C = w.sum() - (w ** 2).sum() / w.sum(); t2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    wr = 1 / (r["se"] ** 2 + t2); m = float((wr * r["logHR"]).sum() / wr.sum())
    return m, float(np.sqrt(1 / wr.sum())), (max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0), t2


OUTS = [("event_cv", "Cardiovascular death", (0.85, 3.4), [1.0, 1.5, 2.0, 3.0]),
        ("event_all", "All-cause death", (0.9, 2.4), [1.0, 1.25, 1.5, 1.75, 2.0]),
        ("event_non", "Non-cardiovascular death", (0.7, 2.2), [0.75, 1.0, 1.5, 2.0])]
fig, axes = plt.subplots(1, 3, figsize=(11.6, 3.8))
for ax, (oc, ol, xlim, ticks) in zip(axes, OUTS):
    res = []
    for c in COH:
        r = fit(pp[pp["cohort"] == c], oc)
        if r:
            res.append({"coh": c, "hr": r[0], "lo": r[1], "hi": r[2],
                        "ev": r[3], "logHR": np.log(r[0]), "se": (np.log(r[2]) - np.log(r[1])) / 3.92})
    if len(res) < 2:
        continue
    r = pd.DataFrame(res)
    m, se, I2, t2 = meta(r)
    labels = [COH_FULL[x] for x in r["coh"]] + ["Pooled (random effects)"]
    hrs = list(r["hr"]) + [float(np.exp(m))]
    los = list(r["lo"]) + [float(np.exp(m - 1.96 * se))]
    his = list(r["hi"]) + [float(np.exp(m + 1.96 * se))]
    ys = np.arange(len(labels))[::-1]
    for i, (h, lo, hi, y) in enumerate(zip(hrs, los, his, ys)):
        pooled = (i == len(hrs) - 1)
        ax.plot([lo, hi], [y, y], color="#1f4e79" if pooled else "#333",
                lw=2.4 if pooled else 1.2, solid_capstyle="butt")
        ax.plot([h], [y], "s" if pooled else "o",
                ms=9 if pooled else 5.8, color="#c0392b" if pooled else "#1f4e79",
                zorder=5)
    ax.axvline(1.0, color="#999", ls="--", lw=0.8, zorder=0)
    ax.set_xscale("log"); ax.set_xlim(*xlim)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t:g}" for t in ticks], fontsize=7.4)
    ax.minorticks_off()
    ax.set_yticks(ys); ax.set_yticklabels(labels, fontsize=7.8)
    ax.set_ylim(-0.75, len(labels) - 0.25)
    ax.set_xlabel("Hazard ratio (95% CI)", fontsize=8.2)
    ax.set_title(f"{ol}\nI² = {I2:.1f}%", fontsize=9, weight="bold")
    # HR 值放在坐标区内部右侧，带白底避免与图形重叠
    for h, lo, hi, y in zip(hrs, los, his, ys):
        ax.text(0.985, y, f"{h:.2f} ({lo:.2f}–{hi:.2f})",
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=7.0, color="#222",
                bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.85), zorder=6)
fig.suptitle("Self-reported heart disease and mortality, by cohort",
             fontsize=10.8, weight="bold", y=1.03)
fig.tight_layout()
fig.subplots_adjust(wspace=0.62)
fig.savefig(os.path.join(FIG, "Figure2_forest.png"), dpi=300)
fig.savefig(os.path.join(FIG, "Figure2_forest.pdf"))
plt.close(fig)
print("Figure 2 完成")

# ================= Figure 3：累积发生率曲线 =================
def aj_curve(sub, event_col, comp_col, tmax=5.0, step=0.25):
    times = np.arange(step, tmax + 1e-9, step)
    s = sub.dropna(subset=["start", "stop", event_col, comp_col])
    ce, S, out = 0.0, 1.0, []
    for t in times:
        at_risk = (s["stop"] >= t - 1e-9).sum()
        d_e = ((s["stop"] <= t) & (s["stop"] > t - step - 1e-9) & (s[event_col] == 1)).sum()
        d_c = ((s["stop"] <= t) & (s["stop"] > t - step - 1e-9) & (s[comp_col] == 1)).sum()
        if at_risk > 0:
            ce += S * d_e / at_risk
            S *= (1 - (d_e + d_c) / at_risk)
        out.append(ce)
    return times, np.array(out)


fig, axes = plt.subplots(1, 4, figsize=(12, 3.1), sharey=True)
for ax, c in zip(axes, COH):
    for lvl, lab, col, ls in [(0, "No heart disease", "#1f4e79", "-"),
                              (1, "Heart disease", "#c0392b", "-")]:
        g = pp[(pp["cohort"] == c) & (pp["x_hearte"] == lvl)]
        if len(g) < 100:
            continue
        t, cif = aj_curve(g, "event_cv", "event_non")
        ax.plot(t, cif * 100, color=col, ls=ls, lw=1.6, label=lab)
    ax.set_title(COH_FULL[c], fontsize=8.6, weight="bold")
    ax.set_xlabel("Years since baseline", fontsize=8)
    ax.set_xlim(0, 5)
    ax.tick_params(labelsize=7.5)
    ax.grid(alpha=0.25, lw=0.5)
axes[0].set_ylabel("Cumulative incidence of CV death (%)", fontsize=8)
axes[0].legend(fontsize=7.2, frameon=False, loc="upper left")
for ax in axes:
    ax.set_ylim(0, 23)
fig.suptitle("Cumulative incidence of cardiovascular death (Aalen–Johansen, competing risk)",
             fontsize=9.8, weight="bold", y=1.06)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "Figure3_CIF.png"), dpi=300)
fig.savefig(os.path.join(FIG, "Figure3_CIF.pdf"))
plt.close(fig)
print("Figure 3 完成")

# ================= Table 2：主分析 + 敏感性汇总 =================
t2 = pd.read_csv(os.path.join(OUT, "S2_敏感性分析.csv"))
t2.to_csv(os.path.join(FIG, "Table2_主分析与敏感性.csv"), index=False, encoding="utf-8-sig")
print("Table 2 完成")
print(f"\n全部图表已输出到：{FIG}")
print(os.listdir(FIG))
