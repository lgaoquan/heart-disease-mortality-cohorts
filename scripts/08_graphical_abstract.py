# -*- coding: utf-8 -*-
"""
Generate the graphical abstract.

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
"""Graphical abstract（EJPC 要求）——重排版"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

FIG = r"DATA_ROOT\图表"
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 300,
                     "savefig.bbox": "tight", "savefig.facecolor": "white"})

fig = plt.figure(figsize=(11.5, 6.6))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
NAVY, RED, GREY, GREEN = "#1f4e79", "#c0392b", "#6b7280", "#2e7d4f"

# ---------- 标题 ----------
ax.text(50, 97.0, "Self-reported heart disease and mortality across four ageing cohorts",
        ha="center", fontsize=13.5, weight="bold", color="#111")
ax.text(50, 92.8, "Pooled analysis of 219,364 participants and 34,571 deaths",
        ha="center", fontsize=10.5, color=GREY)

# ================= 上半区：y 52–88 =================
# ---------- 面板 1 ----------
ax.text(3, 88.0, "1  Four harmonised cohorts", fontsize=10.5, weight="bold", color=NAVY)
cohorts = [("HRS", "USA", 42023, 7426), ("ELSA", "England", 18523, 2706),
           ("SHARE", "Europe", 137684, 21578), ("CHARLS", "China", 21134, 3000)]
for i, (nm, ctry, n, d) in enumerate(cohorts):
    yc = 83.0 - i * 7.6
    ax.add_patch(FancyBboxPatch((3, yc - 3.0), 26, 6.0, boxstyle="round,pad=0.25",
                                fc="#eef3fa", ec=NAVY, lw=1.0))
    ax.text(5.0, yc + 0.9, f"{nm}  ({ctry})", fontsize=9.4, weight="bold", color="#111")
    ax.text(5.0, yc - 1.5, f"n = {n:,}    deaths = {d:,}", fontsize=8.0, color="#333")
ax.add_patch(FancyArrowPatch((31.0, 70.0), (34.8, 70.0), arrowstyle="-|>",
                             mutation_scale=15, color=GREY, lw=1.5))

# ---------- 面板 2 ----------
ax.text(37, 88.0, "2  Exposure and design", fontsize=10.5, weight="bold", color=NAVY)
ax.add_patch(FancyBboxPatch((37, 76.5), 26, 8.6, boxstyle="round,pad=0.3",
                            fc="#fff4e6", ec="#b7791f", lw=1.1))
ax.text(50, 83.4, "Self-reported", ha="center", fontsize=9.2, color="#111")
ax.text(50, 80.6, "physician-diagnosed", ha="center", fontsize=9.2, color="#111")
ax.text(50, 77.8, "heart disease", ha="center", fontsize=9.4, weight="bold", color="#111")
for j, s in enumerate(["time-varying exposure", "5-year follow-up window",
                       "cohort-specific Cox models", "random-effects meta-analysis"]):
    ax.text(50, 72.6 - j * 3.0, s, ha="center", fontsize=8.2, style="italic", color=GREY)
ax.add_patch(FancyArrowPatch((64.0, 70.0), (67.8, 70.0), arrowstyle="-|>",
                             mutation_scale=15, color=GREY, lw=1.5))

# ---------- 面板 3 ----------
ax.text(70, 88.0, "3  Findings", fontsize=10.5, weight="bold", color=NAVY)
ax.text(70, 84.2, "Mortality gradient by cause", fontsize=9.0, weight="bold", color="#111")
def X(v): return 73 + (v - 0.9) / (2.6 - 0.9) * 22
grad = [("Cardiovascular", 1.96, 1.56, 2.46, RED),
        ("All-cause", 1.40, 1.16, 1.68, NAVY),
        ("Non-cardiovascular", 1.16, 0.98, 1.36, GREY)]
for i, (lab, hr, lo, hi, col) in enumerate(grad):
    yc = 78.6 - i * 6.4
    ax.text(72.0, yc + 1.9, lab, fontsize=8.2, color="#111", ha="right")
    ax.plot([X(lo), X(hi)], [yc - 0.8, yc - 0.8], color=col, lw=2.2, solid_capstyle="butt")
    ax.plot([X(hr)], [yc - 0.8], "o", ms=5.8, color=col, zorder=5)
    ax.text(97.6, yc - 0.8, f"{hr:.2f}", fontsize=8.2, va="center", color=col, weight="bold")
ax.plot([X(1.0), X(1.0)], [63.5, 80.5], color="#999", ls="--", lw=0.9)
ax.text(X(1.0), 61.8, "1.0", fontsize=7.4, ha="center", color=GREY)
ax.text(86, 59.2, "Pooled hazard ratio (95% CI)", fontsize=7.6, ha="center", color=GREY)

# ================= 下半区：y 15–52 =================
ax.text(3, 51.5, "4  Absolute excess risk of cardiovascular death over 5 years",
        fontsize=10.5, weight="bold", color=NAVY)
ard = [("CHARLS", 1.99), ("ELSA", 5.70), ("SHARE", 9.67), ("HRS", 16.04)]
maxv, x0, span = 17.5, 10.5, 52
for i, (nm, v) in enumerate(ard):
    yc = 44.0 - i * 6.8
    ax.text(3, yc + 0.6, nm, fontsize=8.8, weight="bold", color="#111")
    w = v / maxv * span
    ax.add_patch(Rectangle((x0, yc - 1.8), w, 3.6, fc="#cfe2f3", ec=NAVY, lw=0.8))
    ax.text(x0 + w + 1.3, yc, f"{v:.2f} pp", fontsize=8.6, va="center",
            weight="bold", color=NAVY)
# 8 倍差距标注放在条下方
ax.annotate("", xy=(x0 + span, 19.6), xytext=(x0, 19.6),
            arrowprops=dict(arrowstyle="<->", color=RED, lw=1.3))
ax.text((x0 + x0 + span) / 2, 16.8,
        "8-fold range in absolute risk difference   vs   1.5-fold range in hazard ratio",
        fontsize=8.6, ha="center", weight="bold", color=RED)

# ---------- 结论条 ----------
ax.add_patch(FancyBboxPatch((3, 1.8), 94, 11.5, boxstyle="round,pad=0.4",
                            fc="#e8f5ec", ec=GREEN, lw=1.3))
ax.text(50, 10.2, "Self-reported heart disease predicts death in every cohort, and the association is cardiovascular-specific.",
        ha="center", fontsize=9.2, weight="bold", color="#1b5e20")
ax.text(50, 6.9, "But its prognostic weight differs substantially between countries, and relative risks conceal far larger differences in absolute risk.",
        ha="center", fontsize=9.0, color="#1b5e20")
ax.text(50, 3.9, "A self-reported diagnosis should not be treated as an equivalent exposure in cross-national research.",
        ha="center", fontsize=8.8, style="italic", color="#2e7d4f")

fig.savefig(os.path.join(FIG, "GraphicalAbstract.png"), dpi=300)
fig.savefig(os.path.join(FIG, "GraphicalAbstract.pdf"))
plt.close(fig)
print("Graphical abstract 完成（重排版）")
