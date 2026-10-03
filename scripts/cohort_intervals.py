# -*- coding: utf-8 -*-
"""
Shared module: build person-interval (start-stop) data for survival analysis.

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
共享模块：向量化构造 person-interval 数据（替代逐人 Python 循环）
修复两个 bug：
  1) cols = EXPO + BASEV 产生重复列名 -> 去重
  2) 逐人循环太慢 -> 全向量化
"""
import numpy as np
import pandas as pd

DEATH = {2, 3, 5, 6}
EXPO = ["diab", "hibp", "hearte", "smokev", "cesd_z", "bmi_cat", "adlfive", "stroke"]
BASEV = ["smokev", "hibp", "hearte", "bmi_cat", "adlfive"]


def prepare(long):
    d = long.copy()
    for c in ["stroke", "diab", "hibp", "hearte", "smokev"]:
        d[c] = d[c].where(d[c].isin([0, 1]))
    for c in ["bmi", "agey", "iwy", "adlfive", "cesd_z"]:
        d.loc[d[c] < -10, c] = np.nan
    d["educl"] = d["educl"].where(d["educl"].isin([1, 2, 3]))
    d["female"] = np.where(d["gender"].isin([1, 2]), (d["gender"] == 2).astype(float), np.nan)
    d["bmi_cat"] = pd.cut(d["bmi"], [0, 18.5, 25, 30, 100], labels=[0, 1, 2, 3]).astype(float)
    d["adlfive"] = d["adlfive"].clip(0, 5)
    d["dead"] = d["iwstat"].isin(DEATH).astype(int)
    d["alive_iw"] = (d["iwstat"] == 1).astype(int)
    return d


def intervals(d, outcome, window=5):
    """返回 person-interval 表：start/stop/event/x_*/{cov}_b/age/female/educl"""
    d = d.sort_values(["cohort", "pid", "wave"]).reset_index(drop=True)
    g = d.groupby(["cohort", "pid"], sort=False)

    # 年份估计
    d["iwy_f"] = g["iwy"].ffill()
    d["w_at_iyw"] = d["wave"].where(d["iwy"].notna())
    d["w_at_iyw"] = d.groupby(["cohort", "pid"], sort=False)["w_at_iyw"].ffill()
    d["year_est"] = d["iwy"].fillna(d["iwy_f"] + 2.0 * (d["wave"] - d["w_at_iyw"]))

    # 基线 = 首个已访波
    g = d.groupby(["cohort", "pid"], sort=False)
    d["is_base"] = (d["alive_iw"] == 1) & (g["alive_iw"].cumsum() == 1)
    b = d[d["is_base"]][["cohort", "pid", "wave", "year_est", "agey", "female", "educl"] + BASEV].copy()
    ren = {"wave": "base_wave", "year_est": "base_year", "agey": "age0",
           "female": "female_b", "educl": "educl_b"}
    for c in BASEV:
        ren[c] = f"{c}_b"
    b = b.rename(columns=ren)          # long 表已有 byear(出生年)/female/educl，基线列全部改名避免冲突

    # 暴露前向填充（时变）
    g = d.groupby(["cohort", "pid"], sort=False)
    for c in EXPO:
        d[f"x_{c}"] = g[c].ffill()

    # 下一个波的信息
    g = d.groupby(["cohort", "pid"], sort=False)
    d["n_year"] = g["year_est"].shift(-1)
    d["n_ev"] = g[outcome].shift(-1)

    d = d.merge(b, on=["cohort", "pid"], how="inner")
    d["start"] = d["year_est"] - d["base_year"]
    d["stop"] = d["n_year"] - d["base_year"]
    d["event"] = d["n_ev"]
    d["age"] = d["age0"] + d["start"]

    m = d["wave"].notna() & (d["wave"] >= d["base_wave"]) & d["start"].notna() & d["stop"].notna()
    m &= (d["stop"] > d["start"]) & (d["start"] >= 0) & (d["start"] < window) & d["event"].notna()
    m &= (d["dead"].shift(0).fillna(0) >= 0)
    out = d[m].copy()
    out["stop"] = np.minimum(out["stop"], window)
    out["event"] = out["event"].astype(int)
    # 基线已死者排除
    base_dead = d[d["is_base"] & (d["dead"] == 1)][["cohort", "pid"]]
    base_dead["_drop"] = 1
    out = out.merge(base_dead, on=["cohort", "pid"], how="left")
    out = out[out["_drop"].isna()].drop(columns=["_drop"])
    # 协变量缺失 -> 单列一类，不删行
    for c in BASEV:
        out[f"{c}_b"] = out[f"{c}_b"].fillna(-1)
    out["educl_b"] = out["educl_b"].fillna(0)
    out["female_b"] = out["female_b"].fillna(-1)
    return out
