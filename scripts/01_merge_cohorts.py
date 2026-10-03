# -*- coding: utf-8 -*-
"""
Merge the four harmonised cohorts into a single long-format table.

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
四库合并：抽取映射变量 -> 统一编码 -> 纵向堆叠为分析用长表
输出：DATA_ROOT\\分析数据\\四库合并长表.parquet  +  覆盖情况报告

设计要点
 1. HRS 侧 = RAND HRS（主） + Harmonized HRS（补 ADL / 定向力 / 教育类别），按 hhidpn 合并
 2. 只用 usecols 读取所需列，避免把 1.8 GB 的 SHARE 全读进内存
 3. 抑郁默认用「单条目二分类」（四库语义一致）；同时保留量表分变量备用
 4. BMI 同时保留连续值与 WHO 分类
 5. 波次年份优先取数据内的访谈年份变量；取不到再报告缺失
"""
import os, re, sys, json
import pandas as pd
import pyreadstat

OUT = r"DATA_ROOT\分析数据"
os.makedirs(OUT, exist_ok=True)

FILES = {
    "HRS":    r"DATA_ROOT\整理后数据库\HRS\数据\RAND_HRS\randhrs1992_2020v1.dta",
    "HRS_HD": r"DATA_ROOT\整理后数据库\HRS\数据\Harmonized_HRS\H_HRS_d.dta",
    "ELSA":   r"DATA_ROOT\整理后数据库\ELSA\数据\Harmonized_ELSA\gh_elsa_h.dta",
    "SHARE":  r"DATA_ROOT\整理后数据库\SHARE\数据\GH_SHARE_g.dta",
    "CHARLS": r"DATA_ROOT\整理后数据库\CHARLS\数据\Harmonized_CHARLS\GH_CHARLS_d2.dta",
}
# 各库取哪个文件
SRC = {"HRS": "HRS", "ELSA": "ELSA", "SHARE": "SHARE", "CHARLS": "CHARLS"}
# 时变概念 -> 各库变量词根（已实测核实）
TV = {
    "stroke":    {"HRS": "stroke",  "ELSA": "stroke",  "SHARE": "stroke",  "CHARLS": "stroke"},
    "diab":      {"HRS": "diab",    "ELSA": "diabe",   "SHARE": "diabe",   "CHARLS": "diabe"},
    "hibp":      {"HRS": "hibp",    "ELSA": "hibpe",   "SHARE": "hibpe",   "CHARLS": "hibpe"},
    "hearte":    {"HRS": "hearte",  "ELSA": "hearte",  "SHARE": "hearte",  "CHARLS": "hearte"},
    "cancre":    {"HRS": "cancre",  "ELSA": "cancre",  "SHARE": "cancre",  "CHARLS": "cancre"},
    "lunge":     {"HRS": "lunge",   "ELSA": "lunge",   "SHARE": "lunge",   "CHARLS": "lunge"},
    "depressed": {"HRS": "depres",  "ELSA": "depres",  "SHARE": "depress", "CHARLS": "depresl"},
    "cesd":      {"HRS": "cesd",    "ELSA": "cesd",    "SHARE": "eurod",   "CHARLS": "cesdl10"},
    "smokev":    {"HRS": "smokev",  "ELSA": "smokev",  "SHARE": "smokev",  "CHARLS": "smokev"},
    "shlt":      {"HRS": "shlt",    "ELSA": "shlt",    "SHARE": "shlt",    "CHARLS": "shlt"},
    "iwstat":    {"HRS": "iwstat",  "ELSA": "iwstat",  "SHARE": "iwstat",  "CHARLS": "iwstat"},
    "hosp":      {"HRS": "hosp",    "ELSA": None,      "SHARE": "hosp1y",  "CHARLS": "hosp1y"},
    "bmi":       {"HRS": "bmi",     "ELSA": "mbmi",    "SHARE": "bmi",     "CHARLS": "mbmi"},
    "adlfive":   {"HRS": "adlfive", "ELSA": "adlfive", "SHARE": "adlfive", "CHARLS": "adlfive"},
    "orient":    {"HRS": "orient",  "ELSA": "orient",  "SHARE": "orient",  "CHARLS": "orient"},
    "cogtot":    {"HRS": "cogtot",  "ELSA": "cogimp",  "SHARE": "cogimp",  "CHARLS": None},
    "agey":      {"HRS": "agey_b",  "ELSA": "agey",    "SHARE": "agey",    "CHARLS": "agey"},
    "iwy":       {"HRS": "iwendy",  "ELSA": "iwy",     "SHARE": "iwy",     "CHARLS": "iwy"},
    # 治疗类（作治疗可及性代理）：HRS 侧从 Harmonized HRS 取
    "rxdiab":    {"HRS": "rxdiabo", "ELSA": "rxdiab",  "SHARE": "rxdiab",  "CHARLS": "rxdiab"},
    "rxhibp":    {"HRS": "rxhibp",  "ELSA": "rxhibp",  "SHARE": "rxhibp",  "CHARLS": "rxhibp"},
    "rxheart":   {"HRS": "rxheart", "ELSA": "rxheart", "SHARE": "rxheart", "CHARLS": "rxheart"},
}
# HRS 从 Harmonized HRS 取的变量
HRS_FROM_HD = {"adlfive", "orient", "rxdiab", "rxhibp", "rxheart"}
# 个体层里需从 Harmonized HRS 取的（HRS 专用）
BASE_FROM_HD = {"educl": "raeducl"}
# 个体层
BASE = {
    "gender": {"HRS": "ragender", "ELSA": "ragender", "SHARE": "ragender", "CHARLS": "ragender"},
    "byear":  {"HRS": "rabyear",  "ELSA": "rabyear",  "SHARE": "rabyear",  "CHARLS": "rabyear"},
    "race":   {"HRS": "raracem",  "ELSA": "raracem",  "SHARE": None,       "CHARLS": None},
    "educl":  {"HRS": "raeducl",  "ELSA": "raeducl",  "SHARE": "raeducl",  "CHARLS": "raeducl"},
    "edyrs":  {"HRS": "raedyrs",  "ELSA": "raedyrs_e","SHARE": "raedyrs",  "CHARLS": None},
}
ID = {"HRS": "hhidpn", "ELSA": "idauniq", "SHARE": "mergeid", "CHARLS": "ID"}

META = {}
for k, p in FILES.items():
    _, m = pyreadstat.read_dta(p, metadataonly=True, output_format="dict")
    META[k] = {v.lower(): v for v in m.column_names}     # lower -> 原始大小写
    print(f"读取元数据 {k:<8} {len(m.column_names):>6} 变量")

def resolve(src, stem, wave):
    """在 src 中找 r<wave><stem> / h<wave><stem> 的原始变量名"""
    if stem is None:
        return None
    for pre in ("r", "h"):
        key = f"{pre}{wave}{stem}".lower()
        if key in META[src]:
            return META[src][key]
    return None

def resolve_base(src, stem):
    if stem is None:
        return None
    return META[src].get(stem.lower())

# ---------- 逐库抽取 ----------
long_frames = []
report = []

for coh in ["HRS", "ELSA", "SHARE", "CHARLS"]:
    src_main = SRC[coh]
    waves = range(1, 16)

    # 需要读取的列
    want_main, want_hd, colmap = set(), set(), {}
    iid = resolve_base(src_main, ID[coh])
    if iid:
        want_main.add(iid)
    for c, stem in BASE.items():
        if coh == "HRS" and c in BASE_FROM_HD:
            v = resolve_base("HRS_HD", BASE_FROM_HD[c])
            if v:
                want_hd.add(v)
                continue
        v = resolve_base(src_main, stem.get(coh))
        if v:
            want_main.add(v); colmap[c] = v
        else:
            report.append({"队列": coh, "类型": "个体层", "概念": c, "变量名": "", "状态": "✘"})

    for c, per in TV.items():
        stem = per.get(coh)
        found = {}
        for w in waves:
            if coh == "HRS" and c in HRS_FROM_HD:
                v = resolve("HRS_HD", stem, w)
                if v:
                    want_hd.add(v); found[w] = v
            else:
                v = resolve(src_main, stem, w)
                if v:
                    want_main.add(v); found[w] = v
        if found:
            colmap[c] = found
            report.append({"队列": coh, "类型": "时变", "概念": c,
                           "变量名": f"{min(found)}-{max(found)}", "状态": f"✔ {len(found)}波"})
        else:
            report.append({"队列": coh, "类型": "时变", "概念": c, "变量名": "", "状态": "✘"})

    print(f"\n{coh}: 主文件读 {len(want_main)} 列" + (f"，Harmonized HRS 读 {len(want_hd)} 列" if want_hd else ""))
    df = pyreadstat.read_dta(FILES[src_main], usecols=sorted(want_main), output_format="pandas")[0]
    df.columns = [c.lower() for c in df.columns]

    if want_hd:
        hd_key = resolve_base("HRS_HD", ID["HRS"])
        want_hd.add(hd_key)
        hd = pyreadstat.read_dta(FILES["HRS_HD"], usecols=sorted(want_hd), output_format="pandas")[0]
        hd.columns = [c.lower() for c in hd.columns]
        df = df.merge(hd, left_on=iid.lower(), right_on=hd_key.lower(), how="left")
        print(f"   已与 Harmonized HRS 合并，{len(df):,} 行")

    # 宽转长
    recs = []
    for w in waves:
        cols, names = [], []
        for c, v in colmap.items():
            if c in BASE:
                continue
            vv = v.get(w) if isinstance(v, dict) else None
            if vv and vv.lower() in df.columns:
                cols.append(vv.lower()); names.append(c)
        if not cols:
            continue
        sub = df[[iid.lower()] + cols].copy()
        sub.columns = ["pid"] + names
        sub.insert(0, "cohort", coh)
        sub.insert(2, "wave", w)
        recs.append(sub)
    if recs:
        long_frames.append(pd.concat(recs, ignore_index=True))

long = pd.concat(long_frames, ignore_index=True)

# 个体层广播
base_rows = []
for coh in ["HRS", "ELSA", "SHARE", "CHARLS"]:
    src = SRC[coh]
    iid = resolve_base(src, ID[coh])
    cols, names = [iid], ["pid"]
    extra = None
    for c, per in BASE.items():
        if coh == "HRS" and c in BASE_FROM_HD:
            v = resolve_base("HRS_HD", BASE_FROM_HD[c])
            if v:
                extra = (v, c)
            continue
        v = resolve_base(src, per.get(coh))
        if v:
            cols.append(v); names.append(c)
    b = pyreadstat.read_dta(FILES[src], usecols=cols, output_format="pandas")[0]
    # 注意：pyreadstat 的 usecols 不保证返回顺序，必须按名重命名，不能按位置赋值
    ren = {iid: "pid"}
    for c, per in BASE.items():
        v = resolve_base(src, per.get(coh))
        if v:
            ren[v] = c
    b = b.rename(columns=ren)[list(ren.values())]
    if extra:
        v, c = extra
        hk = resolve_base("HRS_HD", ID["HRS"])
        h = pyreadstat.read_dta(FILES["HRS_HD"], usecols=[hk, v], output_format="pandas")[0]
        h.columns = ["pid", c]
        b = b.merge(h, on="pid", how="left")
    b.insert(0, "cohort", coh)
    base_rows.append(b)
base = pd.concat(base_rows, ignore_index=True)
base["pid"] = base["pid"].astype(str)
long["pid"] = long["pid"].astype(str)
long = long.merge(base, on=["cohort", "pid"], how="left")

# 抑郁条目二分化：CHARLS 为 1-4 量表，其余为 0/1
if "depressed" in long.columns:
    def to_bin(row):
        v = row["depressed"]
        if pd.isna(v):
            return float("nan")
        if row["cohort"] == "CHARLS":
            return float(v >= 3)          # ≥3 = 偶尔或更多（3-4 天/周）
        return float(v) if v in (0, 1) else float("nan")
    long["dep_bin"] = long.apply(to_bin, axis=1)

# BMI 分类
if "bmi" in long.columns:
    long["bmi_cat"] = pd.cut(long["bmi"], [0, 18.5, 25, 30, 100],
                             labels=["<18.5", "18.5-24.9", "25-29.9", ">=30"])

# 抑郁量表分按队列 z 标准化（量表不同，跨库不可直接比）
if "cesd" in long.columns:
    long["cesd_z"] = long.groupby("cohort")["cesd"].transform(
        lambda s: (s - s.mean()) / s.std())

out_pq = os.path.join(OUT, "四库合并长表.parquet")
long.to_parquet(out_pq, index=False)
rep = pd.DataFrame(report)
rep.to_csv(os.path.join(OUT, "变量覆盖报告.csv"), index=False, encoding="utf-8-sig")

print(f"\n{'='*70}")
print(f"长表：{len(long):,} 行 × {len(long.columns)} 列  ->  {out_pq}")
print(f"{'='*70}")
print(long.groupby("cohort").agg(行数=("pid", "size"), 唯一ID=("pid", "nunique"),
                                 波次数=("wave", "nunique")))
print("\n各列缺失率（前 25）：")
mr = (long.isna().mean() * 100).round(1).sort_values()
print(mr.head(25))
