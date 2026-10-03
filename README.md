# Self-reported heart disease and mortality across four ageing cohorts

Analysis code for:

> **Self-reported heart disease and mortality across four ageing cohorts: a pooled analysis of 219,364 participants and 34,571 deaths**
>
> Ruibin Fu, Liu Liu, Che Jiang, Hanqiang Du, Gaoquan Luo
> General Hospital of the Southern Theater Command, Guangzhou, China

---

## What this code does

The analysis pools individual-level data from four harmonised ageing cohorts and estimates the association between self-reported physician-diagnosed heart disease and mortality, with cross-national heterogeneity as the primary focus.

| Cohort | Country/region | Harmonised version |
|---|---|---|
| HRS | United States | Gateway Harmonized HRS Version D (+ RAND HRS Longitudinal File 2020 V1) |
| ELSA | England | Gateway Harmonized ELSA Version H |
| SHARE | Europe | Gateway Harmonized SHARE Version G |
| CHARLS | China | Harmonized CHARLS Version D2 |

Cause-of-death outcomes use the harmonised End of Life files (HRS Version B, ELSA Version A.2, SHARE Version G, CHARLS Version A).

**Analytical approach.** Cohort-specific Cox models with a time-varying exposure and a uniform 5-year follow-up window, pooled by random-effects meta-analysis (DerSimonian-Laird). Absolute risks are estimated with the Aalen-Johansen estimator, accounting for competing causes of death. Heterogeneity is assessed with I², τ², and leave-one-out re-estimation.

---

## Data access

**Data are not distributed with this repository.** All four cohorts require registration and agreement to their respective terms of use.

| Cohort | Access |
|---|---|
| HRS | https://hrs.isr.umich.edu |
| ELSA | UK Data Service, study 5050 |
| SHARE | https://share-eric.eu |
| CHARLS | http://charls.pku.edu.cn |
| Harmonised versions (all four) | https://g2aging.org |
| Harmonised End of Life files | https://eol.g2aging.org |

The code expects the following files under `$DATA_ROOT`:

```
$DATA_ROOT/
  randhrs1992_2020v1.dta        # RAND HRS Longitudinal File
  H_HRS_d.dta                   # Gateway Harmonized HRS
  gh_elsa_h.dta                 # Gateway Harmonized ELSA
  GH_SHARE_g.dta                # Gateway Harmonized SHARE
  GH_CHARLS_d2.dta              # Harmonized CHARLS
  GH_HRS_EOL_b.dta              # Gateway Harmonized HRS End of Life
  h_elsa_eol_a2.dta             # Harmonized ELSA End of Life
  GH_SHARE_EOL_g.dta            # Gateway Harmonized SHARE End of Life
  H_CHARLS_EOL_a.dta            # Harmonized CHARLS End of Life
```

Outputs are written to `$DATA_ROOT/analysis_output/`.

---

## Requirements

Python 3.13 or later.

```
pip install -r requirements.txt
```

---

## Running the analysis

Set `DATA_ROOT` to the folder holding the data files, then run the scripts in order. Each script reads the outputs of the previous ones.

```bash
export DATA_ROOT=/path/to/data

python scripts/01_merge_cohorts.py            # merge cohorts -> long-format table
python scripts/02_baseline_table.py           # baseline characteristics
python scripts/03_mortality_analysis.py       # main analysis (person-intervals)
python scripts/04_sensitivity_analyses.py     # absolute risk + sensitivity analyses
python scripts/05_exploratory_scan.py         # exploratory scan (supplementary Table S1)
python scripts/06_heterogeneity_diagnostics.py # leave-one-out diagnostics
python scripts/07_figures.py                  # Figure 1-3, Table 1-2
python scripts/08_graphical_abstract.py       # graphical abstract
python scripts/09_supplementary_and_wordcount.py
```

`cohort_intervals.py` is a shared module imported by the other scripts; it is not run directly.

---

## Script reference

| Script | Purpose |
|---|---|
| `01_merge_cohorts.py` | Extracts mapped variables from the five source files, harmonises coding, and stacks them into a long-format table. The HRS cohort requires merging the RAND file with the Gateway Harmonized file on `hhidpn`. |
| `cohort_intervals.py` | Shared module. Builds person-interval (start-stop) data for survival analysis, including interview-year estimation, baseline definition, and carry-forward of time-varying exposures. |
| `02_baseline_table.py` | Baseline characteristics and event counts by cohort. |
| `03_mortality_analysis.py` | Main analysis: cohort-specific time-varying Cox models for all-cause, cardiovascular and non-cardiovascular mortality, random-effects pooling, and leave-one-out diagnostics. |
| `04_sensitivity_analyses.py` | Aalen-Johansen cumulative incidence and absolute risk differences; 3-year window; restriction to HRS and SHARE; sex stratification; treatment-status analysis. |
| `05_exploratory_scan.py` | Fits the same pipeline to all exposure-outcome combinations (supplementary Table S1). |
| `06_heterogeneity_diagnostics.py` | Leave-one-out diagnostics for the highest-heterogeneity combinations. |
| `07_figures.py` | Figure 1 (flow), Figure 2 (forest), Figure 3 (cumulative incidence), Table 1, Table 2. |
| `08_graphical_abstract.py` | Graphical abstract. |
| `09_supplementary_and_wordcount.py` | Supplementary Tables S1-S4 and manuscript word counts. |

---

## Implementation notes

Two properties of the source files required care and are worth flagging for anyone reusing this code.

**`pyreadstat` does not preserve the order of `usecols`.** Columns must be renamed by name, not assigned by position. Assigning positionally silently misaligns columns; during development this placed birth year values in the sex variable.

**Missing-value codes differ across cohorts.** The Gateway harmonised files encode different types of missingness as negative values (-8, -9, -13, -15, -16, -18). Binary variables must be filtered with `isin([0, 1])` and continuous variables with a `< -10` cut-off before any summary statistic is computed. Otherwise means are contaminated; before this was corrected, the smoking prevalence in ELSA was computed as -40.7%.

Two further points. Interview year is named `rWiwendy` in HRS but `rWiwy` in the other three cohorts. The CHARLS depressive symptom item `rWdepresl` is an ordinal 1-4 scale rather than binary, and is dichotomised at ≥3.

---

## Citation

If you use this code, please cite the paper above. Data access is governed by the terms of use of each cohort; users are responsible for complying with them.

## Licence

MIT. See `LICENSE`.
