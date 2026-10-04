# Self-reported heart disease and mortality across four ageing cohorts

This repository contains the analysis scripts used for the manuscript:

> **Self-reported heart disease and mortality across four ageing cohorts: an individual-participant data cohort analysis**

The study uses four harmonised ageing resources: the Health and Retirement Study (HRS), the English Longitudinal Study of Ageing (ELSA), the Survey of Health, Ageing and Retirement in Europe (SHARE), and the China Health and Retirement Longitudinal Study (CHARLS). The current manuscript reports 171,699 participants and 16,128 deaths within the primary five-year analysis window.

## Scope of this release

The repository is a code release. It contains no individual-level data, source cohort files, restricted-use files, derived person-level files, manuscript files, local audit records, credentials, or submission materials. The scripts reproduce the audit-driven reconstruction, survival models, sensitivity analyses, validation checks and figure generation when the user supplies the permitted source data and configures the required local paths.

The revised analysis is exploratory and post hoc rather than preregistered or confirmatory. Numerical results in the manuscript should be interpreted together with its reported heterogeneity, measurement differences and death-timing limitations.

## Data access and non-redistribution

The source resources must be obtained directly from their official custodians or authorised access portals. Registration, permission, citation and non-redistribution conditions apply. The authors do not redistribute individual-level data.

| Resource | Version used in the manuscript | Official access route |
|---|---|---|
| Health and Retirement Study (HRS) | RAND HRS Longitudinal File 2020 Version 1; Harmonized HRS Version D; Harmonized HRS End of Life Version B | [HRS data products](https://hrs.isr.umich.edu/data-products); [Gateway to Global Aging Data](https://g2aging.org/) |
| English Longitudinal Study of Ageing (ELSA) | Harmonized ELSA Version H; Harmonized ELSA End of Life Version A.2 | [UK Data Service, study 5050](https://beta.ukdataservice.ac.uk/datacatalogue/studies/study?id=5050); [Gateway to Global Aging Data](https://g2aging.org/) |
| Survey of Health, Ageing and Retirement in Europe (SHARE) | Harmonized SHARE Version G; Harmonized SHARE End of Life Version G | [SHARE data access](https://share-eric.eu/data/data-access); [Gateway to Global Aging Data](https://g2aging.org/) |
| China Health and Retirement Longitudinal Study (CHARLS) | Harmonized CHARLS Version D2; Harmonized CHARLS End of Life Version A | [CHARLS](http://charls.pku.edu.cn/en); [Gateway to Global Aging Data](https://g2aging.org/) |

Users must check the current custodian terms before obtaining or using any release. The repository does not grant access to the data and does not permit redistribution of source or individual-level files.

## Code availability

The analysis code is publicly available at [github.com/lgaoquan/heart-disease-mortality-cohorts](https://github.com/lgaoquan/heart-disease-mortality-cohorts). The repository does not contain individual-level data or restricted source datasets. Users must obtain the relevant datasets directly from the official custodians and comply with their access, citation and non-redistribution terms.

## Software and execution

The release was developed with Python and R. Package versions used for the reported run are recorded in `requirements.txt`, the scripts and the supplementary statistical methods. The scripts use an environment variable named `ANALYSIS_OUTPUT_ROOT` for a local working directory. `rebuild.py` additionally requires a private source project configured through `SOURCE_PROJECT_ROOT`; this private source project is not included here.

A typical local sequence is:

```bash
export ANALYSIS_OUTPUT_ROOT=/path/to/local/analysis_output
python scripts/rebuild.py
python scripts/make_analysis.py
Rscript scripts/fit_models.R
Rscript scripts/crude_models.R
Rscript scripts/exact_risks.R
python scripts/validate_results.py
python scripts/draw_figures.py
```

The sequence is a guide for authorised users with the required data and source files. It is not a claim that the restricted source data are downloadable from this repository.

## Citation

Please cite the associated manuscript and the source cohort datasets according to each custodian's instructions. See `CITATION.cff` for repository citation metadata.

## License

The code is released under the MIT License. The license does not alter any third-party data-use or non-redistribution condition. See `LICENSE`.
