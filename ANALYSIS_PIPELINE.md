# Analysis pipeline notes

The release sequence is:

1. `rebuild.py` reconstructs interview and endpoint records from the authorised local source project.
2. `make_analysis.py` constructs person-level and counting-process analysis files.
3. `fit_models.R` fits the cohort-specific Cox models, sensitivity analyses, random-effects synthesis and risk estimates.
4. `crude_models.R` and `exact_risks.R` run supplementary checks used by the figure workflow.
5. `validate_results.py` performs structural and numerical consistency checks.
6. `draw_figures.py` exports the publication figures.

The private source project and its restricted data are deliberately excluded from this repository. `build_submission.py`, manuscript-generation helpers, local audit records, intermediate results and submission files are also excluded. The scripts should therefore be treated as a transparent code release for authorised reuse and review, not as a self-contained data package.
