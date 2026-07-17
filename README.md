# Custom analysis code — ER+/TNBC metabolic regulation study

This repository contains the custom R scripts used to generate the figures
and statistical results reported in the manuscript. All scripts call
functions from publicly available CRAN and Bioconductor packages; no new
software tools or algorithms were developed.

## Contents

| Script | Description |
|---|---|
| `01_upset_pathway_overlap.R` | UpSet plot of metabolic pathway overlap between ER-(all), TNBC, and HER2+(ER-) |
| `02_GATA3_YBX1_violin_plots.R` | GATA3/YBX1 single-cell expression, ER+ vs TNBC cancer epithelial cells |
| `03_ER_subtype_ED_LOSO_analysis.R` | Subtype-resolved Euclidean distance, leave-one-subtype-out sensitivity, Jaccard/pathway-overlap analysis |
| `04_RTN_network_reconstruction_metabric.R` | RTN regulatory network reconstruction, METABRIC cohort |
| `05_euclidean_distance_discovery_cohort.R` | Euclidean distance analysis, discovery cohort |
| `06_RTN_network_reconstruction_discovery.R` | RTN regulatory network reconstruction, discovery cohort |
| `07_CNA_pathifier_covariate_regression.R` | Covariate-adjusted regression of Pathifier PDS scores (ER status, proliferation, hypoxia, CNA burden) |
| `08_RTN_regulons_AUCell_scoring.R` | RTN-derived GATA3/YBX1 metabolic regulons scored with AUCell in single-cell data |
| `09_bliss_independence_analysis.R` | Bliss independence analysis of regulator/drug combination effects |
| `10_ER.SurvivalCurves` | Kaplan-Meier survival curves stratified by ER status, covering four endpoints (OS, DSS, RFS, MFS) at both full follow-up and 5-year-truncated time horizons |
| `11_ERα.Pathway.SurvivalCurves` | Pathway-Stratified Survival Analysis  |
| `12_Step1_Paths_ER_Corr` | Pathway Correlation Analysis by ER Status (Discovery vs. Validation, FDR < 0.05)|
| `13_Step2_Paths_ER_Corr` | Opposite-Trend Pathway Correlations Between ER Statuses|
|`14_RTN` |  Gene Regulatory Network Reconstruction and Regulon Activity Profiling |
|`15_RTN_Regulon Count by Fisher test`|Regulon Enrichment for Metabolic Pathway Genes |
|`16_RTN_GSEA1&2` | Focused TF-Regulon Network and Master Regulator Analysis on Metabolic Genes |
|`17_BC_Regact_Corr_PDS` | TF Regulon Activity vs. Metabolic Pathway Correlation (Reproducible, |r| > 0.5) |
|`18_HR_BubblePlot`| Pathway Hazard Ratio Analysis and Cross-Cohort Bubble Plot (ER-stratified) |



Each script begins with a header describing its purpose, expected input
files, and output files.

## Requirements

- R (version used for the analysis: **[fill in — see `sessionInfo()`]**)
- CRAN packages: `data.table`, `Matrix`, `ggplot2`, `ggpubr`, `dplyr`,
  `tidyverse`, `broom`, `purrr`, `ggrepel`, `pheatmap`, `UpSetR`,
  `effsize`, `boot`, `reshape2`, `FSA`, `rstatix`, `multcomp`, `openxlsx`,
  `scales`
- Bioconductor packages: `RTN`, `AUCell`

Install via:

```r
install.packages(c(
  "data.table", "Matrix", "ggplot2", "ggpubr", "dplyr", "tidyverse",
  "broom", "purrr", "ggrepel", "pheatmap", "UpSetR", "effsize",
  "boot", "reshape2", "FSA", "rstatix", "multcomp", "openxlsx", "scales"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("RTN", "AUCell"))
```

## Data

Scripts expect the processed data objects/files described in each header
(e.g. `scRNA_checkpoint.rds`, `PDS_matrix.csv`, `discovery_pathifier.csv`).
Raw and processed data underlying this study are available as described in
the manuscript's Data Availability statement. File paths in these scripts
are placeholders (`path_to_data/...`) and should be edited to point to your
local copy of the data before running.

## Notes

- These scripts were used interactively during the analysis (console-run,
  one argument per line in places) rather than written as a packaged
  pipeline; they are shared here for transparency and reproducibility of
  the specific results reported in the manuscript, not as a general-purpose
  tool.
- `05_euclidean_distance_discovery_cohort.R`: a closing brace missing from
  the original working copy was added when preparing this repository (see
  in-script note at the relevant line). Please cross-check against your
  own tested version before treating this as canonical.

## Citation

If you use this code, please cite the associated manuscript [add citation
on acceptance] and this repository's archived release:
[Zenodo DOI — add after minting].
