###############################################################################
# RTN-DERIVED METABOLIC REGULONS + SINGLE-CELL AUCELL ANALYSIS
#
# PURPOSE
# -------
# 1. Load RTN-derived GATA3 and YBX1 regulons generated from bulk network
#    reconstruction (scripts 04/06) and intersect with a curated metabolic
#    gene set to define core metabolic regulons.
# 2. Score single-cell activity of these regulons using AUCell in cancer
#    epithelial cells, split by ER+/TNBC subtype.
# 3. Compare AUCell scores between subtypes (violin plots + statistics).
#
# INPUT
# -----
# RTN regulon objects (from scripts 04/06)
# Metabolic gene set reference file
# Normalized single-cell expression data (scRNA_checkpoint.rds; see script 02)
#
# OUTPUT
# ------
# AUCell_GATA3_YBX1_MetabolicCoreActivity.csv
# Violin plots of AUCell scores by subtype
#
# PACKAGES
# --------
# data.table, Matrix, AUCell, ggplot2, ggpubr
###############################################################################

# RTN-DERIVED METABOLIC REGULONS + SINGLE-CELL AUCELL ANALYSIS
#
# PURPOSE
# -------
# 1. Load RTN-derived GATA3 and YBX1 regulons generated from bulk
transcriptomic
# data following MI, permutation testing, bootstrap filtering, and
ARACNe DPI.
# 2. Generate metabolic-only regulons by intersecting RTN targets with a
curated
# metabolic gene list.
# 3. Load normalized single-cell RNA-seq data.
# 4. Restrict analysis to cancer epithelial cells.
# 5. Perform differential expression analysis between ER+ and TNBC
cells.
# 6. Construct subtype-specific core metabolic regulons.
# 7. Compute regulon activity using AUCell.
# 8. Compare regulon activity between ER+ and TNBC cells.
# 9. Generate violin-box plots with Wilcoxon statistics.
#
# NOTE
# ----
# This workflow represents the final analysis used for single-cell
validation
# of RTN-derived GATA3 and YBX1 metabolic programs.
###############################################################################
########################
# LOAD PACKAGES
########################
library(data.table)
library(Matrix)
library(AUCell)
library(ggplot2)
library(ggpubr)
library(effsize)
###############################################################################
# PART 1
# GENERATE RTN-DERIVED METABOLIC REGULONS
###############################################################################
########################
# LOAD FULL RTN REGULONS
########################
gata3_rtn <- fread(
"path_to_regulons/RTN_edges_GATA3_eps0e00.csv"
)
ybx1_rtn <- fread(
"path_to_regulons/RTN_edges_YBX1_eps0e00.csv"
)
########################
# LOAD METABOLIC GENE SET
########################
metabolic_genes <- fread(
"path_to_annotations/Genes.csv"
)
metabolic_genes <- unique(
metabolic_genes[[1]]
)
########################
# INTERSECT WITH METABOLIC GENES
########################
gata3_metabolic <- gata3_rtn[
Target %in% metabolic_genes
]
ybx1_metabolic <- ybx1_rtn[
Target %in% metabolic_genes
]
########################
# SAVE METABOLIC REGULONS
########################
fwrite(
gata3_metabolic,
"RTN_edges_GATA3_METABOLIC_ONLY.csv"
)
fwrite(
ybx1_metabolic,
"RTN_edges_YBX1_METABOLIC_ONLY.csv"
)
########################
# CHECK REGULON SIZE
########################
length(unique(gata3_metabolic$Target))
length(unique(ybx1_metabolic$Target))
###############################################################################
# PART 2
# LOAD NORMALIZED SINGLE-CELL DATA
###############################################################################
sc_obj <- readRDS(
"path_to_scRNA_data/scRNA_checkpoint.rds"
)
expr_mat <- sc_obj$expr_mat
cell_meta <- as.data.table(sc_obj$meta_dt)
###############################################################################
# PART 3
# SUBSET CANCER EPITHELIAL CELLS
###############################################################################
cancer_cells <- cell_meta[
celltype_major == "Cancer Epithelial",
NAME
]
mat_cancer <- as(
expr_mat[, cancer_cells, drop = FALSE],
"dgCMatrix"
)
###############################################################################
# PART 4
# DEFINE ER+ AND TNBC CELLS
###############################################################################
er_cells <- intersect(
cell_meta[subtype == "ER+", NAME],
colnames(mat_cancer)
)
tnbc_cells <- intersect(
cell_meta[subtype == "TNBC", NAME],
colnames(mat_cancer)
)
###############################################################################
# PART 5
# DIFFERENTIAL EXPRESSION FUNCTIONS
###############################################################################
wilcox_one_gene <- function(g, mat, grp1, grp2) {
idx <- which(rownames(mat) == g)
if (length(idx) == 0) {
return(list(p = 1, med1 = 0, med2 = 0))
}
x1 <- as.numeric(mat[idx, grp1])
x2 <- as.numeric(mat[idx, grp2])
med1 <- median(x1)
med2 <- median(x2)
if (
(all(x1 == 0) && all(x2 == 0)) ||
(length(unique(x1)) < 2 &&
length(unique(x2)) < 2)
) {
return(list(
p = 1,
med1 = med1,
med2 = med2
))
}
p <- wilcox.test(
x1,
x2
)$p.value
list(
p = p,
med1 = med1,
med2 = med2
)
}
###############################################################################
de_on_targets <- function(genes,
mat,
grp1,
grp2) {
genes <- intersect(
genes,
rownames(mat)
)
res <- lapply(
genes,
function(g) {
d <- wilcox_one_gene(
g,
mat,
grp1,
grp2
)
list(
Gene = g,
p = d$p,
med_ER = d$med1,
med_TNBC = d$med2
)
}
)
dt <- rbindlist(res)
dt[, fdr := p.adjust(
p,
method = "BH"
)]
dt[, med_diff := med_ER - med_TNBC]
dt
}
###############################################################################
# PART 6
# DE ANALYSIS ON METABOLIC REGULON TARGETS
###############################################################################
de_gata3 <- de_on_targets(
unique(gata3_metabolic$Target),
mat_cancer,
er_cells,
tnbc_cells
)
de_ybx1 <- de_on_targets(
unique(ybx1_metabolic$Target),
mat_cancer,
er_cells,
tnbc_cells
)
###############################################################################
# PART 7
# BUILD CORE METABOLIC REGULONS
###############################################################################
N <- 100
########################
# GATA3
# ER+ ENRICHED TARGETS
########################
gata3_core <- de_gata3[
med_diff > 0
][order(-med_diff)][
1:N,
Gene
]
########################
# YBX1
# TNBC ENRICHED TARGETS
########################
ybx1_core <- de_ybx1[
med_diff < 0
][order(med_diff)][
1:N,
Gene
]
########################
# FINAL CORE REGULONS
########################
regulons <- list(
GATA3_core_metabolic = gata3_core,
YBX1_core_metabolic = ybx1_core
)
###############################################################################
# PART 8
# AUCELL RANKING CONSTRUCTION
###############################################################################
mat_ER <- mat_cancer[
,
er_cells,
drop = FALSE
]
mat_TNBC <- mat_cancer[
,
tnbc_cells,
drop = FALSE
]
rank_ER <- AUCell_buildRankings(
mat_ER,
plotStats = FALSE
)
rank_TNBC <- AUCell_buildRankings(
mat_TNBC,
plotStats = FALSE
)
###############################################################################
# PART 9
# CALCULATE AUCELL SCORES
###############################################################################
aucMaxRank <- 2500
auc_ER <- AUCell_calcAUC(
regulons,
rank_ER,
aucMaxRank = aucMaxRank
)
auc_TNBC <- AUCell_calcAUC(
regulons,
rank_TNBC,
aucMaxRank = aucMaxRank
)
###############################################################################
# PART 10
# COMBINE RESULTS
###############################################################################
auc_ER_dt <- as.data.table(
t(getAUC(auc_ER))
)
auc_TNBC_dt <- as.data.table(
t(getAUC(auc_TNBC))
)
auc_ER_dt[, subtype := "ER+"]
auc_TNBC_dt[, subtype := "TNBC"]
auc_final <- rbind(
auc_ER_dt,
auc_TNBC_dt,
fill = TRUE
)
###############################################################################
# PART 11
# STATISTICAL COMPARISON
###############################################################################
compare_regulon <- function(reg) {
er <- auc_final[
subtype == "ER+",
get(reg)
]
tn <- auc_final[
subtype == "TNBC",
get(reg)
]
list(
p_value = wilcox.test(
er,
tn
)$p.value,
cliffs_delta = cliff.delta(
er,
tn
)$estimate,
median_ER = median(er),
median_TNBC = median(tn)
)
}
########################
# RESULTS
########################
compare_regulon(
"GATA3_core_metabolic"
)
compare_regulon(
"YBX1_core_metabolic"
)
###############################################################################
# PART 12
# VIOLIN PLOTS
###############################################################################
########################
# GATA3
########################
ggplot(
auc_final,
aes(
x = subtype,
y = GATA3_core_metabolic,
fill = subtype
)
) +
geom_violin(
trim = FALSE
) +
geom_boxplot(
width = 0.12,
outlier.size = 0.5
) +
stat_compare_means(
method = "wilcox.test"
) +
labs(
x = "",
y = "Regulon activity"
) +
theme_bw()
###############################################################################
# YBX1
###############################################################################
ggplot(
auc_final,
aes(
x = subtype,
y = YBX1_core_metabolic,
fill = subtype
)
) +
geom_violin(
trim = FALSE
) +
geom_boxplot(
width = 0.12,
outlier.size = 0.5
) +
stat_compare_means(
method = "wilcox.test"
) +
labs(
x = "",
y = "Regulon activity"
) +
theme_bw()
###############################################################################
# SAVE RESULTS
###############################################################################
fwrite(
auc_final,
"AUCell_GATA3_YBX1_MetabolicCoreActivity.csv"
)
###############################################################################
# END OF WORKFLOW
