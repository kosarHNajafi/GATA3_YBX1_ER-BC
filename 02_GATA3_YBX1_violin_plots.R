###############################################################################
# GATA3 AND YBX1 GENE EXPRESSION IN ER+ AND TNBC CANCER EPITHELIAL CELLS
#
# PURPOSE
# -------
# Compare GATA3 and YBX1 expression between ER+ and TNBC cancer epithelial
# cells using violin+box plots with Wilcoxon rank-sum statistics.
#
# INPUT
# -----
# scRNA_checkpoint.rds
#   - list with:
#       $expr_mat : genes x cells expression matrix (dgCMatrix)
#       $meta_dt  : cell metadata data.table with columns
#                   NAME, celltype_major, subtype
#
# OUTPUT
# ------
# GATA3_expression_ERplus_vs_TNBC.tiff
# YBX1_expression_ERplus_vs_TNBC.tiff
#
# PACKAGES
# --------
# data.table, Matrix, ggplot2, ggpubr
###############################################################################

# GATA3 AND YBX1 GENE EXPRESSION IN ER+ AND TNBC CANCER EPITHELIAL CELLS
#
# PURPOSE
# -------
# Compare GATA3 and YBX1 expression between ER+ and TNBC cancer
epithelial
# cells using violin-box plots with Wilcoxon statistics.
###############################################################################
########################
# LOAD PACKAGES
########################
library(data.table)
library(Matrix)
library(ggplot2)
library(ggpubr)
###############################################################################
# LOAD SINGLE-CELL DATA
###############################################################################
sc_obj <- readRDS(
"path_to_data/scRNA_checkpoint.rds"
)
expr_mat <- sc_obj$expr_mat
cell_meta <- as.data.table(sc_obj$meta_dt)
###############################################################################
# SUBSET CANCER EPITHELIAL CELLS
###############################################################################
cancer_cells <- cell_meta[
celltype_major == "Cancer Epithelial",
NAME
]
mat_cancer <- expr_mat[
,
cancer_cells,
drop = FALSE
]
###############################################################################
# DEFINE ER+ AND TNBC CELLS
###############################################################################
cell_meta_cancer <- cell_meta[
NAME %in% cancer_cells
]
cell_meta_cancer <- cell_meta_cancer[
subtype %in% c("ER+", "TNBC")
]
###############################################################################
# EXTRACT GATA3 EXPRESSION
###############################################################################
gata3_exp <- data.table(
NAME = colnames(mat_cancer),
Expression = as.numeric(
mat_cancer["GATA3", ]
)
)
gata3_exp <- merge(
gata3_exp,
cell_meta_cancer[, .(NAME, subtype)],
by = "NAME"
)
###############################################################################
# EXTRACT YBX1 EXPRESSION
###############################################################################
ybx1_exp <- data.table(
NAME = colnames(mat_cancer),
Expression = as.numeric(
mat_cancer["YBX1", ]
)
)
ybx1_exp <- merge(
ybx1_exp,
cell_meta_cancer[, .(NAME, subtype)],
by = "NAME"
)
###############################################################################
# GATA3 VIOLIN PLOT
###############################################################################
p_gata3 <- ggplot(
gata3_exp,
aes(
x = subtype,
y = Expression,
fill = subtype
)
) +
geom_violin(
trim = FALSE
) +
geom_boxplot(
width = 0.12,
outlier.size = 0.4
) +
stat_compare_means(
method = "wilcox.test"
) +
labs(
x = "",
y = "GATA3 expression"
) +
theme_bw(base_size = 14) +
theme(
legend.position = "none"
)
print(p_gata3)
###############################################################################
# SAVE GATA3 PLOT
###############################################################################
tiff(
filename = "GATA3_expression_ERplus_vs_TNBC.tiff",
width = 1800,
height = 1600,
res = 300,
compression = "lzw"
)
print(p_gata3)
dev.off()
###############################################################################
# YBX1 VIOLIN PLOT
###############################################################################
p_ybx1 <- ggplot(
ybx1_exp,
aes(
x = subtype,
y = Expression,
fill = subtype
)
) +
geom_violin(
trim = FALSE
) +
geom_boxplot(
width = 0.12,
outlier.size = 0.4
) +
stat_compare_means(
method = "wilcox.test"
) +
labs(
x = "",
y = "YBX1 expression"
) +
theme_bw(base_size = 14) +
theme(
legend.position = "none"
)
print(p_ybx1)
###############################################################################
# SAVE YBX1 PLOT
###############################################################################
tiff(
filename = "YBX1_expression_ERplus_vs_TNBC.tiff",
width = 1800,
height = 1600,
res = 300,
compression = "lzw"
)
print(p_ybx1)
dev.off()
###############################################################################
# SUMMARY STATISTICS
###############################################################################
wilcox.test(
Expression ~ subtype,
data = gata3_exp
)
wilcox.test(
Expression ~ subtype,
data = ybx1_exp
)
