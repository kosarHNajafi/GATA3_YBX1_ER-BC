###############################################################################
# RTN TRANSCRIPTIONAL REGULATORY NETWORK RECONSTRUCTION (DISCOVERY COHORT)
#
# OBJECTIVE
# ---------
# Reconstruction of TF-centered regulatory networks using RTN for the
# discovery cohort (quality control, TF universe definition, permutation
# analysis, bootstrap analysis, and regulon extraction).
#
# INPUT
# -----
# Expression matrix and transcription factor universe (see script body for
# exact file names/objects expected).
#
# OUTPUT
# ------
# Console summary: genes retained, samples, TFs, regulons inferred
#
# PACKAGES
# --------
# RTN
###############################################################################

# RTN TRANSCRIPTIONAL REGULATORY NETWORK RECONSTRUCTION
#
# Objective:
# Reconstruction of TF-centered regulatory networks using RTN
#
# Required Input Files:
# expression_matrix.rds
# tf_list.csv
#
# Generated Output Files:
# tni_after_permutation.rds
# tni_after_bootstrap.rds
# tni_dpi0.rds
# regulon_sizes.csv
#
############################################################
suppressPackageStartupMessages({
library(RTN)
})
############################################################
# USER INPUTS
############################################################
project_dir <- "."
expression_file <- file.path(
project_dir,
"expression_matrix.rds"
)
tf_file <- file.path(
project_dir,
"tf_list.csv"
)
############################################################
# LOAD EXPRESSION MATRIX
############################################################
expr_df <- readRDS(expression_file)
expr_mat <- as.matrix(expr_df)
storage.mode(expr_mat) <- "double"
############################################################
# QUALITY CONTROL
############################################################
# Remove rows containing non-finite values
row_finite <- rowSums(!is.finite(expr_mat)) == 0
expr_mat <- expr_mat[
row_finite,
,
drop = FALSE
]
# Remove genes with zero variance
rv <- apply(
expr_mat,
1,
var,
na.rm = TRUE
)
expr_mat <- expr_mat[
is.finite(rv) & rv > 0,
,
drop = FALSE
]
############################################################
# TRANSCRIPTION FACTOR UNIVERSE
############################################################
tf_df <- read.csv(
tf_file,
check.names = FALSE
)
tf_list <- unique(tf_df$TF)
tf_list <- intersect(
tf_list,
rownames(expr_mat)
)
############################################################
# RTN CONSTRUCTOR
############################################################
rtni <- tni.constructor(
expData = expr_mat,
regulatoryElements = tf_list
)
############################################################
# PERMUTATION ANALYSIS
############################################################
set.seed(1)
rtni <- tni.permutation(
rtni,
estimator = "spearman",
nPermutations = 200,
pAdjustMethod = "BH",
pooledNullDistribution = TRUE,
parChunks = 20,
boxcox = TRUE,
verbose = TRUE
)
saveRDS(
rtni,
file.path(
project_dir,
"tni_after_permutation.rds"
)
)
############################################################
# BOOTSTRAP ANALYSIS
############################################################
rtni <- tni.bootstrap(
rtni,
nBootstraps = 50,
verbose = TRUE
)
saveRDS(
rtni,
file.path(
project_dir,
"tni_after_bootstrap.rds"
)
)
############################################################
# DPI FILTER
############################################################
rtni_dpi0 <- tni.dpi.filter(
rtni,
eps = 0,
verbose = TRUE
)
saveRDS(
rtni_dpi0,
file.path(
project_dir,
"tni_dpi0.rds"
)
)
############################################################
# EXTRACT REGULONS
############################################################
regs <- tni.get(
rtni_dpi0,
"regulons.and.mode"
)
reg_sizes <- sapply(
regs,
function(x) {
if (is.null(x)) {
return(0)
}
nrow(x)
}
)
reg_table <- data.frame(
TF = names(reg_sizes),
RegulonSize = reg_sizes,
stringsAsFactors = FALSE
)
reg_table <- reg_table[
order(
reg_table$RegulonSize,
decreasing = TRUE
),
]
write.csv(
reg_table,
file.path(
project_dir,
"regulon_sizes.csv"
),
row.names = FALSE
)
############################################################
# SUMMARY
############################################################
cat("\n===================================\n")
cat("RTN NETWORK RECONSTRUCTION SUMMARY\n")
cat("===================================\n")
cat("Genes retained:", nrow(expr_mat), "\n")
cat("Samples:", ncol(expr_mat), "\n")
cat("Transcription factors:", length(tf_list), "\n")
cat("Regulons inferred:", length(regs), "\n")
cat("===================================\n")
