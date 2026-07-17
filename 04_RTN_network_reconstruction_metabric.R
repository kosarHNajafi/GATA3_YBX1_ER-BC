###############################################################################
# RTN TRANSCRIPTIONAL REGULATORY NETWORK RECONSTRUCTION (METABRIC COHORT)
#
# DESCRIPTION
# -----------
# Reconstruction of transcription factor regulatory networks using RTN,
# including permutation analysis, bootstrap analysis, and DPI filtering,
# followed by regulon extraction and size summary.
#
# INPUT
# -----
# Expression matrix and transcription factor list (see script body for
# exact file names/objects expected).
#
# OUTPUT
# ------
# regulon_sizes.csv
# Console summary: genes retained, samples, TFs, regulons inferred
#
# PACKAGES
# --------
# RTN
###############################################################################

# RTN TRANSCRIPTIONAL REGULATORY NETWORK RECONSTRUCTION
#
# DESCRIPTION
# Reconstruction of transcription factor regulatory networks
# using RTN.
#
# INPUT FILES
# expression_matrix.rds
# Gene expression matrix
# Rows = genes
# Columns = samples
#
# tf_list.csv
# Single column named "TF"
# Transcription factors to evaluate
#
# OUTPUT FILES
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
# INPUT FILES
############################################################
expression_file <- "expression_matrix.rds"
tf_file <- "tf_list.csv"
############################################################
# LOAD EXPRESSION DATA
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
gene_variance <- apply(
expr_mat,
1,
var,
na.rm = TRUE
)
expr_mat <- expr_mat[
is.finite(gene_variance) &
gene_variance > 0,
,
drop = FALSE
]
############################################################
# LOAD TF LIST
############################################################
tf_df <- read.csv(
tf_file,
check.names = FALSE
)
tf_list <- unique(tf_df$TF)
# Keep TFs present in expression matrix
tf_list <- intersect(
tf_list,
rownames(expr_mat)
)
############################################################
# RTN CONSTRUCTION
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
file = "tni_after_permutation.rds"
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
file = "tni_after_bootstrap.rds"
)
############################################################
# DPI FILTERING
############################################################
rtni_dpi0 <- tni.dpi.filter(
rtni,
eps = 0,
verbose = TRUE
)
saveRDS(
rtni_dpi0,
file = "tni_dpi0.rds"
)
############################################################
# REGULON EXTRACTION
############################################################
regulons <- tni.get(
rtni_dpi0,
"regulons.and.mode"
)
regulon_sizes <- sapply(
regulons,
function(x) {
if (is.null(x)) {
return(0)
}
nrow(x)
}
)
regulon_table <- data.frame(
TF = names(regulon_sizes),
RegulonSize = regulon_sizes,
stringsAsFactors = FALSE
)
regulon_table <- regulon_table[
order(
regulon_table$RegulonSize,
decreasing = TRUE
),
]
write.csv(
regulon_table,
"regulon_sizes.csv",
row.names = FALSE
)
############################################################
# SUMMARY
############################################################
cat("\n")
cat("=====================================\n")
cat("RTN NETWORK RECONSTRUCTION SUMMARY\n")
cat("=====================================\n")
cat("Genes retained:", nrow(expr_mat), "\n")
cat("Samples:", ncol(expr_mat), "\n")
cat("Transcription factors:", length(tf_list), "\n")
cat("Regulons inferred:", length(regulons), "\n")
cat("=====================================\n")
