###############################################################################
# R2.1 FINAL COVARIATE-ADJUSTED ANALYSIS - METABRIC DISCOVERY COHORT
#
# MODEL
# -----
# PDS ~ ER_status + proliferation + hypoxia + CNA_burden
#
# PURPOSE
# -------
# Test whether pathway deregulation scores (Pathifier PDS) are associated
# with ER status after adjusting for proliferation, hypoxia, and CNA burden
# covariates, with BH multiple-testing correction across pathways.
#
# INPUT
# -----
# CNA_burden_METABRIC.csv    : per-sample copy-number alteration burden
# disc_ER_status.csv         : per-sample ER status
# Discovery_hallmarks.csv    : per-sample hallmark (proliferation/hypoxia) scores
# discovery_pathifier.csv    : pathway deregulation score matrix (pathways x samples)
#
# OUTPUT
# ------
# Covariate-adjusted regression results per pathway, with BH-adjusted P values
#
# PACKAGES
# --------
# data.table, dplyr, broom
###############################################################################

# R2.1 FINAL COVARIATE-ADJUSTED ANALYSIS
# METABRIC DISCOVERY COHORT
#
# Model:
# PDS ~ ER_status + proliferation + hypoxia + CNA_burden
#
# Proliferation = mean(E2F targets, G2M checkpoint)
#
# Discovery cohort defined by discovery_pathifier.csv
###############################################################
library(data.table)
library(dplyr)
library(broom)
###############################################################
# FILES
###############################################################
dir <- "path_to_data/"
cna <- fread(file.path(dir, "CNA_burden_METABRIC.csv"))
er <- fread(file.path(dir, "disc_ER_status.csv"))
hk <- fread(file.path(dir, "Discovery_hallmarks.csv"))
pfi <- fread(file.path(dir, "discovery_pathifier.csv"))
###############################################################
# DISCOVERY SAMPLE IDS
###############################################################
disc_ids <- colnames(pfi)[-1]
###############################################################
# KEEP DISCOVERY SAMPLES ONLY
###############################################################
cna <- cna[Sample %in% disc_ids]
###############################################################
# HALLMARK COVARIATES
###############################################################
setnames(
hk,
old = c("Hallmark pathways",
"E2F targets",
"G2M checkpoint",
"Hypoxia"),
new = c("Sample",
"E2F",
"G2M",
"Hypoxia")
)
hk[, Proliferation := rowMeans(
.SD,
na.rm = TRUE
),
.SDcols = c("E2F","G2M")]
###############################################################
# ER STATUS
###############################################################
setnames(er,
old = "ER Status",
new = "ER_status")
er[, ER_status := factor(
ER_status,
levels = c("Negative","Positive")
)]
###############################################################
# MERGE SAMPLE-LEVEL COVARIATES
###############################################################
meta <- merge(
er,
hk[, .(Sample,
Proliferation,
Hypoxia)],
by = "Sample",
all = FALSE
)
meta <- merge(
meta,
cna,
by = "Sample",
all = FALSE
)
###############################################################
# CHECK FINAL SAMPLE COUNT
###############################################################
cat("Final merged sample count:\n")
print(nrow(meta))
cat("\nMissing values per variable:\n")
print(colSums(is.na(meta)))
###############################################################
# CONVERT PATHIFIER TO LONG FORMAT
###############################################################
pfi_long <- melt(
pfi,
id.vars = "Pathway",
variable.name = "Sample",
value.name = "PDS"
)
###############################################################
# MERGE PDS WITH COVARIATES
###############################################################
analysis_df <- merge(
pfi_long,
meta,
by = "Sample",
all = FALSE
)
###############################################################
# CHECK
###############################################################
cat("\nRows in analysis table:\n")
print(nrow(analysis_df))
cat("\nNumber of pathways:\n")
print(length(unique(analysis_df$Pathway)))
###############################################################
# RUN COVARIATE-ADJUSTED MODELS
###############################################################
all_pathways <- unique(analysis_df$Pathway)
results <- lapply(all_pathways, function(pw){
dat <- analysis_df[Pathway == pw]
fit <- lm(
PDS ~ ER_status +
Proliferation +
Hypoxia +
CNA_burden,
data = dat
)
res <- broom::tidy(fit)
er_row <- res %>%
filter(term == "ER_statusPositive")
data.frame(
Pathway = pw,
Beta_ER = er_row$estimate,
SE = er_row$std.error,
t = er_row$statistic,
P = er_row$p.value
)
})
results <- rbindlist(results)
###############################################################
# MULTIPLE TESTING CORRECTION
###############################################################
results$BH_FDR <- p.adjust(
results$P,
method = "BH"
)
results <- results[
order(BH_FDR)
]
###############################################################
# ADD DIRECTION
###############################################################
results[, Direction :=
ifelse(
Beta_ER > 0,
"Higher in ER+",
"Higher in ER-"
)]
###############################################################
# SAVE RESULTS
###############################################################
fwrite(
results,
file.path(
dir,
"R2.1_covariate_adjusted_ER_effects.csv"
)
)
###############################################################
# SUMMARY
###############################################################
cat("\n====================================\n")
cat("Number significant at BH < 0.05:\n")
cat(sum(results$BH_FDR < 0.05), "\n")
cat("\nTop pathways:\n")
print(head(results, 20))
cat("\nResults saved to:\n")
cat(
file.path(
dir,
"R2.1_covariate_adjusted_ER_effects.csv"
),
"\n"
)
###############################################################
# OPTIONAL: OXIDATIVE VS GLYCOLYTIC SUMMARY
###############################################################
sig <- results[results$BH_FDR < 0.05]
cat("\nSignificant pathways:\n")
print(sig[, .(
Pathway,
Beta_ER,
BH_FDR,
Direction
)])
###########################################################################################################################
# Covariate-adjusted pathway deregulation analysis
# METABRIC Discovery Cohort (n = 988)
############################################################
library(tidyverse)
library(broom)
library(purrr)
library(ggplot2)
library(ggrepel)
library(pheatmap)
############################################################
# Load files
############################################################
pathifier_raw <- read.csv(
"discovery_pathifier.csv",
check.names = FALSE
)
er_df <- read.csv(
"disc_ER_status.csv",
check.names = FALSE
)
hall_df <- read.csv(
"Discovery_hallmarks.csv",
check.names = FALSE
)
cna_df <- read.csv(
"CNA_burden_METABRIC.csv",
check.names = FALSE
)
pam50_df <- read.csv(
"PAM50.csv",
check.names = FALSE
)
############################################################
# Standardize headers
############################################################
names(pathifier_raw) <- trimws(names(pathifier_raw))
names(er_df) <- trimws(names(er_df))
names(hall_df) <- trimws(names(hall_df))
names(cna_df) <- trimws(names(cna_df))
names(pam50_df) <- trimws(names(pam50_df))
############################################################
# Discovery sample IDs
############################################################
DISCOVERY_IDS <- colnames(pathifier_raw)[-1]
stopifnot(length(DISCOVERY_IDS) == 988)
############################################################
# Pathifier matrix
############################################################
pds_long <- pathifier_raw %>%
pivot_longer(
cols = -Pathway,
names_to = "Sample",
values_to = "PDS"
)
pds_wide <- pds_long %>%
pivot_wider(
names_from = Pathway,
values_from = PDS,
names_prefix = "PDS_"
)
############################################################
# ER status
############################################################
er_df2 <- er_df %>%
rename(
ER_Status = `ER Status`
) %>%
select(
Sample,
ER_Status
)
############################################################
# Proliferation and hypoxia
############################################################
hall_df2 <- hall_df %>%
rename(
Sample = `Hallmark pathways`,
E2F = `E2F targets`,
G2M = `G2M checkpoint`,
Hypoxia = Hypoxia
) %>%
mutate(
proliferation = (E2F + G2M) / 2,
hypoxia = Hypoxia
) %>%
select(
Sample,
proliferation,
hypoxia
)
############################################################
# CNA burden
############################################################
cna_disc <- cna_df %>%
filter(
Sample %in% DISCOVERY_IDS
) %>%
select(
Sample,
CNA_burden
)
############################################################
# PAM50-derived lineage score
############################################################
pam50_col <- names(pam50_df)[
grepl(
"^Pam50",
names(pam50_df),
ignore.case = TRUE
)
]
lineage_df <- pam50_df %>%
transmute(
Sample,
subtype = .data[[pam50_col]],
lineage_score = case_when(
subtype %in% c("LumA", "LumB") ~ 1,
subtype %in% c(
"Basal",
"TNBC",
"claudin-low"
) ~ -1,
TRUE ~ 0
)
)
############################################################
# Merge metadata
############################################################
meta <- reduce(
list(
er_df2,
hall_df2,
cna_disc,
lineage_df
),
full_join,
by = "Sample"
)
############################################################
# Merge PDS
############################################################
final_block3 <- meta %>%
inner_join(
pds_wide,
by = "Sample"
)
write.csv(
final_block3,
"BLOCK3_merged_discovery.csv",
row.names = FALSE
)
############################################################
# Regression models
############################################################
dat <- final_block3
dat$ER_Status <- factor(
dat$ER_Status,
levels = c(
"Negative",
"Positive"
)
)
pds_cols <- grep(
"^PDS_",
names(dat),
value = TRUE
)
run_model <- function(pathway) {
formula_text <- paste0(
"`",
pathway,
"` ~ ER_Status + proliferation + hypoxia + lineage_score + CNA_burden"
)
fit <- lm(
as.formula(formula_text),
data = dat
)
coef_tab <- summary(fit)$coefficients
tibble(
pathway = pathway,
ER_beta = coef_tab[
"ER_StatusPositive",
"Estimate"
],
ER_pvalue = coef_tab[
"ER_StatusPositive",
"Pr(>|t|)"
],
proliferation_beta = coef_tab[
"proliferation",
"Estimate"
],
hypoxia_beta = coef_tab[
"hypoxia",
"Estimate"
],
lineage_beta = coef_tab[
"lineage_score",
"Estimate"
],
CNA_beta = coef_tab[
"CNA_burden",
"Estimate"
]
)
}
results <- map_dfr(
pds_cols,
run_model
)
results <- results %>%
mutate(
ER_qvalue = p.adjust(
ER_pvalue,
method = "fdr"
)
)
write.csv(
results,
"BLOCK4_regression_summary.csv",
row.names = FALSE
)
############################################################
# Volcano plot
############################################################
key_pathways <- c(
"PDS_Fatty acids oxidation (mitochondrial)",
"PDS_Valine, leucine and isoleucine metabolism ",
"PDS_Glycolysis and gluconeogenesis",
"PDS_Fatty acid metabolism"
)
results$label <- ifelse(
results$pathway %in% key_pathways,
results$pathway,
NA
)
volcano_plot <- ggplot(
results,
aes(
x = ER_beta,
y = -log10(ER_pvalue),
color = ER_qvalue < 0.05
)
) +
geom_point(size = 3) +
geom_vline(
xintercept = 0,
linetype = "dashed"
) +
geom_text_repel(
aes(label = label),
size = 4
) +
scale_color_manual(
values = c(
"grey50",
"red"
)
) +
theme_minimal(base_size = 14)
tiff(
"BLOCK4_volcano.tiff",
width = 8,
height = 6,
units = "in",
res = 300,
compression = "lzw"
)
print(volcano_plot)
dev.off()
############################################################
# Covariate heatmap
############################################################
coef_matrix <- results %>%
select(
pathway,
ER_beta,
proliferation_beta,
hypoxia_beta,
lineage_beta,
CNA_beta
) %>%
column_to_rownames(
"pathway"
)
coef_scaled <- t(
scale(
t(
as.matrix(coef_matrix)
)
)
)
tiff(
"BLOCK4_heatmap.tiff",
width = 8,
height = 14,
units = "in",
res = 300,
compression = "lzw"
)
pheatmap(
coef_scaled,
fontsize = 9,
fontsize_row = 6,
fontsize_col = 10,
clustering_method = "ward.D2",
border_color = NA
)
dev.off()
