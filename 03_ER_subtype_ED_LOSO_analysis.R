###############################################################################
# ER- SUBTYPE RESOLUTION, LEAVE-ONE-SUBTYPE-OUT (LOSO) SENSITIVITY, AND
# PATHWAY-SET UPSET/JACCARD ANALYSIS
#
# PURPOSE
# -------
# 1. Quantify global metabolic deregulation (Euclidean distance, ED) from a
#    normal-sample centroid, resolved by ER+/ER-(all)/TNBC/HER2+(ER-)/
#    Normal(ER-) subtype.
# 2. Test sensitivity of the ER-(all) grouping to exclusion of individual
#    subtypes (leave-one-subtype-out).
# 3. Quantify pathway-set overlap (Jaccard index, intersection counts)
#    between ER-(all), TNBC, and HER2+(ER-) pathway definitions.
#
# INPUT
# -----
# PDS_matrix.csv         : pathway deregulation score matrix (samples x pathways)
# sample_annotation.csv  : sample metadata with columns SampleID, Group
#                           (Group in {ER+, TNBC, HER2+(ER-), Normal(ER-), Normal})
#
# OUTPUT
# ------
# Boxplots of subtype-resolved and LOSO Euclidean distance (printed to device)
# Kruskal-Wallis / Dunn pairwise test results (printed to console)
# Jaccard index and pathway-set overlap counts (printed to console)
#
# PACKAGES
# --------
# tidyverse, ggpubr, FSA, rstatix, multcomp, UpSetR
###############################################################################

# ER- SUBTYPE RESOLUTION, LOSO SENSITIVITY AND UPSET ANALYSIS
############################################################
# ==========================================================
# 1. LOAD LIBRARIES
# ==========================================================
library(tidyverse)
library(ggpubr)
library(FSA)
library(rstatix)
library(multcomp)
library(UpSetR)
# ==========================================================
# 2. INPUT DATA
# ==========================================================
# PDS matrix
# Rows = samples
# Columns = metabolic pathways
pds <- read.csv(
"PDS_matrix.csv",
row.names = 1,
check.names = FALSE
)
# Sample annotation
metadata <- read.csv(
"sample_annotation.csv"
)
# Required columns:
# SampleID
# Group
#
# Example Groups:
# ER+
# TNBC
# HER2+(ER-)
# Normal(ER-)
# Normal
# ==========================================================
# 3. GLOBAL METABOLIC DEREGULATION (EUCLIDEAN DISTANCE)
# ==========================================================
# Determine normal centroid
normal.samples <- metadata$SampleID[
metadata$Group == "Normal"
]
normal.centroid <- colMeans(
pds[normal.samples, ],
na.rm = TRUE
)
# Calculate ED for each sample
ED <- apply(
pds,
1,
function(x)
sqrt(sum((x - normal.centroid)^2))
)
ed.df <- data.frame(
SampleID = names(ED),
ED = ED
) %>%
left_join(metadata, by = "SampleID")
# ==========================================================
# 4. SUBTYPE-RESOLVED ED ANALYSIS
# ==========================================================
ed.subtype <- ed.df %>%
filter(
Group %in% c(
"ER+",
"ER-(all)",
"TNBC",
"HER2+(ER-)",
"Normal(ER-)"
)
)
ggplot(
ed.subtype,
aes(Group, ED, fill = Group)
) +
geom_boxplot(outlier.shape = NA) +
geom_jitter(width = 0.15, alpha = 0.25) +
theme_classic() +
ylab("Euclidean Distance")
# Global test
kruskal.test(
ED ~ Group,
data = ed.subtype
)
# Pairwise comparisons
dunnTest(
ED ~ Group,
data = ed.subtype,
method = "bh"
)
# ==========================================================
# 5. LEAVE-ONE-SUBTYPE-OUT (LOSO) ED ANALYSIS
# ==========================================================
loso.ed <- bind_rows(
ed.df %>%
filter(Group == "ER+") %>%
mutate(Comparison = "ER+"),
ed.df %>%
filter(
Group %in%
c("TNBC",
"HER2+(ER-)",
"Normal(ER-)")
) %>%
mutate(Comparison = "ER-(all)"),
ed.df %>%
filter(
Group %in%
c("TNBC",
"Normal(ER-)")
) %>%
mutate(Comparison = "ER-(all)_minus_HER2"),
ed.df %>%
filter(
Group %in%
c("HER2+(ER-)",
"Normal(ER-)")
) %>%
mutate(Comparison = "ER-(all)_minus_TNBC")
)
ggplot(
loso.ed,
aes(Comparison, ED, fill = Comparison)
) +
geom_boxplot() +
theme_classic() +
ylab("Euclidean Distance")
# Dunnett comparison
fit.ed <- aov(
ED ~ Comparison,
data = loso.ed
)
summary(
glht(
fit.ed,
linfct = mcp(
Comparison = "Dunnett"
)
)
)
# ==========================================================
# 6. DISTANCE-TO-CENTROID ANALYSIS
# ==========================================================
calc.dist2centroid <- function(mat)
{
centroid <- colMeans(mat)
apply(
mat,
1,
function(x)
sqrt(sum((x - centroid)^2))
)
}
groups <- unique(metadata$Group)
disp.list <- list()
for(g in groups)
{
ids <- metadata$SampleID[
metadata$Group == g
]
tmp <- pds[ids, ]
d <- calc.dist2centroid(tmp)
disp.list[[g]] <-
data.frame(
SampleID = names(d),
DistanceToCentroid = d,
Group = g
)
}
disp.df <- bind_rows(disp.list)
# ==========================================================
# 7. LOSO DISTANCE-TO-CENTROID ANALYSIS
# ==========================================================
loso.disp <- bind_rows(
disp.df %>%
filter(Group == "ER+") %>%
mutate(Comparison = "ER+"),
disp.df %>%
filter(
Group %in%
c("TNBC",
"HER2+(ER-)",
"Normal(ER-)")
) %>%
mutate(Comparison = "ER-(all)"),
disp.df %>%
filter(
Group %in%
c("TNBC",
"Normal(ER-)")
) %>%
mutate(Comparison = "ER-(all)_minus_HER2"),
disp.df %>%
filter(
Group %in%
c("HER2+(ER-)",
"Normal(ER-)")
) %>%
mutate(Comparison = "ER-(all)_minus_TNBC")
)
ggplot(
loso.disp,
aes(
Comparison,
DistanceToCentroid,
fill = Comparison
)
) +
geom_boxplot() +
theme_classic() +
ylab("Distance to Centroid")
fit.disp <- aov(
DistanceToCentroid ~ Comparison,
data = loso.disp
)
summary(
glht(
fit.disp,
linfct = mcp(
Comparison = "Dunnett"
)
)
)
# ==========================================================
# 8. UPSET ANALYSIS
# ==========================================================
# Significant pathway lists
# FDR < 0.05 and Delta PDS >= 0.1
ERneg_all <- ERneg_pathways
TNBC <- TNBC_pathways
HER2neg <- HER2_pathways
pathway_list <- list(
"ER-(all)" = ERneg_all,
"TNBC" = TNBC,
"HER2+(ER-)" = HER2neg
)
upset(
fromList(pathway_list),
order.by = "freq",
keep.order = TRUE,
main.bar.color = "purple",
sets.bar.color = c(
"forestgreen",
"purple",
"orange"
),
text.scale = c(
1.3,
1.2,
1,
1,
1.2,
1.2
),
mb.ratio = c(
0.6,
0.4
)
)
# ==========================================================
# 9. JACCARD SIMILARITY
# ==========================================================
jaccard <- function(a, b)
{
length(intersect(a, b)) /
length(union(a, b))
}
J_TNBC <- jaccard(
ERneg_all,
TNBC
)
J_HER2 <- jaccard(
ERneg_all,
HER2neg
)
cat(
"Jaccard ER-(all) vs TNBC =",
round(J_TNBC, 3),
"\n"
)
cat(
"Jaccard ER-(all) vs HER2+(ER-) =",
round(J_HER2, 3),
"\n"
)
# ==========================================================
# 10. PATHWAY-SET OVERLAP COUNTS
# ==========================================================
ER_TNBC_overlap <-
length(
intersect(
ERneg_all,
TNBC
)
)
ER_HER2_overlap <-
length(
intersect(
ERneg_all,
HER2neg
)
)
cat(
"ER-(all) ∩ TNBC =",
ER_TNBC_overlap,
"\n"
)
cat(
"ER-(all) ∩ HER2+(ER-) =",
ER_HER2_overlap,
"\n"
)
