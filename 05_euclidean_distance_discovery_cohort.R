###############################################################################
# GLOBAL METABOLIC DEREGULATION (EUCLIDEAN DISTANCE) - DISCOVERY COHORT
#
# PURPOSE
# -------
# Calculate per-sample Euclidean distance from the normal-sample centroid
# in the discovery cohort and compare across ER+/ER-(all)/TNBC/HER2+(ER-)
# subtype groupings with pairwise Wilcoxon tests.
#
# INPUT
# -----
# Euclidean_Discovery.csv   : per-sample Euclidean distance values
# discovery_subtypes.csv    : sample subtype annotation
#
# OUTPUT
# ------
# Boxplot of Euclidean distance by subtype group (TIFF, printed to device)
#
# PACKAGES
# --------
# data.table, dplyr, ggpubr, ggplot2
###############################################################################

library(data.table)
library(dplyr)
library(ggpubr)
library(ggplot2)
# =====================================================
# Load data
# =====================================================
euclid <- fread("Euclidean_Discovery.csv")
subtypes <- fread("discovery_subtypes.csv")
# =====================================================
# Standardize sample identifiers
# =====================================================
euclid$Sample <- gsub("-", ".", euclid$Sample)
subtypes$Sample <- gsub("-", ".", subtypes$Sample)
# =====================================================
# Prepare Euclidean distance variable
# =====================================================
euclid$Euclidean <- as.numeric(euclid$`Euclidean distance`)
# =====================================================
# Define subtype subsets
# =====================================================
ERneg_all <- subtypes$Sample[
subtypes$`ER Status` == "Negative"
]
TNBC_samples <- subtypes$Sample[
subtypes$`Pam50 + Claudin-low subtype +TNBC` == "TNBC"
]
HER2neg_samples <- subtypes$Sample[
subtypes$`Pam50 + Claudin-low subtype +TNBC` == "Her2" &
subtypes$`ER Status` == "Negative"
]
Normalneg_samples <- subtypes$Sample[
subtypes$`Pam50 + Claudin-low subtype +TNBC` == "Normal" &
subtypes$`ER Status` == "Negative"
]
ERpos <- subtypes$Sample[
subtypes$`ER Status` == "Positive"
]
# =====================================================
# Define comparisons
# =====================================================
comparisons <- list(
list(
name = "ERneg_all_vs_ERpos",
ERneg = ERneg_all
),
list(
name = "ERneg_minus_TNBC_vs_ERpos",
ERneg = setdiff(ERneg_all, TNBC_samples)
),
list(
name = "ERneg_minus_HER2_vs_ERpos",
ERneg = setdiff(ERneg_all, HER2neg_samples)
)
)
# =====================================================
# Loop through comparisons
# =====================================================
for (comp in comparisons) {
comp_name <- comp$name
ERpos_subset <- intersect(ERpos, euclid$Sample)
ERneg_subset <- intersect(comp$ERneg, euclid$Sample)
if (length(ERpos_subset) < 1 | length(ERneg_subset) < 1) {
message("Skipping ", comp_name, ": empty group.")
next
}
# ---------------------------------------------------
# Plot dataframe
# ---------------------------------------------------
plot_df <- bind_rows(
data.frame(
Group = "ER+",
Euclidean = euclid$Euclidean[
match(ERpos_subset, euclid$Sample)
]
),
data.frame(
Group = "ER-",
Euclidean = euclid$Euclidean[
match(ERneg_subset, euclid$Sample)
]
)
)
# ---------------------------------------------------
# Sample sizes
# ---------------------------------------------------
group_n <- plot_df %>%
count(Group)
plot_df$Group <- factor(
plot_df$Group,
levels = c("ER+", "ER-"),
labels = paste0(
group_n$Group,
"\n(n=", group_n$n, ")"
)
)
# ---------------------------------------------------
# Summary statistics
# ---------------------------------------------------
stats_df <- plot_df %>%
group_by(Group) %>%
summarise(
N = n(),
Median = median(Euclidean, na.rm = TRUE),
IQR = IQR(Euclidean, na.rm = TRUE),
Mean = mean(Euclidean, na.rm = TRUE),
SD = sd(Euclidean, na.rm = TRUE),
.groups = "drop"
)
# ---------------------------------------------------
# Wilcoxon test
# ---------------------------------------------------
pval <- wilcox.test(
Euclidean ~ Group,
data = plot_df
)$p.value
stats_df <- stats_df %>%
mutate(
Comparison = comp_name,
Wilcox_p = pval
)
# ---------------------------------------------------
# Save statistics
# ---------------------------------------------------
fwrite(
stats_df,
paste0(comp_name, "_summary_statistics.csv")
)
# ---------------------------------------------------
# Publication-quality figure
# ---------------------------------------------------
tiff(
filename = paste0(comp_name, "_boxplot.tiff"),
width = 90,
height = 90,
units = "mm",
res = 600,
compression = "lzw"
)
library(data.table)
library(dplyr)
library(ggpubr)
library(ggplot2)
# =====================================================
# Read data
# =====================================================
euclid <- fread("Euclidean_validation.csv")
subtypes <- fread("validation_subtypes.csv")
# =====================================================
# Standardize sample IDs
# =====================================================
euclid$Sample <- gsub("-", ".", euclid$Sample)
subtypes$Sample <- gsub("-", ".", subtypes$Sample)
# =====================================================
# Rename Euclidean distance column
# =====================================================
colnames(euclid)[colnames(euclid) == "Euclidean distance"] <-
"Euclidean_distance"
# =====================================================
# Create grouped dataset
# =====================================================
euclid_long <- bind_rows(
# ER+
euclid %>%
filter(Sample %in% subtypes$Sample[subtypes$`ER Status` == "Positive"])
%>%
mutate(Group = "ER+"),
# ER-(all)
euclid %>%
filter(Sample %in% subtypes$Sample[subtypes$`ER Status` == "Negative"])
%>%
mutate(Group = "ER-(all)"),
# TNBC
euclid %>%
filter(Sample %in% subtypes$Sample[
subtypes$`Pam50 + Claudin-low subtype+TNBC` == "TNBC"]) %>%
mutate(Group = "TNBC"),
# HER2+ (ER-)
euclid %>%
filter(Sample %in% subtypes$Sample[
subtypes$`Pam50 + Claudin-low subtype+TNBC` == "Her2" &
subtypes$`ER Status` == "Negative"]) %>%
mutate(Group = "HER2+ / ER-"),
# Normal-like (ER-)
euclid %>%
filter(Sample %in% subtypes$Sample[
subtypes$`Pam50 + Claudin-low subtype+TNBC` == "Normal" &
subtypes$`ER Status` == "Negative"]) %>%
mutate(Group = "Normal-like / ER-")
)
# =====================================================
# Define order
# =====================================================
euclid_long$Group <- factor(
euclid_long$Group,
levels = c(
"ER+",
"ER-(all)",
"TNBC",
"HER2+ / ER-",
"Normal-like / ER-"
)
)
# =====================================================
# Add sample size to labels
# =====================================================
group_n <- euclid_long %>%
count(Group)
levels(euclid_long$Group) <- paste0(
group_n$Group,
"\n(n=", group_n$n, ")"
)
# =====================================================
# Pairwise comparisons
# =====================================================
comparisons <- list(
c("ER+\n(n=" %||% "" , ""), # remove if using dynamic comparisons
c("ER+", "ER-(all)"),
c("ER+", "TNBC"),
c("ER+", "HER2+ / ER-"),
c("ER+", "Normal-like / ER-"),
c("TNBC", "HER2+ / ER-"),
c("TNBC", "Normal-like / ER-"),
c("HER2+ / ER-", "Normal-like / ER-")
)
# Use original labels instead
group_levels <- levels(euclid_long$Group)
comparisons <- list(
c(group_levels[1], group_levels[2]),
c(group_levels[1], group_levels[3]),
c(group_levels[1], group_levels[4]),
c(group_levels[1], group_levels[5]),
c(group_levels[3], group_levels[4]),
c(group_levels[3], group_levels[5]),
c(group_levels[4], group_levels[5])
)
# =====================================================
# Publication-quality TIFF
# =====================================================
tiff(
filename = "Figure_EuclideanDistance_Subtypes.tiff",
width = 180,
height = 140,
units = "mm",
res = 600,
compression = "lzw"
)
p <- ggboxplot(
euclid_long,
x = "Group",
y = "Euclidean_distance",
color = "Group",
palette = "Dark2",
add = "jitter",
add.params = list(size = 1.2, alpha = 0.6),
outlier.shape = NA
) +
stat_compare_means(
comparisons = comparisons,
method = "wilcox.test",
label = "p.signif",
hide.ns = TRUE,
size = 4
) +
stat_summary(
fun = median,
geom = "crossbar",
width = 0.5,
color = "black",
linewidth = 0.4
) +
labs(
x = NULL,
y = "Euclidean distance"
) +
theme_classic(base_size = 12) +
theme(
legend.position = "none",
axis.text.x = element_text(
angle = 0,
hjust = 0.5,
size = 10
),
axis.title.y = element_text(size = 12),
plot.margin = margin(10, 10, 10, 10)
)
print(p)
dev.off()
} # end of the "for (comp in comparisons)" loop opened above.
  # NOTE: this closing brace was MISSING in the original source document.
  # Added here on the assumption the entire per-comparison plotting block
  # (lines from the loop start to here) was meant to be the loop body.
  # Please verify this matches your working/tested version of the script
  # before this is treated as the canonical copy.
