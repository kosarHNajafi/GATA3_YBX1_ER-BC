###############################################################################
# BLISS INDEPENDENCE ANALYSIS
#
# PURPOSE
# -------
# Quantify interactions between regulator knockdown and drug treatments
# using the Bliss independence model, including pairwise and triple
# combination Bliss scores, bootstrap confidence intervals, and
# classification of interaction type (synergy/additivity/antagonism).
#
# INPUT
# -----
# Treatment/viability data (see script body for exact objects expected).
#
# OUTPUT
# ------
# Combined Bliss results table (printed / summarized)
# Bliss heatmap
#
# PACKAGES
# --------
# tidyverse, boot, reshape2
###############################################################################

# BLISS INDEPENDENCE ANALYSIS
#
# PURPOSE
# -------
# Quantify interactions between regulator knockdown and drug treatments
# using the Bliss independence model.
#
# APPLICATION
# -----------
# Used for:
#
# 1. shGATA3 + 2DG + DXR
# 2. shGATA3 + C75 + DXR
# 3. shYBX1 + 2DG + DXR
# 4. shYBX1 + C75 + DXR
# 5. shGATA3 + 2DG + Tam
# 6. shGATA3 + C75 + Tam
#
# Positive Bliss values:
# Synergy
#
# Negative Bliss values:
# Antagonism
#
###############################################################################
########################
# LOAD PACKAGES
########################
library(tidyverse)
library(boot)
library(reshape2)
###############################################################################
# LOAD DATA
###############################################################################
# Example input format:
#
# Condition,Rep1,Rep2,Rep3
# Scrambled,1.00,1.05,0.90
# Scrambled_2DG,0.91,0.92,1.02
# Scrambled_DXR,0.45,0.47,0.64
# Scrambled_2DG_DXR,0.55,0.44,0.52
# shGATA3,0.97,0.99,1.01
# shGATA3_2DG,0.65,0.68,0.76
# shGATA3_DXR,0.86,0.73,0.82
# shGATA3_2DG_DXR,0.60,0.70,0.66
dat <- read.csv(
"viability_data.csv",
check.names = FALSE
)
###############################################################################
# HELPER FUNCTIONS
###############################################################################
get_mean <- function(x){
mean(
as.numeric(x),
na.rm = TRUE
)
}
###############################################################################
# PAIRWISE BLISS
###############################################################################
# Expected = A × B
# Bliss = Expected − Observed
bliss_pair <- function(
effect_A,
effect_B,
observed_AB
){
expected_AB <- effect_A * effect_B
bliss <- expected_AB - observed_AB
return(bliss)
}
###############################################################################
# TRIPLE BLISS
###############################################################################
# Expected = A × B × C
# Bliss = Expected − Observed
bliss_triple <- function(
effect_A,
effect_B,
effect_C,
observed_ABC
){
expected_ABC <-
effect_A *
effect_B *
effect_C
bliss <- expected_ABC -
observed_ABC
return(bliss)
}
###############################################################################
# BOOTSTRAP CONFIDENCE INTERVALS
###############################################################################
bootstrap_triple_bliss <- function(
A,
B,
C,
ABC,
nboot = 10000
){
A <- as.numeric(A)
B <- as.numeric(B)
C <- as.numeric(C)
ABC <- as.numeric(ABC)
bliss_dist <- replicate(
nboot,
{
A_mean <- mean(
sample(
A,
replace = TRUE
)
)
B_mean <- mean(
sample(
B,
replace = TRUE
)
)
C_mean <- mean(
sample(
C,
replace = TRUE
)
)
Obs_mean <- mean(
sample(
ABC,
replace = TRUE
)
)
bliss_triple(
A_mean,
B_mean,
C_mean,
Obs_mean
)
}
)
quantile(
bliss_dist,
probs = c(
0.025,
0.50,
0.975
)
)
}
###############################################################################
# BLISS CLASSIFICATION
###############################################################################
classify_interaction <- function(
bliss,
ci_low,
ci_high
){
if(ci_low > 0){
if(bliss >= 0.10){
return(
"Strong synergy"
)
}
if(bliss >= 0.05){
return(
"Moderate synergy"
)
}
return(
"Weak synergy"
)
}
if(ci_high < 0){
return(
"Antagonism"
)
}
return(
"Additive/Borderline"
)
}
###############################################################################
# EXAMPLE:
# shGATA3 + 2DG + DXR
###############################################################################
shRNA <- c(
0.97,
0.99,
1.01
)
drugA <- c(
0.91,
0.92,
1.02
)
drugB <- c(
0.45,
0.47,
0.64
)
observed <- c(
0.60,
0.70,
0.66
)
bliss_val <- bliss_triple(
effect_A = mean(shRNA),
effect_B = mean(drugA),
effect_C = mean(drugB),
observed_ABC = mean(observed)
)
ci <- bootstrap_triple_bliss(
A = shRNA,
B = drugA,
C = drugB,
ABC = observed,
nboot = 10000
)
result <- data.frame(
Combination =
"shGATA3 + 2DG + DXR",
Bliss = bliss_val,
CI_low = ci[1],
CI_median = ci[2],
CI_high = ci[3],
Classification =
classify_interaction(
bliss_val,
ci[1],
ci[3]
)
)
print(result)
###############################################################################
# COMBINE ALL ANALYSES
###############################################################################
results <- bind_rows(
result
# Add additional
# conditions here
)
###############################################################################
# SAVE TABLE
###############################################################################
write.csv(
results,
"Bliss_results.csv",
row.names = FALSE
)
###############################################################################
# BLISS HEATMAP
###############################################################################
ggplot(
results,
aes(
x = Combination,
y = 1,
fill = Bliss
)
) +
geom_tile(
color = "white"
) +
geom_text(
aes(
label =
round(
Bliss,
3
)
),
size = 4
) +
scale_fill_gradient2(
low = "#2166AC",
mid = "white",
high = "#B2182B",
midpoint = 0
) +
theme_bw() +
theme(
axis.title =
element_blank(),
axis.text.y =
element_blank(),
axis.ticks =
element_blank()
)
ggsave(
filename =
"Bliss_heatmap.tiff",
width = 8,
height = 3,
dpi = 300,
compression = "lzw"
)
###############################################################################
# SUMMARY
###############################################################################
print(results)
