###############################################################################
# UPSET PLOT: PATHWAY OVERLAP BETWEEN ER-(all), TNBC, AND HER2+(ER-)
#
# PURPOSE
# -------
# Visualize shared/unique metabolic pathway membership across ER-(all),
# TNBC, and HER2+(ER-) subtype definitions using an UpSet plot.
#
# INPUT
# -----
# None (pathway set membership is defined in-script from pathway counts
# used elsewhere in the pipeline; see Methods for pathway list source).
#
# OUTPUT
# ------
# Figure_UpSet_ERneg_TNBC_HER2.tiff
#
# PACKAGES
# --------
# UpSetR
###############################################################################

library(UpSetR)
# =====================================================
# Define pathway sets
# =====================================================
ERneg_all <- paste0("Pathway_", 1:29)
# TNBC shares 28 of 29 pathways with ER-
TNBC <- ERneg_all[-29]
# HER2+ (ER-) shares 23 of 29 pathways with ER-
HER2_shared <- ERneg_all[1:23]
HER2_unique <- paste0(
"HER2_unique_pathway_",
1:8
)
HER2neg <- c(
HER2_shared,
HER2_unique
)
# =====================================================
# Create pathway list
# =====================================================
pathway_list <- list(
"ER-(all)" = ERneg_all,
"TNBC" = TNBC,
"HER2+ / ER-" = HER2neg
)
# =====================================================
# Export publication-quality TIFF
# =====================================================
tiff(
filename = "Figure_UpSet_ERneg_TNBC_HER2.tiff",
width = 120,
height = 90,
units = "mm",
res = 600,
compression = "lzw"
)
upset(
fromList(pathway_list),
order.by = "freq",
keep.order = TRUE,
main.bar.color = "#4B0082",
sets.bar.color = c(
"#0072B2", # blue
"#D55E00", # orange
"#009E73" # green
),
text.scale = c(
1.4, # intersection size title
1.2, # intersection size ticks
1.2, # set size title
1.0, # set size ticks
1.2, # set names
1.2 # numbers
),
mb.ratio = c(0.65, 0.35)
)
dev.off()
