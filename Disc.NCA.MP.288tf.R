#No need to biomart because the same gene annotation is used
#---Load libraries-----------------
library(RTN)
#---1.Load Gene Expression Data of 90 metabolic pathways and TFs = a mtrix for counts/assays---------------------------------
#load gene expression txt. Format
disc.mp.tf <- read.delim("C:/NCA.ER/NCA.METABRIC/Final/2.MP.TF/Disc.MP.288tf/Disc.MP.288tf.txt",header = TRUE,row.names = 1)
disc.mp.tf <- disc.mp.tf[order(rownames(disc.mp.tf)),order(colnames(disc.mp.tf))]
colnames(disc.mp.tf) <- gsub("\\.", "_", colnames(disc.mp.tf))
counts.disc.mp.tf <- as.matrix(disc.mp.tf)

str(counts.disc.mp.tf)
#num [1:1712, 1:993] 6.24 5.23 6.78 5.58 5.71 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:1584] "A4GALT" "A4GNT" "AACS" "AADAC" ...
#..$ : chr [1:993] "MB_0005" "MB_0006" "MB_0008" "MB_0014" ...

# Find elements that will turn into NA when coerced to numeric
non_numeric_elements <- counts.disc.mp.tf[is.na(as.numeric(counts.disc.mp.tf))]

length(non_numeric_elements)
#[1] 0

# Display non-numeric elements to understand what they are
print(non_numeric_elements)
#numeric(0)

#---2.Available METABRIC Clinical File as sampleAnnotation.disc =colData----------------------
METABRIC_Manual_Disc <- read.delim("~/NCA.ER/data/METABRIC Clinical.txt", row.names = 1)
METABRIC_Manual_Disc <- METABRIC_Manual_Disc[order(rownames(METABRIC_Manual_Disc)),]
View(METABRIC_Manual_Disc)

# Replace dots with underscores in column names if needed
rownames(METABRIC_Manual_Disc) <- gsub("\\-", "_", rownames(METABRIC_Manual_Disc))
View(METABRIC_Manual_Disc)

#---3.disc_sample_mp.tf_ids Sample/Column Annotation-----------------------------
# Load required library
library(dplyr)

# Assuming your dataset is named METABRIC_Manual_Disc
# Create a new binary dataframe based on transformations
METABRIC_Manual_Disc_binary <- METABRIC_Manual_Disc %>%
  transmute(
    IDs = rownames(METABRIC_Manual_Disc),
    Cohort = Cohort,
    
    OS.time = Overall.Survival..Months.,
    OS.event = ifelse(Overall.Survival.Status == "1:DECEASED", 1, 0),
    DSS.event = ifelse(Patient.s.Vital.Status == "Died of Disease", 1, 0),
    RFS.event = ifelse(Relapse.Free.Status == "1:Recurred", 1, 0),
    RFS.time = Relapse.Free.Status..Months.,
    
    Grade = Neoplasm.Histologic.Grade,
    Size = Tumor.Size,
    LN = Lymph.nodes.examined.positive,
    Age = Age.at.Diagnosis,
    LN = Lymph.nodes.examined.positive,
    Age = Age.at.Diagnosis,
    
    
    # Subtypes for LumA, LumB, Basal, Her2, Normal based on Pam50 subtype
    LumA = ifelse(Pam50...Claudin.low.subtype == "LumA", 1, 0),
    LumB = ifelse(Pam50...Claudin.low.subtype == "LumB", 1, 0),
    Basal = ifelse(Pam50...Claudin.low.subtype == "Basal", 1, 0),
    Her2 = ifelse(Pam50...Claudin.low.subtype == "Her2", 1, 0),
    Normal = ifelse(Pam50...Claudin.low.subtype == "Normal", 1, 0),
    
    # ER and PR status (positive and negative)
    `ER+` = ifelse(ER.Status == "Positive", 1, 0),
    `ER-` = ifelse(ER.Status == "Negative", 1, 0),
    
    
    # Histologic Grade categories G1, G2, G3
    G1 = ifelse(Neoplasm.Histologic.Grade == 1, 1, 0),
    G2 = ifelse(Neoplasm.Histologic.Grade == 2, 1, 0),
    G3 = ifelse(Neoplasm.Histologic.Grade == 3, 1, 0),
    
    # Hormone Therapy (HT)
    HT = ifelse(Hormone.Therapy == "YES", 1, 0),
    
  )

View(METABRIC_Manual_Disc_binary)

# Print head of the colAnnotation for verification
head(METABRIC_Manual_Disc_binary)

disc_sample_mp.tf_ids <- colnames(disc.mp.tf)
View(disc_sample_mp.tf_ids)
length(disc_sample_mp.tf_ids) #[1] 993

sort(disc_sample_mp.tf_ids)
View(disc_sample_mp.tf_ids)
length(disc_sample_mp.tf_ids)

all.equal(disc_sample_mp.tf_ids,colnames(counts.disc.mp.tf))

# Find the common sample IDs between the two datasets
common_samples_disc <- intersect(rownames(METABRIC_Manual_Disc), disc_sample_mp.tf_ids)
common_samples_disc <- sort(common_samples_disc)
View(common_samples_disc)
length(common_samples_disc) #[1] 988

# Subset and reorder both datasets to only include common samples
counts.disc.mp.tf <- counts.disc.mp.tf[, common_samples_disc]
all.equal(rownames(counts.disc.mp.tf), rownames(counts.disc.mp.tf)) #[1] TRUE
View(counts.disc.mp.tf)
dim(counts.disc.mp.tf)  #[1] 1712  988

sampleAnnotation.disc <- METABRIC_Manual_Disc_binary[common_samples_disc, ]
dim(sampleAnnotation.disc) #[1] 988  22

# Verify that they are aligned
all.equal(colnames(counts.disc.mp.tf), rownames(sampleAnnotation.disc))  # [1] TRUE

# If the first is not TRUE, you can match up the samples/columns in
# counts with the samples/rows in sampleAnnotation.disc like this (which is fine
# to run even if the first was TRUE):

#tempindex <- match(colnames(counts.disc.mp.tf), rownames(sample_metadata))
#sampleAnnotation.disc <- sample_metadata[tempindex, ]

#Check again

all.equal(colnames(counts.disc.mp.tf), rownames(sampleAnnotation.disc)) #[1] TRUE

print("Step 3 completed: Discovery SampleAnnotation Prepared")

#---4.Load Ensemble Gene Annotation----------------------------------------------
#Duplicates Removed Manually
gene_annot_mp.tf_disc <- read.delim("~/NCA.ER/Gene_annot_MP.TF_Used_EN113.txt")

dim(gene_annot_mp.tf_disc) #[1] 2912    7
length(gene_annot_mp.tf_disc) #[1] 7
View(gene_annot_mp.tf_disc)
dim(counts.disc.mp.tf) #[1] 1712  988

#Only Validated IDs Taken

# Extract gene IDs from both datasets
ensemble_ids <- gene_annot_mp.tf_disc$external_gene_name
View(ensemble_ids)

ensemble_ids <- sort(ensemble_ids)
length(ensemble_ids) #[1] 2912

disc_mp.tf_gene_ids <- rownames(counts.disc.mp.tf)
View(disc_mp.tf_gene_ids)

common_genes <- intersect(ensemble_ids, disc_mp.tf_gene_ids) # Should match or be close to 1420
common_genes <- sort(common_genes)
length(common_genes) #[1] 1505

# Subset and reorder both datasets to only include common samples
counts.disc.mp.tf <- counts.disc.mp.tf[common_genes, ]
View(counts.disc.mp.tf)

all.equal(common_genes,rownames(counts.disc.mp.tf)) #[1] TRUE

#Set rownames for gene annotation, %n% doesn't let your data go NA 
gene_annot_mp.tf_disc <- gene_annot_mp.tf_disc[gene_annot_mp.tf_disc$external_gene_name %in% common_genes, ]
rownames(gene_annot_mp.tf_disc) <- gene_annot_mp.tf_disc$external_gene_name
gene_annot_mp.tf_disc <- gene_annot_mp.tf_disc[order(rownames(gene_annot_mp.tf_disc)), ]
View(gene_annot_mp.tf_disc)
all.equal(gene_annot_mp.tf_disc$external_gene_name, rownames(counts.disc.mp.tf))
all.equal(gene_annot_mp.tf_disc$external_gene_name,rownames(gene_annot_mp.tf_disc))
all.equal(gene_annot_mp.tf_disc$external_gene_name,common_genes)
all.equal(common_genes,rownames(gene_annot_mp.tf_disc))

#Sort counts.disc.mp.tf and gene_annot_mp.tf_disc based on rownames which is gene names
counts.disc.mp.tf <- counts.disc.mp.tf[order(rownames(counts.disc.mp.tf)),order(colnames(counts.disc.mp.tf)) ]
gene_annot_mp.tf_disc <- gene_annot_mp.tf_disc[order(rownames(gene_annot_mp.tf_disc)), ]

all.equal(rownames(counts.disc.mp.tf), rownames(gene_annot_mp.tf_disc))  # Should return TRUE if aligned
dim(counts.disc.mp.tf) #[1] 1705  988
dim(gene_annot_mp.tf_disc) #[1] 1705    7
View(gene_annot_mp.tf_disc)

#In rtni analysis it needs "SYMBOL" in rowAnnotation
colnames(gene_annot_mp.tf_disc)
#[1] "ensembl_gene_id"    "external_gene_name" "chromosome_name"    "start_position"    
#5] "end_position"       "strand"             "description" 

colnames(gene_annot_mp.tf_disc)[colnames(gene_annot_mp.tf_disc) == "external_gene_name"] <- "SYMBOL"
colnames(gene_annot_mp.tf_disc)[colnames(gene_annot_mp.tf_disc) == "ensembl_gene_id"] <- "ENSEMBL"

# Verify the change
colnames(gene_annot_mp.tf_disc)
#[1] "ENSEMBL"         "SYMBOL"          "chromosome_name" "start_position" 
#[5] "end_position"    "strand"          "description"

# One final check:
stopifnot(rownames(gene_annot_mp.tf_disc) == rownames(counts.disc.mp.tf), # features
          rownames(sampleAnnotation.disc) == colnames(counts.disc.mp.tf)) # samples

save(list = ls(),file = "2.Annotations.RData")

#---5.Create tni_disc.mp.tf list-----------------------
counts.disc.mp.tf <- as.matrix(counts.disc.mp.tf)
View(counts.disc.mp.tf)
str(counts.disc.mp.tf)
#num [1:1705, 1:988] 6.24 5.23 6.78 5.58 5.71 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:1706] "A4GALT" "A4GNT" "AACS" "AADAC" ...
#..$ : chr [1:988] "MB_0005" "MB_0006" "MB_0008" "MB_0014" ...

#As tniData in fletcher2013:
tni_disc.mp.tf <- list(
  expData = counts.disc.mp.tf,                       # expData as an assay
  rowAnnotation = gene_annot_mp.tf_disc,             # Gene annotations (row metadata)
  colAnnotation = sampleAnnotation.disc           # Sample annotations (column metadata)
)

View(tni_disc.mp.tf)

save(list = ls(),file = "1.tniData.RData")

#---6.RegulatoryElements-----------------

# Load TF annotation
data("tfsData")

# Check TF annotation:
# Intersect TFs from Lambert et al. (2018) with gene annotation 
# from the gene expression of 90 metabolic pathway cohort
regulatoryElements <- intersect(tfsData$Lambert2018$SYMBOL, tni_disc.mp.tf$rowAnnotation$SYMBOL)
View(regulatoryElements)
regulatoryElements <- sort(regulatoryElements)
View(regulatoryElements)

#---7.Run the TNI constructor with the extracted matrix for tni_disc.mp.tf--------------
#This dataset consists of a list with 3 objects:
##a named gene expression matrix (tniData$expData),
##a data frame with gene annotations (tniData$rowAnnotation), 
##and a data frame with sample annotations (tniData$colAnnotation).
##alternatively, 'expData' can be a 'SummarizedExperiment' object
rtni_disc_mp.tf <- tni.constructor(expData = tni_disc.mp.tf$expData, 
                                    regulatoryElements = regulatoryElements, 
                                    rowAnnotation = tni_disc.mp.tf$rowAnnotation, 
                                    colAnnotation = tni_disc.mp.tf$colAnnotation)
#-Preprocessing for input data...
#--Mapping 'expData' to 'rowAnnotation'...
#--Mapping 'expData' to 'colAnnotation'...
#--Checking 'regulatoryElements' in 'rowAnnotation'...
#--Checking 'expData'...
#-Preprocessing complete!


all.equal(colnames(tni_disc.mp.tf$expData), rownames(sampleAnnotation.disc))  # Should return TRUE
all.equal(tni_disc.mp.tf$expData, counts.disc.mp.tf)
all.equal(tni_disc.mp.tf$colAnnotation, sampleAnnotation.disc)

all.equal(tni_disc.mp.tf$rowAnnotation,tni_disc.mp.tf$rowAnnotation)
#[1] TRUE
all.equal(tni_disc.mp.tf$colAnnotation,tni_disc.mp.tf$colAnnotation)
#[1] TRUE
all.equal(tni_disc.mp.tf$expData,tni_disc.mp.tf$expData)
#[1] TRUE

save.image("4.rtni.disc.RData")


#1.with pValuCutoff =1e-7 it gives 0 values, even in the summary
#rtni_disc_mp.tf.dpi <- tni.permutation(rtni_disc_mp.tf.dpi, pValueCutoff = 1e-7) #pValueCutoff=1e-7 ?? zero_values

#2.with pValuCutoff= 1e-7 and npermutation = 1000 it gives 0 values, even in the sumary
#rtni_disc_mp.tf.dpi <- tni.permutation(rtni_disc_mp.tf.dpi,nPermutations = 1000, pValueCutoff = 1e-7) #pValueCutoff=1e-7 ?? zero_values

#3.  nPermutations >= 1000 with running snow goes the same zero values
#4.  nPermutations >= 1000 without snow package, and pValueCutoff works!
rtni_disc_mp.tf <- tni.permutation(rtni_disc_mp.tf, 
                                 nPermutations = 1000,
                                 estimator = "spearman",
                                 verbose = TRUE) 

save(rtni_disc_mp.tf,file = "4_5.disc.rtni.afterpermutation.RData")
#Unstable interactions are subsequently removed by bootstrap analysis,
##creates a consensus bootstrap network, referred here as refnet (reference network).
rtni_disc_mp.tf <- tni.bootstrap(rtni_disc_mp.tf)

save(rtni_disc_mp.tf,file = "5.disc.rtni.afterbootstrap.RData")

# Compute the DPI-filtered regulatory network
rtni_disc_mp.tf.NA <- tni.dpi.filter(rtni_disc_mp.tf,eps = NA)
tni.regulon.summary(rtni_disc_mp.tf.NA)
#Regulatory network comprised of 301 regulons. 
#-- DPI-filtered network: 
#  regulatoryElements            Targets              Edges 
#301               1604               8291 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3.0    17.0    25.0    27.5    34.0   120.0 
#-- Reference network: 
#  regulatoryElements            Targets              Edges 
#301               1604             150305 
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#3     353     547     499     661     891 
#---
  
# Save the TNI object for subsequent analyses
save(list = ls(),file = "6.Disc.rtni.NA.RData")
save(rtni_disc_mp.tf.NA,file = "Disc.rtni.NA.forDEG.RData")
#detailed information about a specific regulon
tni.regulon.summary(rtni_disc_mp.tf.NA, regulatoryElements = "GATA3")

regulon.NA <- tni.get(rtni_disc_mp.tf.NA, what = 'regulons.and.mode', idkey = "SYMBOL")
View(regulon.NA)
head(regulon.NA)

# Find the maximum number of genes across all regulons
max_genes <- max(sapply(regulon.NA, function(x) if (is.null(x)) 0 else length(x)))

# Create a list of data frames with equal row lengths
regulon_list <- lapply(names(regulon.NA), function(regulon_name) {
  regulon_data <- regulon.NA[[regulon_name]]
  
  # Handle empty or NULL regulons
  if (is.null(regulon_data) || length(regulon_data) == 0) {
    df <- data.frame(
      Gene = rep(NA, max_genes),
      Value = rep(NA, max_genes),
      stringsAsFactors = FALSE
    )
  } else if (is.vector(regulon_data)) {
    # Ensure matching lengths for Gene and Value
    genes <- names(regulon_data)
    values <- regulon_data
    if (length(genes) == 0) genes <- rep(NA, length(values))
    df <- data.frame(
      Gene = genes,
      Value = values,
      stringsAsFactors = FALSE
    )
  } else if (is.matrix(regulon_data) || is.data.frame(regulon_data)) {
    df <- data.frame(
      Gene = rownames(regulon_data),
      Value = regulon_data[, 1], # Assuming values are in the first column
      stringsAsFactors = FALSE
    )
  } else {
    df <- data.frame(
      Gene = rep(NA, max_genes),
      Value = rep(NA, max_genes),
      stringsAsFactors = FALSE
    )
  }
  
  # Extend to max_genes rows if needed
  if (nrow(df) < max_genes) {
    df <- rbind(df, data.frame(
      Gene = rep(NA, max_genes - nrow(df)),
      Value = rep(NA, max_genes - nrow(df))
    ))
  }
  
  # Rename columns with regulon names
  colnames(df) <- c(paste0(regulon_name, "_Gene"), paste0(regulon_name, "_Value"))
  return(df)
})

# Combine all into one data frame
regulon_df <- do.call(cbind, regulon_list)

# Write to file
write.table(regulon_df, file = "Disc.regulon.NA.txt", sep = "\t", row.names = FALSE, quote = FALSE)
save(regulon.NA,file = "Disc.Regulon.NA.RData")
save.image("7.Disc.before.DEG.RData")
#---8.Linear DEG,Limma-----------

#Since the samples are the same sampleAnnotation.disc_tna will be used
Disc.MP.Genes <- read.delim("~/NCA.ER/data/DISC.MP.Genes.txt",header = TRUE, row.names = 1)
View(Disc.MP.Genes)

#Order by rownames and colnames
Disc.MP.Genes <- Disc.MP.Genes[order(rownames(Disc.MP.Genes)),order(colnames(Disc.MP.Genes))]

# Replace dots with underscores in column names if needed
colnames(Disc.MP.Genes) <- gsub("\\.", "_", colnames(Disc.MP.Genes))

View(Disc.MP.Genes)
dim(Disc.MP.Genes) #[1] 1420  993

counts.mp.disc <- as.matrix(Disc.MP.Genes[,common_samples_disc])
identical(colnames(counts.mp.disc),colnames(counts.disc.mp.tf))
str(counts.mp.disc)
#num [1:1420, 1:988] 6.24 5.23 6.78 5.58 5.71 ...
#- attr(*, "dimnames")=List of 2
#..$ : chr [1:1420] "A4GALT" "A4GNT" "AACS" "AADAC" ...
#..$ : chr [1:988] "MB_0005" "MB_0006" "MB_0008" "MB_0014" ...

sampleAnnotation.disc_tna_mp <- METABRIC_Manual_Disc[common_samples_disc, ]
dim(sampleAnnotation.disc_tna_mp) #[1] 988  36
View(sampleAnnotation.disc_tna_mp)


# Verify that they are aligned
all.equal(colnames(counts.mp.disc), rownames(sampleAnnotation.disc_tna_mp))  # Should return TRUE

# Check that row names of sample_annotations match column names of expression_data
if (!all(rownames(sampleAnnotation.disc_tna_mp) == colnames(counts.mp.disc))) {
  stop("Mismatch between sample annotation rownames and expression data colnames!")
}

# Load required package
library(limma)

# Create the design matrix for the linear model
# Assuming your label column is named "ER_status" with values "ERpos" and "ERneg"
design_mp_pos <- model.matrix(~ 0 + factor(sampleAnnotation.disc_tna_mp$ER.Status))
colnames(design_mp_pos) <- levels(factor(sampleAnnotation.disc_tna_mp$ER.Status))
rownames(design_mp_pos) <- rownames(sampleAnnotation.disc_tna_mp)
View(design_mp_pos)

all.equal(as.vector(design_mp_pos[,"Positive"]),as.vector(sampleAnnotation.disc$`ER+`))
all(design_mp_pos[,"Positive"] == sampleAnnotation.disc$`ER+`)
all(design_mp_pos[,"Negative"] == sampleAnnotation.disc$`ER-`)

# Fit the linear model using limma
fit_mp_pos <- lmFit(counts.mp.disc, design_mp_pos)
View(fit_mp_pos)

# Create contrast matrix to compare ERpos vs ERneg
contrast_matrix_mp_pos <- makeContrasts(Positive_vs_Negative = Positive - Negative, levels = design_mp_pos)

# Apply the contrast matrix
fit2_mp_pos <- contrasts.fit(fit_mp_pos, contrast_matrix_mp_pos)

# Empirical Bayes adjustment
fit2_mp_pos <- eBayes(fit2_mp_pos)

# Extract results (log2 fold changes, p-values, etc.)
phenotype_mp_pos <- topTable(fit2_mp_pos, coef = "Positive_vs_Negative", adjust.method = "BH", number = Inf)

# Order phenotype alphabetically by row names
phenotype_mp_pos <- phenotype_mp_pos[order(rownames(phenotype_mp_pos)), ]

# Save results to a file if needed
write.table(phenotype_mp_pos,file = "DEG.Disc.MP.txt",sep = "\t")
save(phenotype_mp_pos,file = paste0("Disc.DEG.",Sys.Date(),".RData"))

# Output the top results for inspection
View(phenotype_mp_pos)
dim(phenotype_mp_pos)  #[1] 1420    6

d.se.info <- sessionInfo()
save("regulon.NA",file = paste0("Disc.288.regulon.NA_",Sys.Date(),".RData"))
save(list = ls(),file = "8.Disc.MP.288tf.RData")
save.image("9.Disc.MP.288tf.RData")
save(list = c("rtni_disc_mp.tf.NA","regulon.NA","gsea1_disc.mp_h2","gsea1_disc.mp_p","gsea2_disc.mp_h2","gsea2_disc.mp_p","mra_disc.mp_h2","mra_disc.mp_p"), file = paste0("Network.Disc.MP.288TF_",Sys.Date(),".RData"))
write.table(regulatoryElements,file = paste0("Disc.294TF_RegulatoryElements_",Sys.Date(),".txt"),sep = "\t",row.names = FALSE)
