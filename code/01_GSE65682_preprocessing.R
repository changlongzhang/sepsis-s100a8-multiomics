# =========================================================================
# Script Name: 01_GSE65682_preprocessing.R
# Description: Download, format, and annotate the GSE65682 microarray dataset.
# Note for GitHub users: Please update the working directory path before running.
# =========================================================================

# Set working directory (Update this to your local path)
setwd("D:\\桌面\\sepsis-s100a8-multiomics\\data\\01_GSE65682_preprocessing")

library(GEOquery)
library(dplyr)
library(data.table)
library(tibble)

###############################################
### 1. Download GEO data (already log2 RMA, no normalization needed)
###############################################
gset <- getGEO("GSE65682", destdir = ".", AnnotGPL = FALSE, getGPL = FALSE)
expr.mat <- exprs(gset[[1]])   # log2 expression
expr.mat <- as.data.frame(expr.mat)

###############################################
### 2. (Optional) Convert to linear scale, do NOT use TPM
###############################################
linear <- 2^expr.mat    # Convert to linear scale (without normalization)
# If using limma later, log2 expression matrix is directly recommended.
# If for machine learning, linear scale is more reasonable (avoids negative log2 values).

###############################################
### 3. Read custom GPL annotation file
###############################################
gpl <- fread("GPL13667-15572.txt", sep = "\t", header = TRUE)

# Standardize column names
ids <- gpl[, c("ID", "Gene Symbol")]
colnames(ids) <- c("ProbeID", "GeneSymbol")

# Remove content after "//"
ids$GeneSymbol <- gsub("//.*", "", ids$GeneSymbol)

# Clean invalid gene symbols
ids$GeneSymbol[ids$GeneSymbol == "" | ids$GeneSymbol == "---"] <- NA

###############################################
### 4. Merge annotations (keep left_join to prevent probe loss)
###############################################
linear$ProbeID <- rownames(linear)
merged <- left_join(linear, ids, by = "ProbeID")

# Remove probes without gene symbols
merged <- merged %>% filter(!is.na(GeneSymbol))

###############################################
### 5. Aggregate multiple probes to gene (using max value)
###############################################
expr <- merged[, !(colnames(merged) %in% c("ProbeID"))]

gene.mat <- aggregate(. ~ GeneSymbol, data = expr, FUN = max, na.rm = TRUE)

rownames(gene.mat) <- gene.mat$GeneSymbol
gene.mat$GeneSymbol <- NULL

###############################################
### 6. Export final gene matrix
###############################################
write.csv(gene.mat, "GSE65682_gene_matrix.csv")

###############################################
### 7. Export clinical data
###############################################
clinical.all <- pData(gset[[1]])
write.csv(clinical.all, "GSE65682_clinical.csv", row.names = FALSE)


