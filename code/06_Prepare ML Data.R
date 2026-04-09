# =========================================================================
# Script Name: 06_Prepare ML Data.R
# Description: Prepare machine learning training dataset by filtering 
#              intersecting genes, transposing expression matrix, and 
#              merging with numerical clinical group labels.
# Note for GitHub users: Please update the working directory before running.
# =========================================================================

# Set working directory (Update this to your local repository path)
setwd("D:\\桌面\\sepsis-s100a8-multiomics\\data\\06_Prepare ML Data")

# Load required packages
library(openxlsx) 

cat("--- Starting Data Processing ---\n")

#################################################################
## 1. Load All Input Files
#################################################################

# 1.1 Load feature gene list (from "Overlap_Genes.txt")
# Assuming "Overlap_Genes.txt" is a single-column file without a header
genes_to_keep_df <- read.table("Overlap_Genes.txt", header = FALSE, stringsAsFactors = FALSE)
genes_to_keep <- genes_to_keep_df$V1
cat(paste("1. Successfully loaded", length(genes_to_keep), "feature genes.\n"))

# 1.2 Load raw expression matrix (GSE65682_gene.csv)
# Structure: Genes x Samples (Rownames are genes, Colnames are samples)
# row.names = 1: Specifies the first column as rownames (Gene symbols)
# check.names = FALSE: Prevents R from auto-modifying sample names (e.g., converting '-' to '.')
expr_data_raw <- read.csv("GSE65682_gene.csv", row.names = 1, check.names = FALSE)
cat(paste("2. Successfully loaded raw expression matrix:", nrow(expr_data_raw), "genes x", ncol(expr_data_raw), "samples.\n"))

# 1.3 Load sample group metadata (GSE65682分组.csv)
# Structure: Column 1 is 'sample', Column 2 is 'Group'
group_data <- read.csv("GSE65682分组.csv", stringsAsFactors = FALSE)
cat(paste("3. Successfully loaded group metadata:", nrow(group_data), "samples.\n"))

#################################################################
## 2. Data Cleaning and Reshaping
#################################################################

# 2.1 Security Check: Ensure 'sample' and 'Group' columns exist
# This prevents downstream failures due to incorrect CSV formatting
if (!all(c("sample", "Group") %in% colnames(group_data))) {
  stop("Error: The group metadata file must contain 'sample' and 'Group' as headers.")
}
cat("   (Check Passed: 'sample' and 'Group' columns verified.)\n")

# 2.2 Filter expression matrix: Retain only feature genes (Intersecting genes)
common_genes <- intersect(rownames(expr_data_raw), genes_to_keep)
expr_filtered <- expr_data_raw[common_genes, ]
cat(paste("4. Gene filtering complete: Retained", length(common_genes), "intersecting genes.\n"))

# 2.3 Transpose expression matrix (Result: Samples x Genes)
expr_transposed <- as.data.frame(t(expr_filtered))

# 2.4 Add rownames (Sample IDs) as a new column for subsequent merging
expr_transposed$SampleID_col <- rownames(expr_transposed)

# 2.5 Process group data: Convert labels to binary numeric values (1 and 0)
# This is mandatory for machine learning algorithms (e.g., glmnet, pROC)
# Mapping "control" to 1 (Positive/Target), and others (e.g., Normal) to 0
group_data$Group_Numeric <- ifelse(group_data$Group == "control", 1, 0)
cat("5. Converting group labels: control -> 1, Normal -> 0\n")
# Print contingency table of the conversion
print(table(Original_Group = group_data$Group, Numeric_Group = group_data$Group_Numeric))

#################################################################
## 3. Data Merging
#################################################################

# 3.1 Preparation for merging:
# Rename 'sample' column in group_data to 'SampleID_col' to match the transposed matrix
colnames(group_data)[colnames(group_data) == "sample"] <- "SampleID_col"

# 3.2 Merge (Inner Join):
# Retain only samples present in BOTH the expression matrix and the metadata
final_data_merged <- merge(
  expr_transposed,
  # Extract only the Sample ID and the new numeric group column
  group_data[, c("SampleID_col", "Group_Numeric")], 
  by = "SampleID_col"
)
cat(paste("6. Merging complete: Final dataset contains", nrow(final_data_merged), "common samples.\n"))

#################################################################
## 4. Formatting and Exporting
#################################################################

# 4.1 Reorder columns to meet Machine Learning input standards:
# [Col 1: SampleID], [Col 2 to N-1: Genes], [Last Col: Group]
gene_columns <- common_genes
sample_id_column <- "SampleID_col"
group_column <- "Group_Numeric"

# Reorder
final_data_formatted <- final_data_merged[, c(sample_id_column, gene_columns, group_column)]

# 4.2 Standardize column names (to match downstream ML scripts)
colnames(final_data_formatted)[colnames(final_data_formatted) == sample_id_column] <- "SampleID"
colnames(final_data_formatted)[colnames(final_data_formatted) == group_column] <- "Group"

# 4.3 Export to a new CSV file
# row.names = FALSE: Do not save R row numbers
# quote = FALSE: Do not add quotes to strings
output_filename <- "ML_training_data_prepared.csv"
write.csv(final_data_formatted, file = output_filename, row.names = FALSE, quote = FALSE)

cat("--- Data Formatting Completed ---\n")
cat(paste("File saved to:", file.path(getwd(), output_filename), "\n"))
cat("Final data dimensions:", nrow(final_data_formatted), "samples x", ncol(final_data_formatted), "columns\n")

# Print a preview of the structured data for inspection
cat("Data structure preview (First 5 rows, first 3 genes + Group column):\n")
print_cols <- c("SampleID", common_genes[1:3], "Group")
print(final_data_formatted[1:5, print_cols])