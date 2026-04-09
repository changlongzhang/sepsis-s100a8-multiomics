# =========================================================================
# Script Name: 02_DEG.R
# Description: Perform differential expression analysis (Disease vs Normal) 
#              using DESeq2 and generate a customized volcano plot.
# Note for GitHub users: Please update the working directory before running.
# =========================================================================

# ====================== Environment Setup and Package Loading ======================
setwd("D:\\桌面\\sepsis-s100a8-multiomics\\data\\02_DEG") # Update this to your local repository path

# Dynamically check and install required packages
required_packages <- c("dplyr", "ggplot2", "DESeq2", "tibble", "stringr", "ggrepel")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ====================== Data Loading and Preprocessing ======================
cat("Loading gene expression counts data...\n")
expr_data <- read.csv("GSE65682_gene.csv", row.names = 1, header = TRUE)

# --- Fix for DESeq2: Ensure data contains non-negative integers ---
cat("Checking and processing non-integer expression values (rounding) to meet DESeq2 requirements...\n")
expr_data[expr_data < 0] <- 0
expr_data <- round(expr_data)
expr_data <- as.matrix(expr_data)
# --- End of Fix ---

cat("Loading sample metadata...\n")
location_data <- read.csv("GSE65682_Groups.csv", header = TRUE)
colnames(location_data) <- c("Sample", "Group")

# Filter metadata to retain only Disease and Normal groups
location_data <- location_data %>%
  filter(Group %in% c("Disease", "Normal")) %>%
  distinct(Sample, .keep_all = TRUE) 

# Filter for common samples in both expression matrix and metadata
common_samples <- intersect(colnames(expr_data), location_data$Sample)
expr_data <- expr_data[, common_samples]
location_data <- location_data %>% filter(Sample %in% common_samples)

# Reorder sample names to exactly match the columns of the expression matrix
location_data <- location_data %>%
  column_to_rownames("Sample") %>%
  # Ensure factor levels: Normal as reference (first), Disease as treatment (second)
  mutate(Group = factor(Group, levels = c("Normal", "Disease"))) %>%
  .[colnames(expr_data), , drop = FALSE] %>%
  rownames_to_column("Sample")

# Assign to analysis variables
tumor_expr <- expr_data 
tumor_meta <- location_data 

# Calculate group sizes and print summary
disease_count <- sum(tumor_meta$Group == "Disease")
normal_count <- sum(tumor_meta$Group == "Normal")
cat("Total samples:", ncol(tumor_expr), "\n")
cat("Disease group samples:", disease_count, "\n")
cat("Normal group samples:", normal_count, "\n")
if (disease_count < 3 | normal_count < 3) warning("Sample size is too small. Differential analysis results may be unreliable.")

# ====================== Disease vs Normal Differential Analysis ======================
cat("\nStarting Disease vs Normal differential analysis (Log2FC: Disease/Normal)...\n")

col_data_deseq <- tumor_meta %>%
  column_to_rownames("Sample")
# Ensure factor levels for DESeq2 design matrix
col_data_deseq$Group <- factor(col_data_deseq$Group, levels = c("Normal", "Disease")) 

dds <- DESeqDataSetFromMatrix(countData = tumor_expr,
                              colData   = col_data_deseq,
                              design    = ~ Group)
dds <- DESeq(dds)

# Extract results specifically for Disease vs Normal
res <- results(dds, contrast = c("Group", "Disease", "Normal"))

# **************** Format Results ****************
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("Gene") %>%
  mutate(
    Regulation = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Not significant"
    ),
    Significance = case_when(
      padj < 0.001 & abs(log2FoldChange) > 2   ~ "****",
      padj < 0.01  & abs(log2FoldChange) > 1.5 ~ "***",
      padj < 0.05  & abs(log2FoldChange) > 1   ~ "**",
      padj < 0.05                              ~ "*",
      TRUE                                     ~ "NS"
    )
  )

# Export all results
write.csv(res_df,
          "Disease_vs_Normal_all_results_padj.csv",
          row.names = FALSE)

# Filter and export significant DEGs
sig_res <- res_df %>% filter(Regulation != "Not significant")
write.csv(sig_res,
          "Disease_vs_Normal_sig_genes_padj.csv",
          row.names = FALSE)

cat("Number of significant DEGs (padj < 0.05 & |log2FC| > 1):", nrow(sig_res), "\n")
cat("Upregulated in Disease:", sum(sig_res$Regulation == "Up"), "\n")
cat("Downregulated in Disease:", sum(sig_res$Regulation == "Down"), "\n")

# ====================== Customized Volcano Plot ======================
library(ggrepel)
cat("Generating Volcano Plot...\n")

# 1. Extract the target gene (S100A8) for highlighting
target_gene_data <- res_df %>% 
  filter(Gene == "S100A8")

# 2. Build the plot
volcano_plot <-
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  
  # (A) Plot all background points
  geom_point(aes(color = Regulation), alpha = 0.6, size = 1.5) +
  
  # (B) Reference lines (Thresholds)
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  # (C) Custom color palette
  scale_color_manual(values = c(
    "Up"              = "#E31A1C", 
    "Down"            = "#1f78b4", 
    "Not significant" = "grey70"   
  )) +
  
  # (D) Highlight S100A8 with a prominent circle
  geom_point(data = target_gene_data, 
             shape = 21,          # Empty circle
             size = 4,            # Larger size
             stroke = 1.2,        # Thicker border
             color = "black",     
             fill = "transparent" 
  ) +
  
  # (E) Adjust label position for the target gene
  geom_text_repel(
    data = target_gene_data,
    aes(label = Gene),
    size = 6,                  
    fontface = "bold",         
    color = "black",           
    
    # Adjustments to avoid overlap
    box.padding = 0.8,         
    point.padding = 0.5,       
    nudge_x = -1,              # Nudge left
    nudge_y = -5,              # Nudge down to avoid upper margin
    
    min.segment.length = 0,    # Always show line
    segment.color = "black",   
    segment.size = 0.8         
  ) +
  
  # (F) Expand Y-axis limits to accommodate labels
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) + 
  
  # (G) Clean labels and remove titles
  labs(x = "Log2FC",
       y = "−log10(FDR)",
       title = NULL) +         
  
  theme_classic(base_size = 14) +
  theme(legend.title = element_blank(),
        plot.title = element_blank()) 

# Save the finalized plot
ggsave("Disease_vs_Normal_volcano_padj_S100A8_Fixed.pdf",
       volcano_plot, width = 8, height = 6)

cat("Optimized volcano plot saved successfully.\n")