# =========================================================================
# Script Name: 03_WGCNA_Analysis.R
# Description: Weighted Gene Co-expression Network Analysis (WGCNA) for 
#              GSE65682 to identify Disease-associated gene modules.
# Note for GitHub users: Please update the working directory before running.
# =========================================================================

# ====================== Environment Setup and Package Loading ======================
# Set working directory (Update this to your local repository path)
# setwd("D:\\桌面\\sepsis-s100a8-multiomics\data\\03_WGCNA_analysis")

# Current working directory as output directory
output_dir <- getwd()

# Load core analysis packages
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(limma)
  library(WGCNA)
  library(flashClust)
  library(dynamicTreeCut)
  library(RColorBrewer)
  library(tibble)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(factoextra)      # For PCA visualization
  library(ggrepel)         # For PCA labels
  library(ggpubr)          # For boxplot significance annotations
})

# Enable multi-threading for WGCNA (if available)
enableWGCNAThreads()
allowWGCNAThreads()

# Define group constants (Strictly matching the CSV metadata)
GROUP1_NAME <- "Normal"
GROUP2_NAME <- "Disease"

# ====================== Data Loading and Preprocessing ======================
# ==== 1. Load expression data ====
cat("Loading expression data...\n")
# Assuming the first column is Gene ID and data is pre-processed (e.g., log2 or normalized)
expr_data <- read.csv("GSE65682_gene.csv", row.names = 1, header = TRUE)
expr_data <- as.matrix(expr_data)

# ==== 2. Load metadata (Exact matching) ====
cat("Loading group information...\n")
group_data <- read.csv("GSE65682_Groups.csv", header = TRUE)
colnames(group_data) <- c("Sample", "Group") # Ensure columns are named Sample and Group

# Filter valid groups
group_data <- group_data %>%
  # Remove potential leading/trailing spaces
  mutate(Group = trimws(Group)) %>% 
  # Strictly filter rows matching GROUP1_NAME or GROUP2_NAME
  filter(Group %in% c(GROUP1_NAME, GROUP2_NAME)) %>% 
  # Remove duplicate samples
  distinct(Sample, .keep_all = TRUE)

# Validation check
if (nrow(group_data) == 0) {
  stop("Error: No samples remaining after filtering! Please check if group names exactly match 'Normal' and 'Disease' (case-sensitive).")
}

# ==== 3. Intersect common samples ====
cat("Filtering for common samples...\n")
common_samples <- intersect(colnames(expr_data), group_data$Sample)
expr_data <- expr_data[, common_samples]
group_data <- group_data %>% filter(Sample %in% common_samples)

if (ncol(expr_data) == 0) stop("Error: No common samples found between expression data and metadata!")
cat("Valid samples:", ncol(expr_data), "\n")
cat(GROUP1_NAME, "samples:", sum(group_data$Group == GROUP1_NAME), "\n")
cat(GROUP2_NAME, "samples:", sum(group_data$Group == GROUP2_NAME), "\n")

# ==== 4. Data Filtering ====
cat("Filtering low-expressed genes...\n")
# Retain genes expressed (>0) in at least 25% of the samples
keep_genes <- rowSums(expr_data > 0) >= 0.25 * ncol(expr_data)
expr_data <- expr_data[keep_genes, ]
cat("Retained genes:", nrow(expr_data), "\n")

# ==== 5. Data Standardization ====
cat("Data standardization (Z-score)...\n")
# Z-score standardization across samples per gene
expr_data_std <- t(scale(t(expr_data)))

# ==== 6. Sample Outlier Detection (Hierarchical Clustering) ====
cat("Sample outlier detection (Clustering)...\n")
sample_tree <- flashClust(dist(t(expr_data_std)), method = "average")
pdf(file.path(output_dir, "sample_clustering_before_filtering.pdf"), width = 24, height = 18)
plot(sample_tree, main = "", cex = 0.7)
dev.off()

# ==== 7. PCA for Outlier Detection ====
cat("Performing PCA for outlier detection...\n")

# Execute PCA
pca_result <- prcomp(t(expr_data_std), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$Sample <- rownames(pca_df)
pca_df <- pca_df %>%
  left_join(group_data %>% dplyr::select(Sample, Group), by = "Sample")

# Calculate Mahalanobis distance to center
if (nrow(pca_df) > 2) {
  mahalanobis_dist <- mahalanobis(pca_df[, c("PC1", "PC2")],
                                  colMeans(pca_df[, c("PC1", "PC2")]),
                                  cov(pca_df[, c("PC1", "PC2")]))
  
  # Identify outliers (distance > 3 standard deviations)
  outlier_threshold <- mean(mahalanobis_dist) + 3 * sd(mahalanobis_dist)
  outliers <- pca_df$Sample[mahalanobis_dist > outlier_threshold]
} else {
  cat("Insufficient sample size, skipping PCA outlier detection.\n")
  outliers <- character(0)
}

# Plot PCA (Highlighting outliers)
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("Normal" = "#377EB8", "Disease" = "#E41A1C")) +
  labs(
    title = "",
    subtitle = "",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2,2]*100, 1), "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    legend.position = "right"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)

if (length(outliers) > 0) {
  pca_plot <- pca_plot +
    geom_label_repel(
      data = subset(pca_df, Sample %in% outliers),
      aes(label = Sample),
      box.padding = 0.5,
      show.legend = FALSE
    )
}

ggsave(file.path(output_dir, "pca_before_outlier_removal.pdf"), pca_plot, width = 10, height = 8)

# Remove outliers if detected
if (length(outliers) > 0) {
  cat("Outlier(s) detected:", paste(outliers, collapse = ", "), "\n")
  expr_data_std <- expr_data_std[, !colnames(expr_data_std) %in% outliers]
  group_data <- group_data %>% filter(!Sample %in% outliers)
  
  # Update PCA after cleaning
  pca_result_clean <- prcomp(t(expr_data_std), center = TRUE, scale. = TRUE)
  pca_df_clean <- as.data.frame(pca_result_clean$x[, 1:2])
  pca_df_clean$Sample <- rownames(pca_df_clean)
  pca_df_clean <- pca_df_clean %>%
    left_join(group_data %>% dplyr::select(Sample, Group), by = "Sample")
  
  pca_plot_clean <- ggplot(pca_df_clean, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Normal" = "#377EB8", "Disease" = "#E41A1C")) +
    labs(
      title = "",
      subtitle = "",
      x = paste0("PC1 (", round(summary(pca_result_clean)$importance[2,1]*100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result_clean)$importance[2,2]*100, 1), "%)")
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "right"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)
  
  ggsave(file.path(output_dir, "pca_after_outlier_removal.pdf"), pca_plot_clean, width = 10, height = 8)
} else {
  cat("No outlier samples detected.\n")
}

# Update sample info
common_samples_clean <- colnames(expr_data_std)
group_data <- group_data %>% filter(Sample %in% common_samples_clean)
cat("Final sample count:", ncol(expr_data_std), "\n")
cat(GROUP1_NAME, "samples:", sum(group_data$Group == GROUP1_NAME), "\n")
cat(GROUP2_NAME, "samples:", sum(group_data$Group == GROUP2_NAME), "\n")

# ==== 8. Prepare final data for WGCNA ====
expr_data <- expr_data_std

# Plot sample clustering (Final Data)
sample_tree_final <- flashClust(dist(t(expr_data)), method = "average")
pdf(file.path(output_dir, "sample_clustering_final.pdf"), width = 24, height = 18)
plot(sample_tree_final,
     main = "",
     cex = 0.7,
     labels = group_data$Sample)
dev.off()

# Save preprocessed data
save(expr_data, group_data, file = file.path(output_dir, "preprocessed_data.RData"))
cat("Preprocessed data saved as: preprocessed_data.RData\n")

# ====================== WGCNA Network Construction ======================
cat("\n===== Starting Co-expression Network Construction =====\n")

# ==== 1. Soft threshold selection ====
powers <- c(1:20)
sft <- pickSoftThreshold(t(expr_data), powerVector = powers, networkType = "signed", verbose = 5)

# Automatically select soft threshold
softPower <- if (!is.na(sft$powerEstimate)) {
  sft$powerEstimate
} else {
  # Manual selection: First power yielding R^2 > 0.9
  r2_index <- which(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] > 0.9)
  if (length(r2_index) > 0) sft$fitIndices[r2_index[1], 1] else 6
}
cat("Selected soft threshold:", softPower, "\n")

# Soft threshold plots
pdf(file.path(output_dir, "soft_threshold_selection.pdf"), width = 8, height = 6)
par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = "red")
dev.off()

# ==== 2. Construct Co-expression Network ====
cat("Constructing co-expression network...\n")
adjacency <- adjacency(t(expr_data), power = softPower, type = "signed")
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# Gene clustering
geneTree <- flashClust(as.dist(dissTOM), method = "average")

# ==== 3. Module Identification ====
cat("Identifying gene modules...\n")
minModuleSize <- 200
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize,
  verbose = 0
)

# Convert module assignment to colors
dynamicColors <- labels2colors(dynamicMods)

# Module dendrogram
pdf(file.path(output_dir, "gene_dendrogram_modules.pdf"), width = 8, height = 6)
plotDendroAndColors(
  geneTree, dynamicColors, "Dynamic Tree Cut",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = ""
)
dev.off()

# ==== 4. Merge Similar Modules ====
MEList <- moduleEigengenes(t(expr_data), colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- flashClust(as.dist(MEDiss), method = "average")

# Merge threshold
MEDissThres <- 0.25 # Merge modules with correlation > 0.75
merge <- mergeCloseModules(t(expr_data), dynamicColors, cutHeight = MEDissThres)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# Merged module dendrogram
pdf(file.path(output_dir, "merged_modules.pdf"), width = 8, height = 6)
plotDendroAndColors(
  geneTree, cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = ""
)
dev.off()

# Final module assignment
moduleColors <- mergedColors
moduleLabels <- match(moduleColors, standardColors(50)) - 1

# Extract and save genes per module
module_genes <- list()
for (color in unique(moduleColors)) {
  if (color != "grey") {
    module_genes[[color]] <- rownames(expr_data)[moduleColors == color]
    writeLines(module_genes[[color]],
               file.path(output_dir, paste0("module_", color, "_genes.txt")))
  }
}
cat("Identified", length(module_genes), "valid modules.\n")


# ====================== Module-Trait Association ======================
cat("\n===== Analyzing", GROUP1_NAME, "vs", GROUP2_NAME, "Differential Modules =====\n")

# Prepare phenotype/group data
group_info <- group_data %>%
  dplyr::select(Sample, Group) %>%
  mutate(Group = factor(Group, levels = c(GROUP1_NAME, GROUP2_NAME)),
         GroupNum = ifelse(Group == GROUP2_NAME, 1, 0)) %>% # Disease=1, Normal=0
  column_to_rownames("Sample")

# Ensure matching order
mergedMEs <- mergedMEs[rownames(group_info), ]

# ==== Calculate Module-Trait Correlation ====
# Positive correlation indicates enrichment in Disease.
module_group_cor <- as.data.frame(cor(mergedMEs, group_info$GroupNum, use = "p"))
module_group_p <- as.data.frame(corPvalueStudent(as.matrix(module_group_cor), nrow(group_info)))

# Create results dataframe
module_diff_df <- data.frame(
  Module = rownames(module_group_cor),
  Color = gsub("^ME", "", rownames(module_group_cor)),
  Correlation = module_group_cor[, 1],
  Pvalue = module_group_p[, 1],
  row.names = rownames(module_group_cor)
) %>%
  mutate(
    FDR = p.adjust(Pvalue, method = "fdr"),
    Direction = ifelse(Correlation > 0, paste0(GROUP2_NAME, "-enriched"), paste0(GROUP1_NAME, "-enriched")),
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  ) %>%
  arrange(desc(abs(Correlation)))

# Save association results
write.csv(module_diff_df, file.path(output_dir, "module_group_association.csv"), row.names = FALSE)


# ==== Prepare correlation matrix for Heatmap ====
# Note: Disease correlation remains original, Normal is inverted (multiplied by -1)
module_group_cor_for_heatmap <- cbind(-module_diff_df$Correlation, module_diff_df$Correlation)
rownames(module_group_cor_for_heatmap) <- module_diff_df$Module
colnames(module_group_cor_for_heatmap) <- c(GROUP1_NAME, GROUP2_NAME)


# ==== Plot Heatmap (Significant modules only) ====
sig_mod_idx <- which(module_diff_df$FDR < 0.05)

if(length(sig_mod_idx) > 0) {
  sig_mod_names <- module_diff_df$Module[sig_mod_idx]
  sig_colors <- gsub("^ME", "", sig_mod_names)
  
  # Extract correlations
  sig_cor_matrix <- module_group_cor_for_heatmap[sig_mod_names, , drop = FALSE]
  
  # Construct text matrix (Correlation + Significance)
  sig_text <- matrix(
    sprintf("%.2f\n%s", sig_cor_matrix, module_diff_df$Significance[sig_mod_idx]),
    nrow = nrow(sig_cor_matrix),
    ncol = ncol(sig_cor_matrix),
    byrow = FALSE
  )
  
  # Dynamic height based on module count
  pdf_height <- max(6, 0.4 * nrow(sig_cor_matrix))
  
  pdf(file.path(output_dir, "module_group_correlation_heatmap_significant.pdf"),
      width = 8, height = pdf_height)
  
  par(mar = c(5, 8, 4, 2) + 0.1)
  
  labeledHeatmap(
    Matrix = sig_cor_matrix,
    xLabels = c(GROUP1_NAME, GROUP2_NAME),
    yLabels = sig_colors,
    ySymbols = sig_colors,
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = sig_text,
    setStdMargins = FALSE,
    cex.text = 0.8,
    cex.lab.y = 0.8,
    zlim = c(-1, 1),
    main = "",
    xLabelsAngle = 45,
    yLabelsPosition = "left"
  )
  dev.off()
}


# ==== Plot Heatmap (All modules) ====
text_matrix_all <- matrix(
  sprintf("%.2f\n%s", module_group_cor_for_heatmap, module_diff_df$Significance),
  nrow = nrow(module_group_cor_for_heatmap),
  ncol = ncol(module_group_cor_for_heatmap),
  byrow = FALSE
)

pdf_height_all <- max(6, 0.4 * nrow(module_group_cor_for_heatmap))

pdf(file.path(output_dir, "module_group_correlation_heatmap_all.pdf"),
    width = 8, height = pdf_height_all)

par(mar = c(5, 8, 4, 2) + 0.1)

labeledHeatmap(
  Matrix = module_group_cor_for_heatmap,
  xLabels = c(GROUP1_NAME, GROUP2_NAME),
  yLabels = module_diff_df$Color, 
  ySymbols = module_diff_df$Color,
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = text_matrix_all,
  setStdMargins = FALSE,
  cex.text = 0.8,
  cex.lab.y = 0.8,
  zlim = c(-1, 1),
  main = "", 
  xLabelsAngle = 45,
  yLabelsPosition = "left"
)
dev.off()


# ====================== Boxplots for Significant Modules ======================
cat("Plotting boxplots for significant modules...\n")

sig_modules <- module_diff_df %>%
  filter(FDR < 0.05) %>%
  pull(Module)

if (length(sig_modules) > 0) {
  cat("Found", length(sig_modules), "significant differential modules.\n")
  
  for (mod in sig_modules) {
    color <- gsub("^ME", "", mod)
    
    # Prepare plot data
    me_data <- data.frame(
      ME = mergedMEs[, mod],
      Sample = rownames(mergedMEs)
    ) %>%
      inner_join(group_info %>% rownames_to_column("Sample"), by = "Sample")
    
    # Define comparison
    my_comparisons <- list(c(GROUP1_NAME, GROUP2_NAME))
    
    # Render boxplot
    p <- ggplot(me_data, aes(x = Group, y = ME, fill = Group)) +
      geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, color = "black") +
      geom_jitter(width = 0.2, size = 1.5, alpha = 0.5, color = "black") +
      scale_fill_manual(values = c("Normal" = "#377EB8", "Disease" = "#E41A1C")) +
      labs(
        title = "", 
        subtitle = "",
        y = "Module Eigengene Value",
        x = "Group"
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        legend.position = "none",
        axis.line = element_line(color = "black", size = 0.8),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black")
      ) +
      stat_compare_means(comparisons = my_comparisons,
                         method = "wilcox.test",
                         label = "p.signif",
                         label.y.npc = 0.95,
                         size = 6)
    
    # Save boxplot
    ggsave(file.path(output_dir, paste0("module_", color, "_boxplot_half_border.pdf")),
           p, width = 5, height = 6)
    
    # Save gene lists for differential modules
    diff_genes <- rownames(expr_data)[moduleColors == color]
    writeLines(diff_genes, file.path(output_dir, paste0("DIFF_MODULE_", color, "_GENES.txt")))
  }
} else {
  cat("No significant modules (FDR < 0.05) found. Skipping boxplots.\n")
}

cat("\n===== WGCNA Analysis Completed =====\n")

# ====================== Save Workspace ======================
save.image(file.path(output_dir, "wgcna_analysis_complete.RData"))
cat("All data saved to: wgcna_analysis_complete.RData\n")