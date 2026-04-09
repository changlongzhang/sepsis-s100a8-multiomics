# =========================================================================
# Script Name: 09_scRNA-seq & Bulk Immune Infiltration.R
# Description: Comprehensive single-cell RNA-seq pipeline (GSE167363).
#              Includes QC, DoubletFinder, SCTransform, Harmony integration,
#              in silico perturbation (scTenifoldKnk), TF activity (decoupleR),
#              custom CIBERSORT reference building, and bulk deconvolution.
# Note for GitHub users: Please update the working directories before running.
# =========================================================================

rm(list = ls()) # Clear environment
gc()            # Free memory

# --- 1. Set random seed (ensure reproducibility, no random changes) ---
set.seed(1234)

# --- 2. Load dependencies ---
library(Seurat)
library(tidyverse)
library(patchwork)
library(data.table)
library(stringr)
library(ggplot2)
# The following packages must be installed
library(DoubletFinder) # devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(harmony)       # devtools::install_github("immunogenomics/harmony")
library(clustree)

# --- 3. Set paths ---
# Please ensure this is the folder path where files like matrix.mtx.gz are stored
dir_path <- "D:\\桌面\\sepsis-s100a8-multiomics\\data\\09_scRNA-seq & Bulk Immune Infiltration\\Input\\GSE167363_RAW"
setwd(dirname(dir_path)) 
message(paste(">>> Working directory set to:", getwd()))

# =========================================================================
# [Step 1] Read data and build basic objects
# =========================================================================
message(">>> [Step 1] Reading data...")
all_files <- list.files(dir_path)
matrix_files <- all_files[grep("matrix.mtx.gz", all_files)]
samples <- gsub("_matrix.mtx.gz", "", matrix_files)
if(length(samples) == 0) samples <- gsub("matrix.mtx.gz", "", matrix_files)

seurat_list <- list()

for (sample in samples) {
  path_matrix <- file.path(dir_path, paste0(sample, "_matrix.mtx.gz"))
  path_barcodes <- file.path(dir_path, paste0(sample, "_barcodes.tsv.gz"))
  path_features <- file.path(dir_path, paste0(sample, "_genes.tsv.gz"))
  if (!file.exists(path_features)) path_features <- file.path(dir_path, paste0(sample, "_features.tsv.gz"))
  
  if (file.exists(path_matrix) & file.exists(path_barcodes) & file.exists(path_features)) {
    counts <- ReadMtx(mtx = path_matrix, cells = path_barcodes, features = path_features, feature.column = 2)
    
    # Create initial object
    obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
    
    # Add metadata
    # 1. Group
    if (grepl("HC", sample)) {
      obj$Group <- "Control"
    } else {
      obj$Group <- "Sepsis"
    }
    
    # 2. Timepoint
    if (grepl("HC", sample) | grepl("T0", sample)) {
      obj$Timepoint <- "Baseline"
    } else if (grepl("T6", sample)) {
      obj$Timepoint <- "Later"
    } else {
      obj$Timepoint <- "Unknown"
    }
    
    seurat_list[[sample]] <- obj
  }
}

# Merge objects
message(">>> Merging samples...")
sc_combined <- merge(x = seurat_list[[1]], y = seurat_list[2:length(seurat_list)], add.cell.ids = samples, project = "Sepsis")
rm(seurat_list); gc()
sc_combined <- JoinLayers(sc_combined) # Essential for Seurat V5

# =========================================================================
# [Step 2] Initial data filtering (Keep Baseline only)
# =========================================================================
message(">>> [Step 2] Filtering Baseline samples...")
pbmc <- subset(sc_combined, subset = Timepoint == "Baseline")
rm(sc_combined); gc()

# =========================================================================
# [Step 3] Advanced QC
# =========================================================================
message(">>> [Step 3] Executing advanced QC...")
# 3.1 Calculate metrics
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(pbmc))
pbmc[["percent.HB"]] <- PercentageFeatureSet(pbmc, features = HB.genes)

# 3.2 Filtering
# Criteria here: nFeature_RNA 200-6000, mitochondrial < 20%, hemoglobin < 1%
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & 
                 percent.mt < 20 & percent.HB < 1)
message(">>> QC complete. Remaining cell count: ", ncol(pbmc))

# =========================================================================
# [Step 4] Doublet Removal
# =========================================================================
message(">>> [Step 4] Running DoubletFinder (this will take some time)...")

# 4.1 Preprocessing (Only to provide parameters for DoubletFinder)
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:20, seed.use = 1234) # Fixed seed

# 4.2 Find pK
sweep.res.list <- paramSweep(pbmc, PCs = 1:20, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_best <- bcmvn %>% dplyr::arrange(desc(BCmetric)) %>% dplyr::pull(pK) %>% .[1] %>% as.character() %>% as.numeric()
message(paste(">>> Optimal pK value:", pk_best))

# 4.3 Predict and filter
doublet_rate <- ncol(pbmc) * 0.000008 # Estimate doublet rate
nExp_poi <- round(doublet_rate * ncol(pbmc))

pbmc <- doubletFinder(pbmc, PCs = 1:20, pN = 0.25, pK = pk_best, nExp = nExp_poi, reuse.pANN = NULL, sct = FALSE)

# Extract classification column name
meta_cols <- colnames(pbmc@meta.data)
df_class_col <- meta_cols[grep("DF.classifications", meta_cols)]
pbmc$Is_Double <- pbmc@meta.data[, df_class_col]

# Keep only singlets
pbmc <- subset(pbmc, subset = Is_Double == "Singlet")
# Clean metadata
pbmc@meta.data <- pbmc@meta.data[, !grepl("DF.classifications|pANN", colnames(pbmc@meta.data))]
message(">>> Doublet removal complete.")

# =========================================================================
# [Step 5] Cell Cycle Scoring
# =========================================================================
message(">>> [Step 5] Calculating cell cycle scores...")
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
pbmc <- CellCycleScoring(pbmc, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

# =========================================================================
# [Step 6] SCTransform Normalization
# =========================================================================
message(">>> [Step 6] Running SCTransform (regressing out batch and cycle effects)...")
# Regress out mitochondrial, S.Score, G2M.Score
pbmc <- SCTransform(pbmc, vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = TRUE)

# =========================================================================
# [Step 7] Dimensionality Reduction and Batch Correction (Harmony)
# =========================================================================
message(">>> [Step 7] Dimensionality reduction and Harmony batch correction...")
pbmc <- RunPCA(pbmc, verbose = FALSE)
# Harmony integration on orig.ident (different samples)
pbmc <- RunHarmony(pbmc, group.by.vars = "orig.ident", assay.use = "SCT", plot_convergence = FALSE)

# =========================================================================
# [Step 8] Clustering and UMAP (Crucial: Fixed seed)
# =========================================================================
message(">>> [Step 8] Final clustering and UMAP...")
dims_use <- 1:30

# 8.1 Find neighbors (using harmony reduction results)
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = dims_use)

# 8.2 Find clusters (using seed.use to fix results)
pbmc <- FindClusters(pbmc, resolution = 0.5, random.seed = 1234)

# 8.3 Run UMAP (using seed.use to fix results)
pbmc <- RunUMAP(pbmc, reduction = "harmony", dims = dims_use, seed.use = 1234)

# 8.4 Save intermediate results (prevent needing to rerun if it crashes)
message(">>> Saving prepared RDS...")
saveRDS(pbmc, "Sepsis_GSE167363_Prepared.rds")

# Check plotting
p_umap <- DimPlot(pbmc, reduction = "umap", label = T, pt.size = 0.5) + ggtitle("")
ggsave("Plot_01_Initial_UMAP.pdf", p_umap, width = 8, height = 6)

# =========================================================================
# [Part 2 Start] Prepare multi-gene validation dot plot (Analysis 2)
# =========================================================================
message(">>> [Analysis 2] Generating multi-gene validation plot...")

# 1. Switch back to RNA assay and perform standard LogNormalize (for visualization)
DefaultAssay(pbmc) <- "RNA"
pbmc <- JoinLayers(pbmc, assay = "RNA") # Merge V5 layers
pbmc <- NormalizeData(pbmc, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- ScaleData(pbmc, assay = "RNA") # Scale all genes for plotting

# 2. Define multi-gene Marker list (Joint validation of 3-4 genes to prevent mislabeling)
markers_list <- list(
  # --- Myeloid ---
  "Monocytes (Classic)"  = c("CD14", "VCAN", "FCN1", "LYZ"),       # Classical Monocytes (should be abundant)
  "Inflammation/Sepsis"  = c("S100A8", "S100A9", "MNDA"),          # Inflammatory signature (Monocytes/Neutrophils)
  "Non-Classical Mono"   = c("FCGR3A", "LST1", "MS4A7", "CDKN1C"), # CD16+ Non-classical Monocytes
  "Neutrophils"          = c("CSF3R", "FCGR3B", "NAMPT", "CXCR2"), # Neutrophils (if present, low CD14)
  "Dendritic Cells (DC)" = c("FCER1A", "CST3", "HLA-DRA", "HLA-DQA1"), # DC (MHC-II high)
  
  # --- Lymphoid ---
  "T Cells (Pan)"        = c("CD3D", "CD3E", "IL7R", "TRAC"),      # T cell core
  "CD8 T / NK-T"         = c("CD8A", "CD8B"),                      # CD8 T
  "NK Cells"             = c("GNLY", "NKG7", "KLRF1", "KLRD1"),    # NK (CD3 negative, GNLY high)
  "B Cells"              = c("MS4A1", "CD79A", "CD79B", "RALGPS2"),# B cells (CD20)
  
  # --- Others ---
  "Platelets"            = c("PPBP", "PF4", "GP1BB")               # Platelets
)

flat_markers <- unique(unlist(markers_list))

# 3. Draw robust dot plot
p_dot_robust <- DotPlot(pbmc, features = flat_markers, group.by = "seurat_clusters") + 
  coord_flip() + 
  scale_color_gradientn(colors = c("white", "lightgrey", "firebrick")) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 9)) + 
  labs(title = "", y = "Genes", x = "Clusters")

print(p_dot_robust)
ggsave("Plot_02_Annotation_DotPlot.pdf", p_dot_robust, width = 12, height = 14)
message("=========================================================================")


# =========================================================================
# 🟢 Continuing Analysis 2: Fill-in-the-blank Annotation (Based on DotPlot)
# =========================================================================
message(">>> [Step] Applying manual annotations (Based on DotPlot interpretation)...")

# 1. Fill in Cluster numbers (Confirmed based on the plot you sent)
mono_clusters    <- c("9")                # S100A8+, CD14+ (Core focus cluster)
t_clusters       <- c("0", "3", "7", "10") # CD3D+, IL7R+ (Includes CD4 and CD8)
b_clusters       <- c("1", "4", "5", "8")  # MS4A1+, CD79A+
nk_clusters      <- c("2")                # GNLY+, NKG7+, CD3-
platelet_cluster <- c("6")                # PPBP+, PF4+
dc_clusters      <- c("11")               # FCER1A+, HLA-DRA+
# Note: Although Cluster 9 has CSF3R expression, CD14/VCAN/S100A8 features better match inflammatory monocytes, not classified separately as neutrophils for now
neutro_clusters  <- c() 

# 2. Create ID mapping vector
current_ids <- levels(pbmc)
new_ids <- rep("Unknown", length(current_ids))
names(new_ids) <- current_ids

# 3. Mapping logic
if(length(mono_clusters) > 0)    new_ids[names(new_ids) %in% mono_clusters]    <- "Monocytes"
if(length(t_clusters) > 0)       new_ids[names(new_ids) %in% t_clusters]       <- "T Cells"
if(length(b_clusters) > 0)       new_ids[names(new_ids) %in% b_clusters]       <- "B Cells"
if(length(nk_clusters) > 0)      new_ids[names(new_ids) %in% nk_clusters]      <- "NK Cells"
if(length(platelet_cluster) > 0) new_ids[names(new_ids) %in% platelet_cluster] <- "Platelets"
if(length(dc_clusters) > 0)      new_ids[names(new_ids) %in% dc_clusters]      <- "DC"

# 4. Execute renaming
pbmc <- RenameIdents(pbmc, new_ids)
pbmc$Broad_CellType <- Idents(pbmc)

# Print results, you should see: Monocytes, T Cells, B Cells, NK Cells, Platelets, DC
message(">>> Annotation statistics results:")
print(table(pbmc$Broad_CellType))

# 5. Save the fully annotated object (Backup)
saveRDS(pbmc, "Sepsis_GSE167363_Annotated_Full.rds")

# =========================================================================
# 🟢 Analysis 3: Final cleaning (Remove DC, Platelets, bottom debris)
# =========================================================================
message(">>> [Analysis 3] Executing final streamlined cleaning...")

# --- 1. Define "excess" cells to remove ---
# We blacklist Platelets and DC (Cluster 11) together
cells_to_remove <- c("Platelets", "DC")

# Execute removal (invert = TRUE means keeping everything except these)
# Using %in% safely removes them without throwing errors even if some are already deleted
sc_temp <- subset(pbmc, idents = cells_to_remove, invert = TRUE)

message(paste(">>> Interference clusters removed:", paste(cells_to_remove, collapse = ", ")))

# --- 2. Remove bottom debris (Coordinate trimming) ---
# Get UMAP coordinates
umap_coords <- Embeddings(sc_temp, "umap")
sc_temp <- AddMetaData(sc_temp, umap_coords)

# Execute filtering: Keep only UMAP_2 > -10 (Remove bottom impurities)
sc_final_clean <- subset(sc_temp, subset = umap_2 > -8)

message(paste("Total before filtering:", ncol(pbmc)))
message(paste("Total after filtering:", ncol(sc_final_clean)))

# --- 3. Organize factor levels ---
sc_final_clean$Broad_CellType <- droplevels(sc_final_clean$Broad_CellType)
Idents(sc_final_clean) <- sc_final_clean$Broad_CellType

# --- 4. Plot final UMAP (Four main populations version) ---
# Here we only keep 4 main color schemes
my_cols_anno <- c(
  "Monocytes" = "#E31A1C", # Red (Protagonist)
  "T Cells"   = "#33A02C", # Green
  "B Cells"   = "#1F78B4", # Blue
  "NK Cells"  = "#B2DF8A"  # Light Green
)

p_final <- DimPlot(sc_final_clean, reduction = "umap", 
                   group.by = "Broad_CellType",
                   label = TRUE, label.size = 6, repel = TRUE, # Slightly larger font
                   cols = my_cols_anno, pt.size = 0.5) + 
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5))

# Save UMAP
message(">>> Saving final UMAP plot...")
ggsave("Plot_03_Cleaned_UMAP.pdf", p_final, width = 8, height = 6)

# --- 5. Target validation S100A8 (Final Version) ---
message(">>> Plotting S100A8 expression features...")
DefaultAssay(sc_final_clean) <- "RNA"

# Violin plot
p_vln <- VlnPlot(sc_final_clean, 
                 features = c("S100A8", "S100A9", "CD14"), 
                 split.by = "Group", 
                 group.by = "Broad_CellType", 
                 stack = TRUE, 
                 flip = TRUE, 
                 pt.size = 0, 
                 cols = c("Control"="#4DBBD5", "Sepsis"="#E64B35")) + 
  theme(axis.title.x = element_blank(),
        strip.text = element_text(face = "italic")) +
  ggtitle("")

ggsave("Plot_04_S100A8_Violin.pdf", p_vln, width = 8, height = 6)

# --- 6. Save final data ---
saveRDS(sc_final_clean, "Sepsis_GSE167363_Final_Clean.rds")


# 3. (Optional) Draw FeaturePlot heatmap
p_feat <- FeaturePlot(sc_final_clean, features = "S100A8", 
                      cols = c("lightgrey", "firebrick"), 
                      order = TRUE, pt.size = 0.5) + 
  ggtitle("S100A8 Expression")
ggsave("Plot_05_S100A8_FeaturePlot.pdf", p_feat, width = 8, height = 6)


## =========================================================================
## 🟢 Analysis 4: In silico drug grouping (In Silico Perturbation)
## Logic:
## 1. Lock onto Monocytes in the Sepsis group
## 2. Divide into High (top 33%) and Low (bottom 33%) based on S100A8 expression
## 3. Validate grouping significance
## =========================================================================
message(">>> [Analysis 4] Executing virtual grouping (High vs Low)...")

# 0. Load necessary packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)       # For ggboxplot and stat_compare_means
library(clusterProfiler)
library(org.Hs.eg.db) # GO analysis database

# --------------------------------------------------------------------------
# 1. Data Preparation: Extract Sepsis Monocytes
# --------------------------------------------------------------------------
# Ensure sc_final exists
if(!exists("sc_final_clean")) stop("❌ Error: Cannot find sc_final_clean object, please confirm Analysis 3 was executed!")

# Ensure data is pulled from RNA assay
DefaultAssay(sc_final_clean) <- "RNA"

# Extract Monocytes from the Sepsis group
sc_target <- subset(sc_final_clean, idents = "Monocytes", subset = Group == "Sepsis")
sc_target <- JoinLayers(sc_target) # Seurat V5 key step

message(paste(">>> Extracted Sepsis monocytes count:", ncol(sc_target)))

# --------------------------------------------------------------------------
# 2. Calculate tertile thresholds (Tertile Split)
# --------------------------------------------------------------------------
# Get S100A8 expression vector
vals <- LayerData(sc_target, assay = "RNA", layer = "data")["S100A8", ]

# Calculate thresholds based on positive expressing cells (avoids pulling down threshold with 0s, more accurate)
# If positive cells are too few, revert to using all cells
if(sum(vals > 0) > 10) {
  qs <- quantile(vals[vals > 0], probs = c(0.33, 0.66))
} else {
  qs <- quantile(vals, probs = c(0.33, 0.66))
}

message(paste0(">>> Grouping threshold set: Low < ", round(qs[1], 2), " | High > ", round(qs[2], 2)))

# --------------------------------------------------------------------------
# 3. Tagging & Building Analysis Object (sc_mech)
# --------------------------------------------------------------------------
# Initialize tags
sc_target$S100A8_Status <- "Mid"

# Assign High and Low
sc_target$S100A8_Status[vals >= qs[2]] <- "High"
sc_target$S100A8_Status[vals <= qs[1]] <- "Low"

# [Core Step] Keep only High and Low groups (Discard Mid)
sc_mech <- subset(sc_target, subset = S100A8_Status %in% c("High", "Low"))
Idents(sc_mech) <- "S100A8_Status"

# Print grouping counts
print(table(sc_mech$S100A8_Status))

# --------------------------------------------------------------------------
# 4. Plotting Validation (Boxplot + Asterisk Statistics)
# --------------------------------------------------------------------------
# Extract plotting data
df_plot <- data.frame(
  Expr = LayerData(sc_mech, assay = "RNA", layer = "data")["S100A8", ], 
  Group = sc_mech$S100A8_Status
)

# Force factor order: Low on the left, High on the right
df_plot$Group <- factor(df_plot$Group, levels = c("Low", "High"))

# Define comparison groups
my_comparisons <- list(c("Low", "High"))

# Plotting
p_box <- ggboxplot(df_plot, 
                   x = "Group", 
                   y = "Expr", 
                   fill = "Group", 
                   palette = c("#00AFBB", "#E7B800"), # Cyan vs Gold
                   outlier.shape = NA,                # Hide outliers
                   width = 0.5) +                     # Box width
  
  # Scatter layer (Add transparency)
  geom_jitter(width = 0.2, size = 1, alpha = 0.3) +
  
  # Statistics layer (Use asterisks)
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",       # Show asterisks (ns, *, **, ***, ****)
    symnum.args = list(       # Custom asterisk thresholds
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  
  # Labels and themes
  labs(y = "Log-Normalized Expression", x = "", title = "") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black")
  )

# Output
print(p_box)
ggsave("Plot_06_Virtual_Grouping_Boxplot.pdf", p_box, width = 4, height = 5)

message(">>> ✅ Analysis 4 Complete! Grouped object sc_mech generated.")


## =========================================================================
## 🟢 Analysis 5: Mechanism Exploration (High vs Low Differential Comparison)
## Logic: Strictly filter High-specific vs Low-specific pathways, remove shared ambiguous items
## =========================================================================
message(">>> [Analysis 5] Executing differential mechanism exploration (Strict deduplication version)...")

# 1. Differential Analysis (FindMarkers)
# --------------------------------------------------------------------------
message(">>> 1. Calculating DEGs (High vs Low)...")
# Set logfc slightly higher (0.15) to ensure selected genes are more representative
deg_mech <- FindMarkers(sc_mech, ident.1 = "High", ident.2 = "Low", 
                        logfc.threshold = 0.15, min.pct = 0.1, test.use = "wilcox")
write.csv(deg_mech, "S100A8_High_vs_Low_DEG.csv")

# Split up-regulated (High) and down-regulated (Low) genes
genes_high <- rownames(deg_mech)[deg_mech$avg_log2FC > 0.25 & deg_mech$p_val_adj < 0.05]
genes_low  <- rownames(deg_mech)[deg_mech$avg_log2FC < -0.25 & deg_mech$p_val_adj < 0.05]

message(paste0(">>> DEG Statistics: High-specific=", length(genes_high), " | Low-specific=", length(genes_low)))

# 2. GO Enrichment Analysis Function
# --------------------------------------------------------------------------
run_go <- function(genes) {
  if(length(genes) < 5) return(NULL)
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  if(is.null(ego)) return(NULL)
  return(ego@result %>% filter(p.adjust < 0.05))
}

# 3. Calculate enrichment separately
message(">>> 2. Calculating GO enrichment pathways...")
res_high <- run_go(genes_high)
res_low  <- run_go(genes_low)

# 4. [Core Logic] Smart Filtering + Strict Deduplication
# --------------------------------------------------------------------------
process_pathways <- function(df_high, df_low, top_n=10) {
  if(is.null(df_high) || is.null(df_low)) return(NULL)
  
  # A. Keyword filtering (prioritize displaying these related pathways)
  keywords <- "inflammatory|cytokine|immune|NF-kappa|response to bacterium|myeloid|leukocyte|taxis"
  
  df_high <- df_high %>% mutate(is_key = str_detect(Description, regex(keywords, ignore_case = T)))
  df_low  <- df_low  %>% mutate(is_key = str_detect(Description, regex(keywords, ignore_case = T)))
  
  # B. [Strict Deduplication] Remove pathways appearing in both groups (Maintain differential purity)
  common_terms <- intersect(df_high$Description, df_low$Description)
  if(length(common_terms) > 0) {
    message(paste(">>> Removed", length(common_terms), "ambiguous pathways common to both groups, retaining only specific differences."))
    df_high <- df_high %>% filter(!Description %in% common_terms)
    df_low  <- df_low  %>% filter(!Description %in% common_terms)
  }
  
  # C. Select Top N (First by keyword, then by P-value)
  top_high <- df_high %>% arrange(desc(is_key), p.adjust) %>% head(top_n) %>% mutate(Group = "S100A8-High")
  top_low  <- df_low  %>% arrange(desc(is_key), p.adjust) %>% head(top_n) %>% mutate(Group = "S100A8-Low")
  
  return(rbind(top_high, top_low))
}

plot_df <- process_pathways(res_high, res_low, 10) # Take top 10 from each side

# 5. Plotting (Compact bidirectional bar plot)
# --------------------------------------------------------------------------
if(!is.null(plot_df) && nrow(plot_df) > 0) {
  
  # Calculate plotting coordinates
  plot_df <- plot_df %>%
    mutate(
      log10P = -log10(p.adjust),
      # High points right (positive), Low points left (negative)
      plot_val = ifelse(Group == "S100A8-High", log10P, -log10P),
      stars = case_when(p.adjust <= 0.001 ~ "***", p.adjust <= 0.01 ~ "**", TRUE ~ "*"),
      # Truncate overly long text
      ShortName = str_trunc(Description, 45)
    )
  
  # Sort: Make bar lengths orderly
  plot_df$Description <- factor(plot_df$Description, 
                                levels = unique(plot_df$Description[order(plot_df$plot_val)]))
  
  # Dynamically calculate boundaries (Compact layout)
  max_val <- max(abs(plot_df$plot_val))
  max_limit <- max(max_val * 1.15, 2) # Leave 15% margin
  
  fill_cols <- c("S100A8-Low" = "#0074b3", "S100A8-High" = "#982b2b")
  
  p_go <- ggplot(plot_df, aes(x = Description, y = plot_val, fill = Group)) +
    geom_col(width = 0.7, color = "black", size = 0.3) + # Bars with thin black borders, clearer
    coord_flip() +
    
    scale_fill_manual(values = fill_cols, name = NULL) +
    scale_y_continuous(limits = c(-max_limit, max_limit), expand = c(0,0)) +
    
    labs(x = NULL, y = expression(Signed~-log[10]~(P-adjust)), 
         title = "") +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y = element_blank(), # Hide original Y axis
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 11),
      legend.position = "top",
      legend.text = element_text(face = "bold"),
      panel.grid.major.x = element_line(color = "grey92", linetype = "dashed") # Add faint auxiliary lines
    ) +
    
    # Text labels: Low in left whitespace, High in right whitespace
    geom_text(data = subset(plot_df, Group == "S100A8-Low"),
              aes(y = 0.2, label = ShortName), hjust = 0, size = 3.8, fontface = "bold") +
    
    geom_text(data = subset(plot_df, Group == "S100A8-High"),
              aes(y = -0.2, label = ShortName), hjust = 1, size = 3.8, fontface = "bold") +
    
    # Asterisks: Right next to the bars
    geom_text(aes(y = ifelse(plot_val > 0, plot_val + 0.2, plot_val - 0.2), label = stars), 
              vjust = 0.75, size = 4) +
    
    geom_hline(yintercept = 0, size = 0.6) # Central axis
  
  print(p_go)
  ggsave("Plot_07_Mechanism_GO_Barplot.pdf", p_go, width = 8, height = 6)
  
} else {
  warning("No significantly enriched pathways found, or no pathways remaining after deduplication.")
}

message(">>> 🎉 Differential mechanism analysis complete! Displaying only inter-group specific pathways.")


## =========================================================================
## 🟢 Analysis 6: Upstream Transcriptional Regulation (Aggressive cut: K-means forced outlier removal)
## =========================================================================

# 0. Load packages
library(Seurat)
library(tidyverse)
library(ggplot2)
library(decoupleR)
library(OmnipathR)
library(dplyr)
library(patchwork)

message(">>> [Analysis 6] Beginning execution...")

# Ensure sc_mech exists
if(!exists("sc_mech")) stop("❌ Error: Cannot find sc_mech object, please run Analysis 4 first!")

# ==========================================================================
# 1. Calculate TF activity (if run previously, this step is fast)
# ==========================================================================
message(">>> 1. Preparing TF data...")

if(!dir.exists("omnipathr-log")) {
  dir.create("omnipathr-log")
}

# Get network
net <- get_dorothea(organism = "human", levels = c("A", "B", "C"))
mat <- as.matrix(LayerData(sc_mech, assay = "RNA", layer = "data"))

# Run algorithm
tf_acts <- decouple(mat, network = net, .source = "source", .target = "target",
                    statistics = "wmean", args = list(center = TRUE, scale = TRUE))

# Convert format
tf_matrix <- tf_acts %>% 
  dplyr::filter(statistic == "wmean") %>%
  dplyr::select(condition, source, score) %>%
  dplyr::group_by(condition, source) %>%
  dplyr::summarise(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = source, values_from = score) %>% 
  tibble::column_to_rownames("condition") %>% 
  as.matrix() %>%
  t()

tf_matrix[is.na(tf_matrix)] <- 0
sc_mech[["tf"]] <- CreateAssayObject(data = tf_matrix)
DefaultAssay(sc_mech) <- "tf"
sc_mech <- ScaleData(sc_mech)

# ==========================================================================
# 2. Recalculate UMAP
# ==========================================================================
message(">>> 2. Recalculating UMAP coordinates...")
DefaultAssay(sc_mech) <- "RNA"
sc_mech <- FindVariableFeatures(sc_mech, verbose = FALSE)
sc_mech <- ScaleData(sc_mech, verbose = FALSE)
sc_mech <- RunPCA(sc_mech, verbose = FALSE)
sc_mech <- RunUMAP(sc_mech, dims = 1:10, verbose = FALSE)

# ==========================================================================
# 3. 🔴 Core Modification: K-means Clustering Cut Method
# ==========================================================================
message(">>> [Aggressive Cleaning] Using K-means to identify and remove outliers...")

# 1. Get UMAP coordinates
umap_coords <- Embeddings(sc_mech, "umap")

# 2. Apply K-means clustering to coordinates (Set to 2 clusters: Main group vs Outliers)
set.seed(123)
km_res <- kmeans(umap_coords, centers = 2)

# 3. Identify which cluster is the "Main group" (The one with the most cells)
cluster_counts <- table(km_res$cluster)
main_cluster_id <- names(cluster_counts)[which.max(cluster_counts)]
message(paste(">>> Identified main cluster ID:", main_cluster_id, "containing cell count:", max(cluster_counts)))

# 4. Identify cell IDs to retain
cells_to_keep <- names(km_res$cluster)[km_res$cluster == main_cluster_id]

# 5. Execute removal
sc_mech_clean <- subset(sc_mech, cells = cells_to_keep)
message(paste(">>> Removed", ncol(sc_mech) - ncol(sc_mech_clean), "outlier cells."))
sc_mech <- sc_mech_clean # Overwrite

# ==========================================================================
# 4. Plotting (Core Targets)
# ==========================================================================
message(">>> 3. Plotting core targets...")

DefaultAssay(sc_mech) <- "tf"
target_tfs <- c("RELA", "NFKB1", "STAT3", "CEBPB", "JUN", "FOS")
valid_tfs <- intersect(target_tfs, rownames(sc_mech))

if(length(valid_tfs) > 0) {
  p_tf <- FeaturePlot(sc_mech, features = valid_tfs, 
                      cols = c("grey95", "#800080"), 
                      ncol = 3, 
                      pt.size = 2.5, # Enlarge point size to make plot look better
                      order = TRUE,
                      combine = FALSE)
  
  # Retain subtitles and axes
  p_list <- lapply(p_tf, function(x) {
    x + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  })
  
  final_plot <- wrap_plots(p_list, ncol = 3)
  
  print(final_plot)
  ggsave("Plot_08_TF_Activity_FeaturePlot.pdf", final_plot, width = 12, height = 7)
}

message(">>> ✅ Aggressive cleaning complete! Please review Plot_08_TF_Activity_FeaturePlot.pdf")

# ==========================================================================
# 4. Plotting (Core Targets - Clean version without axes)
# ==========================================================================
message(">>> 3. Plotting core targets...")

DefaultAssay(sc_mech) <- "tf"
target_tfs <- c("RELA", "NFKB1", "STAT3", "CEBPB", "JUN", "FOS")
valid_tfs <- intersect(target_tfs, rownames(sc_mech))

if(length(valid_tfs) > 0) {
  p_tf <- FeaturePlot(sc_mech, features = valid_tfs, 
                      cols = c("grey95", "#800080"), 
                      ncol = 3, 
                      pt.size = 2.5, # Enlarge point size
                      order = TRUE,
                      combine = FALSE)
  
  # Modification: Added NoAxes() to remove axes
  p_list <- lapply(p_tf, function(x) {
    x + NoAxes() + theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  })
  
  final_plot <- wrap_plots(p_list, ncol = 3)
  
  print(final_plot)
  ggsave("Plot_09_TF_Activity_FeaturePlot_NoAxes.pdf", final_plot, width = 12, height = 7)
}




## =========================================================================
## 🟢 Analysis 4: In Silico Knockout - Final Adapted Version
## 🎯 Goal: Reproduce Figure 4 I, J, K, L
## 🧬 Scenario: Virtual knockout of S100A8 in Sepsis monocytes
## =========================================================================

# 0. Load necessary packages
library(Seurat)
# If scTenifoldKnk is not installed, run: devtools::install_github("cailab-tamu/scTenifoldKnk")
library(scTenifoldKnk) 
library(Matrix)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)

# Set random seed
set.seed(1234)

## -------------------------------------------------------------------------
## 1. Data Preparation and Deep Cleaning (Methods: Input Pre-processing)
## -------------------------------------------------------------------------
message(">>> [Step 1] Extracting monocytes and removing background noise...")

# Ensure input object exists
if(!exists("sc_final_clean")) stop("❌ Error: Cannot find sc_final_clean object, please run Analysis 3 first!")

# 1.1 Lock onto Monocytes in the Sepsis group
# Note: Use our previously defined Broad_CellType
target_obj <- subset(sc_final_clean, subset = Group == "Sepsis" & Broad_CellType == "Monocytes")

# Seurat V5 key step: Merge layers, ensure counts are not separated
target_obj <- JoinLayers(target_obj)

# 1.2 Extract raw counts matrix (Counts)
counts_matrix <- LayerData(target_obj, assay = "RNA", layer = "counts")
counts_matrix <- as.matrix(counts_matrix)

# 1.3 Key Denoising Step (Noise Removal)
# Remove mitochondrial, ribosomal, and hemoglobin genes
noise_pattern <- "^MT-|^RP[SL]|^HB[AB]" 
noise_genes <- grep(noise_pattern, rownames(counts_matrix), value = TRUE)
clean_genes <- setdiff(rownames(counts_matrix), noise_genes)
sc_clean <- counts_matrix[clean_genes, ]

# 1.4 Filter low expression genes
# Retain genes expressed in at least 5% of cells
n_cells <- ncol(sc_clean)
keep_genes <- rowSums(sc_clean > 0) > (0.05 * n_cells)
sc_input <- sc_clean[keep_genes, ]

# 1.5 Ensure target gene is still present
if(!"S100A8" %in% rownames(sc_input)) {
  stop("❌ Error: S100A8 was lost during filtering! Please try lowering the 0.05 filtering threshold.")
}

# 1.6 Convert to sparse matrix
sc_input_sparse <- Matrix(sc_input, sparse = TRUE)

message(paste(">>> Input matrix ready:", nrow(sc_input_sparse), "Genes x", ncol(sc_input_sparse), "Cells"))

## -------------------------------------------------------------------------
## 2. Run Virtual Knockout Algorithm (Methods: scTenifoldKnk)
## -------------------------------------------------------------------------
message(">>> [Step 2] Building manifold network and simulating knockout (estimated time 5-15 minutes)...")

# Run core function (nCores=1 is most stable)
tryCatch({
  sce_ko <- scTenifoldKnk(
    countMatrix = sc_input_sparse,
    gKO = "S100A8",
    nCores = 1  
  )
}, error = function(e) {
  message("⚠️ Run error. If out of memory, try reducing cell count (subset).")
  stop(e)
})

# Extract and tidy results
ko_results <- sce_ko$diffRegulation
ko_results <- na.omit(ko_results)

# Sort descending by perturbation distance
ko_results <- ko_results %>% arrange(desc(distance))

# Save raw data
write.csv(ko_results, "Result_Analysis4_S100A8_KO.csv")
message(">>> Simulation complete! Plotting...")

## -------------------------------------------------------------------------
## 3. Reproduce Figure 4I: Perturbation Distribution Plot (Half-frame SCI format)
## -------------------------------------------------------------------------
message(">>> [Step 3] Plotting Figure 4I (Half-frame format)...")

# Prepare plotting data (Remove target gene itself S100A8)
plot_data <- ko_results %>% 
  filter(gene != "S100A8") %>%
  arrange(desc(distance)) 

# Plot Z-score vs Distance scatter plot
p_4I <- ggplot(plot_data, aes(x = Z, y = distance)) +
  # Layer 1: Background gray points
  geom_point(color = "#E0E0E0", size = 1, alpha = 0.8) +
  
  # Layer 2: Highlight Top 50 genes
  geom_point(data = head(plot_data, 50), 
             aes(color = distance), size = 2.5) + 
  
  # Layer 3: Color mapping
  scale_color_gradient(low = "#FFD54F", high = "#D32F2F", name = "Distance") +
  
  # Layer 4: Gene labels (Top 15)
  geom_text_repel(data = head(plot_data, 15), 
                  aes(label = gene),
                  size = 4, 
                  fontface = "italic",
                  box.padding = 0.6, 
                  max.overlaps = 50,
                  segment.color = "black") +   
  
  # Layer 5: Theme settings
  theme_classic(base_size = 14) + 
  labs(title = "", subtitle = "", x = "Z-score", y = "Manifold Distance") +
  theme(
    axis.line = element_line(linewidth = 0.8, color = "black"),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold")
  )

print(p_4I)
ggsave("Figure_4I_KO_Distance_Scatter.pdf", p_4I, width = 8, height = 6)

## -------------------------------------------------------------------------
## 4. Reproduce Figure 4J: Top 20 Genes Bar Plot (Removed S100A9 Version)
## -------------------------------------------------------------------------
message(">>> [Step 4] Plotting Figure 4J (Final)...")

# 1. Filter out S100A9 (It is a heterodimer partner, physical binding causes close distance, usually needs removal to see downstream)
plot_data_final <- plot_data %>% 
  filter(gene != "S100A9") %>%  
  head(20)                      

# 2. Plotting
p_4J_final <- ggplot(plot_data_final, aes(x = reorder(gene, distance), y = distance)) +
  geom_bar(stat = "identity", fill = "#3C5488", width = 0.7) + 
  coord_flip() +
  theme_classic(base_size = 14) + 
  labs(title = "", x = "", y = "Manifold Distance") +
  theme(axis.text.y = element_text(size = 11, face = "bold", color = "black"),
        axis.line = element_line(linewidth = 0.8),
        axis.text.x = element_text(color = "black"))

print(p_4J_final)
ggsave("Figure_4J_KO_Distance_Barplot.pdf", p_4J_final, width = 8, height = 6)

## -------------------------------------------------------------------------
## 5. Calculate Enrichment Analysis (GO & KEGG) - Must calculate before plotting
## -------------------------------------------------------------------------
message(">>> [Step 5] Calculating enrichment analysis...")

# Extract significantly affected genes (Top 100) for enrichment
sig_genes_list <- head(plot_data$gene, 100)

# ID Conversion (Symbol -> Entrez)
gene_ids <- bitr(sig_genes_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

if(nrow(gene_ids) > 0) {
  # 5.1 KEGG Calculation
  kegg_res <- enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.5)
  
  # 5.2 GO Calculation
  go_res <- enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.5)
  
} else {
  warning("❌ Gene ID conversion failed, unable to perform enrichment analysis.")
}

## -------------------------------------------------------------------------
## 6. Reproduce Figure 4K, 4L: Minimalist Style Enrichment Analysis Plot
## -------------------------------------------------------------------------
message(">>> [Step 6] Plotting enrichment analysis (Clean Style)...")

# === Define plotting function (Minimalist version) ===
plot_clean_enrichment <- function(enrich_res, analysis_type = "KEGG") {
  
  if (is.null(enrich_res) || nrow(enrich_res@result) == 0) return(NULL)
  
  df <- as.data.frame(enrich_res) %>%
    filter(p.adjust < 0.5) %>%
    head(20) %>%
    mutate(logP = -log10(p.adjust))
  
  # Smart classification logic
  if (analysis_type == "KEGG") {
    df <- df %>%
      mutate(Category = case_when(
        grepl("Differentiation|Antigen|Signaling|Leukocyte|Lymphocyte|Th17|Th1|Cell adhesion", Description, ignore.case = T) ~ "Immune Response",
        grepl("Arthritis|Asthma|Diabetes|Bowel|Graft|Allograft|Autoimmune|Sclerosis|Lupus", Description, ignore.case = T) ~ "Autoimmune & Inflammation",
        grepl("Infection|Influenza|Leishmania|Tuberculosis|Measles|Hepatitis|Viral|Herpes|Epstein|Papilloma", Description, ignore.case = T) ~ "Infectious Diseases",
        TRUE ~ "Others"
      ))
    df$Category <- factor(df$Category, levels = c("Immune Response", "Infectious Diseases", "Autoimmune & Inflammation", "Others"))
    custom_colors <- c("Immune Response"="#E64B35", "Autoimmune & Inflammation"="#00A087", "Infectious Diseases"="#3C5488", "Others"="#8491B4")
    
  } else if (analysis_type == "GO") {
    df <- df %>%
      mutate(Category = case_when(
        grepl("Antigen|MHC|Peptide|Presentation|Leukocyte|Myeloid|Innate", Description, ignore.case = T) & !grepl("T cell|Lymphocyte", Description, ignore.case = T) ~ "Antigen Presentation & Innate",
        grepl("T cell|Cytotoxicity|Killing|Lymphocyte|Cell killing", Description, ignore.case = T) ~ "T Cell Immunity & Cytotoxicity",
        grepl("RNA|Translation|Catabolic|Nucleic|Metabolic|Biosynthetic|Splicing", Description, ignore.case = T) ~ "Gene Expression & Metabolism",
        TRUE ~ "Others"
      ))
    df$Category <- factor(df$Category, levels = c("Antigen Presentation & Innate", "T Cell Immunity & Cytotoxicity", "Gene Expression & Metabolism", "Others"))
    custom_colors <- c("Antigen Presentation & Innate"="#E64B35", "T Cell Immunity & Cytotoxicity"="#3C5488", "Gene Expression & Metabolism"="#00A087", "Others"="#8491B4")
  }
  
  # Sorting
  df <- df %>%
    group_by(Category) %>%
    arrange(desc(logP)) %>%
    ungroup() %>%
    mutate(Description = factor(Description, levels = rev(unique(Description))))
  
  # Plotting
  p <- ggplot(df, aes(x = logP, y = Description)) +
    geom_col(aes(fill = Category), width = 0.65, alpha = 0.9) +
    geom_point(aes(size = Count), color = "black", fill = "black", shape = 21) +
    geom_text(aes(label = Count), color = "white", size = 2.5, fontface = "bold") +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = custom_colors) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "-log10(FDR)", y = NULL, title = "") +
    theme_classic(base_size = 12) +
    theme(
      strip.text.y = element_blank(), # Hide left facet text
      strip.background = element_blank(),
      axis.line = element_line(linewidth = 0.6),
      axis.text = element_text(color = "black"),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
      legend.position = "right",
      legend.title = element_text(face = "bold")
    )
  return(p)
}

# === Execute plotting and saving ===

# KEGG Plot (Figure 4L)
if(exists("kegg_res")) {
  p_4L_clean <- plot_clean_enrichment(kegg_res, analysis_type = "KEGG")
  if(!is.null(p_4L_clean)) {
    print(p_4L_clean)
    ggsave("Figure_4L_KO_KEGG_Enrichment.pdf", p_4L_clean, width = 8, height = 6)
  }
}

# GO Plot (Figure 4K)
if(exists("go_res")) {
  p_4K_clean <- plot_clean_enrichment(go_res, analysis_type = "GO")
  if(!is.null(p_4K_clean)) {
    print(p_4K_clean)
    ggsave("Figure_4K_KO_GO_Enrichment.pdf", p_4K_clean, width = 12, height = 9)
  }
}

message(">>> 🎉 Analysis 4 fully complete! Charts generated.")



## =========================================================================
## Project: Sepsis Multi-omics Validation (Final Version)
## Scenario: Leveraging annotated single-cell data (RData) to empower Bulk immune infiltration analysis
## Output: Figure 5 full set of charts (A-E)
## =========================================================================

# --- 0. Load dependencies ---
library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(corrplot)
library(RColorBrewer)

# Ensure tigeR is installed
# if(!require(tigeR)) devtools::install_github("fayjw/tigeR")
library(tigeR)

# --- 1. Set paths (Strictly following your instructions) ---
# Single-cell data path (Read RData)
sc_dir   <- "D:\\桌面\\sepsis-s100a8-multiomics\\data\\09_scRNA-seq & Bulk Immune Infiltration\\Input"

# Immune infiltration analysis working directory (Read Bulk data & Save results)
work_dir <- "D:\\桌面\\sepsis-s100a8-multiomics\\data\\09_scRNA-seq & Bulk Immune Infiltration\\Input"

# Check directories
if(!dir.exists(work_dir)) dir.create(work_dir, recursive = T)
setwd(work_dir)
message(">>> Working directory set to: ", getwd())

# 1.2 Lock onto core object (Usually called sc_final)
if(exists("sc_final_clean")) {
  sc_obj <- sc_final_clean
  message(">>> Successfully loaded 'sc_final_clean' object.")
} else {
  # Auto-find the largest Seurat object
  objs <- ls()[sapply(ls(), function(x) inherits(get(x), "Seurat"))]
  if(length(objs) > 0) {
    sc_obj <- get(objs[1])
    message(paste(">>> sc_final_clean not found, auto-using object:", objs[1]))
  } else {
    stop("❌ Seurat object not found in RData!")
  }
}
DefaultAssay(sc_obj) <- "RNA"

# 1.3 Check if annotated (Must contain Monocytes)
if(!"Monocytes" %in% levels(sc_obj)) {
  # Try switching Idents to Broad_CellType
  if("Broad_CellType" %in% colnames(sc_obj@meta.data)) {
    Idents(sc_obj) <- "Broad_CellType"
  }
  
  if(!"Monocytes" %in% levels(sc_obj)) {
    print(table(Idents(sc_obj)))
    stop("❌ 'Monocytes' annotation not found in object, please check if RData is the completed result from Part 2.")
  }
}

message(">>> ✅ Data loading verification passed. Current cell types:")
print(table(Idents(sc_obj)))

# 1.4 Execute subpopulation refinement (S100A8 High vs Low)
message(">>> Executing monocyte subpopulation refinement...")

# Extract monocytes
mono_cells <- subset(sc_obj, idents = "Monocytes")

# Calculate threshold (Top 40% as High)
expr_data <- FetchData(mono_cells, vars = "S100A8")
threshold <- quantile(expr_data$S100A8, 0.6)

# Define High/Low
mono_cells$Subtype <- ifelse(expr_data$S100A8 > threshold, "Mono_S100A8_Hi", "Mono_S100A8_Lo")

# Backfill into main object
# First create Ref_CellType column, duplicate original annotations
sc_obj$Ref_CellType <- as.character(Idents(sc_obj))
# Only modify monocyte labels
sc_obj$Ref_CellType[colnames(mono_cells)] <- mono_cells$Subtype

# Format cleaning (Remove spaces)
sc_obj$Ref_CellType <- gsub(" ", "_", sc_obj$Ref_CellType)
sc_obj$Ref_CellType <- gsub("-", "_", sc_obj$Ref_CellType)

# Set as current Idents to build matrix
Idents(sc_obj) <- "Ref_CellType"

message(">>> Subpopulation refinement complete! Building Reference will use the following classifications:")
print(table(Idents(sc_obj)))

## =========================================================================
## =========================================================================
## =========================================================================
## 🟢 Step 2: Build Custom Reference Matrix (Seurat Native Version - Avoid errors)
## =========================================================================
message(">>> [Step 2] Manually building Reference matrix (Seurat Native)...")

# 1. Ensure identities are set correctly
Idents(sc_obj) <- "Ref_CellType"

# 2. Find Marker genes for each subpopulation (Extract features)
# This step ensures we only use the most specific genes as a ruler
message("   -> Calculating Marker genes (this will take some time)...")
all_markers <- FindAllMarkers(
  object = sc_obj, 
  only.pos = TRUE,      # Find only highly expressed genes
  min.pct = 0.25, 
  logfc.threshold = 1.0 # Strict threshold
)

# 3. Filter Top 100 genes
library(dplyr)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)

ref_genes <- unique(top_markers$gene)
message(paste("   -> Filtered", length(ref_genes), "specific genes to build matrix"))

# 4. Calculate average expression (This is exactly what build_CellType_Ref does internally)
# Use Seurat's AverageExpression function
avg_expr <- AverageExpression(
  sc_obj, 
  features = ref_genes, 
  group.by = "Ref_CellType",
  assays = "RNA", 
  slot = "data" 
)$RNA

# 5. Convert to matrix and save
ref_sig_mtr <- as.matrix(avg_expr)
write.csv(ref_sig_mtr, "Sepsis_Custom_Reference_Matrix.csv")

message(">>> ✅ Reference matrix build complete! (Bypassed tigeR errors)")
print(head(ref_sig_mtr[, 1:2]))
## =========================================================================
## 🟢 Step 3: Read Bulk Data (GSE65682)
## =========================================================================
message(">>> [Step 3] Reading Bulk expression matrix...")

bulk_file <- "GSE65682_gene.csv"
if(!file.exists(bulk_file)) stop("❌ Cannot find GSE65682_gene.csv in the immune infiltration folder!")

exp_raw <- read.csv(bulk_file, header = T, check.names = F)

# Handle duplicate genes (Take average)
if(sum(duplicated(exp_raw[,1])) > 0){
  exp_agg <- aggregate(. ~ exp_raw[,1], data = exp_raw, FUN = mean)
  rownames(exp_agg) <- exp_agg[,1]
  bulk_mat <- as.matrix(exp_agg[,-1])
} else {
  rownames(exp_raw) <- exp_raw[,1]
  bulk_mat <- as.matrix(exp_raw[,-1])
}

# Ensure no negative values (CIBERSORT requires linear space)
if(min(bulk_mat) < 0) {
  message("⚠️ Negative values detected, restoring Log data to linear space...")
  bulk_mat <- 2^bulk_mat
}

## =========================================================================
## =========================================================================
## 🟢 Step 4: Execute Deconvolution (Native CIBERSORT Version - Deprecate tigeR)
## Core Objective: Save the matrix from memory as temporary files, pass to CIBERSORT script for computation
## =========================================================================
message(">>> [Step 4] Running native CIBERSORT...")

# 1. Check required script
if(!file.exists("CIBERSORT.R")) {
  stop("❌ Error: Missing 'CIBERSORT.R' file in working directory!\nPlease copy this file to: ", getwd())
}
source("CIBERSORT.R")

# 2. Prepare input files 
# Native CIBERSORT does not accept R variables, only txt file paths
# So we need to first write the matrices calculated in Step 2 and Step 3 to the hard drive

message("   -> Generating temporary input files...")
# Save Reference matrix
write.table(ref_sig_mtr, "Temp_Ref_Sig.txt", sep = "\t", quote = F, col.names = NA)
# Save Bulk expression matrix
write.table(bulk_mat, "Temp_Bulk_Expr.txt", sep = "\t", quote = F, col.names = NA)

# 3. Run algorithm
message("   -> Executing calculation (perm=100)...")

tryCatch({
  # Run CIBERSORT
  # perm=100 for fast results, recommend changing to 1000 for publications
  # QN=TRUE is standard normalization for Bulk data
  cibersort_raw <- CIBERSORT("Temp_Ref_Sig.txt", "Temp_Bulk_Expr.txt", perm = 1000, QN = TRUE)
  
  # 4. Results Cleaning
  # The last three columns of CIBERSORT output are statistics (P-value, Correlation, RMSE)
  # We only want cell proportions for subsequent plotting, so we remove them
  if(ncol(cibersort_raw) > 3) {
    # Keep all columns except the last 3
    res_infiltration <- cibersort_raw[, 1:(ncol(cibersort_raw)-3)]
  } else {
    res_infiltration <- cibersort_raw
  }
  
  # 5. Save final pristine results
  write.csv(res_infiltration, "Final_Deconvolution_Results.csv")
  
  # 6. Delete temporary files (Keep directory clean)
  if(file.exists("Temp_Ref_Sig.txt")) file.remove("Temp_Ref_Sig.txt")
  if(file.exists("Temp_Bulk_Expr.txt")) file.remove("Temp_Bulk_Expr.txt")
  
  message(">>> ✅ Deconvolution successful! Results saved as Final_Deconvolution_Results.csv")
  message(">>> Variable res_infiltration is ready, please proceed running Step 5 for plotting.")
  
}, error = function(e) {
  message("❌ CIBERSORT run failed. Please check input data for null or negative values. Error info:")
  print(e)
})
## =========================================================================
## 🟢 Step 5: Plot Figure 5 Full Set of Charts (A-E)
## =========================================================================
message(">>> [Step 5] Starting chart generation...")

# Read groupings (Must exist)
group_file <- "GSE65682_Groups.csv"
if(!file.exists(group_file)) stop("❌ Cannot find GGSE65682_Groups.csv")
group_info <- read.csv(group_file, header = T, row.names = 1, check.names = F)

# Tidy plotting data
common_sam <- intersect(rownames(res_infiltration), rownames(group_info))
res_plot <- res_infiltration[common_sam, ]
group_plot <- group_info[common_sam, , drop=FALSE]

plot_long <- res_plot %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(Group = group_plot[sample, "Group"]) %>%
  pivot_longer(cols = -c(sample, Group), names_to = "CellType", values_to = "Proportion")


## =========================================================================
## 🚑 Figure 5A Ultimate Fix Version: Resolve PDF corruption issue
## =========================================================================

library(dplyr)
library(pheatmap)
library(Seurat)

message(">>> Executing Figure 5A repair procedure...")

# 0. [Key] Force close all opened plotting devices to prevent file occupation
while(!is.null(dev.list())) dev.off()

# 1. Prepare Marker data (Carry over previous logic)
Idents(sc_obj) <- "Ref_CellType"
if(!exists("all_markers")) {
  message("Calculating Markers (please wait)...")
  all_markers <- FindAllMarkers(sc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.0)
}

# 2. Deduplicate and filter (Maintain your logic)
top_plot_genes_clean <- all_markers %>%
  group_by(gene) %>%
  top_n(n = 1, wt = avg_log2FC) %>% 
  ungroup() %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>% 
  arrange(cluster)

# 3. Extract gene list
plot_genes_list <- top_plot_genes_clean$gene

# 4. Calculate average expression
avg_expr_plot <- AverageExpression(sc_obj, features = plot_genes_list, assays = "RNA")$RNA

# 5. [Key Fix] Ensure matrix and gene list strictly align
# Sometimes AverageExpression loses certain genes, or order changes, must re-align
valid_genes <- intersect(top_plot_genes_clean$gene, rownames(avg_expr_plot))

# Re-filter marker list, keep only those in matrix
top_plot_genes_final <- top_plot_genes_clean %>% 
  filter(gene %in% valid_genes) %>%
  arrange(cluster) # Re-sort by cluster to ensure stair-step pattern

# Extract matrix by final list
plot_mat_clean <- avg_expr_plot[top_plot_genes_final$gene, ]
# Sort by cell name
plot_mat_clean <- plot_mat_clean[, sort(colnames(plot_mat_clean))]

# 6. Standardize & Handle NA values
plot_mat_clean <- t(scale(t(plot_mat_clean)))
# [Fix] Replace NaN generated by scale (genes with 0 variance) with 0 to prevent plot crash
plot_mat_clean[is.na(plot_mat_clean)] <- 0 
# Truncate extreme values
plot_mat_clean[plot_mat_clean > 2.5] <- 2.5
plot_mat_clean[plot_mat_clean < -1.5] <- -1.5

# 7. Build row annotations (Ensure order and row names are 100% identical)
annotation_row <- data.frame(
  CellType = top_plot_genes_final$cluster
)
rownames(annotation_row) <- top_plot_genes_final$gene

# 8. Plotting (Use filename parameter to save directly, safer than pdf() function)
message(">>> Writing to PDF file...")

tryCatch({
  pheatmap(plot_mat_clean, 
           cluster_rows = FALSE, # Key for stair-step
           cluster_cols = FALSE, 
           show_rownames = FALSE, 
           show_colnames = TRUE,
           annotation_names_row = FALSE,
           annotation_row = annotation_row,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
           main = "",
           fontsize_col = 10,
           angle_col = 45,
           # [Key] Specify filename directly here
           filename = "Figure_5A_Reference_Heatmap_Vertical.pdf", 
           width = 6, 
           height = 8)
  
  message(">>> ✅ Success! Please review Figure_5A_Reference_Heatmap_Vertical.pdf")
  
}, error = function(e) {
  message("❌ Plotting failed, error info below:")
  print(e)
  # If it errors out here, dev.off() is already handled in step 0, will not cause file corruption
})


## =========================================================================
## 🚑 Figure 5A Ultimate Full Version (Horizontal Layout + Bottom Horizontal Legend)
## Layout: Rows(Y-axis) = Cell Types | Columns(X-axis) = Marker Genes
## Legend: Located at bottom, arranged horizontally
## =========================================================================

library(dplyr)
library(Seurat)

# Auto-detect and load ComplexHeatmap (Required for plotting)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  message(">>> Installing ComplexHeatmap ...")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}
library(ComplexHeatmap)
library(circlize) # Used for color setup

message(">>> Processing data ...")

# 0. Clear canvas
while(!is.null(dev.list())) dev.off()

# ============================================================
# 1. Prepare Data (Seurat standard pipeline)
# ============================================================
Idents(sc_obj) <- "Ref_CellType"

# Ensure all_markers exists
if(!exists("all_markers")) {
  message("--- Calculating Markers (this may take some time) ---")
  all_markers <- FindAllMarkers(sc_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1.0)
}

# 2. Filter genes (Top 20 per cluster)
top_plot_genes_clean <- all_markers %>%
  group_by(gene) %>%
  top_n(n = 1, wt = avg_log2FC) %>% 
  ungroup() %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) %>% 
  arrange(cluster)

# 3. Calculate average expression
avg_expr_plot <- AverageExpression(sc_obj, features = top_plot_genes_clean$gene, assays = "RNA")$RNA

# 4. Align gene order
valid_genes <- intersect(top_plot_genes_clean$gene, rownames(avg_expr_plot))
top_plot_genes_final <- top_plot_genes_clean %>% 
  filter(gene %in% valid_genes) %>%
  arrange(cluster)

# Extract matrix and sort by gene list
plot_mat <- avg_expr_plot[top_plot_genes_final$gene, ]
# Sort columns by cell type name
plot_mat <- plot_mat[, sort(colnames(plot_mat))]

# ============================================================
# 2. Matrix Transpose and Standardization (Core Step)
# ============================================================
# Logic: Z-score (Scale) on gene dimension first, then transpose, ensures color represents "relative high/low of that gene"
plot_mat_scaled <- t(scale(t(plot_mat)))

# Handle extreme values (Trim outliers to ensure color contrast)
plot_mat_scaled[is.na(plot_mat_scaled)] <- 0
plot_mat_scaled[plot_mat_scaled > 2.5] <- 2.5
plot_mat_scaled[plot_mat_scaled < -1.5] <- -1.5

# Transpose: Becomes [Rows=Cell Types, Cols=Genes]
plot_mat_final <- t(plot_mat_scaled)

# ============================================================
# 3. Build Annotation Bar (Annotation)
# ============================================================
# Here Gene_Group corresponds to which Cluster each gene belongs to
annotation_col <- data.frame(
  Gene_Group = top_plot_genes_final$cluster 
)
rownames(annotation_col) <- top_plot_genes_final$gene

# Set heatmap colors (Blue-White-Red)
col_fun <- colorRamp2(c(-1.5, 0, 2.5), c("navy", "white", "firebrick3"))

# Define top colored annotation bar
# [Key Setting] annotation_legend_param controls legend style
column_ha <- HeatmapAnnotation(
  df = annotation_col,
  which = "column",             # It's a column annotation
  show_annotation_name = FALSE, # Don't show "Gene_Group" text on plot, only in legend
  
  # --- Core Modification: Force horizontal legend ---
  annotation_legend_param = list(
    Gene_Group = list(
      title = "Gene Group",     # Legend title
      direction = "horizontal", # Direction: horizontal
      nrow = 1,                 # Number of rows: 1 (Force into single row)
      title_position = "leftcenter" # Title on left of color blocks
    )
  )
)

# ============================================================
# 4. Plotting and Saving
# ============================================================
message(">>> Plotting Figure 5A (Bottom horizontal legend version) ...")

pdf("Figure_5A_Reference_Heatmap_Horizontal.pdf", width = 10, height = 6) # Wider to accommodate horizontal legend

ht <- Heatmap(plot_mat_final,
              name = "Expression",  # Name of heatmap color bar
              col = col_fun,
              
              # Clustering control (Maintain custom order)
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              
              # Label control
              show_row_names = TRUE,    # Show cell types (Y-axis)
              show_column_names = FALSE,# Don't show gene names (X-axis, avoids overlap)
              row_names_side = "left",  # Cell type names on left
              row_names_gp = gpar(fontsize = 12),
              
              # Bind annotation bar
              top_annotation = column_ha,
              
              # Border
              border = TRUE,
              
              # Set the heatmap's own color bar (Expression) to horizontal as well
              heatmap_legend_param = list(
                direction = "horizontal",
                title_position = "leftcenter"
              )
)

# [Ultimate Step] Draw plot and place legend at bottom
draw(ht, 
     heatmap_legend_side = "bottom",    # Heatmap color bar at bottom
     annotation_legend_side = "bottom", # Grouping legend at bottom
     merge_legend = TRUE)               # Merge these two legends

dev.off()

message(">>> ✅ Complete! Image saved as: Figure_5A_Reference_Heatmap_Horizontal.pdf")
message(">>> Please open and review, the legend should be at the very bottom and arranged in a single row.")


## =========================================================================
## 🚑 Figure 5B Beautification Redraw: Standard publication-grade stacked bar chart
## =========================================================================
## =========================================================================
## 🚑 Figure 5B Ultimate Fix: Custom width stacked bar chart (Patchwork version)
## =========================================================================
library(ggplot2)
library(dplyr)
library(RColorBrewer)
# If you don't have the patchwork package, install first: install.packages("patchwork")
library(patchwork) 

message(">>> Plotting custom ratio Figure 5B ...")

# 1. Prepare colors (Keep consistent with your previous ones)
cell_types <- unique(plot_long$CellType)
nb_cols <- length(cell_types)
cell_cols <- colorRampPalette(brewer.pal(12, "Paired"))(nb_cols)

# 2. Split data
df_disease <- plot_long %>% filter(Group == "Disease") # Or "Sepsis"
df_normal  <- plot_long %>% filter(Group == "Normal")  # Or "Control"

# 3. Plot left figure (Disease)
p_left <- ggplot(df_disease, aes(x = sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = cell_cols) +
  labs(x = "Disease", y = "Estimated Proportion") + # X-axis label set to group name
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none", # No legend on left plot, unify at end
    plot.margin = unit(c(0.2, 0.1, 0.2, 0.2), "cm") # Adjust margins, leave some space on right
  )

# 4. Plot right figure (Normal)
p_right <- ggplot(df_normal, aes(x = sample, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", width = 1) +
  scale_y_continuous(expand = c(0.01, 0)) +
  scale_fill_manual(values = cell_cols) +
  labs(x = "Normal", y = "") + # Y-axis title blank, sharing left side's
  theme_bw() +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_blank(), # Do not show Y-axis tick numbers on right plot
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm") # No space on left, flush against left plot
  )

# 5. Combine plots (Core Step)
# widths = c(3, 1) means the left plot is 3 times the width of the right plot (you can change to 2, 1 or 4, 1)
# guides = "collect" means merge legends
p_final <- p_left + p_right + 
  plot_layout(widths = c(10, 1), guides = "collect") & 
  theme(legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),
        legend.text = element_text(size = 9))

# 6. Save
ggsave("Figure_5B_Immune_Infiltration_StackPlot.pdf", p_final, width = 10, height = 6)
print(p_final)

message(">>> ✅ Figure 5B Fix Complete! Normal group should now be clearly visible.")


## =========================================================================
## 🚑 Figure 5C Advanced Version: Display differences for all cell types
## =========================================================================
library(ggplot2)
library(ggpubr)
library(dplyr)

message(">>> Plotting Figure 5C including all cell types ...")

# 1. Ensure data exists
if(!exists("plot_long")) {
  stop("❌ Error: Cannot find plotting data plot_long. Please ensure you have run the 'Tidy plotting data' section at the start of Step 5.")
}

# 2. Data cleaning (No filtering, but remove samples where Group is NA)
plot_data_all <- plot_long %>% 
  filter(!is.na(Group)) # Remove empty groupings

# 3. (Optional) Sort cells
# To make the plot look better, we can sort by "average abundance" descending, so tall bars are on the left
cell_order <- plot_data_all %>%
  group_by(CellType) %>%
  summarise(mean_prop = mean(Proportion)) %>%
  arrange(desc(mean_prop)) %>%
  pull(CellType)

plot_data_all$CellType <- factor(plot_data_all$CellType, levels = cell_order)

# 4. Plotting
# Note: Because there are many cells, we need to make the plot wider, otherwise X-axis text will crowd
p_5c_all <- ggboxplot(plot_data_all, x = "CellType", y = "Proportion",
                      color = "Group", palette = c("#E31A1C", "#1F78B4"), # Red/Blue color scheme
                      add = "jitter", 
                      short.panel.labs = FALSE) +
  
  # Add statistical testing (Wilcoxon test)
  stat_compare_means(aes(group = Group), label = "p.signif", method = "wilcox.test") +
  
  # Theme adjustments
  theme_bw() +
  labs(x = "", y = "Estimated Proportion", 
       title = "") +
  
  # X-axis label adjustments: Because many cells, rotate 45 degrees to prevent overlap
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
        legend.position = "top")

# 5. Save
# Key point: Set width to 12 or 14 to ensure all cells spread out
ggsave("Figure_5C_Immune_Infiltration_Boxplot.pdf", p_5c_all, width = 12, height = 9)

print(p_5c_all)
message(">>> ✅ All-cell plotting complete! Please review Figure_5C_Immune_Infiltration_Boxplot.pdf")





## =========================================================================
## 🚑 Figure 5D Ultimate Fix Version: Manually build link data, bypass object errors
## =========================================================================

library(linkET)
library(tidyverse)
library(RColorBrewer)

# =========================================================================
message(">>> Reloading missing data from local CSV files...")

# 1. Reload CIBERSORT results (cell_prop_clean)
if(file.exists("Final_Deconvolution_Results.csv")) {
  cell_prop_clean <- read.csv("Final_Deconvolution_Results.csv", row.names = 1, check.names = FALSE)
  message("  - Loaded Final_Deconvolution_Results.csv")
} else {
  stop("❌ Cannot find Final_Deconvolution_Results.csv! Did you successfully run Step 4 (CIBERSORT)?")
}

# 2. Reload Bulk Expression data (exp_mat)
if(file.exists("GSE65682_gene.csv")) {
  exp_raw <- read.csv("GSE65682_gene.csv", header = T, check.names = F)
  # Basic deduplication to create matrix
  if(sum(duplicated(exp_raw[,1])) > 0){
    exp_agg <- aggregate(. ~ exp_raw[,1], data = exp_raw, FUN = mean)
    rownames(exp_agg) <- exp_agg[,1]
    exp_mat <- as.matrix(exp_agg[,-1])
  } else {
    rownames(exp_raw) <- exp_raw[,1]
    exp_mat <- as.matrix(exp_raw[,-1])
  }
  message("  - Loaded GSE65682_gene.csv")
} else {
  stop("❌ Cannot find GSE65682_gene.csv in your working directory!")
}

# 3. Reload Group Info (group_info)
if(file.exists("GSE65682_Groups.csv")) {
  group_info <- read.csv("GSE65682_Groups.csv", header = T, row.names = 1, check.names = F)
  message("  - Loaded GSE65682_Groups.csv")
} else {
  stop("❌ Cannot find GSE65682_Groups.csv in your working directory!")
}

message(">>> ✅ Data reload complete! You can now run the Figure 5D code.")
# =========================================================================


message(">>> [Step 5D-Final-Fix] Preparing merged butterfly plot data...")

# ================= 1. Data Preparation (Same as before) =================

# Check variables
if(!exists("cell_prop_clean") | !exists("exp_mat") | !exists("group_info")) {
  stop("❌ Missing data! Please run previous data reading steps first.")
}

# 1.1 Lock onto samples
if("Disease" %in% group_info$Group) {
  target_group <- "Disease"
} else {
  target_group <- "Sepsis"
}

sam_bf <- rownames(group_info)[group_info$Group == target_group]
sam_bf <- intersect(sam_bf, colnames(exp_mat))
sam_bf <- intersect(sam_bf, rownames(cell_prop_clean))

# 1.2 Extract data
immuno_raw <- cell_prop_clean[sam_bf, ]
gene_expr  <- as.numeric(exp_mat["S100A8", sam_bf]) # Target gene

# ================= 2. Merge Monocytes (Total Monocytes) =================

# FIX: Convert matrix to data frame so $ and [[ ]] operators work correctly
immuno_merged <- as.data.frame(immuno_raw)

cols_hi <- grep("Hi", colnames(immuno_merged), value = T)
cols_lo <- grep("Lo", colnames(immuno_merged), value = T)
cols_hi <- cols_hi[grep("Mono", cols_hi)]
cols_lo <- cols_lo[grep("Mono", cols_lo)]

if(length(cols_hi) > 0 && length(cols_lo) > 0) {
  message(paste(">>> Merging subpopulations:", paste(cols_hi, collapse=", "), "+", paste(cols_lo, collapse=", ")))
  
  # Now that it's a dataframe, these operators will work perfectly
  immuno_merged$Total_Monocytes <- immuno_merged[[cols_hi]] + immuno_merged[[cols_lo]]
  
  # Remove the old columns
  immuno_merged[[cols_hi]] <- NULL
  immuno_merged[[cols_lo]] <- NULL
}

# Remove all-zero columns
immuno_merged <- immuno_merged[, colSums(immuno_merged) > 0]

# ================= 3. Calculate Correlation (Calculate separately to avoid errors) =================

message(">>> 1. Calculating left side heatmap data (Cell-Cell)...")
# correlate here returns a special object, used for direct calling by qcorrplot
cor_cells <- linkET::correlate(immuno_merged, method = "pearson")

message(">>> 2. Manually calculating right side link data (Gene-Cell)...")
# [Key Fix] Manually create a standard data.frame to avoid mutate errors
cells_to_test <- colnames(immuno_merged)
df_link_custom <- data.frame(
  from = "S100A8",       # Origin is gene
  to = cells_to_test,    # Destination is cell
  r = numeric(length(cells_to_test)),
  p = numeric(length(cells_to_test))
)

# Loop to calculate correlation
for(i in 1:nrow(df_link_custom)) {
  cell_vec <- immuno_merged[, df_link_custom$to[i]]
  # Use tryCatch to prevent errors caused by data anomalies
  res <- tryCatch(cor.test(gene_expr, cell_vec), error = function(e) NULL)
  
  if(!is.null(res)) {
    df_link_custom$r[i] <- res$estimate
    df_link_custom$p[i] <- res$p.value
  } else {
    df_link_custom$r[i] <- 0
    df_link_custom$p[i] <- 1
  }
}

# 3.3 Process link styling (df_link_custom is now a standard data.frame, mutate absolutely will not error)
df_link_plot <- df_link_custom %>%
  mutate(significance = case_when(
    p < 0.001 ~ "< 0.001",
    p < 0.01  ~ "< 0.01",
    p < 0.05  ~ "< 0.05",
    TRUE      ~ "ns"
  )) %>%
  filter(p < 0.05) # Keep only significant lines

# Check if there are lines
if(nrow(df_link_plot) == 0) {
  message("⚠️ Warning: No significant correlation, there will be no lines on the plot.")
  # If you want to force display an insignificant line to see the effect, you can comment out the filter above
} else {
  message(paste(">>> Found", nrow(df_link_plot), "significant links."))
}

## =========================================================================
## 🚑 Figure 5D Final Plot Fix: Replace geom_square with geom_tile
## =========================================================================
# =========================================================================
# 🚑 [Bridging code] Align data generated in previous steps to the variable names required by Figure 5D
# =========================================================================
cell_prop_clean <- res_plot
exp_mat <- bulk_mat

# Ensure group_info also exists (already read at the beginning of Step 5)
if(!exists("group_info")) {
  group_info <- read.csv("GSE65682_Groups.csv", header = T, row.names = 1, check.names = F)
}

# ================= 1. Data Preparation =================
library(linkET)
library(tidyverse)
library(RColorBrewer)

message(">>> [Step 5D-Final-Fix] Preparing merged butterfly plot data...")

# Lock onto samples
if("Disease" %in% group_info$Group) {
  target_group <- "Disease"
} else {
  target_group <- "Sepsis"
}

sam_bf <- rownames(group_info)[group_info$Group == target_group]
sam_bf <- intersect(sam_bf, colnames(exp_mat))
sam_bf <- intersect(sam_bf, rownames(cell_prop_clean))

# Extract data
immuno_raw <- cell_prop_clean[sam_bf, ]
gene_expr  <- as.numeric(exp_mat["S100A8", sam_bf]) # Target gene S100A8

message(">>> Data matching successful! Preparing for correlation calculations...")

## =========================================================================
## 💾 Save Figure 5D Statistical Results (Package Conflict Fix Version)
## =========================================================================

library(dplyr) # Confirm load again

message(">>> Exporting detailed statistical data for Figure 5D...")

# 1. Ensure data exists
if(exists("df_link_custom")) {
  
  # 2. Tidy table (Use dplyr:: prefix to prevent conflicts)
  export_table <- df_link_custom %>%
    dplyr::mutate(
      Direction = ifelse(r > 0, "Positive (+)", "Negative (-)"), # Clarify positive/negative correlation
      Correlation_r = round(r, 3),                               # Keep 3 decimal places
      P_value = formatC(p, format = "e", digits = 2),            # Scientific notation
      Significance = case_when(
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE      ~ "ns"
      )
    ) %>%
    # [Key Modification] Add dplyr:: here to force using dplyr's function
    dplyr::select(Gene = from, Cell_Type = to, Correlation_r, Direction, P_value, Significance) %>%
    dplyr::arrange(desc(Correlation_r)) # Sort descending by correlation
  
  # 3. Print preview
  print(head(export_table))
  
  # 4. Save to CSV
  write.csv(export_table, "Figure_5D_Correlation_Stats.csv", row.names = FALSE)
  
  message(">>> ✅ Results saved to: Figure_5D_Correlation_Stats.csv")
  
} else {
  stop("❌ Cannot find df_link_custom variable, please ensure you have completely run previous plotting code.")
}

# --- Fig 5E: Immune Interaction (Sepsis Only) ---
message("Plotting Fig 5E: Immune Interaction Heatmap")
# Assuming Sepsis label in grouping column is called "Disease" (adjust according to your data)
sepsis_sams <- rownames(group_plot)[group_plot$Group == "Disease"] 
if(length(sepsis_sams) == 0) sepsis_sams <- rownames(group_plot)[group_plot$Group == "Sepsis"]

if(length(sepsis_sams) > 5) {
  dat_sepsis <- res_plot[sepsis_sams, ]
  dat_sepsis <- dat_sepsis[, colSums(dat_sepsis) > 0] # Remove all-zero columns
  
  M <- cor(dat_sepsis)
  res_m <- cor.mtest(dat_sepsis)
  
  pdf("Figure_5E_Immune_Interaction_Corrplot.pdf", width = 7, height = 7)
  corrplot(M, p.mat = res_m$p, type = "upper", order = "hclust",
           sig.level = 0.05, insig = "blank",
           col = colorRampPalette(c("blue", "white", "red"))(200),
           tl.col = "black", title = "", mar=c(0,0,2,0))
  dev.off()
}

message(">>> 🎉 All complete! Please review the generated files.")
