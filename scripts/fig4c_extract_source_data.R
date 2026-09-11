suppressPackageStartupMessages(library(Seurat))

project <- "D:/桌面/sepsis"
object_path <- file.path(
  project, "sepsis-s100a8-multiomics", "data",
  "09_scRNA-seq & Bulk Immune Infiltration", "output",
  "Sepsis_GSE167363_Prepared.rds"
)
output_path <- "C:/codex_sepsis_exchange/Plot_02_Annotation_DotPlot_source_data.csv"

markers <- c(
  "CD14", "VCAN", "FCN1", "LYZ", "S100A8", "S100A9", "MNDA",
  "FCGR3A", "LST1", "MS4A7", "CDKN1C", "CSF3R", "FCGR3B",
  "NAMPT", "CXCR2", "FCER1A", "CST3", "HLA-DRA", "HLA-DQA1",
  "CD3D", "CD3E", "IL7R", "TRAC", "CD8A", "CD8B", "GNLY",
  "NKG7", "KLRF1", "KLRD1", "MS4A1", "CD79A", "CD79B",
  "RALGPS2", "PPBP", "PF4", "GP1BB"
)

pbmc <- readRDS(object_path)
DefaultAssay(pbmc) <- "RNA"
dot_data <- DotPlot(pbmc, features = markers, group.by = "seurat_clusters")$data
write.csv(dot_data, output_path, row.names = FALSE)
