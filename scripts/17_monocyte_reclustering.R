# 输入：GSE167363 基线 Seurat 对象中的已注释 monocytes。
# 输出：不使用 S100A8 定义的无监督重聚类对象、分辨率审计、marker 与 UMAP。
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")
source(file.path("functions", "scrna_functions.R"), encoding = "UTF-8")

with_script_log("17_monocyte_reclustering", {
  suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(cluster)})
  in_file <- file.path(PROJECT_ROOT, "单细胞测序", "Sepsis_GSE167363_Final_Clean.rds")
  assert_file(in_file)
  full <- readRDS(in_file)
  cells <- rownames(full@meta.data)[as.character(full$Broad_CellType) == "Monocytes"]
  mono <- subset(full, cells = cells)
  mono$donor_id <- parse_gse167363_donor(mono$orig.ident)
  mono$group <- as.character(mono$Group)
  DefaultAssay(mono) <- "RNA"
  mono[["SCT"]] <- NULL
  mono <- NormalizeData(mono, verbose = FALSE)
  mono <- FindVariableFeatures(mono, selection.method = "vst", nfeatures = min(2000, nrow(mono)), verbose = FALSE)
  mono <- ScaleData(mono, features = rownames(mono), verbose = FALSE)
  mono <- RunPCA(mono, npcs = 30, verbose = FALSE)
  mono <- FindNeighbors(mono, reduction = "pca", dims = 1:20, verbose = FALSE)
  resolutions <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.8)
  mono <- FindClusters(mono, resolution = resolutions, random.seed = MASTER_SEED, verbose = FALSE)
  mono <- RunUMAP(mono, reduction = "pca", dims = 1:20, seed.use = MASTER_SEED, verbose = FALSE)

  pca <- Embeddings(mono, "pca")[, 1:10, drop = FALSE]
  audits <- list()
  for (res in resolutions) {
    col <- paste0("RNA_snn_res.", res)
    cl <- as.character(mono@meta.data[[col]])
    ct <- table(cl); donor_tab <- table(cl, mono$donor_id)
    donor_coverage <- rowSums(donor_tab > 0)
    max_donor_fraction <- apply(donor_tab, 1, max) / rowSums(donor_tab)
    sil <- if (length(unique(cl)) > 1L) mean(cluster::silhouette(as.integer(factor(cl)), stats::dist(pca))[, "sil_width"]) else NA_real_
    audits[[col]] <- data.frame(resolution = res, n_clusters = length(ct), min_cluster_cells = min(ct), median_cluster_cells = median(ct),
      median_donor_coverage = median(donor_coverage), max_single_donor_fraction = max(max_donor_fraction), mean_silhouette = sil,
      acceptable = min(ct) >= 20 && median(donor_coverage) >= 3 && max(max_donor_fraction) <= 0.8)
  }
  audit <- do.call(rbind, audits)
  eligible <- audit[audit$acceptable, ]
  if (nrow(eligible)) selected <- eligible$resolution[which.max(eligible$mean_silhouette)] else selected <- audit$resolution[which.max(audit$mean_silhouette)]
  selected_col <- paste0("RNA_snn_res.", selected)
  mono$state_cluster <- factor(mono@meta.data[[selected_col]])
  Idents(mono) <- "state_cluster"
  audit$selected <- audit$resolution == selected
  write_csv_utf8(audit, file.path("06_monocyte_state", "monocyte_resolution_audit.csv"))

  comp <- as.data.frame(table(cluster = mono$state_cluster, donor_id = mono$donor_id, group = mono$group), stringsAsFactors = FALSE)
  comp <- comp[comp$Freq > 0, ]; names(comp)[4] <- "n_cells"
  totals <- aggregate(n_cells ~ cluster, comp, sum); names(totals)[2] <- "cluster_total"
  comp <- merge(comp, totals, by = "cluster"); comp$donor_fraction <- comp$n_cells / comp$cluster_total
  write_csv_utf8(comp, file.path("06_monocyte_state", "monocyte_cluster_donor_composition.csv"))

  markers <- FindAllMarkers(mono, assay = "RNA", slot = "data", only.pos = TRUE, min.pct = 0.15, logfc.threshold = 0.25, test.use = "wilcox", verbose = FALSE)
  write_csv_utf8(markers, file.path("06_monocyte_state", "monocyte_cluster_markers_celllevel_descriptive.csv"))
  saveRDS(mono, file.path("06_monocyte_state", "GSE167363_monocyte_reclustered.rds"), compress = "xz")

  um_xy <- as.data.frame(Embeddings(mono, "umap"))
  names(um_xy)[1:2] <- c("UMAP_1", "UMAP_2")
  um <- data.frame(um_xy[, 1:2, drop = FALSE], cluster = mono$state_cluster, donor_id = mono$donor_id, group = mono$group, cell_id = colnames(mono))
  cluster_cols<-setNames(grDevices::hcl.colors(length(unique(um$cluster)),"Zissou 1"),sort(unique(um$cluster)))
  p1 <- ggplot(um, aes(UMAP_1, UMAP_2, colour = cluster)) + geom_point(size = 0.42, alpha = 0.72,stroke=0) +
    scale_colour_manual(values=cluster_cols)+labs(x = "UMAP 1", y = "UMAP 2", colour = "Cluster") + theme_umap()+
    guides(colour=guide_legend(override.aes=list(size=2,alpha=1)))
  save_nature_plot(p1, file.path("09_figures", "Fig_Monocyte_UnsupervisedClusters"), 89, 72, source_data = um)
  donor_cols<-setNames(grDevices::hcl.colors(length(unique(um$donor_id)),"Dark 3"),sort(unique(um$donor_id)))
  p2 <- ggplot(um, aes(UMAP_1, UMAP_2, colour = donor_id)) + geom_point(size = 0.35, alpha = 0.55,stroke=0) +
    scale_colour_manual(values=donor_cols)+labs(x = "UMAP 1", y = "UMAP 2", colour = "Donor") + theme_umap()+
    guides(colour=guide_legend(override.aes=list(size=1.8,alpha=1),ncol=2))
  save_nature_plot(p2, file.path("09_figures", "Fig_Monocyte_UMAP_ByDonor"), 89, 72, source_data = um)

  if (requireNamespace("clustree", quietly = TRUE)) {
    tryCatch({
      cp <- clustree::clustree(mono@meta.data, prefix = "RNA_snn_res.") +
        labs(colour="Resolution",size="Cells",edge_alpha="Incoming proportion") +
        theme_nature() + theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_blank(),legend.position="right")
      save_nature_plot(cp, file.path("09_figures", "Fig_Monocyte_Clustree"), 183, 100)
    }, error = function(e) {
      write_csv_utf8(data.frame(component = "clustree", status = "Unavailable because of package compatibility", detail = conditionMessage(e)),
                     file.path("06_monocyte_state", "clustree_compatibility_note.csv"))
    })
  }
})
