# 输入：本地 GSE167363 基线 Seurat 对象。
# 输出：核查后的 donor/sample 元数据、细胞数、细胞类型数及供者级质控图。
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")
source(file.path("functions", "scrna_functions.R"), encoding = "UTF-8")

with_script_log("14_scrna_metadata_qc", {
  suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
  in_file <- file.path(PROJECT_ROOT, "单细胞测序", "Sepsis_GSE167363_Final_Clean.rds")
  assert_file(in_file, "GSE167363 Seurat")
  obj <- readRDS(in_file)
  required <- c("orig.ident", "Group", "Timepoint", "Broad_CellType")
  if (!all(required %in% names(obj@meta.data))) stop("Seurat 元数据缺列: ", paste(setdiff(required, names(obj@meta.data)), collapse = ", "))
  md <- obj@meta.data
  md$cell_id <- rownames(md)
  md$sample_id <- as.character(md$orig.ident)
  md$donor_id <- parse_gse167363_donor(md$sample_id)
  md$derived_time <- parse_gse167363_time(md$sample_id)
  md$group <- as.character(md$Group)
  md$cell_type <- as.character(md$Broad_CellType)
  if (any(is.na(md$donor_id) | md$donor_id == "")) stop("存在无法解析的 donor ID。")
  donor_sample <- unique(md[, c("donor_id", "sample_id", "group", "Timepoint", "derived_time")])
  donor_sample$donor_sample_consistent <- donor_sample$Timepoint == donor_sample$derived_time
  if (any(duplicated(donor_sample[, c("donor_id", "Timepoint")]))) stop("同一 donor/timepoint 对应多个样本，需人工核查。")
  if (any(!donor_sample$donor_sample_consistent)) stop("样本名与 Timepoint 元数据不一致。")

  donor_counts <- as.data.frame(table(md$donor_id, md$group), stringsAsFactors = FALSE)
  names(donor_counts) <- c("donor_id", "group", "n_cells")
  donor_counts <- donor_counts[donor_counts$n_cells > 0, ]
  donor_ct <- as.data.frame(table(md$donor_id, md$group, md$cell_type), stringsAsFactors = FALSE)
  names(donor_ct) <- c("donor_id", "group", "cell_type", "n_cells")
  donor_ct <- donor_ct[donor_ct$n_cells > 0, ]
  mono <- donor_ct[donor_ct$cell_type == "Monocytes", c("donor_id", "n_cells")]
  names(mono)[2] <- "n_monocytes"
  donor_meta <- merge(donor_sample, mono, by = "donor_id", all.x = TRUE)
  donor_meta$n_monocytes[is.na(donor_meta$n_monocytes)] <- 0L
  donor_meta$eligible_20 <- donor_meta$n_monocytes >= 20
  donor_meta$eligible_50 <- donor_meta$n_monocytes >= 50
  donor_meta$orig_ident_is_donor <- donor_meta$sample_id == donor_meta$donor_id
  donor_meta$note <- ifelse(donor_meta$orig_ident_is_donor, "orig.ident equals donor", "orig.ident is sample; donor parsed and verified from sample name")

  write_csv_utf8(donor_meta, file.path("05_single_cell_pseudobulk", "donor_metadata_clean.csv"))
  write_csv_utf8(donor_counts, file.path("05_single_cell_pseudobulk", "donor_cell_counts.csv"))
  write_csv_utf8(donor_ct, file.path("05_single_cell_pseudobulk", "donor_celltype_counts.csv"))

  p1 <- ggplot(donor_counts, aes(x = n_cells, y = reorder(donor_id, n_cells), colour = group)) +
    geom_segment(aes(x = 0, xend = n_cells, yend = reorder(donor_id, n_cells)), colour = "#D8D8D8", linewidth = .55) +
    geom_point(size = 2.3) + scale_colour_manual(values = c(Control = palette_nature[["Normal"]], Sepsis = palette_nature[["Sepsis"]])) +
    labs(x = "Cells per donor", y = NULL, colour = NULL) + theme_nature() + theme(legend.position="top",legend.justification="left")
  save_nature_plot(p1, file.path("09_figures", "Fig_Donor_CellNumber"), 89, 70, source_data = donor_counts)

  umap_cols <- intersect(c("umap_1", "umap_2"), names(md))
  if (length(umap_cols) != 2L) stop("对象缺少已保存的 UMAP 坐标。")
  um <- md[, c("cell_id", "donor_id", "group", umap_cols)]
  donor_cols <- setNames(grDevices::hcl.colors(length(unique(um$donor_id)), "Dark 3"), sort(unique(um$donor_id)))
  p2 <- ggplot(um, aes(x = umap_1, y = umap_2, colour = donor_id)) + geom_point(size = 0.16, alpha = 0.5, stroke=0) +
    scale_colour_manual(values=donor_cols)+labs(x = "UMAP 1", y = "UMAP 2", colour = "Donor") + theme_umap() +
    guides(colour=guide_legend(override.aes=list(size=1.7,alpha=1),ncol=2))+
    theme(legend.key.height = grid::unit(2.5, "mm"),legend.position="right")
  save_nature_plot(p2, file.path("09_figures", "Fig_Donor_UMAP"), 89, 78, source_data = um)

  totals <- aggregate(n_cells ~ donor_id, donor_ct, sum)
  comp <- merge(donor_ct, totals, by = "donor_id"); comp$proportion <- comp$n_cells.x / comp$n_cells.y
  comp$group <- donor_meta$group[match(comp$donor_id, donor_meta$donor_id)]
  cell_cols <- setNames(grDevices::hcl.colors(length(unique(comp$cell_type)), "Oslo", rev=TRUE), sort(unique(comp$cell_type)))
  p3 <- ggplot(comp, aes(x = donor_id, y = proportion, fill = cell_type)) + geom_col(width = 0.88,colour="white",linewidth=.08) +
    facet_grid(~ group, scales = "free_x", space = "free_x") + scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values=cell_cols)+labs(x = NULL, y = "Cell-type proportion", fill = "Cell type") + theme_nature() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position="right")
  save_nature_plot(p3, file.path("09_figures", "Fig_Donor_Composition"), 183, 75, source_data = comp)

  audit <- data.frame(metric = c("donors_total", "control_donors", "sepsis_donors", "cells_total", "monocytes_total", "donors_mono_ge20", "donors_mono_ge50", "orig_ident_equals_donor"),
                      value = c(nrow(donor_meta), sum(donor_meta$group == "Control"), sum(donor_meta$group == "Sepsis"), nrow(md), sum(md$cell_type == "Monocytes"), sum(donor_meta$eligible_20), sum(donor_meta$eligible_50), all(donor_meta$orig_ident_is_donor)))
  write_csv_utf8(audit, file.path("05_single_cell_pseudobulk", "scrna_metadata_audit_summary.csv"))
})
