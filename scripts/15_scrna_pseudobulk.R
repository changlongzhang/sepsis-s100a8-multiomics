# 输入：GSE167363 原始 RNA counts 与脚本14核查后的 donor 元数据。
# 输出：donor×monocyte pseudo-bulk edgeR QL 结果；20-cell主分析和50-cell严格敏感性分析。
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")
source(file.path("functions", "scrna_functions.R"), encoding = "UTF-8")

with_script_log("15_scrna_pseudobulk", {
  suppressPackageStartupMessages({library(Seurat); library(edgeR); library(ggplot2); library(Matrix)})
  in_file <- file.path(PROJECT_ROOT, "单细胞测序", "Sepsis_GSE167363_Final_Clean.rds")
  donor_file <- file.path("05_single_cell_pseudobulk", "donor_metadata_clean.csv")
  assert_file(in_file); assert_file(donor_file)
  obj <- readRDS(in_file); md <- obj@meta.data
  md$donor_id <- parse_gse167363_donor(md$orig.ident); md$group <- as.character(md$Group)
  mono_cells <- rownames(md)[as.character(md$Broad_CellType) == "Monocytes"]
  if (length(mono_cells) < 100L) stop("已注释单核细胞过少，无法进行计划分析。")
  counts <- get_counts_layer(obj, "RNA")[, mono_cells, drop = FALSE]
  cell_md <- md[mono_cells, , drop = FALSE]
  agg <- aggregate_counts_by_donor(counts, cell_md$donor_id)
  ncell <- table(cell_md$donor_id)
  donor_group <- unique(cell_md[, c("donor_id", "group")])
  donor_group <- donor_group[match(colnames(agg), donor_group$donor_id), ]
  donor_group$n_monocytes <- as.integer(ncell[donor_group$donor_id])
  if (anyNA(donor_group$group)) stop("聚合 counts 与 donor metadata 未完全匹配。")
  saveRDS(agg, file.path("05_single_cell_pseudobulk", "monocyte_donor_pseudobulk_counts.rds"))
  write_csv_utf8(donor_group, file.path("05_single_cell_pseudobulk", "monocyte_pseudobulk_sample_metadata.csv"))

  thresholds <- c(primary_20 = 20L, strict_50 = 50L)
  focus <- c("S100A8","S100A9","S100A12","IL1B","FCN1","CTSS","NFKBIA","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA","IL1R2","RETN","LILRB1","LILRB3","SOCS3","ARG1","MPO","ELANE")
  all_focus <- list(); feasibility <- list()
  for (nm in names(thresholds)) {
    keep <- donor_group$n_monocytes >= thresholds[[nm]]
    sm <- donor_group[keep, , drop = FALSE]
    feasibility[[nm]] <- data.frame(analysis = nm, min_cells = thresholds[[nm]], n_donors = nrow(sm), n_control = sum(sm$group == "Control"), n_sepsis = sum(sm$group == "Sepsis"), estimable = length(unique(sm$group)) == 2 && all(table(sm$group) >= 1))
    if (!feasibility[[nm]]$estimable || nrow(sm) < 3L) next
    res <- run_edger_pseudobulk(agg[, sm$donor_id, drop = FALSE], sm, min_cpm_samples = 2L)
    write_csv_utf8(res$table, file.path("05_single_cell_pseudobulk", paste0("pseudobulk_edgeR_", nm, "_all_genes.csv")))
    ft <- merge(data.frame(gene = focus), res$table, by = "gene", all.x = TRUE)
    ft$analysis <- nm; ft$min_cells <- thresholds[[nm]]
    all_focus[[nm]] <- ft
    cpm_mat <- edgeR::cpm(res$dge, log = TRUE, prior.count = 2)
    cpm_out <- data.frame(gene = rownames(cpm_mat), cpm_mat, check.names = FALSE)
    write_csv_utf8(cpm_out, file.path("05_single_cell_pseudobulk", paste0("pseudobulk_logCPM_", nm, ".csv")))
    if (nm == "primary_20" && TARGET_GENE %in% rownames(cpm_mat)) {
      pd <- data.frame(donor_id = colnames(cpm_mat), group = sm$group[match(colnames(cpm_mat), sm$donor_id)], logCPM = as.numeric(cpm_mat[TARGET_GENE, ]), n_monocytes = sm$n_monocytes[match(colnames(cpm_mat), sm$donor_id)])
      write_csv_utf8(pd, file.path("05_single_cell_pseudobulk", "S100A8_donor_logCPM.csv"))
      pd$donor_label<-sub("_.*$","",pd$donor_id)
      p <- ggplot(pd, aes(x = group, y = logCPM, colour = group)) +
        geom_boxplot(width=.42,outlier.shape=NA,colour="#5F6368",fill=NA,linewidth=.4) +
        geom_point(size = 2.2, position = position_jitter(width = 0.035, height = 0)) +
        ggrepel::geom_text_repel(aes(label=donor_label),size=2.1,show.legend=FALSE,max.overlaps=Inf,box.padding=.18,point.padding=.15,segment.colour="#B0B0B0",segment.size=.25) +
        scale_colour_manual(values = c(Control = palette_nature[["Normal"]], Sepsis = palette_nature[["Sepsis"]])) +
        labs(x = NULL, y = expression(S100A8~"donor pseudo-bulk logCPM"), colour = NULL) + theme_nature() + theme(legend.position = "none")
      save_nature_plot(p, file.path("09_figures", "Fig_Pseudobulk_S100A8_Donor"), 89, 70, source_data = pd)
    }
  }
  write_csv_utf8(do.call(rbind, feasibility), file.path("05_single_cell_pseudobulk", "pseudobulk_threshold_feasibility.csv"))
  if (length(all_focus)) write_csv_utf8(do.call(rbind, all_focus), file.path("05_single_cell_pseudobulk", "Table_Pseudobulk_FocusGenes.csv"))
})
