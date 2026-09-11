# 输入：脚本17无监督重聚类对象。
# 输出：多基因状态评分、数据驱动状态注释以及供者级状态评分。
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")

with_script_log("18_monocyte_state_annotation", {
  suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
  in_file <- file.path("06_monocyte_state", "GSE167363_monocyte_reclustered.rds")
  assert_file(in_file); mono <- readRDS(in_file); DefaultAssay(mono) <- "RNA"
  sets <- list(
    inflammatory = c("S100A8","S100A9","S100A12","IL1B","FCN1","CTSS","NFKBIA","CEBPB"),
    antigen_presentation = c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA"),
    immunosuppression = c("IL1R2","RETN","LILRB1","LILRB3","SOCS3"),
    immature_emergency = c("MPO","ELANE","AZU1","CTSG","LTF","OLFM4","LCN2","MMP8","MMP9"),
    monocyte_core = c("LST1","FCER1G","CTSS","LILRB1","CTSD","TYROBP"),
    classical = c("VCAN","FCN1","CTSS","LILRA1","SLC11A1"),
    nonclassical = c("FCGR3A","MS4A7","LST1","IFITM3"),
    lymphoid_contamination = c("CD3D","TRAC","CD79A","MS4A1","NKG7","GNLY")
  )
  sets <- lapply(sets, intersect, y = rownames(mono))
  dat <- as.matrix(SeuratObject::LayerData(mono[["RNA"]], layer = "data"))
  score <- sapply(sets, function(g) if (length(g)) colMeans(dat[g, , drop = FALSE]) else rep(NA_real_, ncol(dat)))
  score <- scale(score)
  for (nm in colnames(score)) mono[[paste0(nm, "_score")]] <- score[, nm]
  mono$combined_inflammatory_HLADRlow_score <- mono$inflammatory_score - mono$antigen_presentation_score

  score_cols <- c(paste0(names(sets), "_score"), "combined_inflammatory_HLADRlow_score")
  cell_scores <- data.frame(cell_id = colnames(mono), donor_id = mono$donor_id, group = mono$group,
                            cluster = as.character(mono$state_cluster), mono@meta.data[, score_cols, drop = FALSE], check.names = FALSE)
  cluster_scores <- aggregate(cell_scores[, score_cols], list(cluster = cell_scores$cluster), mean, na.rm = TRUE)
  cluster_scores$n_cells <- as.integer(table(cell_scores$cluster)[cluster_scores$cluster])

  # 先识别明显的淋巴细胞污染，再依据未使用 S100A8 分组的多基因程序命名。
  cluster_scores$state_label <- NA_character_
  contam <- cluster_scores$lymphoid_contamination_score > cluster_scores$monocyte_core_score + 0.25
  cluster_scores$state_label[contam] <- "Lymphoid-like contamination"
  genuine <- which(!contam)
  if (length(genuine)) {
    em <- genuine[which.max(cluster_scores$immature_emergency_score[genuine])]
    if (cluster_scores$immature_emergency_score[em] > 0) cluster_scores$state_label[em] <- "Emergency myeloid / immature neutrophil-like"
  }
  left <- which(is.na(cluster_scores$state_label))
  for (i in left) {
    infl_high <- cluster_scores$inflammatory_score[i] > median(cluster_scores$inflammatory_score[genuine])
    hla_low <- cluster_scores$antigen_presentation_score[i] < median(cluster_scores$antigen_presentation_score[genuine])
    suppress_high <- cluster_scores$immunosuppression_score[i] > median(cluster_scores$immunosuppression_score[genuine])
    cluster_scores$state_label[i] <- if (infl_high && hla_low && suppress_high) "Inflammatory / HLA-DR-low monocyte" else if (cluster_scores$nonclassical_score[i] > cluster_scores$classical_score[i]) "FCGR3A+ non-classical monocyte" else "VCAN+ classical monocyte"
  }
  label_map <- setNames(cluster_scores$state_label, cluster_scores$cluster)
  mono$state_label <- unname(label_map[as.character(mono$state_cluster)])
  cell_scores$state_label <- mono$state_label

  donor_scores <- aggregate(cell_scores[, score_cols], list(donor_id = cell_scores$donor_id, group = cell_scores$group), mean, na.rm = TRUE)
  write_csv_utf8(cell_scores, file.path("06_monocyte_state", "monocyte_cell_state_scores.csv"))
  write_csv_utf8(cluster_scores, file.path("06_monocyte_state", "monocyte_cluster_state_annotation.csv"))
  write_csv_utf8(donor_scores, file.path("06_monocyte_state", "monocyte_donor_state_scores.csv"))
  saveRDS(mono, file.path("06_monocyte_state", "GSE167363_monocyte_state_annotated.rds"), compress = "xz")

  long <- do.call(rbind,lapply(score_cols,function(sc)data.frame(cluster=cluster_scores$cluster,state_label=cluster_scores$state_label,program=sub("_score$","",sc),mean_z_score=cluster_scores[[sc]],stringsAsFactors=FALSE)))
  rownames(long) <- NULL
  program_order<-c("classical","nonclassical","monocyte_core","immature_emergency","inflammatory","antigen_presentation","immunosuppression","lymphoid_contamination","combined_inflammatory_HLADRlow")
  long$program<-factor(long$program,levels=program_order)
  long$cluster_label<-paste0("Cluster ",long$cluster,"  ·  ",long$state_label)
  lim<-max(abs(long$mean_z_score),na.rm=TRUE)
  p <- ggplot(long, aes(x = program, y = cluster_label, fill = mean_z_score)) + geom_tile(colour = "white", linewidth = 0.55) +
    scale_fill_gradient2(low = "#315F7D", mid = "#F7F7F7", high = "#C65A46", midpoint = 0,limits=c(-lim,lim),oob=scales::squish) +
    scale_x_discrete(labels=c(classical="Classical",nonclassical="Non-classical",monocyte_core="Monocyte core",immature_emergency="Immature",inflammatory="Inflammatory",antigen_presentation="Antigen",immunosuppression="Suppression",lymphoid_contamination="Lymphoid",combined_inflammatory_HLADRlow="Inflam. − antigen"))+
    labs(x = NULL, y = NULL, fill = "Program\nz-score") + theme_heatmap() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),legend.position="right")
  save_nature_plot(p, file.path("09_figures", "Fig_Monocyte_StatePrograms"), 183, 85, source_data = long)
})
