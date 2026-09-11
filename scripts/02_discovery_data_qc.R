# 输入：本地GSE65682 series matrix、GPL注释、分组表和临床表
# 输出：清洗后的log2基因矩阵RDS、样本QC、metadata、类别分布和PCA图
# 随机种子：20260718；不复用旧的线性尺度基因矩阵作为主分析输入。
rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "geo_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")
source(file.path("functions", "validation_functions.R"), encoding = "UTF-8")

with_script_log("02_discovery_data_qc", {
  stopifnot(requireNamespace("GEOquery", quietly=TRUE), requireNamespace("data.table", quietly=TRUE),
            requireNamespace("ggplot2", quietly=TRUE), requireNamespace("svglite", quietly=TRUE),
            requireNamespace("ragg", quietly=TRUE))
  series <- file.path(PROJECT_ROOT, "GSE65682", "GSE65682_series_matrix.txt.gz")
  gpl <- file.path(PROJECT_ROOT, "GSE65682", "GPL13667-15572.txt")
  group_file <- file.path(PROJECT_ROOT, "GSE65682", "GSE65682分组.csv")
  clinical_file <- file.path(PROJECT_ROOT, "GSE65682", "GSE65682_clinical.csv")
  lapply(c(series,gpl,group_file,clinical_file), assert_file)

  eset <- read_local_eset(series)
  cat("CHECKPOINT: series matrix loaded\n")
  probe_expr <- Biobase::exprs(eset)
  pdata <- Biobase::pData(eset)
  qx <- stats::quantile(probe_expr, c(0,.01,.25,.5,.75,.99,1), na.rm=TRUE)
  if (qx[6] > 20 || qx[7]-qx[1] > 50) stop("GSE65682 series matrix不像log2尺度；停止以避免错误转换", call.=FALSE)
  ann <- read_affy_annotation(gpl)
  cat("CHECKPOINT: annotation loaded", nrow(ann), "rows\n")
  gene_expr <- aggregate_probes_mean(probe_expr, ann)
  cat("CHECKPOINT: probes aggregated", nrow(gene_expr), "genes\n")
  if (!TARGET_GENE %in% rownames(gene_expr)) stop("聚合后未找到S100A8", call.=FALSE)

  groups <- data.table::fread(group_file, data.table=FALSE)
  clinical <- data.table::fread(clinical_file, data.table=FALSE, check.names=FALSE)
  names(groups)[tolower(names(groups)) == "sample"] <- "sample_id"
  names(groups)[tolower(names(groups)) == "group"] <- "group"
  groups$group <- ifelse(tolower(groups$group) %in% c("disease","sepsis","control"), "Sepsis",
                         ifelse(tolower(groups$group) %in% c("normal","healthy"), "Normal", NA))
  names(clinical)[names(clinical) == "geo_accession"] <- "sample_id"
  meta <- merge(groups, clinical, by="sample_id", all.x=TRUE, sort=FALSE)
  cat("CHECKPOINT: metadata merged", nrow(meta), "rows\n")
  meta <- meta[match(colnames(gene_expr), meta$sample_id), ]
  meta$group <- factor(meta$group, levels=c("Normal","Sepsis"))
  meta$all_characteristics <- combine_characteristics(meta)
  meta$label_supported_by_title <- (meta$group == "Normal" & grepl("healthy", meta$title, ignore.case=TRUE)) |
    (meta$group == "Sepsis" & grepl("intensive-care unit patient", meta$title, ignore.case=TRUE))
  validate_expression_metadata(gene_expr, meta, "sample_id")
  cat("CHECKPOINT: expression-metadata validation passed\n")
  if (anyNA(meta$group)) stop("存在缺失分组", call.=FALSE)

  s100_probes <- ann$probe_id[which(!is.na(ann$gene_symbol) & ann$gene_symbol == TARGET_GENE &
                                      ann$probe_id %in% rownames(probe_expr))]
  sample_qc <- data.frame(
    sample_id=colnames(gene_expr), group=meta$group,
    array_median=apply(probe_expr,2,median,na.rm=TRUE), array_IQR=apply(probe_expr,2,IQR,na.rm=TRUE),
    missing_fraction=colMeans(!is.finite(probe_expr)), S100A8=as.numeric(gene_expr[TARGET_GENE,]),
    duplicate_sample_id=duplicated(colnames(gene_expr)), label_supported_by_title=meta$label_supported_by_title
  )
  write_csv_utf8(sample_qc, file.path("02_class_imbalance","discovery_sample_qc.csv"))
  cat("CHECKPOINT: sample QC written\n")
  write_csv_utf8(meta, file.path("02_class_imbalance","sample_metadata_clean.csv"))
  cat("CHECKPOINT: metadata written\n")
  saveRDS(gene_expr, file.path("02_class_imbalance","GSE65682_gene_log2_mean_probes.rds"), compress=FALSE)
  cat("CHECKPOINT: gene matrix RDS written\n")
  saveRDS(list(probe_expr=probe_expr[s100_probes,,drop=FALSE], annotation=ann[ann$probe_id %in% s100_probes,]),
          file.path("02_class_imbalance","S100A8_probe_audit.rds"))
  cat("CHECKPOINT: S100A8 probe audit written\n")

  class_df <- as.data.frame(table(meta$group)); names(class_df) <- c("Group","n")
  p_class <- ggplot2::ggplot(class_df, ggplot2::aes(Group,n,fill=Group)) +
    ggplot2::geom_col(width=.65,colour="black",linewidth=.35) +
    ggplot2::geom_text(ggplot2::aes(label=n),vjust=-.25,size=2.7) +
    ggplot2::scale_fill_manual(values=palette_nature[c("Normal","Sepsis")]) +
    ggplot2::labs(x=NULL,y="Number of participants") + theme_nature() + ggplot2::theme(legend.position="none")
  cat("CHECKPOINT: class plot object created\n")
  save_nature_plot(p_class,file.path("02_class_imbalance","class_distribution"),89,70,source_data=class_df)
  cat("CHECKPOINT: class plot exported\n")

  vars <- matrixStats::rowVars(gene_expr)
  top <- order(vars,decreasing=TRUE)[seq_len(min(2000,length(vars)))]
  pca <- stats::prcomp(t(gene_expr[top,,drop=FALSE]),center=TRUE,scale.=TRUE)
  cat("CHECKPOINT: PCA completed\n")
  pca_df <- data.frame(sample_id=rownames(pca$x),group=meta$group,PC1=pca$x[,1],PC2=pca$x[,2])
  ve <- 100*pca$sdev^2/sum(pca$sdev^2)
  p_pca <- ggplot2::ggplot(pca_df,ggplot2::aes(PC1,PC2,colour=group))+
    ggplot2::geom_point(size=1.5,alpha=.75)+ggplot2::scale_colour_manual(values=palette_nature[c("Normal","Sepsis")])+
    ggplot2::labs(x=sprintf("PC1 (%.1f%%)",ve[1]),y=sprintf("PC2 (%.1f%%)",ve[2]),colour=NULL)+theme_nature()
  save_nature_plot(p_pca,file.path("02_class_imbalance","PCA_before_balancing"),89,75,source_data=pca_df)

  qc_summary <- data.frame(
    item=c("probe_count","gene_count","sample_count","sepsis_count","normal_count","duplicate_sample_ids","missing_group","log2_scale_confirmed","S100A8_probe_count","S100A8_gene_present"),
    value=c(nrow(probe_expr),nrow(gene_expr),ncol(gene_expr),sum(meta$group=="Sepsis"),sum(meta$group=="Normal"),anyDuplicated(colnames(gene_expr)),sum(is.na(meta$group)),TRUE,length(s100_probes),TARGET_GENE %in% rownames(gene_expr))
  )
  write_csv_utf8(qc_summary,file.path("02_class_imbalance","discovery_data_summary.csv"))
  cat("完成：",nrow(gene_expr),"genes x",ncol(gene_expr),"samples; S100A8 probes:",paste(s100_probes,collapse=", "),"\n")
})
