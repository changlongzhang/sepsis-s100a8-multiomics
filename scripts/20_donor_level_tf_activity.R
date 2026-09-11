# 输入：GSE167363 单核细胞表达与 DoRothEA 调控网络。
# 输出：RELA/NFKB1/STAT3/CEBPB/JUN/FOS 的供者级 ULM 活性和精确置换统计。
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")

with_script_log("20_donor_level_tf_activity", {
  suppressPackageStartupMessages({library(Seurat); library(decoupleR); library(dorothea); library(ggplot2)})
  in_file <- file.path("06_monocyte_state", "GSE167363_monocyte_state_annotated.rds")
  assert_file(in_file); mono <- readRDS(in_file); DefaultAssay(mono) <- "RNA"
  data(dorothea_hs, package = "dorothea")
  tfs <- c("RELA","NFKB1","STAT3","CEBPB","JUN","FOS")
  net <- dorothea_hs[dorothea_hs$tf %in% tfs & dorothea_hs$confidence %in% c("A","B","C"), c("tf","target","mor")]
  mat <- SeuratObject::LayerData(mono[["RNA"]], layer = "data")
  net <- net[net$target %in% rownames(mat), ]
  act <- decoupleR::run_ulm(mat = mat, network = net, .source = tf, .target = target, .mor = mor, minsize = 5)
  names(act)[names(act) == "condition"] <- "cell_id"; names(act)[names(act) == "source"] <- "tf"
  score_col <- intersect(c("score","statistic"), names(act))[1]
  if (is.na(score_col)) stop("decoupleR 输出缺少 activity score。")
  act$donor_id <- mono$donor_id[match(act$cell_id, colnames(mono))]
  act$group <- mono$group[match(act$cell_id, colnames(mono))]
  donor <- aggregate(act[[score_col]], list(donor_id=act$donor_id,group=act$group,tf=act$tf), mean, na.rm=TRUE)
  names(donor)[4] <- "activity"
  stats <- do.call(rbind, lapply(split(donor, donor$tf), function(d) {
    obs <- mean(d$activity[d$group=="Sepsis"]) - mean(d$activity[d$group=="Control"])
    cmb <- combn(seq_len(nrow(d)), sum(d$group=="Control")); null <- apply(cmb,2,function(ix) mean(d$activity[-ix])-mean(d$activity[ix]))
    b <- replicate(2000, mean(sample(d$activity[d$group=="Sepsis"],replace=TRUE))-mean(sample(d$activity[d$group=="Control"],replace=TRUE)))
    data.frame(tf=d$tf[1],n_control=sum(d$group=="Control"),n_sepsis=sum(d$group=="Sepsis"),effect=obs,ci_low=quantile(b,.025),ci_high=quantile(b,.975),permutation_p=mean(abs(null)>=abs(obs)-1e-12))
  })); stats$fdr <- p.adjust(stats$permutation_p,"BH")
  write_csv_utf8(donor, file.path("05_single_cell_pseudobulk", "donor_level_TF_activity.csv"))
  write_csv_utf8(stats, file.path("05_single_cell_pseudobulk", "Table_DonorLevel_TFActivity.csv"))
  stats$tf<-factor(stats$tf,levels=stats$tf[order(stats$effect)])
  p <- ggplot(stats,aes(effect,tf))+
    geom_vline(xintercept=0,linetype=2,colour="#8C8C8C",linewidth=.35)+
    geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=0,linewidth=.55,colour="#5F6368")+
    geom_point(aes(fill=fdr<.05),shape=21,size=2.5,colour="white",stroke=.35)+
    scale_fill_manual(values=c(`TRUE`="#C65A46",`FALSE`="#9B9B9B"),guide="none")+
    labs(x="Sepsis − control difference\nin donor-mean TF activity (95% bootstrap CI)",y=NULL)+theme_forest()
  save_nature_plot(p,file.path("09_figures","Fig_DonorLevel_TFActivity"),120,72,source_data=stats)
})
