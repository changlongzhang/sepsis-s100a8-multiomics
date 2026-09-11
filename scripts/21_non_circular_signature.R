# 输入：无监督状态、cluster markers 与原始 counts。
# 输出：删除 S100A8（及同时删除 S100A8/S100A9）的跨供者稳定独立签名、UCell评分和LODO稳健性结果。
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")

with_script_log("21_non_circular_signature", {
  suppressPackageStartupMessages({library(Seurat); library(UCell); library(ggplot2)})
  obj_file <- file.path("06_monocyte_state", "GSE167363_monocyte_state_annotated.rds")
  marker_file <- file.path("06_monocyte_state", "monocyte_cluster_markers_celllevel_descriptive.csv")
  assert_file(obj_file); assert_file(marker_file)
  mono <- readRDS(obj_file); DefaultAssay(mono) <- "RNA"; markers <- read.csv(marker_file, check.names=FALSE)
  ann <- read.csv(file.path("06_monocyte_state","monocyte_cluster_state_annotation.csv"),check.names=FALSE)
  genuine <- ann[ann$state_label != "Lymphoid-like contamination",]
  target_cluster <- genuine$cluster[which.max(genuine$combined_inflammatory_HLADRlow_score)]
  target_label <- genuine$state_label[genuine$cluster==target_cluster]
  target_cells <- colnames(mono)[as.character(mono$state_cluster)==as.character(target_cluster)]
  counts <- SeuratObject::LayerData(mono[["RNA"]],layer="counts")
  exclude_pattern <- "^(MT-|RPL|RPS|HB[ABDGEMQZ])"
  cand <- markers[as.character(markers$cluster)==as.character(target_cluster) & markers$p_val_adj<0.05 & markers$avg_log2FC>0.5 & markers$pct.1>=0.25,]
  cand <- cand[!grepl(exclude_pattern,cand$gene),]
  det <- sapply(unique(mono$donor_id),function(d){cc<-target_cells[mono$donor_id[match(target_cells,colnames(mono))]==d]; if(!length(cc)) rep(0,nrow(cand)) else Matrix::rowMeans(counts[cand$gene,cc,drop=FALSE]>0)})
  cand$donors_detected_ge10pct <- rowSums(det>=0.10)
  cand <- cand[cand$donors_detected_ge10pct>=ceiling(length(unique(mono$donor_id))/2),]
  cand <- cand[order(-cand$avg_log2FC),]
  sig_a8 <- head(setdiff(cand$gene,"S100A8"),30)
  sig_a8a9 <- head(setdiff(cand$gene,c("S100A8","S100A9")),30)
  antigen <- intersect(c("HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA"),rownames(mono))
  signature_table <- rbind(data.frame(signature="Independent inflammatory, no S100A8",gene=sig_a8),data.frame(signature="Independent inflammatory, no S100A8/S100A9",gene=sig_a8a9),data.frame(signature="Antigen presentation",gene=antigen))
  signature_table$target_state <- target_label
  write_csv_utf8(signature_table,file.path("07_non_circular_signature","independent_signature_genes.csv"))
  write_csv_utf8(cand,file.path("07_non_circular_signature","independent_signature_marker_audit.csv"))

  mono <- UCell::AddModuleScore_UCell(mono,features=list(Independent_NoA8=sig_a8,Independent_NoA8A9=sig_a8a9,Antigen=antigen),assay="RNA",slot="counts",name="_UCell",ncores=1)
  score_cols <- c("Independent_NoA8_UCell","Independent_NoA8A9_UCell","Antigen_UCell")
  mono$Combined_NoA8_HLADRlow <- mono$Independent_NoA8_UCell - mono$Antigen_UCell
  mono$Combined_NoA8A9_HLADRlow <- mono$Independent_NoA8A9_UCell - mono$Antigen_UCell
  score_cols <- c(score_cols,"Combined_NoA8_HLADRlow","Combined_NoA8A9_HLADRlow")
  donor <- aggregate(mono@meta.data[,score_cols],list(donor_id=mono$donor_id,group=mono$group),mean)
  write_csv_utf8(donor,file.path("07_non_circular_signature","donor_independent_signature_scores.csv"))

  # Leave-one-donor-out 是内部稳健性分析：cluster 本身已在全数据中确定，不能称为独立验证。
  dat <- SeuratObject::LayerData(mono[["RNA"]],layer="data")
  lodo <- list()
  for (held in unique(mono$donor_id)) {
    train <- colnames(mono)[mono$donor_id!=held & mono$state_label!="Lymphoid-like contamination"]
    target <- train[as.character(mono$state_cluster[match(train,colnames(mono))])==as.character(target_cluster)]
    other <- setdiff(train,target)
    mu1<-Matrix::rowMeans(dat[,target,drop=FALSE]); mu0<-Matrix::rowMeans(dat[,other,drop=FALSE]); lfc<-mu1-mu0
    dfrac<-Matrix::rowMeans(counts[,target,drop=FALSE]>0)
    genes<-names(sort(lfc[lfc>0.25 & dfrac>=0.1 & !grepl(exclude_pattern,names(lfc))],decreasing=TRUE))
    genes<-head(setdiff(genes,c("S100A8","S100A9")),30)
    hc<-colnames(mono)[mono$donor_id==held]
    sc<-if(length(genes)) colMeans(as.matrix(dat[genes,hc,drop=FALSE])) else rep(NA_real_,length(hc))
    lodo[[held]]<-data.frame(heldout_donor=held,group=unique(mono$group[mono$donor_id==held]),n_signature_genes=length(genes),heldout_mean_score=mean(sc),genes=paste(genes,collapse=";"),validation_scope="LODO internal robustness; clustering defined in all donors")
  }
  lodo<-do.call(rbind,lodo); write_csv_utf8(lodo,file.path("07_non_circular_signature","leave_one_donor_out_signature_robustness.csv"))
  saveRDS(mono,file.path("07_non_circular_signature","GSE167363_monocyte_independent_scores.rds"),compress="xz")
  long<-reshape(donor,varying=score_cols,v.names="score",timevar="signature",times=score_cols,direction="long");rownames(long)<-NULL
  long$signature<-factor(long$signature,levels=score_cols,labels=c("Independent program (no S100A8)","Independent program (no S100A8/S100A9)","Antigen presentation","Combined program (no S100A8)","Combined program (no S100A8/S100A9)"))
  levels(long$signature)<-c("Inflam. (−S100A8)","Inflam. (−S100A8/A9)","Antigen presentation","Combined (−S100A8)","Combined (−S100A8/A9)")
  long$donor_short<-sub("_.*$","",long$donor_id)
  p<-ggplot(long,aes(group,score,fill=group))+
    geom_boxplot(width=.42,outlier.shape=NA,colour="#5F6368",fill=NA,linewidth=.35)+
    geom_point(shape=21,size=2,position=position_jitter(width=.045),colour="white",stroke=.25)+
    facet_wrap(~signature,scales="free_y",ncol=3)+
    scale_fill_manual(values=c(Control=palette_nature[["Normal"]],Sepsis=palette_nature[["Sepsis"]]))+
    labs(x=NULL,y="Donor-mean UCell score",fill=NULL)+theme_nature()+
    theme(legend.position="top",legend.justification="left",strip.text=element_text(hjust=0))
  save_nature_plot(p,file.path("09_figures","Fig_IndependentSignature_DonorScores"),183,108,source_data=long)
})
