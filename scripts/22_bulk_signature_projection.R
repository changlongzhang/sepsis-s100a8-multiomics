# 输入：删除 S100A8/S100A9 的独立单细胞签名和各 bulk 队列表达矩阵。
# 输出：ssGSEA 投射、组间差异、S100A8相关及控制疾病组/单核细胞背景后的偏相关。
rm(list=ls()); options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","geo_functions.R"),encoding="UTF-8")
source(file.path("functions","theme_nature.R"),encoding="UTF-8")

map_probe_matrix <- function(mat,pkg,keytype="PROBEID") {
  if (pkg == "hgu133plus2.db") {
    smap <- get("hgu133plus2SYMBOL", envir = asNamespace(pkg), inherits = FALSE)
    vals <- mget(rownames(mat), smap, ifnotfound = NA)
    sym <- vapply(vals, function(z) if (length(z) && !all(is.na(z))) as.character(z[[1]]) else NA_character_, character(1))
  } else {
    db <- get(pkg,envir=asNamespace(pkg),inherits=FALSE)
    sym <- AnnotationDbi::mapIds(db,keys=rownames(mat),column="SYMBOL",keytype=keytype,multiVals="first")
  }
  ann <- data.frame(probe_id=names(sym),gene_symbol=unname(sym),stringsAsFactors=FALSE)
  aggregate_probes_mean(mat,ann)
}
ssgsea_scores <- function(mat,sets) {
  sets <- lapply(sets,intersect,y=rownames(mat)); sets <- sets[lengths(sets)>=3]
  as.matrix(GSVA::gsva(GSVA::ssgseaParam(mat,sets,minSize=3,normalize=TRUE),verbose=FALSE))
}
partial_spearman <- function(x,y,covars) {
  rx <- residuals(lm(rank(x,na.last="keep")~.,data=covars)); ry <- residuals(lm(rank(y,na.last="keep")~.,data=covars))
  unname(cor.test(rx,ry,method="pearson")$estimate)
}

with_script_log("22_bulk_signature_projection", {
  suppressPackageStartupMessages({library(GSVA);library(AnnotationDbi);library(ggplot2);library(pROC)})
  sig <- read.csv(file.path("07_non_circular_signature","independent_signature_genes.csv"),check.names=FALSE)
  sets <- split(sig$gene,sig$signature)
  sets$Monocyte_core <- c("LST1","FCER1G","TYROBP","CTSS","LILRB1","CTSD")
  sample_group <- read.csv(file.path("03_clinical_validation","clinical_validation_S100A8_samples.csv"),check.names=FALSE)
  mats <- list(); groups <- list()
  mats$GSE65682 <- readRDS(file.path("02_class_imbalance","GSE65682_gene_log2_mean_probes.rds"))
  md656 <- read.csv(file.path("02_class_imbalance","sample_metadata_clean.csv"),check.names=FALSE)
  groups$GSE65682 <- setNames(ifelse(md656$group=="Normal","Healthy",md656$group),md656$sample_id)
  for(gse in EXTERNAL_GSE) {
    f<-file.path(PROJECT_ROOT,"sepsis-s100a8-multiomics","data","08_External Validation","Input",paste0(gse,"_series_matrix.txt.gz")); es<-read_local_eset(f)
    pkg<-if(Biobase::annotation(es)=="GPL570")"hgu133plus2.db" else "illuminaHumanv4.db"
    mats[[gse]]<-map_probe_matrix(Biobase::exprs(es),pkg)
    d<-sample_group[sample_group$dataset==gse,]; groups[[gse]]<-setNames(d$group,d$sample_id)
  }
  g134<-readRDS(file.path("01_data_inventory","downloaded_data","GSE134347_eset.rds")); mats$GSE134347<-map_probe_matrix(Biobase::exprs(g134),"hta20transcriptcluster.db")
  d<-sample_group[sample_group$dataset=="GSE134347",];groups$GSE134347<-setNames(d$group,d$sample_id)
  cnt<-openxlsx::read.xlsx(file.path("01_data_inventory","downloaded_data","GSE131411_rawcounts_CS_SS.xlsx")); ens<-sub("\\..*$","",as.character(cnt[[1]]));cm<-as.matrix(cnt[,-1]);rownames(cm)<-ens;storage.mode(cm)<-"numeric";cm<-cm[,grepl("T1$",colnames(cm)),drop=FALSE]
  dge<-edgeR::DGEList(cm);keep<-edgeR::filterByExpr(dge);dge<-edgeR::calcNormFactors(dge[keep,,keep.lib.sizes=FALSE]);lcpm<-edgeR::cpm(dge,log=TRUE,prior.count=2)
  mats$GSE131411<-map_probe_matrix(lcpm,"org.Hs.eg.db",keytype="ENSEMBL");d<-sample_group[sample_group$dataset=="GSE131411",];groups$GSE131411<-setNames(d$group,d$sample_id)

  samples<-list(); comparisons<-list(); cors<-list()
  primary_name<-"Independent inflammatory, no S100A8/S100A9"; antigen_name<-"Antigen presentation"
  for(ds in names(mats)) {
    mat<-mats[[ds]]; sc<-ssgsea_scores(mat,sets); ids<-intersect(colnames(mat),names(groups[[ds]])); grp<-unname(groups[[ds]][ids]); keep<-!is.na(grp)&nzchar(grp)&grp!="Unresolved";ids<-ids[keep];grp<-grp[keep]
    if(!length(ids)||!TARGET_GENE%in%rownames(mat)||!all(c(primary_name,antigen_name,"Monocyte_core")%in%rownames(sc)))next
    dd<-data.frame(dataset=ds,sample_id=ids,group=grp,S100A8=as.numeric(mat[TARGET_GENE,ids]),independent_score=as.numeric(sc[primary_name,ids]),antigen_score=as.numeric(sc[antigen_name,ids]),monocyte_core_score=as.numeric(sc["Monocyte_core",ids]))
    dd$combined_score<-dd$independent_score-dd$antigen_score;samples[[ds]]<-dd
    for(cp in setdiff(unique(dd$group),"Sepsis")) {
      z<-dd[dd$group%in%c("Sepsis",cp),]; if(min(table(z$group))<5)next
      roc<-pROC::roc(factor(z$group,levels=c(cp,"Sepsis")),z$combined_score,quiet=TRUE,direction="<")
      comparisons[[paste(ds,cp)]]<-data.frame(dataset=ds,comparator=cp,n_sepsis=sum(z$group=="Sepsis"),n_comparator=sum(z$group==cp),mean_difference=mean(z$combined_score[z$group=="Sepsis"])-mean(z$combined_score[z$group==cp]),auc=as.numeric(pROC::auc(roc)),wilcoxon_p=wilcox.test(combined_score~group,data=z,exact=FALSE)$p.value)
    }
    ct<-cor.test(dd$S100A8,dd$independent_score,method="spearman",exact=FALSE)
    cov1<-data.frame(group=factor(dd$group));cov2<-data.frame(group=factor(dd$group),monocyte_core=dd$monocyte_core_score)
    cors[[ds]]<-data.frame(dataset=ds,n=nrow(dd),spearman_rho=unname(ct$estimate),p_value=ct$p.value,partial_rho_group=partial_spearman(dd$S100A8,dd$independent_score,cov1),partial_rho_group_monocyte=partial_spearman(dd$S100A8,dd$independent_score,cov2))
  }
  sample_out<-do.call(rbind,samples);comp_out<-do.call(rbind,comparisons);cor_out<-do.call(rbind,cors)
  write_csv_utf8(sample_out,file.path("08_bulk_singlecell_integration","bulk_independent_signature_scores.csv"))
  write_csv_utf8(comp_out,file.path("08_bulk_singlecell_integration","Table_BulkSignature_Comparisons.csv"))
  write_csv_utf8(cor_out,file.path("08_bulk_singlecell_integration","Table_BulkS100A8_SignatureCorrelation.csv"))
  pd<-sample_out[sample_out$dataset%in%c("GSE65682","GSE134347","GSE28750","GSE131411"),]
  bulk_cols<-c(Healthy=palette_nature[["Normal"]],Sepsis=palette_nature[["Sepsis"]],"Non-infectious ICU"="#D39B44","Postoperative sterile inflammation"="#4F8A76","Cardiogenic shock"="#777777")
  p<-ggplot(pd,aes(group,combined_score,fill=group))+
    geom_violin(width=.78,trim=TRUE,scale="width",alpha=.28,colour=NA)+
    geom_boxplot(outlier.shape=NA,width=.22,fill="white",colour="#454545",linewidth=.35)+
    geom_point(size=.48,alpha=.32,position=position_jitter(width=.09),colour="#303030")+
    facet_wrap(~dataset,scales="free_x",nrow=1)+scale_fill_manual(values=bulk_cols,na.value="#A0A0A0")+
    labs(x=NULL,y="Independent myeloid/HLA-DR-low ssGSEA score",fill=NULL)+theme_nature()+
    theme(axis.text.x=element_text(angle=28,hjust=1),legend.position="none",strip.text=element_text(hjust=0))
  save_nature_plot(p,file.path("09_figures","Fig_Bulk_IndependentSignature"),183,95,source_data=pd)
})
