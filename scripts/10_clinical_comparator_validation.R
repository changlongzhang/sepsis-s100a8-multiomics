# 输入：4个本地外部队列、GSE134347、GSE131411原始counts
# 输出：逐样本S100A8表达和各临床比较的效应量/ROC/PR/分类指标及Bootstrap CI
rm(list=ls()); options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","geo_functions.R"),encoding="UTF-8")
source(file.path("functions","stat_functions.R"),encoding="UTF-8")
source(file.path("functions","ml_functions.R"),encoding="UTF-8")

metric_vector <- function(group,expr) {
  group<-factor(group,levels=c("Comparator","Sepsis")); ro<-pROC::roc(group,expr,levels=c("Comparator","Sepsis"),direction="<",quiet=TRUE)
  co<-pROC::coords(ro,"best",best.method="youden",ret=c("threshold","sensitivity","specificity"),transpose=FALSE)
  thr<-as.numeric(unlist(co["threshold"])[1]); pred<-expr>=thr; truth<-group=="Sepsis"
  tp<-sum(pred&truth); tn<-sum(!pred&!truth); fp<-sum(pred&!truth); fn<-sum(!pred&truth)
  sens<-tp/(tp+fn); spec<-tn/(tn+fp); ppv<-tp/(tp+fp); npv<-tn/(tn+fn)
  lrpos<-sens/(1-spec); lrneg<-(1-sens)/spec; mccd<-sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); mcc<-if(mccd>0)(tp*tn-fp*fn)/mccd else NA
  c(AUC=as.numeric(pROC::auc(ro)),PR_AUC=pr_auc_score(group,expr),mean_difference=mean(expr[truth])-mean(expr[!truth]),
    SMD=standardized_mean_difference(expr[truth],expr[!truth]),threshold=thr,sensitivity=sens,specificity=spec,
    LR_positive=lrpos,LR_negative=lrneg,PPV=ppv,NPV=npv,balanced_accuracy=mean(c(sens,spec)),MCC=mcc)
}

comparison_summary <- function(df,B=1000,seed=MASTER_SEED) {
  point<-metric_vector(df$binary_group,df$expression); ids0<-which(df$binary_group=="Comparator"); ids1<-which(df$binary_group=="Sepsis")
  set.seed(seed); boot<-replicate(B,{ii<-c(sample(ids0,length(ids0),TRUE),sample(ids1,length(ids1),TRUE)); metric_vector(df$binary_group[ii],df$expression[ii])})
  lo<-apply(boot,1,quantile,.025,na.rm=TRUE); hi<-apply(boot,1,quantile,.975,na.rm=TRUE)
  data.frame(metric=names(point),estimate=as.numeric(point),ci_low=lo[names(point)],ci_high=hi[names(point)],row.names=NULL)
}

extract_existing <- function(gse) {
  f<-file.path(PROJECT_ROOT,"sepsis-s100a8-multiomics","data","08_External Validation","Input",paste0(gse,"_series_matrix.txt.gz")); es<-read_local_eset(f)
  mat<-Biobase::exprs(es); pd<-Biobase::pData(es); txt<-tolower(paste(as.character(pd$title),as.character(pd$source_name_ch1),combine_characteristics(pd)))
  if(Biobase::annotation(es)=="GPL570") probes<-intersect(c("202917_s_at","214370_at"),rownames(mat)) else probes<-intersect("ILMN_1729801",rownames(mat))
  if(!length(probes)) stop(gse,"未找到经注释包验证的S100A8探针")
  # Keep the external-cohort aggregation rule aligned with the manuscript and
  # the original discovery workflow: one value per gene/sample is the maximum
  # across annotation-confirmed probes.  The previous revision script used a
  # probe mean here, which created conflicting AUCs between Fig. 3A, the text,
  # and the source-data tables.
  val<-if(length(probes)>1) apply(mat[probes,,drop=FALSE],2,max,na.rm=TRUE) else as.numeric(mat[probes,])
  group<-rep("Unresolved",ncol(mat))
  group[grepl("sepsis|septic shock",txt)]<-"Sepsis"
  group[grepl("healthy|normal control|normal children",txt)]<-"Healthy"
  group[grepl("post[_ -]?surgical|postoperative",txt)]<-"Postoperative sterile inflammation"
  group[grepl("uninfected type 2 diabetes",txt)]<-"Uninfected chronic disease control"
  data.frame(dataset=gse,sample_id=colnames(mat),group=group,expression=val,platform=Biobase::annotation(es),probe=paste(probes,collapse=";"))
}

with_script_log("10_clinical_comparator_validation", {
  needed<-c("pROC","PRROC","edgeR","openxlsx","hta20transcriptcluster.db")
  miss<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]; if(length(miss))stop("缺少包:",paste(miss,collapse=", "))
  samples<-do.call(rbind,lapply(EXTERNAL_GSE,extract_existing))

  g134<-readRDS(file.path("01_data_inventory","downloaded_data","GSE134347_eset.rds")); pd134<-Biobase::pData(g134)
  probe134<-"TC01003261.hg.1"; if(!probe134%in%rownames(Biobase::exprs(g134)))stop("GSE134347 S100A8 transcript cluster missing")
  g134_group<-c(healthy="Healthy",noninfectious="Non-infectious ICU",sepsis="Sepsis")[pd134[["disease state:ch1"]]]
  samples<-rbind(samples,data.frame(dataset="GSE134347",sample_id=colnames(g134),group=unname(g134_group),
    expression=as.numeric(Biobase::exprs(g134)[probe134,]),platform="GPL17586",probe=probe134))

  count_file<-file.path("01_data_inventory","downloaded_data","GSE131411_rawcounts_CS_SS.xlsx")
  cnt<-openxlsx::read.xlsx(count_file,sheet=1,colNames=TRUE); genes<-as.character(cnt[[1]]); cnt<-as.matrix(cnt[,-1]); rownames(cnt)<-genes; storage.mode(cnt)<-"numeric"
  keep_cols<-grepl("T1$",colnames(cnt)); cnt<-cnt[,keep_cols,drop=FALSE]
  dge<-edgeR::DGEList(counts=cnt); keep<-edgeR::filterByExpr(dge); dge<-edgeR::calcNormFactors(dge[keep,,keep.lib.sizes=FALSE]); logcpm<-edgeR::cpm(dge,log=TRUE,prior.count=2)
  ens<-"ENSG00000143546"; if(!ens%in%rownames(logcpm))stop("GSE131411 counts中未找到S100A8 Ensembl ID")
  md<-read.csv(file.path("01_data_inventory","downloaded_data","GSE131411_GPL16791_metadata.csv"),check.names=FALSE)
  key<-sub("_.*$","",md$title); grp<-ifelse(grepl("SepticShock",md$title),"Sepsis",ifelse(grepl("CardiogenicShock",md$title),"Cardiogenic shock","Other")); names(grp)<-key
  g<-unname(grp[colnames(logcpm)])
  samples<-rbind(samples,data.frame(dataset="GSE131411",sample_id=colnames(logcpm),group=g,expression=as.numeric(logcpm[ens,]),platform="RNA-seq",probe=ens))
  write_csv_utf8(samples,file.path("03_clinical_validation","clinical_validation_S100A8_samples.csv"))

  specs<-list(
    GSE26440_Healthy=c("GSE26440","Healthy"),GSE28750_Healthy=c("GSE28750","Healthy"),
    GSE28750_Postoperative=c("GSE28750","Postoperative sterile inflammation"),GSE69528_Healthy=c("GSE69528","Healthy"),
    GSE69528_ChronicDisease=c("GSE69528","Uninfected chronic disease control"),GSE9692_Healthy=c("GSE9692","Healthy"),
    GSE134347_Healthy=c("GSE134347","Healthy"),GSE134347_NoninfectiousICU=c("GSE134347","Non-infectious ICU"),
    GSE131411_CardiogenicShock=c("GSE131411","Cardiogenic shock"))
  summaries<-list()
  for(nm in names(specs)) {
    sp<-specs[[nm]]; d<-samples[samples$dataset==sp[1] & samples$group%in%c("Sepsis",sp[2]),]
    d$binary_group<-ifelse(d$group=="Sepsis","Sepsis","Comparator")
    if(min(table(d$binary_group))<5){cat("跳过样本量不足:",nm,"\n");next}
    ss<-comparison_summary(d,B=1000,seed=MASTER_SEED+match(nm,names(specs)))
    ss$comparison_id<-nm; ss$dataset<-sp[1]; ss$comparator<-sp[2]; ss$n_sepsis<-sum(d$binary_group=="Sepsis"); ss$n_comparator<-sum(d$binary_group=="Comparator")
    summaries[[nm]]<-ss; cat(nm,"n=",nrow(d),"AUC=",ss$estimate[ss$metric=="AUC"],"\n")
  }
  out<-do.call(rbind,summaries); write_csv_utf8(out,file.path("03_clinical_validation","Table_ClinicalValidation_long.csv"))
})
