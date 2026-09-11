# 输入：含临床相似对照的逐样本 S100A8 表达。
# 输出：分层5折交叉验证概率的探索性DCA、校准与Brier；不把病例抽样比例当真实患病率。
rm(list=ls());options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","theme_nature.R"),encoding="UTF-8")

cv_prob <- function(d,k=5,seed=MASTER_SEED) {
  set.seed(seed);fold<-integer(nrow(d))
  for(g in unique(d$outcome)) {ix<-which(d$outcome==g);fold[ix]<-sample(rep(seq_len(k),length.out=length(ix)))}
  p<-rep(NA_real_,nrow(d))
  for(i in seq_len(k)) {tr<-fold!=i;te<-fold==i;fit<-glm(outcome~S100A8_z,data=d[tr,],family=binomial());p[te]<-predict(fit,newdata=d[te,],type="response")}
  pmin(pmax(p,1e-6),1-1e-6)
}
net_benefit <- function(y,p,thresholds) sapply(thresholds,function(pt){pred<-p>=pt;mean(pred&y==1)-mean(pred&y==0)*pt/(1-pt)})

with_script_log("11_diagnostic_metrics_dca", {
  suppressPackageStartupMessages({library(ggplot2);library(pROC)})
  f<-file.path("03_clinical_validation","clinical_validation_S100A8_samples.csv");assert_file(f);x<-read.csv(f,check.names=FALSE)
  specs<-list(GSE28750_Postoperative=c("GSE28750","Postoperative sterile inflammation"),GSE134347_NoninfectiousICU=c("GSE134347","Non-infectious ICU"),GSE131411_CardiogenicShock=c("GSE131411","Cardiogenic shock"))
  th<-seq(.05,.80,by=.01);curves<-list();perf<-list();preds<-list()
  for(nm in names(specs)) {
    sp<-specs[[nm]];d<-x[x$dataset==sp[1]&x$group%in%c("Sepsis",sp[2]),];d<-d[nzchar(d$group),];d$outcome<-as.integer(d$group=="Sepsis");d$S100A8_z<-as.numeric(scale(d$expression));d$probability<-cv_prob(d,5,MASTER_SEED+match(nm,names(specs)));preds[[nm]]<-cbind(comparison_id=nm,d)
    prev<-mean(d$outcome);nb_model<-net_benefit(d$outcome,d$probability,th);nb_all<-prev-(1-prev)*th/(1-th)
    curves[[nm]]<-rbind(data.frame(comparison_id=nm,threshold=th,strategy="S100A8 (cross-validated)",net_benefit=nb_model),data.frame(comparison_id=nm,threshold=th,strategy="Treat all",net_benefit=nb_all),data.frame(comparison_id=nm,threshold=th,strategy="Treat none",net_benefit=0))
    ro<-pROC::roc(d$outcome,d$probability,quiet=TRUE,direction="<");cal<-glm(outcome~qlogis(probability),data=d,family=binomial())
    perf[[nm]]<-data.frame(comparison_id=nm,dataset=sp[1],comparator=sp[2],n=nrow(d),sample_case_fraction=prev,cv_auc=as.numeric(auc(ro)),brier=mean((d$outcome-d$probability)^2),calibration_intercept=coef(cal)[1],calibration_slope=coef(cal)[2],traditional_biomarkers="Not available in public metadata")
  }
  curves<-do.call(rbind,curves);perf<-do.call(rbind,perf);preds<-do.call(rbind,preds)
  write_csv_utf8(curves,file.path("03_clinical_validation","decision_curve_source_data.csv"));write_csv_utf8(perf,file.path("03_clinical_validation","Table_DCA_Calibration.csv"));write_csv_utf8(preds,file.path("03_clinical_validation","clinical_comparator_cv_predictions.csv"))
  dca_labels<-c(GSE28750_Postoperative="Postoperative sterile inflammation",GSE134347_NoninfectiousICU="Non-infectious ICU",GSE131411_CardiogenicShock="Cardiogenic shock")
  p<-ggplot(curves,aes(threshold,net_benefit,colour=strategy,linetype=strategy))+
    geom_hline(yintercept=0,colour="#BEBEBE",linewidth=.3)+geom_line(linewidth=.72)+
    facet_wrap(~comparison_id,scales="free_y",labeller=as_labeller(dca_labels),nrow=1)+
    scale_colour_manual(values=c("S100A8 (cross-validated)"="#C65A46","Treat all"="#6E7073","Treat none"="#B8B8B8"))+
    scale_linetype_manual(values=c("S100A8 (cross-validated)"=1,"Treat all"=2,"Treat none"=3))+
    labs(x="Threshold probability",y="Net benefit",colour=NULL,linetype=NULL)+theme_curve()+
    theme(legend.position="top",legend.justification="left",strip.text=element_text(hjust=0),aspect.ratio=.72)
  save_nature_plot(p,file.path("09_figures","Fig_DecisionCurve"),183,85,source_data=curves)
})
