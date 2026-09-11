# 输入：临床验证的标准化效应量及脚本12预后Meta结果。
# 输出：健康对照与临床相似对照分层的随机效应汇总；不同临床对照的合并仅为探索性异质性描述。
rm(list=ls());options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","theme_nature.R"),encoding="UTF-8")

with_script_log("13_clinical_meta_analysis", {
  suppressPackageStartupMessages({library(metafor);library(ggplot2)})
  x<-read.csv(file.path("03_clinical_validation","Table_ClinicalValidation_long.csv"),check.names=FALSE)
  x<-x[x$metric=="SMD"&is.finite(x$estimate)&is.finite(x$ci_low)&is.finite(x$ci_high),]
  x$se<-(x$ci_high-x$ci_low)/(2*1.96);x$family<-ifelse(x$comparator=="Healthy","Healthy case-control","Clinically similar comparator")
  fits<-list();summary<-list()
  for(fam in unique(x$family)) {d<-x[x$family==fam&x$se>0,];if(nrow(d)<2)next;fit<-metafor::rma(yi=d$estimate,sei=d$se,method="REML");fits[[fam]]<-fit;summary[[fam]]<-data.frame(comparator_family=fam,k=fit$k,pooled_SMD=as.numeric(fit$b),ci_low=fit$ci.lb,ci_high=fit$ci.ub,p_value=fit$pval,I2=fit$I2,interpretation=if(fam=="Clinically similar comparator")"Exploratory only because comparator phenotypes differ" else "Case-control evidence")}
  out<-do.call(rbind,summary);write_csv_utf8(out,file.path("03_clinical_validation","clinical_effect_random_effects_meta.csv"));write_csv_utf8(x,file.path("03_clinical_validation","clinical_meta_input.csv"))
  pd<-x;pd$label<-paste(pd$dataset,pd$comparator,sep=" vs ")
  pd$label<-factor(pd$label,levels=rev(unique(pd$label[order(pd$family,pd$estimate)])))
  p<-ggplot(pd,aes(estimate,label,colour=family))+
    geom_vline(xintercept=0,linetype=2,colour="#8C8C8C",linewidth=.35)+
    geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=0,linewidth=.55)+geom_point(size=2,shape=18)+
    facet_grid(family~.,scales="free_y",space="free_y")+
    scale_colour_manual(values=c("Healthy case-control"="#315F7D","Clinically similar comparator"="#C65A46"))+
    labs(x="Standardized mean difference (95% CI)",y=NULL,colour=NULL)+theme_forest()+
    theme(legend.position="none",strip.text.y=element_text(angle=0,hjust=0),strip.placement="outside")
  save_nature_plot(p,file.path("09_figures","Fig_ClinicalComparator_EffectForest"),183,115,source_data=pd)
})
