# 输入：GSE65682、GSE26440、GSE9692的S100A8及公开结局元数据
# 输出：死亡logistic/Cox关联、可用协变量调整、随机效应Meta分析和缺失严重度变量审计
rm(list=ls()); options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","geo_functions.R"),encoding="UTF-8")

fit_logistic <- function(d,formula_text,dataset,model_name) {
  f<-stats::as.formula(formula_text); fit<-stats::glm(f,data=d,family=stats::binomial())
  cf<-summary(fit)$coefficients; term<-grep("S100A8_z",rownames(cf),value=TRUE)[1]
  est<-cf[term,"Estimate"]; se<-cf[term,"Std. Error"]
  data.frame(dataset=dataset,outcome="mortality",model=model_name,n=stats::nobs(fit),events=sum(stats::model.response(stats::model.frame(fit))==1),
    effect_type="OR per 1 SD",estimate=exp(est),ci_low=exp(est-1.96*se),ci_high=exp(est+1.96*se),p_value=cf[term,"Pr(>|z|)"],
    log_effect=est,se=se,converged=fit$converged)
}

with_script_log("12_prognosis_severity_analysis", {
  # GSE65682
  expr<-readRDS(file.path("02_class_imbalance","GSE65682_gene_log2_mean_probes.rds")); meta<-read.csv(file.path("02_class_imbalance","sample_metadata_clean.csv"),check.names=FALSE)
  meta<-meta[match(colnames(expr),meta$sample_id),]; d656<-meta[meta$group=="Sepsis",]
  d656$S100A8<-as.numeric(expr[TARGET_GENE,d656$sample_id]); d656$S100A8_z<-as.numeric(scale(d656$S100A8))
  d656$mortality<-suppressWarnings(as.numeric(d656[["mortality_event_28days:ch1"]])); d656$time28<-suppressWarnings(as.numeric(d656[["time_to_event_28days:ch1"]]))
  d656$age<-suppressWarnings(as.numeric(d656[["age:ch1"]])); d656$sex<-factor(d656[["gender:ch1"]])

  # GSE26440和GSE9692，仅保留明确septic shock样本。
  sample_expr<-read.csv(file.path("03_clinical_validation","clinical_validation_S100A8_samples.csv"),check.names=FALSE)
  outcome_sets<-list(GSE65682=d656)
  for(gse in c("GSE26440","GSE9692")) {
    f<-file.path(PROJECT_ROOT,"sepsis-s100a8-multiomics","data","08_External Validation","Input",paste0(gse,"_series_matrix.txt.gz")); es<-read_local_eset(f); pd<-Biobase::pData(es)
    txt<-combine_characteristics(pd); is_case<-grepl("septic shock",txt,ignore.case=TRUE)
    out<-ifelse(grepl("non.?survivor|non survivor",txt,ignore.case=TRUE),1,ifelse(grepl("survivor",txt,ignore.case=TRUE),0,NA))
    dd<-data.frame(sample_id=rownames(pd),mortality=out,age=NA_real_,sex=NA_character_,stringsAsFactors=FALSE)
    if("age (years):ch1"%in%names(pd))dd$age<-suppressWarnings(as.numeric(pd[["age (years):ch1"]]))
    if("gender:ch1"%in%names(pd))dd$sex<-as.character(pd[["gender:ch1"]])
    dd<-dd[is_case,]; ee<-sample_expr[sample_expr$dataset==gse,c("sample_id","expression")]
    dd<-merge(dd,ee,by="sample_id"); names(dd)[names(dd)=="expression"]<-"S100A8"; dd$S100A8_z<-as.numeric(scale(dd$S100A8)); dd$sex<-factor(dd$sex)
    outcome_sets[[gse]]<-dd
  }
  outcome_samples<-do.call(rbind,lapply(names(outcome_sets),function(nm){d<-outcome_sets[[nm]]; d$dataset<-nm; d[,c("dataset","sample_id","S100A8","S100A8_z","mortality","age","sex")]}))
  write_csv_utf8(outcome_samples,file.path("04_prognosis_severity","prognosis_samples_clean.csv"))

  res<-list(); res[[1]]<-fit_logistic(subset(d656,!is.na(mortality)),"mortality~S100A8_z","GSE65682","unadjusted")
  adj656<-subset(d656,complete.cases(mortality,S100A8_z,age,sex)); if(nrow(adj656)>100)res[[2]]<-fit_logistic(adj656,"mortality~S100A8_z+age+sex","GSE65682","age_sex_adjusted")
  for(gse in c("GSE26440","GSE9692")) {
    d<-outcome_sets[[gse]]; d<-d[!is.na(d$mortality),]
    res[[length(res)+1L]]<-fit_logistic(d,"mortality~S100A8_z",gse,"unadjusted")
    if(sum(complete.cases(d[,c("mortality","S100A8_z","age")]))>=20 && length(unique(na.omit(d$age)))>2)
      res[[length(res)+1L]]<-fit_logistic(d,"mortality~S100A8_z+age",gse,"age_adjusted")
  }
  assoc<-do.call(rbind,res)
  assoc$ph_test_p<-NA_real_

  cox_d<-subset(d656,complete.cases(mortality,time28,S100A8_z) & time28>0)
  if(nrow(cox_d)>20) {
    cfit<-survival::coxph(survival::Surv(time28,mortality)~S100A8_z,data=cox_d); cc<-summary(cfit)$coefficients[1,]; ci<-summary(cfit)$conf.int[1,]
    ph<-survival::cox.zph(cfit)$table[1,"p"]
    assoc<-rbind(assoc,data.frame(dataset="GSE65682",outcome="28-day time-to-death",model="unadjusted Cox",n=nrow(cox_d),events=sum(cox_d$mortality),
      effect_type="HR per 1 SD",estimate=ci["exp(coef)"],ci_low=ci["lower .95"],ci_high=ci["upper .95"],p_value=cc["Pr(>|z|)"],log_effect=cc["coef"],se=cc["se(coef)"],converged=TRUE,ph_test_p=ph))
  }
  write_csv_utf8(assoc,file.path("04_prognosis_severity","Table_Prognosis_Associations.csv"))

  meta_rows<-assoc[assoc$model=="unadjusted" & assoc$effect_type=="OR per 1 SD",]
  if(nrow(meta_rows)>=2) {
    mf<-metafor::rma(yi=meta_rows$log_effect,sei=meta_rows$se,method="REML")
    meta_out<-data.frame(k=mf$k,pooled_OR=exp(as.numeric(mf$b)),ci_low=exp(mf$ci.lb),ci_high=exp(mf$ci.ub),p_value=mf$pval,tau2=mf$tau2,I2=mf$I2)
    write_csv_utf8(meta_out,file.path("04_prognosis_severity","prognosis_random_effects_meta.csv"))
  }
  avail<-data.frame(variable=c("28-day mortality","time-to-death","septic shock","SOFA","APACHE II","organ failure","mechanical ventilation","vasopressor use","ICU length of stay","treatment response"),
    availability=c("GSE65682/GSE26440/GSE9692","GSE65682","Case definition only","Not available in analysed public metadata","Not available in analysed public metadata","Not consistently available","Not available","Not available","Not available","Not available"),
    status=c("Analysed","Analysed","Not an analyzable within-sepsis endpoint",rep("Data unavailable",7)))
  write_csv_utf8(avail,file.path("04_prognosis_severity","severity_variable_availability.csv"))
})
