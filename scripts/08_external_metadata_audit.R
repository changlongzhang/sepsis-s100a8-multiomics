# 输入：四个本地外部队列series matrix
# 输出：逐样本原始元数据、保守重分组建议和人工复核标记
# 原则：不把infection/pneumonia/bacterial等词自动等同于Sepsis。
rm(list=ls()); options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","geo_functions.R"),encoding="UTF-8")

classify_conservative <- function(text) {
  z<-tolower(text); group<-rep("Unresolved",length(z)); comp<-rep("Ambiguous",length(z)); reason<-rep("",length(z))
  # 从最具体的临床相似对照开始；只有明确措辞才纳入。
  is_healthy<-grepl("healthy|normal volunteer|healthy control|normal control|normal children",z)
  is_sirs<-grepl("\\bsirs\\b|systemic inflammatory response syndrome",z) & !grepl("sepsis|septic",z)
  is_trauma<-grepl("trauma|injury|burn",z) & !grepl("sepsis|septic",z)
  is_nonsep_icu<-grepl("non[- ]?septic|icu control|critical illness control",z)
  is_viral<-grepl("influenza|covid|sars-cov|viral infection|virus",z) & !grepl("sepsis|septic",z)
  is_pneu_nonsep<-grepl("pneumonia",z) & grepl("without sepsis|non[- ]?septic|no sepsis",z)
  is_inf_nonsep<-grepl("infection",z) & grepl("without (organ dysfunction|sepsis)|uncomplicated|non[- ]?septic",z)
  is_postsurg<-grepl("post[_ -]?surgical|postoperative",z) & !grepl("sepsis|septic",z)
  is_uninfected_disease<-grepl("uninfected type 2 diabetes|uninfected disease control",z)
  is_sepsis<-grepl("\\bsepsis\\b|septic shock|septic patient|severe sepsis",z)
  group[is_healthy]<-"Healthy"; comp[is_healthy]<-"Healthy control"
  group[is_sirs]<-"SIRS"; comp[is_sirs]<-"Sterile/undifferentiated inflammation"
  group[is_trauma]<-"Trauma"; comp[is_trauma]<-"Sterile inflammation"
  group[is_nonsep_icu]<-"Non-septic ICU"; comp[is_nonsep_icu]<-"Critical illness control"
  group[is_viral]<-"Viral infection"; comp[is_viral]<-"Infectious mimic"
  group[is_pneu_nonsep]<-"Pneumonia without sepsis"; comp[is_pneu_nonsep]<-"Infectious mimic"
  group[is_inf_nonsep]<-"Infection without organ dysfunction"; comp[is_inf_nonsep]<-"Infectious mimic"
  group[is_postsurg]<-"Postoperative sterile inflammation"; comp[is_postsurg]<-"Sterile inflammation"
  group[is_uninfected_disease]<-"Uninfected chronic disease control"; comp[is_uninfected_disease]<-"Non-infectious disease control"
  group[is_sepsis]<-ifelse(grepl("shock",z[is_sepsis]),"Septic shock","Sepsis"); comp[is_sepsis]<-"Sepsis case"
  reason[group=="Unresolved"]<-"Metadata lacks an explicit conservative group definition"
  list(group=group,comparator=comp,reason=reason)
}

with_script_log("08_external_metadata_audit", {
  rows<-list()
  for(gse in EXTERNAL_GSE) {
    f<-file.path(PROJECT_ROOT,"sepsis-s100a8-multiomics","data","08_External Validation","Input",paste0(gse,"_series_matrix.txt.gz"))
    assert_file(f,paste0(gse," series matrix")); eset<-read_local_eset(f); pd<-Biobase::pData(eset)
    raw_chars<-combine_characteristics(pd)
    cols<-grep("title|source_name|characteristics|description",names(pd),ignore.case=TRUE,value=TRUE)
    raw_text<-apply(pd[,cols,drop=FALSE],1,function(z) paste(na.omit(as.character(z)),collapse="; "))
    cl<-classify_conservative(raw_text)
    rows[[gse]]<-data.frame(dataset=gse,sample_id=rownames(pd),raw_title=as.character(pd$title),
      raw_source=as.character(pd$source_name_ch1),raw_characteristics=raw_chars,raw_metadata_text=raw_text,
      original_group="Re-audit required",revised_group=cl$group,comparator_type=cl$comparator,
      inclusion=cl$group!="Unresolved",exclusion_reason=cl$reason,
      reviewer_relevance=ifelse(cl$group%in%c("SIRS","Trauma","Non-septic ICU","Infection without organ dysfunction","Pneumonia without sepsis","Viral infection"),"High","Contextual"),
      manual_review_required=cl$group=="Unresolved",stringsAsFactors=FALSE)
    cat(gse,":",nrow(pd),"samples; groups:",paste(names(table(cl$group)),table(cl$group),collapse="; "),"\n")
  }
  out<-do.call(rbind,rows); write_csv_utf8(out,file.path("03_clinical_validation","sample_group_mapping_all_cohorts.csv"))
  # 便于人工核查的去重描述表
  patterns<-unique(out[,c("dataset","raw_title","raw_source","raw_characteristics","revised_group","comparator_type")])
  write_csv_utf8(patterns,file.path("03_clinical_validation","external_metadata_unique_patterns.csv"))
})
