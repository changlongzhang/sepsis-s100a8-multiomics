# 输入：公开GEO；候选队列来自预设关键词和人工审查
# 输出：筛选清单、GSE134347/GSE131411本地series对象、metadata、平台注释和下载校验清单
rm(list=ls()); options(stringsAsFactors=FALSE,warn=1,timeout=1200)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")

with_script_log("09_clinical_dataset_search", {
  stopifnot(requireNamespace("GEOquery",quietly=TRUE))
  screening<-data.frame(
    dataset=c("GSE134347","GSE131411","GSE32707","GSE63042","GSE110487","GSE272769","GSE205672","GSE279448","GSE12624","GSE5760","GSE296830","GSE54514"),
    organism="Homo sapiens",
    tissue=c("Whole blood","Whole blood","Whole blood","Whole blood","Whole blood","Whole blood","PBMC/monocytes","PBMC","Whole blood","Whole blood","Blood","Whole blood"),
    design=c("Sepsis vs non-infectious ICU vs healthy","Septic shock vs cardiogenic shock, longitudinal","Sepsis/Sepsis+ARDS vs SIRS vs ventilated ICU control, repeated day 0/day 7","SIRS vs sepsis with survivor/non-survivor outcome","Septic shock early SOFA responder vs non-responder","Sepsis 30-day mortality and septic shock, n=161","Independent sepsis PBMC/monocyte cohort with ex vivo stimulation","Sepsis vs ICU controls vs healthy; adult/pediatric","Trauma at admission; later sepsis vs no sepsis","Polytrauma prognostic/genetic study","Six patients longitudinal severe sepsis/recovery","Sepsis survivor vs non-survivor longitudinal"),
    decision=c("Include","Include","Include","Include","Include","Include","Include","Reserve","Exclude from diagnostic specificity","Exclude","Exclude","Include"),
    reason=c("Direct clinically similar ICU comparator; n=298","Direct alternative-shock comparator; use T1 and one sample per donor","Direct same-study SIRS and ventilated ICU comparators","Direct SIRS comparator plus mortality outcome","SOFA-defined early treatment-response endpoint","Large independent mortality and septic-shock validation cohort","Independent monocyte-state and stimulation validation","Strong comparator but PBMC and recent dataset; reserve if processed matrix is accessible","Baseline precedes later sepsis and contains technical replicates; not a concurrent diagnostic comparison","Not a direct sepsis-vs-mimic expression design","n=6 with repeated samples; below primary group-size criterion","Longitudinal outcome cohort"),
    source_url=paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",c("GSE134347","GSE131411","GSE32707","GSE63042","GSE110487","GSE272769","GSE205672","GSE279448","GSE12624","GSE5760","GSE296830","GSE54514")),
    stringsAsFactors=FALSE)
  write_csv_utf8(screening,file.path("03_clinical_validation","clinical_dataset_screening.csv"))

  dest<-file.path("01_data_inventory","downloaded_data"); assert_dir(dest)
  targets<-unique(c("GSE134347","GSE131411",CLINICAL_MIMIC_GSE,PROGNOSIS_EXTENSION_GSE,MONOCYTE_VALIDATION_GSE))
  manifest<-list()
  for(gse in targets) {
    cat("下载/读取",gse,"\n")
    res<-NULL; last_error<-NULL
    for(attempt in 1:3) {
      res<-tryCatch(GEOquery::getGEO(gse,destdir=dest,GSEMatrix=TRUE,getGPL=FALSE),error=function(e){last_error<<-conditionMessage(e);NULL})
      if(!is.null(res)) break
      cat("重试",attempt,"失败:",last_error,"\n")
    }
    if(is.null(res)) stop(gse,"下载失败，三次重试后仍不可用: ",last_error)
    if(!is.list(res)) res<-list(res)
    for(i in seq_along(res)) {
      eset<-res[[i]]; suffix<-if(length(res)>1) paste0("_",Biobase::annotation(eset)) else ""
      rds<-file.path(dest,paste0(gse,suffix,"_eset.rds")); if(!file.exists(rds)) saveRDS(eset,rds,compress=FALSE)
      pdfile<-file.path(dest,paste0(gse,suffix,"_metadata.csv")); write_csv_utf8(cbind(sample_id=rownames(Biobase::pData(eset)),Biobase::pData(eset)),pdfile)
      platform<-Biobase::annotation(eset)
      manifest[[length(manifest)+1L]]<-data.frame(accession=gse,platform=platform,samples=ncol(Biobase::exprs(eset)),features=nrow(Biobase::exprs(eset)),rds=rds,metadata=pdfile)
      if(nzchar(platform)) {
        gpl_csv<-file.path(dest,paste0(platform,"_annotation.csv"))
        if(!file.exists(gpl_csv)) {
          gp<-tryCatch(GEOquery::getGEO(platform,destdir=dest),error=function(e)NULL)
          if(!is.null(gp)) write_csv_utf8(GEOquery::Table(gp),gpl_csv)
        }
      }
    }
  }
  man<-do.call(rbind,manifest); man$download_date<-as.character(Sys.Date())
  man$rds_md5<-unname(tools::md5sum(man$rds)); man$metadata_md5<-unname(tools::md5sum(man$metadata))
  old_manifest<-file.path(dest,"download_manifest.csv")
  if(file.exists(old_manifest)) {
    old<-read.csv(old_manifest,check.names=FALSE,stringsAsFactors=FALSE)
    all_names<-union(names(old),names(man)); for(nm in setdiff(all_names,names(old)))old[[nm]]<-NA; for(nm in setdiff(all_names,names(man)))man[[nm]]<-NA
    man<-rbind(old[,all_names,drop=FALSE],man[,all_names,drop=FALSE]); man<-man[!duplicated(paste(man$accession,man$platform,sep="|"),fromLast=TRUE),]
  }
  write_csv_utf8(man,old_manifest)
})
