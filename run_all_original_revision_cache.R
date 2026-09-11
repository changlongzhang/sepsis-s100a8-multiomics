# 从空白环境运行全部返修分析。重计算轻量步骤；对已通过完整性检查的高耗时步骤使用本地缓存。
rm(list=ls());options(stringsAsFactors=FALSE,warn=1)
root<-normalizePath(".",winslash="/",mustWork=TRUE);dir.create("11_logs",recursive=TRUE,showWarnings=FALSE)
logfile<-file.path("11_logs","run_all_console.log");sink(logfile,split=TRUE);on.exit(sink(),add=TRUE)
cat("run_all start:",format(Sys.time()),"\nproject:",root,"\n")
rscript<-file.path(R.home("bin"),"Rscript.exe");if(!file.exists(rscript))rscript<-file.path(R.home("bin"),"Rscript")
count_csv<-function(p) if(file.exists(p))nrow(tryCatch(read.csv(p),error=function(e)data.frame())) else 0L
cache_ok<-list(
  `03_class_weighted_models.R`=count_csv("02_class_imbalance/weighted_cv_fold_metrics.csv")==400,
  `04_repeated_downsampling.R`=count_csv("02_class_imbalance/Table_Downsampling_S100A8_1000.csv")==1000&&count_csv("02_class_imbalance/Table_Downsampling_ModelPerformance.csv")==800,
  `05_stratified_bootstrap.R`=count_csv("02_class_imbalance/Table_Bootstrap_Results.csv")==1000,
  `06_balanced_cross_validation.R`=count_csv("02_class_imbalance/balance_strategy_fold_metrics_unweighted_down.csv")==800,
  `09_clinical_dataset_search.R`=file.exists("01_data_inventory/downloaded_data/GSE134347_eset.rds")&&file.exists("01_data_inventory/downloaded_data/GSE131411_rawcounts_CS_SS.xlsx"),
  `10_clinical_comparator_validation.R`=count_csv("03_clinical_validation/Table_ClinicalValidation_long.csv")>=100,
  `14_scrna_metadata_qc.R`=count_csv("05_single_cell_pseudobulk/donor_metadata_clean.csv")==7,
  `15_scrna_pseudobulk.R`=file.exists("05_single_cell_pseudobulk/Table_Pseudobulk_FocusGenes.csv"),
  `16_scrna_mixed_effects.R`=file.exists("05_single_cell_pseudobulk/Table_MixedEffects_FocusGenes.csv"),
  `17_monocyte_reclustering.R`=file.exists("06_monocyte_state/GSE167363_monocyte_reclustered.rds"),
  `18_monocyte_state_annotation.R`=file.exists("06_monocyte_state/GSE167363_monocyte_state_annotated.rds"),
  `19_donor_level_abundance.R`=file.exists("06_monocyte_state/Table_DonorLevel_StateAbundance.csv"),
  `20_donor_level_tf_activity.R`=file.exists("05_single_cell_pseudobulk/Table_DonorLevel_TFActivity.csv"),
  `21_non_circular_signature.R`=file.exists("07_non_circular_signature/independent_signature_genes.csv"),
  `22_bulk_signature_projection.R`=file.exists("08_bulk_singlecell_integration/Table_BulkSignature_Comparisons.csv")
)
scripts<-sprintf("%02d_%s.R",0:26,c("project_audit","package_environment","discovery_data_qc","class_weighted_models","repeated_downsampling","stratified_bootstrap","balanced_cross_validation","class_imbalance_summary","external_metadata_audit","clinical_dataset_search","clinical_comparator_validation","diagnostic_metrics_dca","prognosis_severity_analysis","clinical_meta_analysis","scrna_metadata_qc","scrna_pseudobulk","scrna_mixed_effects","monocyte_reclustering","monocyte_state_annotation","donor_level_abundance","donor_level_tf_activity","non_circular_signature","bulk_signature_projection","integrated_evidence_summary","generate_revision_tables","generate_manuscript_text","generate_reviewer_response"))
Sys.setenv(CV_REPEATS="20")
for(s in scripts){p<-file.path("scripts",s);if(!file.exists(p))stop("Missing script: ",p)
  if(!is.null(cache_ok[[s]])&&isTRUE(cache_ok[[s]])){cat("CACHED VERIFIED:",s,"\n");next}
  cat("RUN:",s,"at",format(Sys.time()),"\n");out<-system2(rscript,p,stdout=TRUE,stderr=TRUE);cat(paste(out,collapse="\n"),"\n");status<-attr(out,"status");if(!is.null(status)&&status!=0)stop("Script failed: ",s)
}
cat("run_all completed:",format(Sys.time()),"\n")
