# Portable cached numerical runner for the final revision.
rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(Sys.getenv("SEPSIS_PROJECT_ROOT", unset = ".."), winslash = "/", mustWork = TRUE)
Sys.setenv(SEPSIS_PROJECT_ROOT = root, CV_REPEATS = "20")
source(file.path("config", "project_config.R"), encoding = "UTF-8")
dir.create("11_logs", recursive = TRUE, showWarnings = FALSE)
logfile <- file.path("11_logs", "run_analysis_console.log")
sink(logfile, split = TRUE)
on.exit(sink(), add = TRUE)

rscript <- file.path(R.home("bin"), "Rscript.exe")
if (!file.exists(rscript)) rscript <- file.path(R.home("bin"), "Rscript")
count_csv <- function(path) if (file.exists(path)) nrow(tryCatch(read.csv(path), error = function(e) data.frame())) else 0L
cache_ok <- list(
  `03_class_weighted_models.R` = count_csv("02_class_imbalance/weighted_cv_fold_metrics.csv") == 400,
  `04_repeated_downsampling.R` = count_csv("02_class_imbalance/Table_Downsampling_S100A8_1000.csv") == 1000,
  `05_stratified_bootstrap.R` = count_csv("02_class_imbalance/Table_Bootstrap_Results.csv") == 1000,
  `06_balanced_cross_validation.R` = count_csv("02_class_imbalance/balance_strategy_fold_metrics_unweighted_down.csv") == 800,
  `09_clinical_dataset_search.R` = file.exists("03_clinical_validation/clinical_dataset_screening.csv"),
  `10_clinical_comparator_validation.R` = count_csv("03_clinical_validation/Table_ClinicalValidation_long.csv") >= 100,
  `14_scrna_metadata_qc.R` = count_csv("05_single_cell_pseudobulk/donor_metadata_clean.csv") == 7,
  `15_scrna_pseudobulk.R` = file.exists("05_single_cell_pseudobulk/Table_Pseudobulk_FocusGenes.csv"),
  `17_monocyte_reclustering.R` = file.exists("06_monocyte_state/GSE167363_monocyte_reclustered.rds"),
  `18_monocyte_state_annotation.R` = file.exists("06_monocyte_state/GSE167363_monocyte_state_annotated.rds"),
  `19_donor_level_abundance.R` = file.exists("06_monocyte_state/Table_DonorLevel_StateAbundance.csv"),
  `20_donor_level_tf_activity.R` = file.exists("05_single_cell_pseudobulk/Table_DonorLevel_TFActivity.csv"),
  `21_non_circular_signature.R` = file.exists("07_non_circular_signature/independent_signature_genes.csv"),
  `22_bulk_signature_projection.R` = file.exists("08_bulk_singlecell_integration/Table_BulkSignature_Comparisons.csv")
)

scripts <- list.files("scripts", pattern = "^(0[1-9]|1[0-9]|2[0-4])_.*\\.R$", full.names = TRUE)
scripts <- sort(scripts)
for (path in scripts) {
  name <- basename(path)
  if (!is.null(cache_ok[[name]]) && isTRUE(cache_ok[[name]])) {
    cat("CACHED VERIFIED:", name, "\n")
    next
  }
  cat("RUN:", name, "at", format(Sys.time()), "\n")
  output <- system2(rscript, path, stdout = TRUE, stderr = TRUE)
  cat(paste(output, collapse = "\n"), "\n")
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) stop("Script failed: ", name)
}
cat("Completed:", format(Sys.time()), "\n")
