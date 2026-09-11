# Data conversion only: no graphics are produced in this script.
load(file.path("sepsis-s100a8-multiomics", "data", "08_External Validation", "output", "08_External Validation.RData"))
roc_rows <- do.call(rbind, lapply(names(final_results), function(nm) {
  z <- final_results[[nm]]
  data.frame(
    dataset = nm,
    specificity = as.numeric(z$roc$specificities),
    sensitivity = as.numeric(z$roc$sensitivities),
    AUC = as.numeric(z$auc),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  roc_rows,
  file.path("review_revision", "09_figures", "Fig3A_ExternalROC_source_data.csv"),
  row.names = FALSE
)
