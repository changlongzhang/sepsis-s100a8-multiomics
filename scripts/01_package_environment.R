# 输入：本机R 4.5.1包库
# 输出：01_data_inventory/package_versions.csv、sessionInfo.txt
# 目的：记录实际分析环境，不静默伪装缺失依赖。
rm(list = ls())
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")

with_script_log("01_package_environment", {
  packages <- c("data.table","GEOquery","limma","edgeR","DESeq2","caret","glmnet","mboost","xgboost","ranger","pROC","PRROC","boot","rsample","recipes","yardstick","dcurves","rmda","survival","metafor","Seurat","SingleCellExperiment","muscat","scuttle","scran","variancePartition","dreamlet","glmmTMB","lme4","speckle","decoupleR","UCell","GSVA","ComplexHeatmap","patchwork","ggplot2","ggrepel","svglite","ragg","openxlsx","officer","here","renv")
  installed <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  versions <- vapply(packages, function(x) if (requireNamespace(x, quietly = TRUE)) as.character(utils::packageVersion(x)) else NA_character_, character(1))
  out <- data.frame(package = packages, installed = installed, version = versions)
  write_csv_utf8(out, file.path("01_data_inventory", "package_versions.csv"))
  capture.output(sessionInfo(), file = file.path("01_data_inventory", "sessionInfo.txt"))
  cat("已安装", sum(installed), "/", length(packages), "个优先包。缺失包:", paste(packages[!installed], collapse=", "), "\n")
})

