# 输入：02_discovery_data_qc.R生成的GSE65682 log2基因矩阵和metadata
# 输出：类别加权、5折分层交叉验证重复20次的折级性能、预测和特征稳定性
# 所有监督特征筛选、标准化和调参均在训练折内完成。
rm(list=ls())
options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","stat_functions.R"),encoding="UTF-8")
source(file.path("functions","ml_functions.R"),encoding="UTF-8")

with_script_log("03_class_weighted_models", {
  needed <- c("caret","limma","glmnet","mboost","ranger","xgboost","pROC","PRROC")
  miss <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
  if(length(miss)) stop("缺少模型包: ",paste(miss,collapse=", "))
  expr_file <- file.path("02_class_imbalance","GSE65682_gene_log2_mean_probes.rds")
  meta_file <- file.path("02_class_imbalance","sample_metadata_clean.csv")
  assert_file(expr_file); assert_file(meta_file)
  expr <- readRDS(expr_file)
  meta <- read.csv(meta_file,check.names=FALSE)
  meta <- meta[match(colnames(expr),meta$sample_id),]
  x <- t(expr); rownames(x) <- meta$sample_id
  y <- as.character(meta$group)
  repeats <- as.integer(Sys.getenv("CV_REPEATS","20")); if(!is.finite(repeats)||repeats<1) repeats<-20L
  folds <- stratified_outer_folds(y,k=5,repeats=repeats,seed=MASTER_SEED)
  all_metrics <- list(); all_pred <- list(); all_feat <- list(); all_screen <- list()
  for(i in seq_along(folds)) {
    cat(sprintf("Weighted outer fold %d/%d\n",i,length(folds)))
    rr <- run_outer_fold(x,y,folds[[i]],strategy="weighted",seed=MASTER_SEED+i,top_n=100)
    for(j in seq_along(rr)) {
      key <- length(all_metrics)+1L
      all_metrics[[key]] <- cbind(fold=names(folds)[i],strategy="weighted",rr[[j]]$metrics)
      all_pred[[key]] <- cbind(fold=names(folds)[i],strategy="weighted",rr[[j]]$predictions)
      if(nrow(rr[[j]]$features)) all_feat[[key]] <- cbind(fold=names(folds)[i],strategy="weighted",rr[[j]]$features)
      all_screen[[key]] <- data.frame(fold=names(folds)[i],strategy="weighted",algorithm=rr[[j]]$metrics$algorithm,
                                      gene=rr[[j]]$screened_genes)
    }
  }
  metrics <- do.call(rbind,all_metrics); pred <- do.call(rbind,all_pred)
  feat <- if(length(all_feat)) do.call(rbind,all_feat) else data.frame()
  screen <- do.call(rbind,all_screen)
  write_csv_utf8(metrics,file.path("02_class_imbalance","weighted_cv_fold_metrics.csv"))
  write_csv_utf8(pred,file.path("02_class_imbalance","weighted_cv_predictions.csv"))
  write_csv_utf8(feat,file.path("02_class_imbalance","weighted_cv_feature_importance.csv"))
  write_csv_utf8(screen,file.path("02_class_imbalance","weighted_cv_training_screen.csv"))
  saveRDS(folds,file.path("02_class_imbalance","outer_folds_5x20.rds"))
  cat("完成",length(folds),"个外层折、",nrow(metrics),"个算法-折结果。\n")
})

