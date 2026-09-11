# 输入：固定的5折×20次外层划分
# 输出：原始不平衡、类别加权、训练折内下采样三策略比较；SMOTE不作为主要证据
rm(list=ls()); options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","stat_functions.R"),encoding="UTF-8")
source(file.path("functions","ml_functions.R"),encoding="UTF-8")

with_script_log("06_balanced_cross_validation", {
  expr<-readRDS(file.path("02_class_imbalance","GSE65682_gene_log2_mean_probes.rds")); x<-t(expr)
  meta<-read.csv(file.path("02_class_imbalance","sample_metadata_clean.csv"),check.names=FALSE)
  meta<-meta[match(rownames(x),meta$sample_id),]; rownames(x)<-meta$sample_id; y<-as.character(meta$group)
  repeats<-as.integer(Sys.getenv("CV_REPEATS","20")); folds<-stratified_outer_folds(y,5,repeats,MASTER_SEED)
  strategies<-c("unweighted","downsample"); metrics<-list(); preds<-list(); feats<-list()
  for(s in strategies) for(i in seq_along(folds)) {
    cat(s,"outer fold",i,"/",length(folds),"\n")
    rr<-run_outer_fold(x,y,folds[[i]],strategy=s,seed=MASTER_SEED+i+match(s,strategies)*100000L,top_n=100)
    for(j in seq_along(rr)) {
      metrics[[length(metrics)+1L]]<-cbind(fold=names(folds)[i],strategy=s,rr[[j]]$metrics)
      preds[[length(preds)+1L]]<-cbind(fold=names(folds)[i],strategy=s,rr[[j]]$predictions)
      if(nrow(rr[[j]]$features)) feats[[length(feats)+1L]]<-cbind(fold=names(folds)[i],strategy=s,rr[[j]]$features)
    }
  }
  write_csv_utf8(do.call(rbind,metrics),file.path("02_class_imbalance","balance_strategy_fold_metrics_unweighted_down.csv"))
  write_csv_utf8(do.call(rbind,preds),file.path("02_class_imbalance","balance_strategy_predictions_unweighted_down.csv"))
  write_csv_utf8(do.call(rbind,feats),file.path("02_class_imbalance","balance_strategy_features_unweighted_down.csv"))
})
