# 输入：GSE65682清洗后log2基因矩阵和metadata
# 输出：1000次1:1下采样的S100A8指标，以及至少200次四模型稳定性结果
rm(list=ls()); options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","stat_functions.R"),encoding="UTF-8")
source(file.path("functions","ml_functions.R"),encoding="UTF-8")

with_script_log("04_repeated_downsampling", {
  expr<-readRDS(file.path("02_class_imbalance","GSE65682_gene_log2_mean_probes.rds"))
  meta<-read.csv(file.path("02_class_imbalance","sample_metadata_clean.csv"),check.names=FALSE)
  meta<-meta[match(colnames(expr),meta$sample_id),]; y<-as.character(meta$group); x<-t(expr); rownames(x)<-meta$sample_id
  n_single<-as.integer(Sys.getenv("DOWNSAMPLE_SINGLE_N","1000")); n_model<-as.integer(Sys.getenv("DOWNSAMPLE_MODEL_N","200"))
  normal_ids<-which(y=="Normal"); sepsis_ids<-which(y=="Sepsis")
  single<-vector("list",n_single)
  for(b in seq_len(n_single)) {
    seed<-MASTER_SEED+b; set.seed(seed); ids<-c(normal_ids,sample(sepsis_ids,length(normal_ids)))
    yy<-factor(y[ids],levels=c("Normal","Sepsis")); xx<-as.numeric(x[ids,TARGET_GENE])
    ro<-pROC::roc(yy,xx,levels=c("Normal","Sepsis"),direction="<",quiet=TRUE)
    co<-pROC::coords(ro,"best",best.method="youden",ret=c("threshold","sensitivity","specificity"),transpose=FALSE)
    rng<-range(xx,na.rm=TRUE); prob_scaled<-(xx-rng[1])/(rng[2]-rng[1]); threshold_scaled<-(as.numeric(co["threshold"])-rng[1])/(rng[2]-rng[1])
    met<-binary_metrics(yy,prob_scaled,threshold=threshold_scaled)
    single[[b]]<-data.frame(iteration=b,seed=seed,AUC=as.numeric(pROC::auc(ro)),PR_AUC=pr_auc_score(yy,xx),
                            mean_difference=mean(xx[yy=="Sepsis"])-mean(xx[yy=="Normal"]),
                            standardized_mean_difference=standardized_mean_difference(xx[yy=="Sepsis"],xx[yy=="Normal"]),
                            sensitivity=as.numeric(co["sensitivity"]),specificity=as.numeric(co["specificity"]),
                            balanced_accuracy=mean(as.numeric(co[c("sensitivity","specificity")])),
                            MCC=met["MCC"],F1=met["F1"])
  }
  single_df<-do.call(rbind,single); write_csv_utf8(single_df,file.path("02_class_imbalance","Table_Downsampling_S100A8_1000.csv"))

  algorithms<-c("LASSO","glmBoost","Ranger","XGBoost"); model_rows<-list(); feat_rows<-list()
  for(b in seq_len(n_model)) {
    seed<-MASTER_SEED+10000L+b; set.seed(seed); ids<-c(normal_ids,sample(sepsis_ids,length(normal_ids)))
    yy<-y[ids]; train_local<-as.vector(caret::createDataPartition(factor(yy),p=.7,list=FALSE)); test_local<-setdiff(seq_along(ids),train_local)
    fs<-select_training_features(x[ids[train_local],,drop=FALSE],yy[train_local],top_n=100)
    st<-scale_train_test(x[ids[train_local],fs$genes,drop=FALSE],x[ids[test_local],fs$genes,drop=FALSE])
    orig<-colnames(st$train); safe<-make.names(orig,unique=TRUE); map<-setNames(orig,safe); colnames(st$train)<-safe; colnames(st$test)<-safe
    for(ai in seq_along(algorithms)) {
      alg<-algorithms[ai]; fit<-fit_fast_stability_model(st$train,yy[train_local],st$test,alg,seed+ai*1000L)
      imp<-fit$importance; if(length(imp)) names(imp)<-unname(map[names(imp)])
      mets<-evaluate_probabilities(yy[test_local],fit$probability)
      model_rows[[length(model_rows)+1L]]<-data.frame(iteration=b,seed=seed,algorithm=alg,t(as.data.frame(mets)),check.names=FALSE)
      ranks<-seq_along(imp); feat_rows[[length(feat_rows)+1L]]<-data.frame(iteration=b,seed=seed,algorithm=alg,gene=names(imp),
        importance=as.numeric(imp),rank=ranks,selected=if(alg%in%c("Ranger","XGBoost")) ranks<=20 else TRUE)
    }
    if(b%%20==0) cat("Downsampling model iteration",b,"/",n_model,"\n")
  }
  write_csv_utf8(do.call(rbind,model_rows),file.path("02_class_imbalance","Table_Downsampling_ModelPerformance.csv"))
  write_csv_utf8(do.call(rbind,feat_rows),file.path("02_class_imbalance","Table_Downsampling_ModelFeatures.csv"))
})
