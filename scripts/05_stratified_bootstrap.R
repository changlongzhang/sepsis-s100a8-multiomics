# 输入：GSE65682 log2矩阵和metadata
# 输出：1000次类别内有放回Bootstrap的S100A8效应/AUC及四模型特征稳定性
rm(list=ls()); options(stringsAsFactors=FALSE,warn=1)
source(file.path("config","project_config.R"),encoding="UTF-8")
source(file.path("functions","io_functions.R"),encoding="UTF-8")
source(file.path("functions","stat_functions.R"),encoding="UTF-8")
source(file.path("functions","ml_functions.R"),encoding="UTF-8")

with_script_log("05_stratified_bootstrap", {
  expr<-readRDS(file.path("02_class_imbalance","GSE65682_gene_log2_mean_probes.rds")); x<-t(expr)
  meta<-read.csv(file.path("02_class_imbalance","sample_metadata_clean.csv"),check.names=FALSE); meta<-meta[match(rownames(x),meta$sample_id),]
  y<-as.character(meta$group); normal<-which(y=="Normal"); sepsis<-which(y=="Sepsis")
  B<-as.integer(Sys.getenv("BOOTSTRAP_N","1000")); model_B<-as.integer(Sys.getenv("BOOTSTRAP_MODEL_N","1000"))
  simple<-vector("list",B)
  for(b in seq_len(B)) {
    seed<-MASTER_SEED+20000L+b; set.seed(seed); ids<-c(sample(normal,length(normal),TRUE),sample(sepsis,length(sepsis),TRUE))
    yy<-factor(y[ids],levels=c("Normal","Sepsis")); xx<-x[ids,TARGET_GENE]
    ro<-pROC::roc(yy,xx,levels=c("Normal","Sepsis"),direction="<",quiet=TRUE)
    simple[[b]]<-data.frame(iteration=b,seed=seed,AUC=as.numeric(pROC::auc(ro)),
      mean_difference=mean(xx[yy=="Sepsis"])-mean(xx[yy=="Normal"]),
      standardized_mean_difference=standardized_mean_difference(xx[yy=="Sepsis"],xx[yy=="Normal"]))
  }
  write_csv_utf8(do.call(rbind,simple),file.path("02_class_imbalance","Table_Bootstrap_Results.csv"))

  algorithms<-c("LASSO","glmBoost","Ranger","XGBoost"); features<-list()
  for(b in seq_len(model_B)) {
    seed<-MASTER_SEED+30000L+b; set.seed(seed); ids<-c(sample(normal,length(normal),TRUE),sample(sepsis,length(sepsis),TRUE))
    fs<-select_training_features(x[ids,,drop=FALSE],y[ids],top_n=100)
    st<-scale_train_test(x[ids,fs$genes,drop=FALSE],x[ids,fs$genes,drop=FALSE]); orig<-colnames(st$train); safe<-make.names(orig,unique=TRUE); map<-setNames(orig,safe); colnames(st$train)<-safe; colnames(st$test)<-safe
    for(ai in seq_along(algorithms)) {
      alg<-algorithms[ai]; fit<-fit_fast_stability_model(st$train,y[ids],st$test,alg,seed+ai*1000L); imp<-fit$importance
      if(length(imp)) names(imp)<-unname(map[names(imp)]); ranks<-seq_along(imp)
      features[[length(features)+1L]]<-data.frame(iteration=b,seed=seed,algorithm=alg,gene=names(imp),importance=as.numeric(imp),rank=ranks,
        selected=if(alg%in%c("Ranger","XGBoost")) ranks<=20 else TRUE)
    }
    if(b%%50==0) cat("Bootstrap model iteration",b,"/",model_B,"\n")
  }
  write_csv_utf8(do.call(rbind,features),file.path("02_class_imbalance","Table_Bootstrap_ModelFeatures.csv"))
})

