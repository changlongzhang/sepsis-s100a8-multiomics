# 类别不平衡与嵌套交叉验证公共函数
select_training_features <- function(x_train, y_train, top_n = 100L) {
  # x_train: samples x genes；只使用训练折标签进行limma排序。
  y_train <- factor(y_train, levels = c("Normal", "Sepsis"))
  design <- stats::model.matrix(~ y_train)
  fit <- limma::eBayes(limma::lmFit(t(x_train), design), robust = TRUE)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
  genes <- rownames(tt)[seq_len(min(top_n, nrow(tt)))]
  list(genes = genes, table = tt)
}

scale_train_test <- function(x_train, x_test) {
  mu <- colMeans(x_train, na.rm = TRUE)
  sig <- apply(x_train, 2, stats::sd, na.rm = TRUE)
  keep <- is.finite(sig) & sig > 1e-8 & is.finite(mu)
  x_train <- x_train[, keep, drop = FALSE]
  x_test <- x_test[, keep, drop = FALSE]
  list(train = sweep(sweep(x_train, 2, mu[keep], "-"), 2, sig[keep], "/"),
       test = sweep(sweep(x_test, 2, mu[keep], "-"), 2, sig[keep], "/"),
       center = mu[keep], scale = sig[keep])
}

inverse_frequency_weights <- function(y) {
  tab <- table(y)
  w <- length(y)/(length(tab)*as.numeric(tab[y]))
  as.numeric(w/mean(w))
}

model_grid <- function(algorithm, p) {
  if (algorithm == "LASSO") return(expand.grid(alpha = 1, lambda = 10^seq(-3, 0, length.out = 5)))
  if (algorithm == "glmBoost") return(expand.grid(mstop = c(50,100,200), prune = "no"))
  if (algorithm == "Ranger") return(expand.grid(mtry = unique(pmax(1, pmin(p, round(c(sqrt(p), p/5))))), splitrule = "gini", min.node.size = c(1,5)))
  if (algorithm == "XGBoost") return(expand.grid(nrounds=c(50,100), max_depth=c(2,3), eta=c(0.05), gamma=0,
                                                   colsample_bytree=0.8, min_child_weight=1, subsample=0.8))
  stop("未知算法: ", algorithm, call. = FALSE)
}

caret_method <- function(algorithm) switch(algorithm, LASSO="glmnet", glmBoost="glmboost", Ranger="ranger", XGBoost="xgbTree")

fit_inner_model <- function(x_train, y_train, x_test, algorithm, weighted = FALSE, seed = 1L) {
  set.seed(seed)
  y_train <- factor(y_train, levels = c("Sepsis", "Normal"))
  ctrl <- caret::trainControl(method="cv", number=3, classProbs=TRUE,
                              summaryFunction=caret::twoClassSummary, savePredictions="none",
                              allowParallel=FALSE)
  args <- list(x=x_train, y=y_train, method=caret_method(algorithm), metric="ROC",
               trControl=ctrl, tuneGrid=model_grid(algorithm,ncol(x_train)))
  if (weighted) args$weights <- inverse_frequency_weights(y_train)
  if (algorithm == "Ranger") args$importance <- "permutation"
  if (algorithm == "XGBoost") args$verbose <- 0
  model <- do.call(caret::train, args)
  prob <- as.numeric(stats::predict(model, newdata=x_test, type="prob")[,"Sepsis"])
  imp <- extract_model_importance(model, algorithm)
  list(model=model, probability=prob, importance=imp)
}

extract_model_importance <- function(model, algorithm) {
  if (algorithm == "LASSO") {
    cc <- as.matrix(stats::coef(model$finalModel, s=model$bestTune$lambda))[,1]
    cc <- cc[names(cc) != "(Intercept)" & cc != 0]
    out <- abs(cc)
  } else if (algorithm == "glmBoost") {
    cc <- stats::coef(model$finalModel)
    cc <- unlist(cc); cc <- cc[names(cc) != "(Intercept)" & cc != 0]
    out <- abs(cc)
  } else if (algorithm == "Ranger") {
    out <- model$finalModel$variable.importance
  } else {
    xi <- xgboost::xgb.importance(model=model$finalModel)
    out <- xi$Gain; names(out) <- xi$Feature
  }
  out <- sort(out[is.finite(out)], decreasing=TRUE)
  out
}

pr_auc_score <- function(truth, prob, positive="Sepsis") {
  y <- as.numeric(truth == positive)
  if (length(unique(y)) < 2) return(NA_real_)
  as.numeric(PRROC::pr.curve(scores.class0=prob[y==1], scores.class1=prob[y==0], curve=FALSE)$auc.integral)
}

evaluate_probabilities <- function(truth, prob, positive="Sepsis") {
  truth <- factor(truth, levels=c("Normal","Sepsis"))
  roc_obj <- pROC::roc(truth, prob, levels=c("Normal","Sepsis"), direction="<", quiet=TRUE)
  p_clip <- pmin(pmax(prob,1e-6),1-1e-6)
  slope <- tryCatch(stats::coef(stats::glm(as.numeric(truth==positive) ~ stats::qlogis(p_clip), family=stats::binomial()))[2],
                    error=function(e) NA_real_)
  c(ROC_AUC=as.numeric(pROC::auc(roc_obj)), PR_AUC=pr_auc_score(truth,prob,positive),
    binary_metrics(truth,prob,0.5,positive), calibration_slope=unname(slope))
}

stratified_outer_folds <- function(y, k=5L, repeats=20L, seed=MASTER_SEED) {
  set.seed(seed)
  caret::createMultiFolds(factor(y), k=k, times=repeats)
}

run_outer_fold <- function(x, y, train_idx, algorithms=c("LASSO","glmBoost","Ranger","XGBoost"),
                           strategy=c("weighted","unweighted","downsample"), seed=1L, top_n=100L) {
  strategy <- match.arg(strategy)
  test_idx <- setdiff(seq_len(nrow(x)), train_idx)
  tr_idx <- train_idx
  if (strategy == "downsample") {
    set.seed(seed)
    ids0 <- tr_idx[y[tr_idx]=="Normal"]; ids1 <- tr_idx[y[tr_idx]=="Sepsis"]
    nmin <- min(length(ids0),length(ids1))
    tr_idx <- c(sample(ids0,nmin),sample(ids1,nmin))
  }
  fs <- select_training_features(x[tr_idx,,drop=FALSE], y[tr_idx], top_n=top_n)
  genes <- fs$genes
  st <- scale_train_test(x[tr_idx,genes,drop=FALSE],x[test_idx,genes,drop=FALSE])
  original_names <- colnames(st$train)
  safe_names <- make.names(original_names, unique=TRUE)
  name_map <- stats::setNames(original_names, safe_names)
  colnames(st$train) <- safe_names; colnames(st$test) <- safe_names
  results <- lapply(seq_along(algorithms), function(ai) {
    alg <- algorithms[ai]
    fit <- fit_inner_model(st$train,y[tr_idx],st$test,alg,weighted=(strategy=="weighted"),seed=seed+ai*1000L)
    met <- evaluate_probabilities(y[test_idx],fit$probability)
    imp <- fit$importance
    if(length(imp)) names(imp) <- unname(name_map[names(imp)])
    ranks <- if(length(imp)) seq_along(imp) else integer()
    feat <- data.frame(algorithm=alg,gene=names(imp),importance=as.numeric(imp),rank=ranks,
                       selected=if(alg %in% c("Ranger","XGBoost")) ranks<=20 else TRUE,
                       stringsAsFactors=FALSE)
    list(metrics=data.frame(algorithm=alg,t(as.data.frame(met)),check.names=FALSE),
         predictions=data.frame(sample_id=rownames(x)[test_idx],truth=y[test_idx],probability=fit$probability,algorithm=alg),
         features=feat,screened_genes=genes,best_tune=fit$model$bestTune)
  })
  results
}

fit_fast_stability_model <- function(x_train, y_train, x_test, algorithm, seed=1L) {
  set.seed(seed)
  y01 <- as.numeric(y_train=="Sepsis")
  if(algorithm=="LASSO") {
    fit <- glmnet::cv.glmnet(as.matrix(x_train),y01,family="binomial",alpha=1,nfolds=3,
                             type.measure="deviance",standardize=FALSE)
    prob <- as.numeric(stats::predict(fit,newx=as.matrix(x_test),s="lambda.1se",type="response"))
    cc <- as.matrix(stats::coef(fit,s="lambda.1se"))[,1]; cc<-cc[names(cc)!="(Intercept)" & cc!=0]
    imp <- sort(abs(cc),decreasing=TRUE)
  } else if(algorithm=="glmBoost") {
    fit <- mboost::glmboost(x=as.matrix(x_train),y=factor(y_train,levels=c("Normal","Sepsis")),family=mboost::Binomial(),
                            control=mboost::boost_control(mstop=100,nu=.1))
    prob <- as.numeric(stats::predict(fit,newdata=as.matrix(x_test),type="response"))
    cc <- unlist(stats::coef(fit)); cc<-cc[names(cc)!="(Intercept)" & cc!=0]
    imp <- sort(abs(cc),decreasing=TRUE)
  } else if(algorithm=="Ranger") {
    dat <- data.frame(.y=factor(y_train,levels=c("Normal","Sepsis")),x_train,check.names=FALSE)
    fit <- ranger::ranger(.y~.,data=dat,probability=TRUE,num.trees=300,
                          mtry=max(1,round(sqrt(ncol(x_train)))),importance="permutation",seed=seed)
    prob <- as.numeric(stats::predict(fit,data=as.data.frame(x_test,check.names=FALSE))$predictions[,"Sepsis"])
    imp <- sort(fit$variable.importance,decreasing=TRUE)
  } else if(algorithm=="XGBoost") {
    dtr <- xgboost::xgb.DMatrix(as.matrix(x_train),label=y01)
    fit <- xgboost::xgb.train(params=list(objective="binary:logistic",eval_metric="logloss",max_depth=2,
                                          eta=.08,subsample=.8,colsample_bytree=.8,nthread=1),
                              data=dtr,nrounds=60,verbose=0)
    prob <- as.numeric(stats::predict(fit,as.matrix(x_test)))
    xi <- xgboost::xgb.importance(model=fit); imp<-xi$Gain; names(imp)<-xi$Feature; imp<-sort(imp,decreasing=TRUE)
  } else stop("未知算法",call.=FALSE)
  list(probability=prob,importance=imp)
}
