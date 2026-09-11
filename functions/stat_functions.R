# 统计量与诊断性能公共函数
standardized_mean_difference <- function(x1, x0) {
  n1 <- sum(is.finite(x1)); n0 <- sum(is.finite(x0))
  s_pool <- sqrt(((n1 - 1) * stats::var(x1, na.rm = TRUE) + (n0 - 1) * stats::var(x0, na.rm = TRUE)) / (n1 + n0 - 2))
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / s_pool
}

binary_metrics <- function(truth, prob, threshold = 0.5, positive = "Sepsis") {
  truth <- factor(truth)
  pred <- ifelse(prob >= threshold, positive, setdiff(levels(truth), positive)[1])
  tp <- sum(pred == positive & truth == positive); tn <- sum(pred != positive & truth != positive)
  fp <- sum(pred == positive & truth != positive); fn <- sum(pred != positive & truth == positive)
  sens <- if ((tp + fn) > 0) tp/(tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0) tn/(tn + fp) else NA_real_
  ppv <- if ((tp + fp) > 0) tp/(tp + fp) else NA_real_
  npv <- if ((tn + fn) > 0) tn/(tn + fn) else NA_real_
  f1 <- if (is.finite(ppv + sens) && (ppv + sens) > 0) 2*ppv*sens/(ppv+sens) else NA_real_
  mcc_d <- sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
  mcc <- if (mcc_d > 0) (tp*tn-fp*fn)/mcc_d else NA_real_
  c(sensitivity=sens, specificity=spec, balanced_accuracy=mean(c(sens,spec),na.rm=TRUE),
    PPV=ppv, NPV=npv, F1=f1, MCC=mcc, Brier=mean((as.numeric(truth==positive)-prob)^2))
}

