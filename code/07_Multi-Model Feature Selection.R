# =========================================================================
# Script Name: 07_Multi-Model Feature Selection.R
# Description: Ensemble Machine Learning Feature Selection Pipeline (v5.0).
#              Integrates LASSO, Random Forest, glmBOOST, and XGBoost.
#              Features targeted SHAP interpretability (Top 6 genes) and 
#              specific biomarker (MAFG/S100A8) ROC validation.
# Note for GitHub users: Please update the working directory before running.
# =========================================================================

rm(list = ls())

# ====================== Environment Setup ======================
# Set working directory (Update this to your local repository path)
setwd("D:\\桌面\\sepsis-s100a8-multiomics\\data\\07_Multi-Model Feature Selection\\Input")

# Load and install required packages
required_packages <- c(
  "data.table", "dplyr", "glmnet", "caret", "randomForest", "mboost",
  "pROC", "ggplot2", "UpSetR", "pheatmap", "RColorBrewer",
  "ggpubr", "reshape2", "doParallel", "patchwork", "cowplot", 
  "ggsci", "reshape", "xgboost", "fastshap", "ggbeeswarm", 
  "shapviz", "kernelshap"
)

invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}))

# Set global ggplot theme after ensuring cowplot is loaded
theme_set(theme_cowplot())

# Set global seed for reproducibility
set.seed(123)

# ====================== Data Loading & Preprocessing ======================
data <- read.csv("Train_Data.csv", header = TRUE, row.names = 1)
if (!"Group" %in% colnames(data)) {
  colnames(data)[ncol(data)] <- "Group"
}
data <- as.data.frame(data)
data$Group <- factor(as.character(data$Group))
# Ensure 'Normal' is set as the baseline reference level
data$Group <- relevel(data$Group, ref = "Normal") 

# Extract expression matrix and remove Near-Zero Variance (NZV) features
expr <- data[, setdiff(colnames(data), "Group"), drop = FALSE]
nzv <- caret::nearZeroVar(expr)
if (length(nzv) > 0) expr <- expr[, -nzv, drop = FALSE]

# Impute missing values with median if any exist
if (anyNA(expr)) {
  for (j in seq_len(ncol(expr))) expr[is.na(expr[, j]), j] <- median(expr[, j], na.rm = TRUE)
}

# Scale expression data
expr_scaled <- scale(expr)

# Create output directory for final plots
dir.create("FinalPlots", showWarnings = FALSE)

# Define common caret trainControl for cross-validation
control <- trainControl(method = "repeatedcv", number = 10, repeats = 3,
                        summaryFunction = twoClassSummary, classProbs = TRUE,
                        savePredictions = "final", verboseIter = FALSE,
                        allowParallel = TRUE) 

# ==============================================================================
# 1) LASSO Regression (glmnet) 
# ==============================================================================
cat("Running LASSO (glmnet)...\n")

x_matrix <- as.matrix(expr_scaled)
y_vector <- as.numeric(data$Group == "Disease")

set.seed(123)
fit <- glmnet(x_matrix, y_vector, family = "binomial", alpha = 1, maxit = 10000)
cvfit <- cv.glmnet(x_matrix, y_vector, family = "binomial", alpha = 1, maxit = 10000)

best_lambda <- cvfit$lambda.min
coefs <- coef(fit, s = best_lambda)
lasso_genes <- setdiff(rownames(coefs)[which(coefs != 0)], "(Intercept)")

# Extract LASSO 1se genes
coefs_1se <- coef(fit, s = cvfit$lambda.1se)
lasso_genes_1se <- setdiff(rownames(coefs_1se)[which(coefs_1se != 0)], "(Intercept)")

# --- 1.1 LASSO Coefficient Path Plot ---
x_coef <- coef(fit)
tmp <- as.data.frame(as.matrix(x_coef))
tmp$coef <- rownames(tmp)
tmp <- reshape::melt(tmp, id = "coef")
tmp$variable <- as.numeric(gsub("s", "", tmp$variable))
tmp$coef <- gsub('_','-',tmp$coef)
tmp$lambda <- fit$lambda[tmp$variable+1]

min_y_path <- min(tmp$value)

p_lasso_path <- ggplot(tmp, aes(log(lambda), value, color = coef)) +
  geom_vline(xintercept = log(cvfit$lambda.min), size=1, color='red', alpha=0.8, linetype=2) +
  geom_vline(xintercept = log(cvfit$lambda.1se), size=1, color='blue', alpha=0.8, linetype=2) +
  geom_line(size=0.8) + 
  xlab(expression(paste("Log(", lambda, ")"))) + 
  ylab('Coefficients (Regularized)') +
  theme_bw(base_size = 14) + 
  scale_color_discrete() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=14, color='black'),
        axis.text = element_text(size=12, color='black'),
        legend.title = element_blank(),
        legend.position = 'none') +
  annotate("text", x = log(cvfit$lambda.min), y = min_y_path * 0.95, label = "lambda_min", 
           color = 'red', hjust = 1.1, vjust = 1, size = 4) +
  annotate("text", x = log(cvfit$lambda.1se), y = min_y_path * 0.95, label = "lambda_1se", 
           color = 'blue', hjust = -0.1, vjust = 1, size = 4)

ggsave("FinalPlots/1_LASSO_Path_Classic_Labeled_Text.pdf", plot = p_lasso_path, width = 16, height = 12)

# --- 1.2 LASSO CV Error Plot ---
xx <- data.frame(lambda=cvfit[["lambda"]],
                 cvm=cvfit[["cvm"]],
                 cvsd=cvfit[["cvsd"]],
                 cvup=cvfit[["cvup"]],
                 cvlo=cvfit[["cvlo"]],
                 nozezo=cvfit[["nzero"]])
xx$ll <- log(xx$lambda)

min_y_cv <- min(xx$cvlo)

p_lasso_cv <- ggplot(xx, aes(ll, cvm)) + 
  geom_errorbar(aes(x=ll, ymin=cvlo, ymax=cvup), width=0.05, size=0.8, alpha=0.5, color="gray50") +
  geom_vline(xintercept = log(cvfit$lambda.min), size=1, color='red', alpha=0.8, linetype=2) +
  geom_vline(xintercept = log(cvfit$lambda.1se), size=1, color='blue', alpha=0.8, linetype=2) +
  geom_point(size=2, color="black") +
  xlab(expression(paste("Log(", lambda, ")"))) +
  ylab('Cross-Validated Error (Deviance)') +
  theme_bw(base_size = 14) + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=14, color='black'),
        axis.text = element_text(size=12, color='black'),
        legend.position = 'none') + 
  scale_x_continuous(sec.axis = sec_axis(~ ., breaks = xx$ll, labels = xx$nozezo, name = "Number of Nonzero Variables")) +
  annotate("text", x = log(cvfit$lambda.min), y = min_y_cv * 0.999, label = "lambda_min", 
           color = 'red', hjust = 1.1, vjust = 1, size = 4) +
  annotate("text", x = log(cvfit$lambda.1se), y = min_y_cv * 0.999, label = "lambda_1se", 
           color = 'blue', hjust = -0.1, vjust = 1, size = 4)

ggsave("FinalPlots/2_LASSO_CV_Classic_Labeled_Text.pdf", plot = p_lasso_cv, width = 16, height = 12)

# ==============================================================================
# 2) Random Forest (RF) 
# ==============================================================================
cat("Running Random Forest...\n")
set.seed(123)

rf_train <- train(x = expr_scaled, y = data$Group, method = "rf",
                  metric = "ROC", trControl = control, importance = TRUE, tuneLength = 5)

rf_imp_raw <- importance(rf_train$finalModel) 
Importance_Metric <- if("MeanDecreaseGini" %in% colnames(rf_imp_raw)) rf_imp_raw[, "MeanDecreaseGini", drop = FALSE] else rf_imp_raw[, "MeanDecreaseAccuracy", drop = FALSE]
rf_imp_df <- data.frame(Gene = rownames(Importance_Metric), Importance = Importance_Metric[, 1], row.names = NULL) %>% arrange(desc(Importance))
rf_genes <- head(rf_imp_df$Gene, 15)

# --- 2.1 RF OOB Error Curve Plot ---
oob_error_data <- as.data.frame(rf_train$finalModel$err.rate)
oob_error_data$Ntree <- 1:nrow(oob_error_data)
oob_melted <- reshape2::melt(oob_error_data, id.vars = "Ntree", variable.name = "ErrorType", value.name = "ErrorRate")

p_rf_oob <- ggplot(oob_melted, aes(x = Ntree, y = ErrorRate, color = ErrorType)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("OOB" = "#E64B35FF", "Normal" = "#00A087FF", "Disease" = "#3C5488FF"),
                     labels = c("OOB Error", "Normal Error", "Disease Error")) +
  labs(x = "Number of Trees", y = "Error Rate", color = "Error Type") +
  theme_classic(base_size = 14) +
  theme(legend.position = c(0.8, 0.8), 
        legend.background = element_rect(fill = "white", color = "black"),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5))

ggsave("FinalPlots/3B_RF_OOB_Error_Curve.pdf", plot = p_rf_oob, width = 8, height = 6)

# ==============================================================================
# 3) glmBOOST (mboost)
# ==============================================================================
cat("Running glmBOOST (mboost)...\n")
set.seed(123)

glmboost_train <- train(x = expr_scaled, y = data$Group, method = "glmboost",
                        metric = "ROC", trControl = control, tuneLength = 5)

final_model <- glmboost_train$finalModel
glmboost_coefs_raw <- coef(final_model) 
glmboost_coef_df <- data.frame(Gene = names(glmboost_coefs_raw), Coefficient = as.vector(glmboost_coefs_raw)) %>%
  filter(Gene != "(Intercept)", Coefficient != 0) %>% arrange(desc(abs(Coefficient))) 
glmboost_genes <- glmboost_coef_df$Gene 

if (length(glmboost_genes) == 0) {
  glmboost_imp <- varImp(glmboost_train, scale = FALSE)
  glmboost_coef_df <- data.frame(Gene = rownames(glmboost_imp$importance), Coefficient = glmboost_imp$importance$Overall) %>% arrange(desc(Coefficient))
  glmboost_genes <- head(glmboost_coef_df$Gene, 50)
}

# --- 3.1 glmBOOST CV Risk Curve ---
cv_res <- cvrisk(final_model) 
mstop_opt <- mstop(cv_res)
risk_vector <- as.numeric(cv_res)

risk_df <- data.frame(
  Iteration = 1:length(risk_vector),
  Risk = risk_vector
)

p_boost_cv <- ggplot(risk_df, aes(x = Iteration, y = Risk)) +
  geom_line(color = "#336699", linewidth = 1) +
  geom_vline(xintercept = mstop_opt, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(x = "Boosting Iterations (m)", y = "Cross-Validated Risk") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5))

ggsave("FinalPlots/3A_GLMBOOST_CV_Risk_Curve.pdf", plot = p_boost_cv, width = 8, height = 6)

# --- 3.2 glmBOOST Coefficient Path Plot ---
coef_path_raw <- mboost:::coef.mboost(final_model, which = "", by = 1, off2int = TRUE)
coef_path_matrix <- do.call(cbind, coef_path_raw)
colnames(coef_path_matrix) <- paste0("m", 1:ncol(coef_path_matrix))
coef_path_df <- as.data.frame(coef_path_matrix)
coef_path_df$Gene <- rownames(coef_path_df)

coef_path_melted <- reshape2::melt(coef_path_df, id.vars = "Gene", variable.name = "Iteration_Str", value.name = "Coefficient")
coef_path_melted$Iteration <- as.numeric(gsub("m", "", coef_path_melted$Iteration_Str))

coef_path_plot_data <- coef_path_melted %>%
  filter(Gene != "(Intercept)", Iteration <= mstop_opt, Coefficient != 0)

p_boost_path <- ggplot(coef_path_plot_data, aes(x = Iteration, y = Coefficient, group = Gene, color = Gene)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_vline(xintercept = mstop_opt, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(x = "Boosting Iterations (m)", y = "Coefficient Value") +
  theme_classic(base_size = 14) +
  theme(plot.title = element_blank(),
        legend.position = "none", 
        axis.line = element_line(colour = "black", linewidth = 0.5),
        axis.ticks = element_line(colour = "black", linewidth = 0.5))

ggsave("FinalPlots/3B_GLMBOOST_Coefficient_Path.pdf", plot = p_boost_path, width = 8, height = 6)


# ==============================================================================
# 4) XGBoost (xgbTree) 
# ==============================================================================
cat("Running XGBoost...\n")
set.seed(123)

xgb_train <- train(x = as.matrix(expr_scaled), y = data$Group, method = "xgbTree",
                   metric = "ROC", trControl = control, tuneLength = 5)

importance_matrix <- xgb.importance(feature_names = colnames(expr_scaled), model = xgb_train$finalModel)
xgb_imp_df <- data.frame(Gene = importance_matrix$Feature, Gain = importance_matrix$Gain) %>% arrange(desc(Gain))
xgb_genes <- head(xgb_imp_df$Gene, 50)


# ==============================================================================
# 5) Integration and Comparison (ROC & UpSet)
# ==============================================================================
cat("Starting Integration Analysis...\n")
lasso_caret <- train(x = as.matrix(expr_scaled), y = data$Group, method = "glmnet", trControl = control, metric = "ROC", tuneLength = 5)

roc_lasso <- pROC::roc(response = lasso_caret$pred$obs, predictor = lasso_caret$pred$Disease, levels = rev(levels(lasso_caret$pred$obs)))
roc_rf <- pROC::roc(response = rf_train$pred$obs, predictor = rf_train$pred$Disease, levels = rev(levels(rf_train$pred$obs)))
roc_xgb <- pROC::roc(response = xgb_train$pred$obs, predictor = xgb_train$pred$Disease, levels = rev(levels(xgb_train$pred$obs)))
roc_glmboost <- pROC::roc(response = glmboost_train$pred$obs, predictor = glmboost_train$pred$Disease, levels = rev(levels(glmboost_train$pred$obs)))

pdf("FinalPlots/5A_MultiModel_CV_ROC_Comparison.pdf", width = 7, height = 6)
plot(roc_lasso, col = "#E64B35FF", lwd = 2, main = "", legacy.axes = TRUE)
plot(roc_rf, add = TRUE, col = "#00A087FF", lwd = 2) 
plot(roc_glmboost, add = TRUE, col = "#3C5488FF", lwd = 2) 
plot(roc_xgb, add = TRUE, col = "#ff7f0e", lwd = 2)
legend("bottomright", 
       legend = c(paste0("LASSO (AUC=", round(auc(roc_lasso),3),")"),
                  paste0("RF (AUC=", round(auc(roc_rf),3),")"),
                  paste0("glmBOOST (AUC=", round(auc(roc_glmboost),3),")"),
                  paste0("XGBoost (AUC=", round(auc(roc_xgb),3),")")),
       col = c("#E64B35FF","#00A087FF","#3C5488FF","#ff7f0e"), lwd = 2, cex = 0.8)
dev.off()

gene_lists_all <- list(LASSO = lasso_genes, RF = rf_genes, GLMBOOST = glmboost_genes, XGBOOST = xgb_genes)
gene_lists_all <- gene_lists_all[sapply(gene_lists_all, length) > 0] 

pdf("FinalPlots/5B_UpSet_Feature_Overlap.pdf", width = 10, height = 7)
UpSetR::upset(UpSetR::fromList(gene_lists_all), sets = names(gene_lists_all), order.by = "freq", main.bar.color = "#e31a1c", sets.bar.color = "#1f78b4", text.scale = 1.5)      
dev.off()

# ==============================================================================
# 6) Unified Style Plotting Module (Gradient Style, No Titles)
# ==============================================================================
cat("Regenerating all Bar Plots with Unified Gradient Style...\n")
col_gradient_low <- "#99CCCC"  
col_gradient_high <- "#3C5488FF" 
my_bar_theme <- theme_classic(base_size = 14) + theme(plot.title = element_blank(), axis.text.y = element_text(face = "italic", color = "black"), axis.text.x = element_text(color = "black"), axis.line = element_line(colour = "black", linewidth = 0.6), legend.position = "none")

# LASSO
lasso_coefs_raw <- coef(fit, s = cvfit$lambda.min) 
lasso_imp_df <- data.frame(Gene = rownames(lasso_coefs_raw), Coefficient = as.vector(lasso_coefs_raw)) %>%
  filter(Gene != "(Intercept)", Coefficient != 0)
lasso_imp_df1 <- lasso_imp_df %>% arrange(desc(Coefficient)) %>% head(5)
p_lasso_bar <- ggplot(lasso_imp_df1, aes(x = reorder(Gene, Coefficient), y = Coefficient, fill = Coefficient)) +
  geom_col(width = 0.7) + coord_flip() + scale_fill_gradient(low = col_gradient_low, high = col_gradient_high) + labs(x = "", y = "LASSO Coefficient") + my_bar_theme
ggsave("FinalPlots/1_LASSO_Coefficients_Barplot.pdf", plot = p_lasso_bar, width = 6, height = 5)

# Random Forest
rf_top5 <- head(rf_imp_df, 5)
p_rf_bar <- ggplot(rf_top5, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance)) +
  geom_col(width = 0.8) + coord_flip() + scale_fill_gradient(low = col_gradient_low, high = col_gradient_high) + labs(x = "", y = "Feature Importance (Gini)") + my_bar_theme
ggsave("FinalPlots/3A_RF_Importance_Top5.pdf", plot = p_rf_bar, width = 6, height = 5)

# glmBOOST
if (nrow(glmboost_coef_df) > 0) {
  boost_plot_data <- head(glmboost_coef_df, 5)
  p_boost_bar <- ggplot(boost_plot_data, aes(x = reorder(Gene, Coefficient), y = Coefficient, fill = Coefficient)) +
    geom_col(width = 0.7) + coord_flip() + scale_fill_gradient(low = col_gradient_low, high = col_gradient_high) + labs(x = "", y = "glmBoost Coefficient") + my_bar_theme
  ggsave("FinalPlots/3C_GLMBOOST_Final_Coefficients.pdf", plot = p_boost_bar, width = 6, height = 5)
}

# XGBoost
xgb_top5 <- head(xgb_imp_df, 5)
p_xgb_bar <- ggplot(xgb_top5, aes(x = reorder(Gene, Gain), y = Gain, fill = Gain)) +
  geom_col(width = 0.8) + coord_flip() + scale_fill_gradient(low = col_gradient_low, high = col_gradient_high) + labs(x = "", y = "Feature Importance (Gain)") + my_bar_theme
ggsave("FinalPlots/4A_XGBoost_Importance_Top5.pdf", plot = p_xgb_bar, width = 6, height = 5)


# ==============================================================================
# 7) Specific Core Gene ROC Validation (MAFG, S100A8)
# ==============================================================================
cat("Generating Single Gene ROC Curves for MAFG and S100A8...\n")

if("MAFG" %in% colnames(data) & "S100A8" %in% colnames(data)) {
  roc_mafg <- pROC::roc(response = data$Group, predictor = data$MAFG, levels = c("Normal", "Disease"), direction = "<")
  roc_s100a8 <- pROC::roc(response = data$Group, predictor = data$S100A8, levels = c("Normal", "Disease"), direction = "<")
  
  pdf("FinalPlots/7_Single_Gene_ROC_MAFG_S100A8.pdf", width = 6, height = 6)
  plot(roc_mafg, col = "red", lwd = 1.5, legacy.axes = TRUE, xlab = "1 - Specificity", ylab = "Sensitivity", main = "")
  plot(roc_s100a8, add = TRUE, col = "green4", lwd = 1.5)
  legend("bottomright", 
         legend = c(paste0("MAFG (AUC=", sprintf("%.3f", roc_mafg$auc), ")"),
                    paste0("S100A8 (AUC=", sprintf("%.3f", roc_s100a8$auc), ")")),
         col = c("red", "green4"), text.col = c("black", "red"), lwd = 1.5, bty = "n", cex = 1.1) 
  dev.off()
} else {
  message("Warning: MAFG or S100A8 not found in dataset.")
}


# ==============================================================================
# 8) SHAP Interpretability (Fixed 6 Genes)
# ==============================================================================
cat("Running SHAP Interpretability Analysis...\n")

X_train <- as.matrix(expr_scaled)
pred_fun <- function(object, newdata) predict(object, newdata = newdata, type = "prob")[, "Disease"]

set.seed(123)
bg_size <- min(520, nrow(X_train))
shap_values <- kernelshap(object = xgb_train, X = X_train, bg_X = X_train[1:bg_size, ], pred_fun = pred_fun)

feature_importance <- colMeans(abs(as.matrix(shap_values$S)))
top6_genes <- head(names(sort(feature_importance, decreasing = TRUE)), 6)
shap_top6 <- shapviz(shap_values$S[, top6_genes], X = X_train[, top6_genes])

viz_theme <- theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14),
        axis.title = element_text(size = 12),
        panel.background = element_rect(fill = "white", color = NA), 
        plot.background = element_rect(fill = "white", color = NA))

pdf("FinalPlots/8A_SHAP_Feature_Importance_Barplot.pdf", width = 8, height = 6)
print(sv_importance(shap_top6, kind = "bar", max_display = 6, show_numbers = TRUE) + viz_theme + labs(title = "", x = "Mean |SHAP value|", y = "Feature"))
dev.off()

pdf("FinalPlots/8B_SHAP_BeeSwarm_Plot.pdf", width = 9, height = 7)
print(sv_importance(shap_top6, kind = "bee", max_display = 6) + viz_theme + labs(title = "", x = "SHAP value", y = "Feature"))
dev.off()

pdf("FinalPlots/8C_SHAP_Feature_Dependence.pdf", width = 10, height = 8)
print(sv_dependence(shap_top6, v = top6_genes) + viz_theme)
dev.off()

pdf("FinalPlots/8D_SHAP_Sample_Waterfall.pdf", width = 9, height = 6)
print(sv_waterfall(shap_top6, row_id = 2) + theme(panel.background = element_rect(fill = "white", color = NA)) + labs(title = "", subtitle = "", x = "Feature contribution", y = "Prediction"))
dev.off()

cat("Pipeline Completed! All plots saved in 'FinalPlots/' folder.\n")


