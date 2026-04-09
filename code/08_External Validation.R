# =========================================================================
# Script Name: 08_External Validation.R
# Description: Targeted external validation pipeline for the identified biomarker.
#              Includes robust probe mapping, smart metadata parsing, and 
#              comprehensive performance evaluation (ROC, PR, Calibration, Brier).
# Note for GitHub users: Please update the working directory before running.
# =========================================================================

# 1. Environment Setup and Dependency Check
required_packages <- c("GEOquery", "limma", "dplyr", "ggplot2", "ggpubr", "pROC", "RColorBrewer", "precrec", "rms")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg, character.only = TRUE)
}

library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(pROC)
library(RColorBrewer)
library(precrec)
library(rms)

# Set working directory
# setwd("D:/sepsis-s100a8-multiomics/data/08_External_Validation")

# Global Configuration
biomarker_symbol <- "S100A8"
independent_cohorts <- c("GSE9692", "GSE26440", "GSE28750", "GSE69528")

# ==============================================================================
# 2. Core Processing Functions (Methodological Robustness)
# ==============================================================================

# 2.1 Super-Robust Probe Mapping Function
# Addresses the inconsistency in platform annotations (GPL)
find_probe_robust <- function(gset, target = "S100A8") {
  platform <- annotation(gset)
  
  # A. Hard-coded mapping for common platforms to ensure 100% accuracy
  if (platform %in% c("GPL570", "GPL571", "GPL96", "GPL1352", "GPL97", "GPL13667")) {
    return(c("202917_s_at", "214370_at", "202917_at"))
  } else if (platform %in% c("GPL6947", "GPL10558", "GPL6883", "GPL6884")) {
    return(c("ILMN_1693114", "ILMN_2038777", "ILMN_1762590"))
  }
  
  # B. Heuristic search in feature data (fData)
  fd <- fData(gset)
  if (is.null(fd) || nrow(fd) == 0) return(NULL)
  
  candidate_probes <- c()
  check_cols <- grep("symbol|gene|title|description|name", colnames(fd), ignore.case = TRUE)
  if(length(check_cols) == 0) check_cols <- 1:min(ncol(fd), 5)
  
  for (col_idx in check_cols) {
    col_val <- as.character(fd[, col_idx])
    matches <- grep(paste0("^", target, "$"), col_val, ignore.case = TRUE)
    if (length(matches) > 0) { candidate_probes <- rownames(fd)[matches]; break }
    
    matches_fuzzy <- grep(target, col_val, ignore.case = TRUE)
    if (length(matches_fuzzy) > 0) {
      is_real <- grepl(paste0("\\b", target, "\\b"), col_val[matches_fuzzy], ignore.case = TRUE)
      if (any(is_real)) { candidate_probes <- rownames(fd)[matches_fuzzy[is_real]]; break }
    }
  }
  return(candidate_probes)
}

# 2.2 Intelligent Metadata-based Sample Grouping
assign_groups_intelligent <- function(pdata) {
  cols_to_use <- intersect(c("title", "source_name_ch1", "characteristics_ch1", "description"), colnames(pdata))
  if(length(cols_to_use) == 0) return(rep(NA, nrow(pdata)))
  
  meta_text <- tolower(apply(pdata[, cols_to_use, drop=FALSE], 1, function(x) paste(x, collapse = " ")))
  group_vec <- rep(NA, nrow(pdata))
  
  # Keyword-based classification
  k_norm <- c("healthy", "control", "normal", "volunteer", "donor", "non-septic")
  k_sep  <- c("sepsis", "septic", "shock", "infection", "pneumonia", "bacterial", "cap")
  
  group_vec[grepl(paste(k_norm, collapse = "|"), meta_text)] <- "Normal"
  group_vec[grepl(paste(k_sep, collapse = "|"), meta_text) & is.na(group_vec)] <- "Sepsis"
  
  return(group_vec)
}

# 2.3 Main Pipeline for Individual GSE Processing
process_external_cohort <- function(gse_id) {
  filename <- paste0(gse_id, "_series_matrix.txt.gz")
  message(paste0("\n>>> Analyzing Cohort [", gse_id, "]..."))
  
  # A. Load Data (Support both local files and online download)
  gset <- tryCatch({
    if (file.exists(filename)) {
      message("  - Local file detected. Loading...")
      getGEO(filename = filename, getGPL = FALSE) 
    } else {
      message("  - Local file missing. Attempting online download...")
      res <- getGEO(GEO = gse_id, destdir = getwd(), getGPL = FALSE)
      if(inherits(res, "list")) res[[1]] else res
    }
  }, error = function(e) return(NULL))
  
  if(is.null(gset)) { message("  - Error: Failed to load ExpressionSet."); return(NULL) }
  
  # B. Matrix Extraction and Adaptive Log2 Transformation
  expr_mat <- exprs(gset)
  qx <- as.numeric(quantile(expr_mat, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  if ((qx[5] > 100) || (qx[6]-qx[1] > 50 && qx[2] > 0)) {
    expr_mat[expr_mat <= 0] <- NaN
    expr_mat <- log2(expr_mat + 1)
    message("  - Automated Log2 transformation applied.")
  }
  
  # C. Probe Identification and Value Extraction
  probes <- find_probe_robust(gset, biomarker_symbol)
  valid_probes <- probes[probes %in% rownames(expr_mat)]
  if (length(valid_probes) == 0) { message("  - Skipped: No matching probes found."); return(NULL) }
  
  # Use max value if multiple probes match the gene
  sub_mat <- expr_mat[valid_probes, , drop = FALSE]
  val <- if(nrow(sub_mat) > 1) apply(sub_mat, 2, max, na.rm=T) else as.numeric(sub_mat)
  
  # D. Grouping and Data Cleaning
  group_vec <- assign_groups_intelligent(pData(gset))
  df <- data.frame(Expression = val, Group = group_vec) %>% 
    filter(!is.na(Group)) %>%
    filter(is.finite(Expression))
  
  if (sum(df$Group=="Normal") < 3 || sum(df$Group=="Sepsis") < 3) {
    message("  - Skipped: Insufficient sample size or unidentified groups."); return(NULL)
  }
  
  df$Group <- factor(df$Group, levels = c("Normal", "Sepsis"))
  
  # E. ROC Performance Calculation
  roc_obj <- tryCatch({
    pROC::roc(response = df$Group, predictor = df$Expression, 
              levels = c("Normal", "Sepsis"), direction = "<", quiet = TRUE)
  }, error = function(e) return(NULL))
  
  if(is.null(roc_obj) || !inherits(roc_obj, "roc")) {
    message("  - Error: ROC calculation failed."); return(NULL)
  }
  
  auc_val <- as.numeric(pROC::auc(roc_obj))
  message(paste0("  - Success! AUC: ", round(auc_val, 3)))
  
  return(list(gse = gse_id, data = df, roc = roc_obj, auc = auc_val))
}

# ==============================================================================
# 3. Execution Phase
# ==============================================================================
validation_results <- list()

for (gse in independent_cohorts) {
  res <- process_external_cohort(gse)
  if (!is.null(res)) validation_results[[gse]] <- res
}

# ==============================================================================
# 4. Visualization and Statistical Output
# ==============================================================================
if (length(validation_results) > 0) {
  message("\n>>> Generating final publication-quality figures...")
  plot_colors <- colorRampPalette(brewer.pal(8, "Set1"))(length(validation_results))
  names(plot_colors) <- names(validation_results)
  
  # 4.1 Boxplots with Wilcoxon Rank Sum Test
  for (item in validation_results) {
    p <- ggboxplot(item$data, x = "Group", y = "Expression", color = "Group", 
                   palette = c("#3C5488FF", "#E64B35FF"), add = "jitter") +
      stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
      labs(title = item$gse, x = "", y = paste(biomarker_symbol, "Expression")) + 
      theme_classic()
    ggsave(paste0("Validation_Boxplot_", item$gse, ".pdf"), p, width = 4, height = 5)
  }
  
  # 4.2 Combined ROC Curve
  pdf("Validation_Combined_ROC.pdf", width = 6, height = 6)
  plot(0,0, type="n", xlim=c(1,0), ylim=c(0,1), xlab="Specificity", ylab="Sensitivity", main="External Validation: ROC Curve")
  abline(a=1, b=-1, lty=2, col="gray")
  for (i in seq_along(validation_results)) {
    plot(validation_results[[i]]$roc, add=TRUE, col=plot_colors[i], lwd=2)
  }
  legend("bottomright", legend=paste0(names(validation_results), " (AUC=", round(sapply(validation_results, function(x) x$auc), 3), ")"), 
         col=plot_colors, lwd=2, bty="n", cex=0.75)
  dev.off()
  
  # 4.3 Combined Precision-Recall (PR) Curve
  pdf("Validation_Combined_PRC.pdf", width = 6, height = 6)
  plot(0,0, type="n", xlim=c(0,1), ylim=c(0,1), xlab="Recall", ylab="Precision", main="External Validation: PR Curve")
  for (i in seq_along(validation_results)) {
    y_true <- ifelse(validation_results[[i]]$data$Group == "Sepsis", 1, 0)
    scores <- validation_results[[i]]$data$Expression
    precrec_obj <- evalmod(scores = scores, labels = y_true)
    pr_data <- subset(fortify(precrec_obj), curvetype == "PRC")
    lines(pr_data$x, pr_data$y, col=plot_colors[i], lwd=2)
  }
  legend("bottomleft", legend=names(validation_results), col=plot_colors, lwd=2, bty="n", cex=0.75)
  dev.off()
  
  # 4.4 Brier Score Analysis
  # Brier score measures the accuracy of probabilistic predictions
  brier_vals <- sapply(validation_results, function(x) {
    y_binary <- ifelse(x$data$Group == "Sepsis", 1, 0)
    # Min-max normalization of expression to simulate probability (0-1)
    prob_scaled <- (x$data$Expression - min(x$data$Expression))/(max(x$data$Expression)-min(x$data$Expression))
    mean((prob_scaled - y_binary)^2)
  })
  brier_df <- data.frame(GSE = names(validation_results), Brier = brier_vals)
  
  p_brier <- ggplot(brier_df, aes(x=GSE, y=Brier, fill=GSE)) +
    geom_bar(stat="identity", color="black", width = 0.7) + 
    scale_fill_manual(values=plot_colors) +
    geom_text(aes(label=round(Brier, 3)), vjust=-0.3) +
    theme_classic() + labs(title="Model Calibration: Brier Scores", x="", y="Brier Score") +
    theme(legend.position = "none")
  ggsave("Validation_Brier_Scores.pdf", p_brier, width = 6, height = 4)
  
  message(">>> Validation complete. High-resolution PDF outputs saved to working directory.")
} else {
  message(">>> [Warning] No cohorts were successfully processed. Check data connectivity or local files.")
}