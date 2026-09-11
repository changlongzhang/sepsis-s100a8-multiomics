options(stringsAsFactors = FALSE)
.libPaths(c("C:/Users/ZCL/AppData/Local/R/win-library/4.5", .libPaths()))

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(caret)
  library(glmnet)
  library(randomForest)
  library(mboost)
  library(xgboost)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(grid)
})

root <- "C:/Users/ZCL/sepsis_workspace_link"
panel_dir <- file.path(root, "review_revision", "10_figure_sources", "Fig1_S1_corrected_panels")
formal_dir <- file.path(root, "\u6587\u7ae0", "\u8fd4\u4fee\u6295\u7a3f\u65b0\u4e3b\u56fe")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

font_family <- "Arial"
ink <- "#20252A"
mid <- "#7A858D"
light <- "#D9DEE2"
pale <- "#F2F4F5"
navy <- "#315B78"
teal <- "#5B958C"
coral <- "#D56A54"
red <- "#D66A5E"
blue <- "#4C78A8"
green <- "#59A14F"

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(colour = ink),
      axis.text = element_text(size = 8, colour = ink),
      axis.title = element_text(size = 8.5, colour = ink),
      axis.line = element_line(linewidth = 0.35, colour = ink),
      axis.ticks = element_line(linewidth = 0.3, colour = ink),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.height = unit(3.3, "mm"),
      plot.tag = element_text(size = 12, face = "bold", family = font_family, colour = ink),
      plot.tag.position = c(0, 1),
      plot.margin = ggplot2::margin(4, 5, 4, 5)
    )
}

save_gg_panel <- function(plot, stem, width, height, dpi = 300) {
  ggsave(file.path(panel_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf, family = font_family, bg = "white")
  ggsave(file.path(panel_dir, paste0(stem, ".svg")), plot, width = width, height = height,
         device = svglite, bg = "white")
  ggsave(file.path(panel_dir, paste0(stem, ".tif")), plot, width = width, height = height,
         device = ragg::agg_tiff, dpi = dpi, compression = "lzw", bg = "white")
}

cat("STAGE\tLOAD_AND_LIMMA\n")
expr_linear <- read.csv(
  file.path(root, "sepsis-s100a8-multiomics", "data", "02_DEG", "Input", "GSE65682_gene.csv"),
  row.names = 1, check.names = FALSE
)
stopifnot(all(as.matrix(expr_linear) > 0, na.rm = TRUE))
expr <- log2(as.matrix(expr_linear))
meta <- read.csv(
  file.path(root, "sepsis-s100a8-multiomics", "data", "02_DEG", "Input", "GSE65682_Groups.csv"),
  check.names = FALSE
)
colnames(meta)[1:2] <- c("sample_id", "group")
mi <- match(colnames(expr), meta$sample_id)
stopifnot(!anyNA(mi))
group <- factor(meta$group[mi], levels = c("Normal", "Disease"))

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
fit_deg <- eBayes(
  contrasts.fit(lmFit(expr, design), makeContrasts(Disease - Normal, levels = design)),
  robust = TRUE
)
deg <- topTable(fit_deg, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("Gene") %>%
  transmute(
    Gene,
    baseMean = rowMeans(expr[Gene, , drop = FALSE], na.rm = TRUE),
    log2FoldChange = logFC,
    lfcSE = ifelse(t == 0, NA_real_, abs(logFC / t)),
    stat = t,
    pvalue = P.Value,
    padj = adj.P.Val,
    Regulation = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 1 ~ "Up",
      !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE ~ "Not significant"
    )
  )
sig <- deg %>% filter(Regulation != "Not significant")

deg_out <- file.path(root, "sepsis-s100a8-multiomics", "data", "02_DEG", "output")
write.csv(deg, file.path(deg_out, "Disease_vs_Normal_all_results_padj.csv"), row.names = FALSE)
write.csv(sig, file.path(deg_out, "Disease_vs_Normal_sig_genes_padj.csv"), row.names = FALSE)
write.csv(deg, file.path(panel_dir, "Fig1A_source_data.csv"), row.names = FALSE)

black <- unique(trimws(readLines(file.path(
  root, "sepsis-s100a8-multiomics", "data", "03_WGCNA_analysis", "output", "DIFF_MODULE_black_GENES.txt"
), warn = FALSE)))
black <- black[nzchar(black)]
overlap <- intersect(black, sig$Gene)
stopifnot(length(overlap) == 200, "S100A8" %in% overlap)

writeLines(sig$Gene, file.path(root, "sepsis-s100a8-multiomics", "data", "04_Extract overlapping genes", "Input", "DEGs.txt"))
writeLines(overlap, file.path(root, "sepsis-s100a8-multiomics", "data", "05_Enrichment Analysis", "Input", "Overlap_Genes.txt"))
writeLines(overlap, file.path(root, "sepsis-s100a8-multiomics", "data", "06_Prepare ML Data", "Input", "Overlap_Genes.txt"))
writeLines(overlap, file.path(panel_dir, "corrected_overlap_200_genes.txt"))

cat("DEG_N\t", nrow(sig), "\n", sep = "")
cat("OVERLAP_N\t", length(overlap), "\n", sep = "")
cat("S100A8_LOG2FC\t", deg$log2FoldChange[deg$Gene == "S100A8"], "\n", sep = "")
cat("S100A8_FDR\t", deg$padj[deg$Gene == "S100A8"], "\n", sep = "")

cat("STAGE\tFIG1_PANELS\n")
volcano <- ggplot(deg, aes(log2FoldChange, -log10(pmax(padj, 1e-300)))) +
  geom_point(aes(colour = Regulation), size = 1.05, alpha = 0.65) +
  geom_vline(xintercept = c(-1, 1), colour = mid, linewidth = 0.35, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), colour = mid, linewidth = 0.35, linetype = "dashed") +
  geom_point(
    data = filter(deg, Gene == "S100A8"), shape = 21, size = 2.7,
    stroke = 0.65, colour = "black", fill = NA
  ) +
  geom_text_repel(
    data = filter(deg, Gene == "S100A8"), aes(label = Gene),
    size = 8 / ggplot2::.pt, family = font_family, fontface = "bold",
    nudge_x = -0.7, nudge_y = -8, box.padding = 0.5, point.padding = 0.3,
    min.segment.length = 0, segment.size = 0.35, seed = 123
  ) +
  scale_colour_manual(
    values = c("Down" = blue, "Not significant" = "#C8CDD0", "Up" = red),
    breaks = c("Up", "Down", "Not significant"), labels = c("Up", "Down", "NS")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Log2FC", y = expression(-log[10](FDR)), colour = NULL) +
  theme_pub() +
  theme(legend.position = c(0.82, 0.46), legend.background = element_blank(),
        legend.key.width = unit(3, "mm"), plot.margin = ggplot2::margin(3, 3, 3, 4))
save_gg_panel(volcano, "Fig1A_volcano_corrected", 3.42, 2.25)

only_black <- length(setdiff(black, sig$Gene))
only_deg <- length(setdiff(sig$Gene, black))
theta <- seq(0, 2 * pi, length.out = 400)
venn_df <- bind_rows(
  data.frame(x = 0.40 + 0.27 * cos(theta), y = 0.50 + 0.31 * sin(theta), set = "WGCNA"),
  data.frame(x = 0.60 + 0.27 * cos(theta), y = 0.50 + 0.31 * sin(theta), set = "DEGs")
)
venn <- ggplot() +
  geom_polygon(data = venn_df, aes(x, y, group = set, fill = set), alpha = 0.80,
               colour = "black", linewidth = 0.65) +
  annotate("text", x = 0.27, y = 0.50, label = only_black, family = font_family, fontface = "bold", size = 16 / ggplot2::.pt) +
  annotate("text", x = 0.50, y = 0.50, label = length(overlap), family = font_family, fontface = "bold", size = 16 / ggplot2::.pt) +
  annotate("text", x = 0.73, y = 0.50, label = only_deg, family = font_family, fontface = "bold", size = 16 / ggplot2::.pt) +
  annotate("text", x = 0.24, y = 0.85, label = "WGCNA", family = font_family, size = 11 / ggplot2::.pt) +
  annotate("text", x = 0.76, y = 0.85, label = "DEGs", family = font_family, size = 11 / ggplot2::.pt) +
  scale_fill_manual(values = c("WGCNA" = "#D25E5E", "DEGs" = "#4E8AC8")) +
  coord_fixed(xlim = c(0.03, 0.97), ylim = c(0.08, 0.93), clip = "off") +
  theme_void(base_family = font_family) + theme(legend.position = "none", plot.margin = ggplot2::margin(0, 0, 0, 0))
save_gg_panel(venn, "Fig1D_overlap_200", 6, 6)

cat("STAGE\tENRICHMENT\n")
mapped <- suppressMessages(bitr(overlap, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
universe_symbols <- rownames(expr)
universe_mapped <- suppressMessages(bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
gene_ids <- unique(mapped$ENTREZID)
universe_ids <- unique(universe_mapped$ENTREZID)

go_result <- enrichGO(
  gene = gene_ids, universe = universe_ids, OrgDb = org.Hs.eg.db,
  ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
  minGSSize = 10, readable = TRUE
)
kegg_result <- enrichKEGG(
  gene = gene_ids, universe = universe_ids, organism = "hsa", keyType = "kegg",
  pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10
)
kegg_result <- setReadable(kegg_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
go_all <- as.data.frame(go_result)
kegg_all <- as.data.frame(kegg_result)
write.csv(go_all, file.path(panel_dir, "Fig1E_GO_all_results.csv"), row.names = FALSE)
write.csv(kegg_all, file.path(panel_dir, "Fig1F_KEGG_all_results.csv"), row.names = FALSE)

go_top <- go_all %>%
  group_by(ONTOLOGY) %>% arrange(pvalue, .by_group = TRUE) %>% slice_head(n = 5) %>% ungroup() %>%
  mutate(
    ONTOLOGY = factor(ONTOLOGY, levels = c("MF", "CC", "BP")),
    Description_plot = str_wrap(Description, width = 45)
  ) %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description_plot = factor(Description_plot, levels = unique(Description_plot)))

go_plot <- ggplot(go_top, aes(Count, Description_plot, fill = ONTOLOGY)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = Count), hjust = -0.25, size = 8 / ggplot2::.pt, family = font_family) +
  scale_fill_manual(values = c("BP" = red, "CC" = blue, "MF" = green)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.17))) +
  labs(x = "Gene count", y = NULL, fill = NULL) +
  theme_pub() +
  theme(
    axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, colour = "#E7EBEE"),
    legend.position = c(0.91, 0.15), legend.background = element_blank(),
    plot.margin = ggplot2::margin(4, 4, 4, 4)
  )
save_gg_panel(go_plot, "Fig1E_GO_corrected", 4.17, 2.42)

kegg_top <- kegg_all %>% arrange(pvalue) %>% slice_head(n = 5) %>%
  mutate(
    Category = case_when(
      grepl("Complement|Fc epsilon|TNF|IL-17|immune|phagocytosis", Description, ignore.case = TRUE) ~ "Immune",
      grepl("metabolism|biosynthesis|steroidogenesis", Description, ignore.case = TRUE) ~ "Metabolism",
      TRUE ~ "Related"
    ),
    Description_plot = str_wrap(Description, width = 22)
  ) %>%
  arrange(Category, Count) %>%
  mutate(Description_plot = factor(Description_plot, levels = unique(Description_plot)))

kegg_plot <- ggplot(kegg_top, aes(Count, Description_plot, fill = Category)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = Count), hjust = -0.25, size = 8 / ggplot2::.pt, family = font_family) +
  scale_fill_manual(values = c("Immune" = red, "Metabolism" = green, "Related" = blue)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(x = "Gene count", y = NULL, fill = NULL,
       caption = "No pathway passed FDR < 0.05") +
  theme_pub() +
  theme(
    axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(linewidth = 0.25, colour = "#E7EBEE"),
    legend.position = "none", plot.caption = element_text(size = 8, hjust = 0, colour = mid),
    plot.margin = ggplot2::margin(4, 4, 3, 4)
  )
save_gg_panel(kegg_plot, "Fig1F_KEGG_corrected", 2.72, 2.42)

save(go_result, kegg_result, go_all, kegg_all, overlap, deg, sig,
     file = file.path(panel_dir, "Fig1_corrected_analysis.RData"))

cat("GO_FDR_SIG_N\t", sum(go_all$p.adjust < 0.05, na.rm = TRUE), "\n", sep = "")
cat("KEGG_FDR_SIG_N\t", sum(kegg_all$p.adjust < 0.05, na.rm = TRUE), "\n", sep = "")

cat("STAGE\tS1_MODELS\n")
genes <- intersect(overlap, rownames(expr))
x <- as.data.frame(t(expr[genes, , drop = FALSE]), check.names = FALSE)
nzv <- nearZeroVar(x)
if (length(nzv)) x <- x[, -nzv, drop = FALSE]
for (j in seq_len(ncol(x))) if (anyNA(x[[j]])) x[[j]][is.na(x[[j]])] <- median(x[[j]], na.rm = TRUE)
x_scaled <- scale(x)
y_num <- as.numeric(group == "Disease")

set.seed(123)
lasso_fit <- glmnet(x_scaled, y_num, family = "binomial", alpha = 1, maxit = 10000)
set.seed(123)
lasso_cv <- cv.glmnet(x_scaled, y_num, family = "binomial", alpha = 1, nfolds = 10, maxit = 10000)
lasso_coef <- as.matrix(coef(lasso_fit, s = lasso_cv$lambda.min))
lasso_imp <- data.frame(Gene = rownames(lasso_coef), Value = abs(lasso_coef[, 1])) %>%
  filter(Gene != "(Intercept)", Value > 0) %>% arrange(desc(Value))
lasso_genes <- lasso_imp$Gene

control <- trainControl(
  method = "repeatedcv", number = 10, repeats = 3,
  summaryFunction = twoClassSummary, classProbs = TRUE,
  savePredictions = "final", verboseIter = FALSE, allowParallel = TRUE
)

workers <- min(6L, max(1L, parallel::detectCores(logical = TRUE) - 2L))
cl <- parallel::makePSOCKcluster(workers)
doParallel::registerDoParallel(cl)

set.seed(123)
rf_train <- train(
  x = x_scaled, y = group, method = "rf", metric = "ROC",
  trControl = control, importance = TRUE, tuneLength = 5
)
rf_raw <- randomForest::importance(rf_train$finalModel)
rf_metric <- if ("MeanDecreaseGini" %in% colnames(rf_raw)) "MeanDecreaseGini" else colnames(rf_raw)[1]
rf_imp <- data.frame(Gene = rownames(rf_raw), Value = rf_raw[, rf_metric]) %>% arrange(desc(Value))
rf_genes <- head(rf_imp$Gene, 15)

set.seed(123)
boost_train <- train(
  x = x_scaled, y = group, method = "glmboost", metric = "ROC",
  trControl = control, tuneLength = 5
)
boost_coef_raw <- coef(boost_train$finalModel)
boost_imp <- data.frame(Gene = names(boost_coef_raw), Value = abs(as.numeric(boost_coef_raw))) %>%
  filter(Gene != "(Intercept)", Value > 0) %>% arrange(desc(Value))
if (!nrow(boost_imp)) {
  bvi <- varImp(boost_train, scale = FALSE)$importance
  boost_imp <- data.frame(Gene = rownames(bvi), Value = bvi$Overall) %>% arrange(desc(Value))
}
boost_genes <- boost_imp$Gene

log_file <- file.path(panel_dir, "S1_xgboost_training.log")
zz <- file(log_file, open = "wt")
sink(zz, type = "output")
sink(zz, type = "message")
old_warn <- getOption("warn")
options(warn = -1)
set.seed(123)
xgb_train <- train(
  x = as.matrix(x_scaled), y = group, method = "xgbTree", metric = "ROC",
  trControl = control, tuneLength = 5
)
options(warn = old_warn)
sink(type = "message")
sink(type = "output")
close(zz)

xgb_raw <- xgb.importance(feature_names = colnames(x_scaled), model = xgb_train$finalModel)
xgb_imp <- data.frame(Gene = xgb_raw$Feature, Value = xgb_raw$Gain) %>% arrange(desc(Value))
xgb_genes <- head(xgb_imp$Gene, 50)

parallel::stopCluster(cl)
foreach::registerDoSEQ()
cl <- NULL

contrib <- predict(xgb_train$finalModel, as.matrix(x_scaled), predcontrib = TRUE)
feature_cols <- setdiff(colnames(contrib), "BIAS")
shap_all <- as.data.frame(contrib[, feature_cols, drop = FALSE], check.names = FALSE)
baseline <- mean(contrib[, "BIAS"])
shap_importance <- sort(colMeans(abs(as.matrix(shap_all))), decreasing = TRUE)
top6 <- head(names(shap_importance), 6)
shap_top <- shap_all[, top6, drop = FALSE]
x_top <- as.data.frame(x_scaled[, top6, drop = FALSE], check.names = FALSE)

rank_df <- bind_rows(
  head(lasso_imp, 5) %>% mutate(Model = "LASSO"),
  head(boost_imp, 5) %>% mutate(Model = "glmBoost"),
  head(xgb_imp, 5) %>% mutate(Model = "XGBoost"),
  head(rf_imp, 5) %>% mutate(Model = "Ranger")
) %>% select(Gene, Model, Value)

write.csv(data.frame(Sample = rownames(shap_top), shap_top), file.path(panel_dir, "S1_shap_values_corrected.csv"), row.names = FALSE)
write.csv(data.frame(Sample = rownames(x_top), x_top), file.path(panel_dir, "S1_feature_values_corrected.csv"), row.names = FALSE)
write.csv(rank_df, file.path(panel_dir, "S1_model_rankings_corrected.csv"), row.names = FALSE)
write.csv(data.frame(Baseline_log_odds = baseline), file.path(panel_dir, "S1_shap_baseline_corrected.csv"), row.names = FALSE)
write.csv(data.frame(
  lambda = lasso_cv$lambda, ll = log(lasso_cv$lambda), cvm = lasso_cv$cvm,
  cvsd = lasso_cv$cvsd, cvup = lasso_cv$cvup, cvlo = lasso_cv$cvlo,
  nzero = lasso_cv$nzero
), file.path(panel_dir, "S1_lasso_cv_corrected.csv"), row.names = FALSE)

path_df <- as.data.frame(as.matrix(lasso_fit$beta)) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Step", values_to = "Coefficient") %>%
  mutate(StepIndex = as.integer(sub("s", "", Step)), Lambda = lasso_fit$lambda[StepIndex + 1L], LogLambda = log(Lambda))
write.csv(path_df, file.path(panel_dir, "S1_lasso_path_corrected.csv"), row.names = FALSE)

cat("S1_INPUT_N\t", ncol(x_scaled), "\n", sep = "")
cat("S1_TOP6_SHAP\t", paste(top6, collapse = ","), "\n", sep = "")
cat("S1_LASSO\t", paste(lasso_genes, collapse = ","), "\n", sep = "")
cat("S1_ALL4\t", paste(Reduce(intersect, list(lasso_genes, rf_genes, boost_genes, xgb_genes)), collapse = ","), "\n", sep = "")

cat("STAGE\tS1_PLOTS\n")
cv_df <- data.frame(
  LogLambda = log(lasso_cv$lambda), Mean = lasso_cv$cvm,
  Lower = lasso_cv$cvlo, Upper = lasso_cv$cvup
)
pA <- ggplot(cv_df, aes(LogLambda, Mean)) +
  geom_linerange(aes(ymin = Lower, ymax = Upper), colour = light, linewidth = 0.35) +
  geom_point(size = 1.35, shape = 21, fill = navy, colour = "white", stroke = 0.15) +
  geom_vline(xintercept = log(lasso_cv$lambda.min), colour = coral, linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = log(lasso_cv$lambda.1se), colour = teal, linetype = "dashed", linewidth = 0.5) +
  annotate("text", x = log(lasso_cv$lambda.min), y = min(cv_df$Lower), label = expression(lambda[min]),
           hjust = 1.05, vjust = -0.2, colour = coral, size = 8 / ggplot2::.pt, family = font_family) +
  annotate("text", x = log(lasso_cv$lambda.1se), y = min(cv_df$Lower), label = expression(lambda[1*SE]),
           hjust = -0.05, vjust = -0.2, colour = teal, size = 8 / ggplot2::.pt, family = font_family) +
  labs(x = expression(log(lambda)), y = "CV deviance", tag = "A") + theme_pub()

path_colors <- c("S100A8" = navy, "MAFG" = coral)
pB <- ggplot(path_df, aes(LogLambda, Coefficient, group = Gene)) +
  geom_line(data = filter(path_df, !Gene %in% c("S100A8", "MAFG")), colour = light, linewidth = 0.32, alpha = 0.85) +
  geom_line(data = filter(path_df, Gene %in% c("S100A8", "MAFG")), aes(colour = Gene), linewidth = 0.75) +
  geom_vline(xintercept = log(lasso_cv$lambda.min), colour = coral, linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = log(lasso_cv$lambda.1se), colour = teal, linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = path_colors) +
  labs(x = expression(log(lambda)), y = "Regularized coefficient", colour = NULL, tag = "B") +
  theme_pub() + theme(legend.position = c(0.82, 0.82), legend.background = element_blank())

rank_plot_df <- rank_df %>%
  group_by(Model) %>% arrange(desc(Value), .by_group = TRUE) %>% mutate(Rank = row_number()) %>% ungroup()
model_levels <- c("LASSO", "glmBoost", "XGBoost", "Ranger")
gene_levels <- rank_plot_df %>% group_by(Gene) %>% summarise(MinRank = min(Rank), MeanRank = mean(Rank), .groups = "drop") %>%
  arrange(MinRank, MeanRank) %>% pull(Gene)
tile_df <- tidyr::expand_grid(Gene = gene_levels, Model = model_levels) %>%
  left_join(rank_plot_df, by = c("Gene", "Model")) %>%
  mutate(
    Gene = factor(Gene, levels = rev(gene_levels)),
    Model = factor(Model, levels = model_levels),
    FillScore = ifelse(is.na(Rank), 0, 6 - Rank)
  )
pC <- ggplot(tile_df, aes(Model, Gene)) +
  geom_tile(aes(fill = FillScore), colour = "white", linewidth = 0.6, width = 0.86, height = 0.86) +
  geom_text(aes(label = ifelse(is.na(Rank), "", Rank), colour = ifelse(!is.na(Rank) & Rank <= 2, "white", ink)),
            size = 8 / ggplot2::.pt, fontface = "bold", family = font_family, show.legend = FALSE) +
  scale_fill_gradient(low = pale, high = teal, limits = c(0, 5), guide = "none") +
  scale_colour_identity() +
  scale_x_discrete(labels = c("Ranger" = "Random forest")) +
  labs(x = "Within-model rank (1 = highest)", y = NULL, tag = "C") +
  theme_minimal(base_size = 8, base_family = font_family) +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8, face = "italic", colour = ink),
        panel.grid = element_blank(), plot.tag = element_text(size = 12, face = "bold"),
        plot.tag.position = c(0, 1), plot.margin = ggplot2::margin(4, 5, 4, 5))

shap_long <- shap_top %>% mutate(Row = row_number()) %>% pivot_longer(-Row, names_to = "Gene", values_to = "SHAP") %>%
  left_join(x_top %>% mutate(Row = row_number()) %>% pivot_longer(-Row, names_to = "Gene", values_to = "FeatureValue"),
            by = c("Row", "Gene")) %>%
  mutate(Gene = factor(Gene, levels = rev(top6)))
set.seed(123)
pD <- ggplot(shap_long, aes(SHAP, Gene, colour = FeatureValue)) +
  geom_vline(xintercept = 0, colour = mid, linewidth = 0.35) +
  geom_jitter(height = 0.15, width = 0, size = 1.0, alpha = 0.82) +
  scale_colour_gradient2(low = "#3569A8", mid = "#F3F0E8", high = "#CB4C64", midpoint = 0,
                         limits = c(-2.5, 2.5), oob = scales::squish) +
  labs(x = "SHAP value (log-odds)", y = "Feature", colour = "Feature value", tag = "D") +
  theme_pub() + theme(axis.text.y = element_text(face = "italic"), legend.position = "right")

row_id <- min(2L, nrow(shap_top))
wf <- data.frame(Gene = top6, Value = as.numeric(shap_top[row_id, top6]), Feature = as.numeric(x_top[row_id, top6])) %>%
  arrange(abs(Value)) %>% mutate(
    Start = baseline + lag(cumsum(Value), default = 0), End = Start + Value,
    Y = row_number(), Label = sprintf("%+.3f", Value),
    GeneLabel = sprintf("%s  (%.2f)", Gene, Feature)
  )
final_margin <- baseline + sum(wf$Value)
pE <- ggplot(wf) +
  geom_rect(aes(xmin = pmin(Start, End), xmax = pmax(Start, End), ymin = Y - 0.31, ymax = Y + 0.31, fill = Value >= 0),
            colour = "white", linewidth = 0.3) +
  geom_text(aes(x = (Start + End) / 2, y = Y, label = Label), colour = "white", size = 8 / ggplot2::.pt, family = font_family) +
  geom_vline(xintercept = baseline, colour = mid, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = final_margin, colour = ink, linetype = "dashed", linewidth = 0.4) +
  scale_fill_manual(values = c("TRUE" = coral, "FALSE" = navy), guide = "none") +
  scale_y_continuous(breaks = wf$Y, labels = wf$GeneLabel, expand = expansion(mult = c(0.03, 0.12))) +
  labs(x = "SHAP contribution (log-odds)", y = NULL, tag = "E") +
  annotate("text", x = baseline, y = max(wf$Y) + 0.65, label = sprintf("Baseline = %.3f", baseline),
           hjust = 0, size = 8 / ggplot2::.pt, colour = mid, family = font_family) +
  annotate("text", x = final_margin, y = max(wf$Y) + 0.65, label = sprintf("Prediction = %.3f", final_margin),
           hjust = 1, size = 8 / ggplot2::.pt, colour = ink, family = font_family) +
  theme_pub() + theme(axis.text.y = element_text(face = "italic"))

dep_plots <- lapply(seq_along(top6), function(i) {
  gene <- top6[i]
  partner <- if (gene == "S100A8") setdiff(top6, gene)[1] else "S100A8"
  if (!partner %in% top6) partner <- setdiff(top6, gene)[1]
  dd <- data.frame(X = x_top[[gene]], SHAP = shap_top[[gene]], Color = x_top[[partner]])
  ggplot(dd, aes(X, SHAP, colour = Color)) +
    geom_hline(yintercept = 0, colour = light, linewidth = 0.3) +
    geom_point(size = 1.0, alpha = 0.82) +
    scale_colour_gradient2(low = "#3569A8", mid = "#F3F0E8", high = "#CB4C64", midpoint = 0,
                           limits = c(-2.5, 2.5), oob = scales::squish, guide = "none") +
    labs(x = gene, y = "SHAP value", tag = if (i == 1) "F" else NULL) +
    annotate("text", x = Inf, y = -Inf, label = paste0("Color: ", partner), hjust = 1.05, vjust = -0.5,
             size = 7.5 / ggplot2::.pt, colour = mid, family = font_family) +
    theme_pub() + theme(axis.title.x = element_text(face = "italic"), plot.margin = ggplot2::margin(3, 4, 4, 4))
})
pF <- wrap_plots(dep_plots, ncol = 3, nrow = 2)

save_gg_panel(pA, "S1A_lasso_cv_corrected", 3.3, 1.85)
save_gg_panel(pB, "S1B_lasso_path_corrected", 3.3, 1.85)
save_gg_panel(pC, "S1C_model_ranks_corrected", 3.3, 1.9)
save_gg_panel(pD, "S1D_shap_beeswarm_corrected", 3.3, 1.9)
save_gg_panel(pE, "S1E_shap_waterfall_corrected", 6.7, 1.35)
save_gg_panel(pF, "S1F_shap_dependence_corrected", 6.7, 2.3)

s1 <- ((pA | pB) / (pC | pD) / pE / pF) +
  plot_layout(heights = c(1.90, 1.95, 1.40, 2.40)) &
  theme(plot.margin = ggplot2::margin(4, 5, 4, 5))

ggsave(file.path(formal_dir, "S1_fig_R_master.pdf"), s1, width = 7.5, height = 8.75,
       device = cairo_pdf, family = font_family, bg = "white")
ggsave(file.path(panel_dir, "S1_fig_corrected_preview.tif"), s1, width = 7.5, height = 8.75,
       device = ragg::agg_tiff, dpi = 300, compression = "lzw", bg = "white")

save(
  lasso_fit, lasso_cv, lasso_imp, rf_train, rf_imp, boost_train, boost_imp,
  xgb_train, xgb_imp, shap_all, shap_top, x_top, rank_df, top6, baseline,
  overlap, file = file.path(panel_dir, "S1_corrected_analysis.RData")
)

cat("DONE\n")
