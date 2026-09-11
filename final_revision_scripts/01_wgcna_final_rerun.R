options(stringsAsFactors = FALSE, warn = 1)
.libPaths(c("C:/Users/ZCL/AppData/Local/R/win-library/4.5", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(WGCNA)
  library(dynamicTreeCut)
  library(matrixStats)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

# Bioconductor packages loaded after WGCNA can mask WGCNA::cor with an S4
# generic that does not accept the weight arguments used internally by
# blockwiseModules. Bind the WGCNA implementation explicitly; this is a
# compatibility fix and does not alter the requested Pearson correlation.
cor <- WGCNA::cor

set.seed(20260908)
allowWGCNAThreads(nThreads = max(1L, min(6L, parallel::detectCores(logical = TRUE) - 2L)))

root <- Sys.getenv("SEPSIS_PROJECT_ROOT", unset = normalizePath("..", winslash = "/", mustWork = TRUE))
out_dir <- file.path(root, "work", "final_revision_20260908", "wgcna_rerun")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(out_dir, "WGCNA_final_rerun.log"), open = "wt", encoding = "UTF-8")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message", append = TRUE)
on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(sink(type = "output"), silent = TRUE)
  try(close(log_con), silent = TRUE)
}, add = TRUE)

cat("WGCNA FINAL RERUN\n")
cat("Started:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
cat("R:", R.version.string, "\n")
cat("WGCNA:", as.character(packageVersion("WGCNA")), "\n")
cat("Seed: 20260908\n")

input_matrix <- file.path(root, "sepsis-s100a8-multiomics", "data", "02_DEG", "Input", "GSE65682_gene.csv")
input_meta <- file.path(root, "sepsis-s100a8-multiomics", "data", "02_DEG", "Input", "GSE65682_Groups.csv")
deg_file <- file.path(root, "sepsis-s100a8-multiomics", "data", "02_DEG", "output", "Disease_vs_Normal_all_results_padj.csv")
stopifnot(file.exists(input_matrix), file.exists(input_meta), file.exists(deg_file))

cat("Reading maximum-probe gene matrix without modifying source file...\n")
raw_dt <- fread(input_matrix, data.table = FALSE, check.names = FALSE)
genes <- raw_dt[[1L]]
linear_mat <- as.matrix(raw_dt[, -1L, drop = FALSE])
storage.mode(linear_mat) <- "double"
rownames(linear_mat) <- genes
rm(raw_dt, genes)
if (any(!is.finite(linear_mat)) || any(linear_mat <= 0)) {
  stop("The archived maximum-probe matrix contains non-finite or non-positive values; log2 conversion is not valid.")
}
expr_log2 <- log2(linear_mat)
rm(linear_mat)
gc()

meta <- fread(input_meta, data.table = FALSE, check.names = FALSE)
names(meta)[1:2] <- c("sample_id", "group_original")
meta$group <- ifelse(meta$group_original %in% c("Disease", "Sepsis"), "Sepsis",
                     ifelse(meta$group_original %in% c("Normal", "Healthy"), "Normal", NA_character_))
idx <- match(colnames(expr_log2), meta$sample_id)
if (anyNA(idx)) stop("Expression samples are missing from phenotype metadata.")
meta <- meta[idx, , drop = FALSE]
if (!identical(colnames(expr_log2), meta$sample_id)) stop("Expression and phenotype sample order is not identical.")
if (anyNA(meta$group)) stop("Unsupported or missing phenotype labels.")
if (ncol(expr_log2) != 802L || sum(meta$group == "Sepsis") != 760L || sum(meta$group == "Normal") != 42L) {
  stop("Discovery-cohort count check failed; expected 802 total, 760 sepsis, and 42 normal controls.")
}

input_summary <- data.frame(
  item = c("input_scale", "gene_aggregation", "genes_before_filter", "samples_before_qc",
           "sepsis_samples", "normal_samples", "nonfinite_values", "duplicate_genes", "duplicate_samples"),
  value = c("log2 RMA", "maximum probe per gene", nrow(expr_log2), ncol(expr_log2),
            sum(meta$group == "Sepsis"), sum(meta$group == "Normal"),
            sum(!is.finite(expr_log2)), sum(duplicated(rownames(expr_log2))),
            sum(duplicated(colnames(expr_log2))))
)
write.csv(input_summary, file.path(out_dir, "input_matrix_summary.csv"), row.names = FALSE)

# Prespecified, result-independent low-information filter: retain genes with MAD
# strictly above the median MAD across all genes. No disease labels are used.
gene_mad <- rowMads(expr_log2, na.rm = TRUE)
gene_var <- rowVars(expr_log2, na.rm = TRUE)
mad_cut <- median(gene_mad[is.finite(gene_mad)], na.rm = TRUE)
keep_low_info <- is.finite(gene_mad) & is.finite(gene_var) & gene_var > 0 & gene_mad > mad_cut
filter_audit <- data.frame(
  gene = rownames(expr_log2), MAD = gene_mad, variance = gene_var,
  retained = keep_low_info, reason = ifelse(!is.finite(gene_mad) | !is.finite(gene_var), "nonfinite_dispersion",
                                    ifelse(gene_var <= 0, "zero_variance",
                                    ifelse(gene_mad <= mad_cut, "MAD_at_or_below_median", "retained")))
)
write.csv(filter_audit, file.path(out_dir, "gene_filter_audit.csv"), row.names = FALSE)
expr_filtered <- expr_log2[keep_low_info, , drop = FALSE]
cat("MAD threshold:", format(mad_cut, digits = 10), "\n")
cat("Genes retained:", nrow(expr_filtered), "of", nrow(expr_log2), "\n")
cat("S100A8 retained:", "S100A8" %in% rownames(expr_filtered), "\n")

datExpr0 <- t(expr_filtered)
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  bad_samples <- rownames(datExpr0)[!gsg$goodSamples]
  bad_genes <- colnames(datExpr0)[!gsg$goodGenes]
  writeLines(bad_samples, file.path(out_dir, "goodSamplesGenes_failed_samples.txt"))
  writeLines(bad_genes, file.path(out_dir, "goodSamplesGenes_failed_genes.txt"))
  if (length(bad_samples)) stop("goodSamplesGenes identified failed samples; automatic sample deletion is prohibited in this rerun.")
  datExpr0 <- datExpr0[, gsg$goodGenes, drop = FALSE]
}
datExpr <- datExpr0
rm(datExpr0, expr_filtered)
gc()

# Sample QC is diagnostic only. No samples are removed based on group separation,
# hierarchical clustering, PCA, or connectivity.
sample_tree <- hclust(dist(datExpr), method = "average")
pdf(file.path(out_dir, "sample_clustering_all_802.pdf"), width = 24, height = 12, family = "sans")
plot(sample_tree, labels = FALSE, main = "All 802 discovery samples retained", xlab = "", sub = "")
dev.off()

pca_genes <- order(apply(datExpr, 2, var), decreasing = TRUE)[seq_len(min(2000L, ncol(datExpr)))]
pca <- prcomp(datExpr[, pca_genes, drop = FALSE], center = TRUE, scale. = TRUE)
pca_df <- data.frame(sample_id = rownames(datExpr), group = meta$group,
                     PC1 = pca$x[, 1], PC2 = pca$x[, 2], stringsAsFactors = FALSE)
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)
sample_cor <- cor(t(datExpr), use = "pairwise.complete.obs", method = "pearson")
sample_connectivity <- rowSums((1 + sample_cor) / 2, na.rm = TRUE) - 1
connectivity_z <- as.numeric(scale(sample_connectivity))
pca_df$sample_connectivity <- sample_connectivity[pca_df$sample_id]
pca_df$connectivity_z <- connectivity_z[match(pca_df$sample_id, names(sample_connectivity))]
pca_df$excluded <- FALSE
write.csv(pca_df, file.path(out_dir, "sample_qc_all_802.csv"), row.names = FALSE)
rm(sample_cor)
gc()

p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = group)) +
  geom_point(size = 1.0, alpha = 0.70) +
  scale_colour_manual(values = c(Normal = "#4C78A8", Sepsis = "#D66A5E")) +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), colour = NULL) +
  theme_classic(base_size = 8, base_family = "sans")
ggsave(file.path(out_dir, "PCA_all_802.pdf"), p_pca, width = 5.5, height = 4.2, device = cairo_pdf, family = "sans")

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", corFnc = "cor",
                         corOptions = list(use = "p"), verbose = 5)
softPower <- if (!is.na(sft$powerEstimate)) {
  as.numeric(sft$powerEstimate)
} else {
  signed_r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  candidate <- which(signed_r2 > 0.90)
  if (length(candidate)) sft$fitIndices[candidate[1], 1] else 6
}
sft_out <- as.data.frame(sft$fitIndices)
write.csv(sft_out, file.path(out_dir, "soft_threshold_diagnostics.csv"), row.names = FALSE)
cat("Selected soft threshold:", softPower, "\n")

pdf(file.path(out_dir, "soft_threshold_selection.pdf"), width = 8, height = 4, family = "sans")
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
signed_r2 <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
plot(sft$fitIndices[, 1], signed_r2, type = "n", xlab = "Soft-threshold power", ylab = "Signed scale-free fit, R²")
text(sft$fitIndices[, 1], signed_r2, labels = sft$fitIndices[, 1], cex = 0.7, col = "#D66A5E")
abline(h = 0.90, lty = 2, col = "#7A858D")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft-threshold power", ylab = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = sft$fitIndices[, 1], cex = 0.7, col = "#4C78A8")
dev.off()

cat("Building one-block signed network and signed TOM...\n")
net <- blockwiseModules(
  datExpr,
  power = softPower,
  maxBlockSize = 6000,
  networkType = "signed",
  TOMType = "signed",
  corType = "pearson",
  deepSplit = 2,
  minModuleSize = 200,
  mergeCutHeight = 0.25,
  pamRespectsDendro = FALSE,
  numericLabels = FALSE,
  randomSeed = 20260908,
  saveTOMs = FALSE,
  verbose = 3
)
if (length(unique(net$blocks)) != 1L) stop("The filtered network unexpectedly required more than one block.")
moduleColors <- net$colors
names(moduleColors) <- colnames(datExpr)
MEs <- orderMEs(net$MEs)
trait <- as.numeric(meta$group == "Sepsis")
names(trait) <- meta$sample_id
if (!identical(rownames(MEs), meta$sample_id)) stop("Module eigengene sample order mismatch.")

module_cor <- as.numeric(cor(MEs, trait, use = "p", method = "pearson"))
module_p <- as.numeric(corPvalueStudent(module_cor, nSamples = nrow(datExpr)))
module_trait <- data.frame(
  module = sub("^ME", "", colnames(MEs)),
  eigengene = colnames(MEs),
  module_size = as.integer(table(moduleColors)[sub("^ME", "", colnames(MEs))]),
  correlation_sepsis = module_cor,
  p_value = module_p,
  FDR = p.adjust(module_p, method = "BH"),
  direction = ifelse(module_cor > 0, "Sepsis-enriched", "Normal-enriched"),
  stringsAsFactors = FALSE
)
module_trait <- module_trait[order(-abs(module_trait$correlation_sepsis)), ]
write.csv(module_trait, file.path(out_dir, "module_trait_all_results.csv"), row.names = FALSE)
write.csv(data.frame(gene = names(moduleColors), module = unname(moduleColors)),
          file.path(out_dir, "module_assignments_all_genes.csv"), row.names = FALSE)
write.csv(data.frame(module = names(table(moduleColors)), module_size = as.integer(table(moduleColors))),
          file.path(out_dir, "module_sizes_all.csv"), row.names = FALSE)

positive_sig <- subset(module_trait, correlation_sepsis > 0 & FDR < 0.05)
if (nrow(positive_sig)) {
  primary_module <- positive_sig$module[which.max(positive_sig$correlation_sepsis)]
} else {
  primary_module <- module_trait$module[which.max(module_trait$correlation_sepsis)]
}

all_kme <- signedKME(datExpr, MEs, outputColumnName = "kME", corFnc = "cor", corOptions = "use='p'")
gs <- as.numeric(cor(datExpr, trait, use = "p", method = "pearson"))
gs_p <- as.numeric(corPvalueStudent(gs, nSamples = nrow(datExpr)))
own_me_col <- paste0("kME", moduleColors)
own_mm <- vapply(seq_along(moduleColors), function(i) all_kme[i, own_me_col[i]], numeric(1))
gene_stats <- data.frame(
  gene = colnames(datExpr), module = moduleColors, module_membership = own_mm,
  gene_significance = gs, gene_significance_p = gs_p,
  gene_significance_FDR = p.adjust(gs_p, method = "BH"), stringsAsFactors = FALSE
)
write.csv(gene_stats, file.path(out_dir, "gene_module_membership_and_significance.csv"), row.names = FALSE)

s100 <- gene_stats[gene_stats$gene == "S100A8", , drop = FALSE]
if (!nrow(s100)) {
  s100 <- data.frame(gene = "S100A8", module = NA_character_, module_membership = NA_real_,
                     gene_significance = NA_real_, gene_significance_p = NA_real_, gene_significance_FDR = NA_real_)
}
write.csv(s100, file.path(out_dir, "S100A8_WGCNA_statistics.csv"), row.names = FALSE)

deg <- fread(deg_file, data.table = FALSE, check.names = FALSE)
if (!all(c("Gene", "padj", "log2FoldChange") %in% names(deg))) stop("Corrected DEG file lacks required columns.")
deg$Regulation_final <- ifelse(!is.na(deg$padj) & deg$padj < 0.05 & deg$log2FoldChange > 1, "Up",
                               ifelse(!is.na(deg$padj) & deg$padj < 0.05 & deg$log2FoldChange < -1, "Down", "Not significant"))
sig_deg <- deg$Gene[deg$Regulation_final != "Not significant"]
primary_genes <- names(moduleColors)[moduleColors == primary_module]
overlap <- intersect(primary_genes, sig_deg)
writeLines(primary_genes, file.path(out_dir, paste0("primary_sepsis_module_", primary_module, "_genes.txt")))
writeLines(overlap, file.path(out_dir, "primary_sepsis_module_corrected_DEG_overlap.txt"))
write.csv(deg, file.path(out_dir, "corrected_DEG_results_used.csv"), row.names = FALSE)

universe_symbols <- intersect(colnames(datExpr), deg$Gene)
mapped_overlap <- suppressMessages(bitr(overlap, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
mapped_universe <- suppressMessages(bitr(universe_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db))
overlap_ids <- unique(mapped_overlap$ENTREZID)
universe_ids <- unique(mapped_universe$ENTREZID)

go_all <- data.frame()
kegg_all <- data.frame()
if (length(overlap_ids) >= 10L) {
  go_obj <- enrichGO(gene = overlap_ids, universe = universe_ids, OrgDb = org.Hs.eg.db,
                     ont = "ALL", pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                     minGSSize = 10, readable = TRUE)
  kegg_obj <- enrichKEGG(gene = overlap_ids, universe = universe_ids, organism = "hsa",
                         keyType = "kegg", pvalueCutoff = 1, pAdjustMethod = "BH", minGSSize = 10)
  if (!is.null(kegg_obj) && nrow(as.data.frame(kegg_obj))) {
    kegg_obj <- setReadable(kegg_obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  }
  go_all <- as.data.frame(go_obj)
  kegg_all <- as.data.frame(kegg_obj)
}
write.csv(go_all, file.path(out_dir, "overlap_GO_all_results.csv"), row.names = FALSE)
write.csv(kegg_all, file.path(out_dir, "overlap_KEGG_all_results.csv"), row.names = FALSE)

parameters <- data.frame(
  parameter = c("input_scale", "aggregation", "gene_filter", "missing_filter", "sample_exclusion",
                "sample_n", "gene_n", "soft_threshold_rule", "soft_threshold_power", "networkType",
                "TOMType", "correlation", "minModuleSize", "deepSplit", "pamRespectsDendro",
                "mergeCutHeight", "module_eigengene", "module_trait_test", "module_test_correction",
                "random_seed", "maxBlockSize", "blocks"),
  value = c("log2 RMA", "maximum probe per gene", "MAD above across-gene median; labels not used",
            "finite values, nonzero variance, goodSamplesGenes", "none; clustering/PCA/connectivity diagnostic only",
            nrow(datExpr), ncol(datExpr), "WGCNA estimate; else first signed R2 > 0.90; else 6", softPower,
            "signed", "signed", "Pearson", 200, 2, FALSE, 0.25, "WGCNA moduleEigengenes via blockwiseModules",
            "Pearson correlation with Sepsis=1, Normal=0", "BH across all module-trait tests",
            20260908, 6000, length(unique(net$blocks)))
)
write.csv(parameters, file.path(out_dir, "WGCNA_parameters.csv"), row.names = FALSE)

summary_out <- data.frame(
  metric = c("samples", "sepsis", "normal", "genes_before_filter", "genes_in_network", "modules_non_grey",
             "primary_sepsis_module", "primary_module_size", "primary_module_correlation", "primary_module_p",
             "primary_module_FDR", "S100A8_module", "S100A8_MM", "S100A8_GS", "corrected_DEG_n",
             "primary_module_DEG_overlap_n", "GO_FDR_significant_n", "KEGG_FDR_significant_n"),
  value = c(nrow(datExpr), sum(trait == 1), sum(trait == 0), nrow(expr_log2), ncol(datExpr),
            length(setdiff(unique(moduleColors), "grey")), primary_module,
            module_trait$module_size[module_trait$module == primary_module],
            module_trait$correlation_sepsis[module_trait$module == primary_module],
            module_trait$p_value[module_trait$module == primary_module],
            module_trait$FDR[module_trait$module == primary_module],
            s100$module[1], s100$module_membership[1], s100$gene_significance[1], length(sig_deg), length(overlap),
            if (nrow(go_all)) sum(go_all$p.adjust < 0.05, na.rm = TRUE) else 0,
            if (nrow(kegg_all)) sum(kegg_all$p.adjust < 0.05, na.rm = TRUE) else 0)
)
write.csv(summary_out, file.path(out_dir, "WGCNA_final_summary.csv"), row.names = FALSE)

saveRDS(expr_log2, file.path(out_dir, "GSE65682_log2_RMA_genelevel_maxprobe_final.rds"), compress = "xz")
saveRDS(list(moduleColors = moduleColors, MEs = MEs, moduleTrait = module_trait,
             geneStats = gene_stats, parameters = parameters, softThreshold = sft_out,
             primaryModule = primary_module, overlap = overlap, GO = go_all, KEGG = kegg_all,
             sampleQC = pca_df), file.path(out_dir, "WGCNA_final_results.rds"), compress = "xz")

cat("Primary sepsis-associated module:", primary_module, "\n")
cat("Primary module size:", length(primary_genes), "\n")
cat("Primary module-DEG overlap:", length(overlap), "\n")
cat("S100A8 module/MM/GS:", s100$module[1], s100$module_membership[1], s100$gene_significance[1], "\n")
cat("GO FDR significant:", if (nrow(go_all)) sum(go_all$p.adjust < 0.05, na.rm = TRUE) else 0, "\n")
cat("KEGG FDR significant:", if (nrow(kegg_all)) sum(kegg_all$p.adjust < 0.05, na.rm = TRUE) else 0, "\n")
cat("Completed:", format(Sys.time(), tz = "Asia/Shanghai"), "\n")
