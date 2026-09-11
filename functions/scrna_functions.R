# GSE167363 供者级单细胞分析公共函数。
# 原则：供者是独立统计单位；原始 counts 用于 pseudo-bulk 差异分析。

parse_gse167363_donor <- function(sample_id) {
  x <- sub("^GSM[0-9]+_", "", as.character(sample_id))
  sub("_(T0|T6)$", "", x)
}

parse_gse167363_time <- function(sample_id) {
  x <- as.character(sample_id)
  ifelse(grepl("_T6$", x), "T6", "Baseline")
}

get_counts_layer <- function(object, assay = "RNA") {
  if (!assay %in% Seurat::Assays(object)) stop("Seurat 对象缺少 assay: ", assay, call. = FALSE)
  SeuratObject::LayerData(object[[assay]], layer = "counts")
}

aggregate_counts_by_donor <- function(counts, donor) {
  donor <- factor(donor)
  mm <- Matrix::sparse.model.matrix(~ 0 + donor)
  colnames(mm) <- sub("^donor", "", colnames(mm))
  counts %*% mm
}

run_edger_pseudobulk <- function(counts, sample_meta, min_cpm_samples = 2L) {
  stopifnot(identical(colnames(counts), sample_meta$donor_id))
  sample_meta$group <- stats::relevel(factor(sample_meta$group), ref = "Control")
  if (length(unique(sample_meta$group)) < 2L || min(table(sample_meta$group)) < 1L) {
    stop("pseudo-bulk 设计缺少可比较的 Control/Sepsis 供者。", call. = FALSE)
  }
  y <- edgeR::DGEList(counts = counts, samples = sample_meta)
  keep <- rowSums(edgeR::cpm(y) > 1) >= min(min_cpm_samples, ncol(y))
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  design <- stats::model.matrix(~ group, data = y$samples)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = "groupSepsis")
  tab <- edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table
  tab$gene <- rownames(tab)
  tab <- tab[, c("gene", setdiff(names(tab), "gene"))]
  list(table = tab, dge = y, fit = fit)
}

donor_summary_ci <- function(values, group, conf = 0.95) {
  z <- stats::qnorm(1 - (1 - conf) / 2)
  spl <- split(values, group)
  do.call(rbind, lapply(names(spl), function(g) {
    x <- spl[[g]]; n <- sum(is.finite(x)); m <- mean(x, na.rm = TRUE)
    se <- stats::sd(x, na.rm = TRUE) / sqrt(n)
    data.frame(group = g, n = n, mean = m, ci_low = m - z * se, ci_high = m + z * se)
  }))
}
