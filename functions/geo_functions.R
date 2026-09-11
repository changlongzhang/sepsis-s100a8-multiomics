# GEO本地文件读取、探针注释和基因聚合函数
read_local_eset <- function(series_matrix) {
  assert_file(series_matrix, "GEO series matrix")
  x <- GEOquery::getGEO(filename = series_matrix, getGPL = FALSE)
  if (inherits(x, "list")) x <- x[[1]]
  if (!inherits(x, "ExpressionSet")) stop("GEO文件未解析为ExpressionSet: ", series_matrix, call. = FALSE)
  x
}

read_affy_annotation <- function(path) {
  assert_file(path, "GPL annotation")
  header_line <- grep("^ID\\t", readLines(path, warn = FALSE, n = 200))[1]
  if (is.na(header_line)) stop("GPL文件中未找到ID表头: ", path, call. = FALSE)
  ann <- data.table::fread(path, skip = header_line - 1L, check.names = FALSE, data.table = FALSE)
  symbol_col <- grep("^Gene Symbol$|gene.?symbol", names(ann), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(symbol_col) || !"ID" %in% names(ann)) stop("GPL注释缺少ID或Gene Symbol列", call. = FALSE)
  symbol <- trimws(as.character(ann[[symbol_col]]))
  symbol <- sub("\\s*///.*$", "", symbol)
  symbol <- sub("\\s*//.*$", "", symbol)
  symbol[symbol %in% c("", "---", "NA")] <- NA_character_
  data.frame(probe_id = as.character(ann$ID), gene_symbol = symbol, stringsAsFactors = FALSE)
}

aggregate_probes_mean <- function(expr_probe, annotation) {
  stopifnot(is.matrix(expr_probe) || is.data.frame(expr_probe))
  expr_probe <- as.matrix(expr_probe)
  idx <- match(rownames(expr_probe), annotation$probe_id)
  sym <- annotation$gene_symbol[idx]
  keep <- !is.na(sym) & nzchar(sym)
  expr_probe <- expr_probe[keep, , drop = FALSE]
  sym <- sym[keep]
  # limma::avereps在log2尺度上按基因对多探针取均值；快速、确定且不使用分组标签。
  gene_mat <- limma::avereps(expr_probe, ID = sym)
  storage.mode(gene_mat) <- "double"
  gene_mat
}

combine_characteristics <- function(pdata) {
  cc <- grep("^characteristics_ch1", names(pdata), value = TRUE)
  apply(pdata[, cc, drop = FALSE], 1, function(z) paste(na.omit(as.character(z)), collapse = "; "))
}
