# 数据质量与泄漏检查
validate_expression_metadata <- function(expr, meta, sample_col = "sample_id") {
  problems <- character()
  if (anyDuplicated(colnames(expr))) problems <- c(problems, "表达矩阵存在重复样本列")
  if (anyDuplicated(meta[[sample_col]])) problems <- c(problems, "metadata存在重复样本ID")
  if (!setequal(colnames(expr), meta[[sample_col]])) problems <- c(problems, "表达矩阵与metadata样本集合不一致")
  if (anyNA(expr)) problems <- c(problems, "表达矩阵含NA")
  if (length(problems)) stop(paste(problems, collapse = "; "), call. = FALSE)
  invisible(TRUE)
}

