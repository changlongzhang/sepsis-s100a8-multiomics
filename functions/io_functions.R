# 通用输入输出与日志函数
assert_file <- function(path, label = "input") {
  if (!file.exists(path)) stop(sprintf("缺少%s文件: %s", label, path), call. = FALSE)
  invisible(normalizePath(path, winslash = "/", mustWork = TRUE))
}

assert_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("无法创建目录: ", path, call. = FALSE)
  }
  invisible(path)
}

with_script_log <- function(script_name, expr, seed = MASTER_SEED) {
  assert_dir("11_logs")
  log_file <- file.path("11_logs", paste0(script_name, ".log"))
  start <- Sys.time()
  warnings_seen <- character()
  error_seen <- NA_character_
  status <- "Success"
  outputs_verified <- FALSE
  sink(log_file, split = TRUE)
  on.exit(sink(), add = TRUE)
  cat("脚本:", script_name, "\n开始时间:", format(start), "\n随机种子:", seed, "\n")
  set.seed(seed)
  value <- tryCatch(
    withCallingHandlers(force(expr), warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    }),
    error = function(e) {
      status <<- "Failed"
      error_seen <<- conditionMessage(e)
      cat("ERROR:", error_seen, "\n")
      NULL
    }
  )
  end <- Sys.time()
  cat("结束时间:", format(end), "\n耗时(秒):", as.numeric(difftime(end, start, units = "secs")), "\n")
  print(sessionInfo())
  status_row <- data.frame(
    script = script_name, status = status, start_time = format(start), end_time = format(end),
    runtime = as.numeric(difftime(end, start, units = "secs")),
    warnings = paste(unique(warnings_seen), collapse = " | "), errors = error_seen,
    outputs_verified = outputs_verified, stringsAsFactors = FALSE
  )
  assert_dir("11_logs/status_parts")
  write.csv(status_row, file.path("11_logs/status_parts", paste0(script_name, ".csv")), row.names = FALSE, fileEncoding = "UTF-8")
  if (status == "Failed") stop(sprintf("脚本%s失败，详见%s: %s", script_name, log_file, error_seen), call. = FALSE)
  invisible(value)
}

write_csv_utf8 <- function(x, path) {
  assert_dir(dirname(path))
  write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8", na = "")
  if (!file.exists(path) || file.info(path)$size == 0) stop("输出文件写入失败: ", path, call. = FALSE)
  invisible(path)
}

