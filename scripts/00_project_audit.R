# 输入：原始项目 D:/桌面/sepsis（只读）
# 输出：00_project_audit/file_inventory.csv、重复文件摘要、现有R脚本快照与审计摘要
# 目的：建立可追踪的原始项目基线；不删除、不移动、不覆盖任何原始文件。

rm(list = ls())
options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
set.seed(MASTER_SEED)

start_time <- Sys.time()
log_file <- file.path("11_logs", "00_project_audit.log")
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
sink(log_file, split = TRUE)
on.exit({
  cat("\n结束时间:", format(Sys.time()), "\n耗时(秒):", as.numeric(difftime(Sys.time(), start_time, units = "secs")), "\n")
  print(sessionInfo())
  sink()
}, add = TRUE)

cat("开始时间:", format(start_time), "\n随机种子:", MASTER_SEED, "\n")

all_files <- list.files(PROJECT_ROOT, recursive = TRUE, full.names = TRUE,
                        all.files = TRUE, no.. = TRUE, include.dirs = FALSE)
# 排除本次返修目录，确保清单只描述原始项目基线。
all_files <- all_files[!startsWith(normalizePath(all_files, winslash = "/", mustWork = FALSE),
                                   paste0(REVISION_ROOT, "/"))]
info <- file.info(all_files)

infer_dataset <- function(x) {
  hits <- regmatches(x, gregexpr("GSE[0-9]+", x, ignore.case = TRUE))
  vapply(hits, function(z) if (length(z) && !identical(z, character(0))) paste(unique(toupper(z)), collapse = ";") else NA_character_, character(1))
}

infer_content <- function(path, ext) {
  x <- tolower(path)
  if (grepl("series_matrix|soft", x)) return("GEO expression/metadata")
  if (grepl("matrix\\.mtx|barcodes|features", x)) return("10x raw counts component")
  if (ext %in% c("rds", "rdata")) return("R serialized object")
  if (ext %in% c("r", "rmd", "qmd")) return("analysis code")
  if (ext %in% c("docx", "pdf", "tex")) return("manuscript/review/document")
  if (grepl("clinical|生存|survival|group|分组|metadata|phenotype", x)) return("sample/clinical metadata")
  if (grepl("gpl|platform|annot", x)) return("platform annotation")
  if (ext %in% c("csv", "tsv", "txt", "xlsx")) return("tabular data/result")
  if (ext %in% c("png", "tif", "tiff", "svg", "ai")) return("figure/image")
  "other"
}

infer_use <- function(dataset, content, path) {
  if (grepl("GSE65682", dataset %||% "")) return("discovery/class-imbalance/prognosis")
  if (grepl("GSE26440|GSE28750|GSE69528|GSE9692", dataset %||% "")) return("external clinical validation")
  if (grepl("GSE167363", dataset %||% "")) return("single-cell donor-level analysis")
  if (content == "analysis code") return("existing workflow audit/reuse")
  if (content == "manuscript/review/document") return("manuscript/reviewer context")
  "audit/reference"
}

`%||%` <- function(x, y) if (is.na(x) || !nzchar(x)) y else x
ext <- sub("^\\.", "", tolower(tools::file_ext(all_files)))
dataset <- infer_dataset(all_files)
content <- mapply(infer_content, all_files, ext, USE.NAMES = FALSE)
inventory <- data.frame(
  full_path = normalizePath(all_files, winslash = "/", mustWork = FALSE),
  file_name = basename(all_files),
  extension = ext,
  file_size = as.numeric(info$size),
  modified_time = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
  inferred_dataset = dataset,
  inferred_content = content,
  intended_use = mapply(infer_use, dataset, content, all_files, USE.NAMES = FALSE),
  status = ifelse(info$size > 0, "available_local", "empty_file_review"),
  stringsAsFactors = FALSE
)
inventory <- inventory[order(inventory$full_path), ]
write.csv(inventory, file.path("00_project_audit", "file_inventory.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 仅对同尺寸文件进一步计算MD5，识别明显重复副本。
size_tab <- table(inventory$file_size)
dup_sizes <- as.numeric(names(size_tab[size_tab > 1 & as.numeric(names(size_tab)) > 0 & as.numeric(names(size_tab)) <= 500 * 1024^2]))
dup_idx <- inventory$file_size %in% dup_sizes
dup <- inventory[dup_idx, c("full_path", "file_name", "file_size", "inferred_dataset")]
dup$md5 <- if (nrow(dup)) unname(tools::md5sum(dup$full_path)) else character(0)
dup$duplicate_count <- if (nrow(dup)) ave(seq_len(nrow(dup)), dup$md5, FUN = length) else integer(0)
dup <- dup[dup$duplicate_count > 1, ]
write.csv(dup, file.path("00_project_audit", "duplicate_file_candidates.csv"), row.names = FALSE, fileEncoding = "UTF-8")

# 大于500 MB的同尺寸文件只列为“可能重复”，避免首次审计对数GB对象反复做全文件哈希。
large_dup_sizes <- as.numeric(names(size_tab[size_tab > 1 & as.numeric(names(size_tab)) > 500 * 1024^2]))
large_dup <- inventory[inventory$file_size %in% large_dup_sizes,
                       c("full_path", "file_name", "file_size", "inferred_dataset")]
write.csv(large_dup, file.path("00_project_audit", "large_same_size_candidates_unhashed.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# 快照原有分析脚本，供审计和复用；原文件保持不变。
code_files <- inventory$full_path[inventory$inferred_content == "analysis code"]
snapshot_root <- file.path("00_project_audit", "existing_code_snapshot")
dir.create(snapshot_root, recursive = TRUE, showWarnings = FALSE)
if (length(code_files)) {
  safe_name <- gsub("[^A-Za-z0-9._-]", "_", substring(code_files, nchar(PROJECT_ROOT) + 2))
  copied <- file.copy(code_files, file.path(snapshot_root, safe_name), overwrite = FALSE)
  code_index <- data.frame(original_path = code_files, snapshot_file = file.path(snapshot_root, safe_name), copied = copied)
  write.csv(code_index, file.path("00_project_audit", "existing_code_index.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

summary_lines <- c(
  "# Project audit summary",
  "",
  paste0("- Audit time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("- Original files inventoried: ", nrow(inventory)),
  paste0("- R/Rmd/Qmd scripts: ", sum(inventory$inferred_content == "analysis code")),
  paste0("- GEO series/raw-count files: ", sum(grepl("GEO expression|10x raw", inventory$inferred_content))),
  paste0("- RDS/RData objects: ", sum(inventory$inferred_content == "R serialized object")),
  paste0("- Manuscript/review documents: ", sum(inventory$inferred_content == "manuscript/review/document")),
  paste0("- Exact duplicate candidates: ", nrow(dup)),
  "",
  "No original files were modified, moved, or deleted."
)
writeLines(summary_lines, file.path("00_project_audit", "audit_summary.md"), useBytes = TRUE)
cat(paste(summary_lines, collapse = "\n"), "\n")
