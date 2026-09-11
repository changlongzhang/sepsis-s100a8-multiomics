options(stringsAsFactors = FALSE)
root <- "C:/Users/ZCL/sepsis_workspace_link"
env <- new.env(parent = emptyenv())
load(file.path(root, "sepsis-s100a8-multiomics", "data", "03_WGCNA_analysis", "output", "03_WGCNA_Analysis.RData"), envir = env)

stopifnot(exists("mergedMEs", envir = env), exists("group_info", envir = env))
me <- get("mergedMEs", envir = env)
group_info <- get("group_info", envir = env)
stopifnot("MEblack" %in% colnames(me))

out <- data.frame(
  Sample = rownames(me),
  ME = as.numeric(me[, "MEblack"]),
  Group = as.character(group_info[rownames(me), "Group"]),
  stringsAsFactors = FALSE
)
stopifnot(!anyNA(out$ME), !anyNA(out$Group))
out$Group <- factor(out$Group, levels = c("Normal", "Disease"))

dest <- file.path(root, "review_revision", "10_figure_sources", "Fig1_S1_corrected_panels", "Fig1C_source_data.csv")
write.csv(out, dest, row.names = FALSE)
cat("ROWS\t", nrow(out), "\n", sep = "")
cat("GROUPS\t", paste(names(table(out$Group)), as.integer(table(out$Group)), collapse = ";"), "\n", sep = "")
cat("OUTPUT\t", dest, "\n", sep = "")
