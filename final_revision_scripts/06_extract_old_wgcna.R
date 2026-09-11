options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(WGCNA))

old_file <- "C:/Users/ZCL/sepsis_workspace_link/WGCNA/wgcna_analysis_complete.RData"
out_dir <- "C:/Users/ZCL/sepsis_workspace_link/work/final_revision_20260908/old_wgcna_extract"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(old_file)
object_inventory <- data.frame(
  object = ls(),
  class = vapply(ls(), function(nm) paste(class(get(nm)), collapse = ";"), character(1)),
  stringsAsFactors = FALSE
)
write.csv(object_inventory, file.path(out_dir, "old_RData_object_inventory.csv"), row.names = FALSE)

get_first <- function(candidates) {
  hit <- candidates[candidates %in% ls(envir = .GlobalEnv)]
  if (!length(hit)) return(NULL)
  get(hit[[1]], envir = .GlobalEnv)
}

datExpr_old <- get_first(c("datExpr", "datExpr0", "expr_data", "expr_data_std"))
module_colors_old <- get_first(c("mergedColors", "moduleColors", "dynamicColors"))
trait_old <- get_first(c("traitData", "datTraits", "Traits", "trait", "group_info", "group_data"))

if (is.null(datExpr_old)) stop("No datExpr/datExpr0 object found in old RData")
if (is.null(module_colors_old)) stop("No module color vector found in old RData")
if (length(module_colors_old) == nrow(datExpr_old)) datExpr_old <- t(datExpr_old)
if (length(module_colors_old) != ncol(datExpr_old)) stop("Module colors do not align with datExpr columns")

module_colors_old <- as.character(module_colors_old)
names(module_colors_old) <- colnames(datExpr_old)

if (!"S100A8" %in% colnames(datExpr_old)) stop("S100A8 absent from old datExpr")
s100_module <- unname(module_colors_old[["S100A8"]])

MEs_old <- orderMEs(moduleEigengenes(datExpr_old, colors = module_colors_old)$eigengenes)
me_name <- paste0("ME", s100_module)
if (!me_name %in% colnames(MEs_old)) stop("S100A8 module eigengene not found")
s100_mm <- as.numeric(WGCNA::cor(datExpr_old[, "S100A8"], MEs_old[, me_name], use = "p"))

resolve_sepsis_trait <- function(x, sample_names) {
  if (is.null(x)) return(NULL)
  if (is.vector(x) && length(x) == length(sample_names)) {
    vals <- x
  } else {
    x <- as.data.frame(x)
    if (!is.null(rownames(x)) && all(sample_names %in% rownames(x))) x <- x[sample_names, , drop = FALSE]
    numeric_cols <- names(x)[vapply(x, is.numeric, logical(1))]
    binary_cols <- numeric_cols[vapply(x[numeric_cols], function(z) length(unique(na.omit(z))) == 2, logical(1))]
    preferred <- grep("sepsis|disease|status|group", names(x), ignore.case = TRUE, value = TRUE)
    preferred <- intersect(preferred, binary_cols)
    use_col <- if (length(preferred)) preferred[[1]] else if (length(binary_cols)) binary_cols[[1]] else NA_character_
    if (is.na(use_col)) return(NULL)
    vals <- x[[use_col]]
  }
  as.numeric(vals)
}

sepsis_trait <- resolve_sepsis_trait(trait_old, rownames(datExpr_old))
s100_gs <- if (is.null(sepsis_trait)) NA_real_ else as.numeric(WGCNA::cor(datExpr_old[, "S100A8"], sepsis_trait, use = "p"))

stats <- data.frame(
  gene = "S100A8",
  module = s100_module,
  module_membership = s100_mm,
  gene_significance = s100_gs,
  sample_n = nrow(datExpr_old),
  gene_n = ncol(datExpr_old),
  stringsAsFactors = FALSE
)
write.csv(stats, file.path(out_dir, "old_S100A8_WGCNA_statistics.csv"), row.names = FALSE)
write.csv(as.data.frame(table(module = module_colors_old)), file.path(out_dir, "old_module_sizes.csv"), row.names = FALSE)
saveRDS(list(datExpr = datExpr_old, colors = module_colors_old, MEs = MEs_old, trait = sepsis_trait),
        file.path(out_dir, "old_wgcna_core_objects.rds"), compress = TRUE)
print(stats)
