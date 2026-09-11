options(stringsAsFactors = FALSE)
.libPaths(c("C:/Users/ZCL/AppData/Local/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages(library(data.table))

root <- Sys.getenv("SEPSIS_PROJECT_ROOT", unset = normalizePath("..", winslash = "/", mustWork = TRUE))
out_dir <- file.path(root, "work", "final_revision_20260908", "wgcna_rerun")
res <- readRDS(file.path(out_dir, "WGCNA_final_results.rds"))
raw <- fread(file.path(root, "sepsis-s100a8-multiomics", "data", "02_DEG", "Input", "GSE65682_gene.csv"), data.table = FALSE, check.names = FALSE)
gene_col <- names(raw)[1]
idx <- which(raw[[gene_col]] == "S100A8")
stopifnot(length(idx) == 1L)
s100 <- log2(as.numeric(raw[idx, -1, drop = TRUE]))
names(s100) <- names(raw)[-1]
stopifnot(identical(rownames(res$MEs), names(s100)))

r <- vapply(res$MEs, function(me) stats::cor(s100, me, method = "pearson", use = "complete.obs"), numeric(1))
p <- 2 * stats::pt(-abs(r * sqrt((length(s100) - 2) / pmax(1e-15, 1 - r^2))), df = length(s100) - 2)
projection <- data.frame(
  module = sub("^ME", "", names(r)),
  correlation_with_S100A8 = unname(r),
  p_value = unname(p),
  FDR = p.adjust(p, method = "BH"),
  stringsAsFactors = FALSE
)
projection <- projection[order(-abs(projection$correlation_with_S100A8)), ]
write.csv(projection, file.path(out_dir, "S100A8_posthoc_module_eigengene_projection.csv"), row.names = FALSE)

trait <- as.numeric(fread(file.path(root, "sepsis-s100a8-multiomics", "data", "02_DEG", "Input", "GSE65682_Groups.csv"), data.table = FALSE)[[2]] == "Disease")
gs <- stats::cor(s100, trait, method = "pearson")
gs_p <- 2 * stats::pt(-abs(gs * sqrt((length(s100) - 2) / (1 - gs^2))), df = length(s100) - 2)
filter_audit <- fread(file.path(out_dir, "gene_filter_audit.csv"), data.table = FALSE)
s100_filter <- filter_audit[filter_audit$gene == "S100A8", ]
audit <- data.frame(
  gene = "S100A8", retained_in_network = FALSE, module_membership = NA_real_,
  gene_significance_all_802 = gs, gene_significance_p = gs_p,
  MAD = s100_filter$MAD, MAD_threshold = median(filter_audit$MAD[is.finite(filter_audit$MAD)]),
  posthoc_projection_only = TRUE,
  interpretation = "S100A8 was not assigned to a module; eigengene correlations are post hoc projections and are not module membership.",
  stringsAsFactors = FALSE
)
write.csv(audit, file.path(out_dir, "S100A8_WGCNA_statistics.csv"), row.names = FALSE)
print(projection)
print(audit)
