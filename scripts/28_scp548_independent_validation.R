# Purpose: independent donor-level validation in SCP548 (Reyes et al.).
# Only cells with an exact public-portal barcode match are used. Donor, not
# cell, is the statistical unit. Outputs are individual figures.
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")

with_script_log("28_scp548_independent_validation", {
  suppressPackageStartupMessages({
    library(data.table); library(SingleCellExperiment); library(SummarizedExperiment)
    library(Matrix); library(edgeR); library(ggplot2)
  })
  ddir <- file.path("01_data_inventory", "downloaded_data")
  sce_file <- file.path(ddir, "20220219_pbmc130k_Reyes2020.rds")
  ann_files <- c("SCP548_Cohort_cell_values.tsv", "SCP548_donor_id_cell_values.tsv",
                 "SCP548_Cell_State_cell_values.tsv", "SCP548_Cell_Type_cell_values.tsv")
  assert_file(sce_file); invisible(lapply(file.path(ddir, ann_files), assert_file))

  ann <- Reduce(function(x, y) merge(x, y, by = "NAME", all = TRUE),
                lapply(file.path(ddir, ann_files), data.table::fread))
  ann$rds_cell <- gsub("-", ".", ann$NAME, fixed = TRUE)
  sce <- readRDS(sce_file)
  idx <- match(colnames(sce), ann$rds_cell)
  matched <- !is.na(idx)
  if (sum(matched) < 20000) stop("Too few SCP548 cells map exactly to the public annotations")
  cell_md <- as.data.frame(ann[idx[matched], ])
  rownames(cell_md) <- colnames(sce)[matched]
  counts <- SummarizedExperiment::assay(sce, "counts")[, matched, drop = FALSE]
  genes <- as.character(SummarizedExperiment::rowData(sce)$gene)
  if (length(genes) != nrow(counts)) stop("SCP548 gene metadata mismatch")
  rownames(counts) <- make.unique(genes)

  donor_cohort <- unique(cell_md[, c("donor_id", "Cohort")])
  if (anyDuplicated(donor_cohort$donor_id)) stop("A donor maps to multiple SCP548 cohorts")
  donor_counts <- counts %*% Matrix::sparse.model.matrix(~ 0 + factor(cell_md$donor_id))
  colnames(donor_counts) <- sub("^factor\\(cell_md\\$donor_id\\)", "", colnames(donor_counts))
  donor_cohort <- donor_cohort[match(colnames(donor_counts), donor_cohort$donor_id), ]
  n_cells <- table(cell_md$donor_id)
  donor_cohort$n_monocytes <- as.integer(n_cells[donor_cohort$donor_id])
  saveRDS(donor_counts, file.path("05_single_cell_pseudobulk", "SCP548_monocyte_donor_pseudobulk_counts.rds"))
  write_csv_utf8(donor_cohort, file.path("05_single_cell_pseudobulk", "SCP548_donor_metadata.csv"))

  y <- edgeR::DGEList(donor_counts)
  keep <- rowSums(edgeR::cpm(y) > 1) >= 3
  y <- edgeR::calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  gene_value <- function(g) {
    hit <- which(toupper(rownames(logcpm)) == toupper(g))
    if (!length(hit)) return(rep(NA_real_, ncol(logcpm)))
    colMeans(logcpm[hit, , drop = FALSE])
  }
  score <- function(gs) {
    m <- do.call(rbind, lapply(intersect(gs, rownames(logcpm)), gene_value))
    if (is.null(dim(m))) m <- matrix(m, nrow = 1)
    colMeans(t(scale(t(m))), na.rm = TRUE)
  }
  signature <- read.csv(file.path("07_non_circular_signature", "independent_signature_genes.csv"), check.names = FALSE)
  sig_genes <- unique(signature$gene[signature$signature == "Independent inflammatory, no S100A8/S100A9"])
  hla_genes <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA")
  nfkb_genes <- c("NFKBIA", "NFKB1", "RELA", "TNFAIP3", "IL1B", "CXCL8")
  donor <- data.frame(donor_id = colnames(logcpm),
                      cohort = donor_cohort$Cohort[match(colnames(logcpm), donor_cohort$donor_id)],
                      n_monocytes = donor_cohort$n_monocytes[match(colnames(logcpm), donor_cohort$donor_id)],
                      S100A8 = gene_value("S100A8"),
                      independent_state_score = score(sig_genes),
                      antigen_presentation_score = score(hla_genes),
                      NFkB_score = score(nfkb_genes))

  state_tab <- as.data.frame.matrix(table(cell_md$donor_id, cell_md$Cell_State))
  state_tab$donor_id <- rownames(state_tab)
  state_cols <- intersect(c("MS1", "MS2", "MS3", "MS4"), names(state_tab))
  state_tab$total <- rowSums(state_tab[, state_cols, drop = FALSE])
  for (s in state_cols) state_tab[[paste0(s, "_fraction")]] <- state_tab[[s]] / state_tab$total
  donor <- merge(donor, state_tab[, c("donor_id", paste0(state_cols, "_fraction"))], by = "donor_id", all.x = TRUE)
  write_csv_utf8(donor, file.path("05_single_cell_pseudobulk", "SCP548_donor_level_scores.csv"))

  comparisons <- list(c("Control", "Leuk-UTI"), c("Leuk-UTI", "URO"),
                      c("Leuk-UTI", "Int-URO"), c("ICU-NoSEP", "ICU-SEP"),
                      c("Control", "Bac-SEP"))
  measures <- c("S100A8", "independent_state_score", "antigen_presentation_score", "NFkB_score", "MS1_fraction")
  tests <- list(); k <- 0L
  for (m in measures) for (cmp in comparisons) {
    dd <- donor[donor$cohort %in% cmp & is.finite(donor[[m]]), ]
    if (all(table(factor(dd$cohort, levels = cmp)) >= 3)) {
      wt <- wilcox.test(dd[[m]][dd$cohort == cmp[[2]]], dd[[m]][dd$cohort == cmp[[1]]], exact = FALSE)
      k <- k + 1L
      tests[[k]] <- data.frame(measure = m, reference = cmp[[1]], comparison = cmp[[2]],
                               n_reference = sum(dd$cohort == cmp[[1]]), n_comparison = sum(dd$cohort == cmp[[2]]),
                               median_difference = median(dd[[m]][dd$cohort == cmp[[2]]]) - median(dd[[m]][dd$cohort == cmp[[1]]]),
                               p_value = wt$p.value)
    }
  }
  tests <- do.call(rbind, tests); tests$fdr <- p.adjust(tests$p_value, method = "BH")
  write_csv_utf8(tests, file.path("05_single_cell_pseudobulk", "Table_SCP548_DonorComparisons.csv"))

  cohort_order <- c("Control", "Leuk-UTI", "URO", "Int-URO", "ICU-NoSEP", "ICU-SEP", "Bac-SEP")
  donor$cohort <- factor(donor$cohort, levels = cohort_order)
  cohort_cols <- c(Control = palette_nature[["Normal"]], `Leuk-UTI` = "#729CB3", URO = "#D39B44",
                   `Int-URO` = "#B87945", `ICU-NoSEP` = "#7C8792", `ICU-SEP` = palette_nature[["Sepsis"]],
                   `Bac-SEP` = "#9F3F32")
  make_donor_plot <- function(yvar, ylab, stem, height = 76) {
    p <- ggplot(donor, aes(cohort, .data[[yvar]], colour = cohort)) +
      geom_boxplot(width = .5, outlier.shape = NA, colour = "#5F6368", fill = NA, linewidth = .38) +
      geom_point(size = 1.65, alpha = .88, position = position_jitter(width = .065, height = 0)) +
      scale_colour_manual(values = cohort_cols, drop = FALSE) +
      labs(x = NULL, y = ylab, colour = NULL) + theme_nature() +
      theme(axis.text.x = element_text(angle = 32, hjust = 1), legend.position = "none")
    save_nature_plot(p, file.path("09_figures", stem), 120, height, source_data = donor[, c("donor_id", "cohort", "n_monocytes", yvar)])
  }
  make_donor_plot("S100A8", expression(S100A8~"donor pseudo-bulk logCPM"), "Fig_SCP548_DonorPseudobulk_S100A8")
  make_donor_plot("MS1_fraction", "MS1 fraction per donor", "Fig_SCP548_MS1_Abundance")
  make_donor_plot("antigen_presentation_score", "Antigen-presentation score (scaled)", "Fig_SCP548_AntigenPresentation")
  make_donor_plot("independent_state_score", "Independent emergency-myeloid score", "Fig_SCP548_ClinicalGradient")

  ct <- cor.test(donor$S100A8, donor$antigen_presentation_score, method = "spearman", exact = FALSE)
  corr <- data.frame(dataset = "SCP548", n_donors = nrow(donor), rho = unname(ct$estimate), p_value = ct$p.value)
  write_csv_utf8(corr, file.path("08_bulk_singlecell_integration", "SCP548_S100A8_HLAII_correlation.csv"))
  p_cor <- ggplot(donor, aes(S100A8, antigen_presentation_score, colour = cohort)) +
    geom_point(size = 1.9, alpha = .9) + geom_smooth(method = "lm", se = TRUE, colour = "#5F6368", fill = "#D8D8D8", linewidth = .55) +
    scale_colour_manual(values = cohort_cols) +
    labs(x = expression(S100A8~"donor logCPM"), y = "Antigen-presentation score", colour = NULL) + theme_nature() +
    theme(legend.position = "right")
  save_nature_plot(p_cor, file.path("09_figures", "Fig_SCP548_S100A8_HLAII_Correlation"), 105, 76, source_data = donor)

  audit <- data.frame(resource = c("Public portal annotations", "Secondary public count object", "Exact matched analysis set"),
                      cells = c(nrow(ann), ncol(sce), sum(matched)),
                      donors = c(length(unique(ann$donor_id)), NA, length(unique(cell_md$donor_id))),
                      rule = c("All portal monocytes", "Public processed RDS", "Exact transformed barcode match only"))
  write_csv_utf8(audit, file.path("01_data_inventory", "SCP548_mapping_audit.csv"))
})
