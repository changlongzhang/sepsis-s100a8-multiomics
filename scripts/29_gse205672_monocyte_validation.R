# Purpose: independent purified-monocyte validation in GSE205672.
# Only unstimulated samples are included; one donor/sample is one unit.
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")

with_script_log("29_gse205672_monocyte_validation", {
  suppressPackageStartupMessages({library(GEOquery); library(Biobase); library(data.table); library(ggplot2)})
  ddir <- file.path("01_data_inventory", "downloaded_data")
  es <- GEOquery::getGEO(filename = file.path(ddir, "GSE205672_series_matrix.txt.gz"), getGPL = FALSE)
  pd <- Biobase::pData(es)
  x <- data.table::fread(file.path(ddir, "GSE205672_AllSamples.GeneExpression.FPKM.txt.gz"))
  sample_cols <- grep("_FPKM$", names(x), value = TRUE)
  sample_key <- sub("_FPKM$", "", sample_cols)
  title <- as.character(pd$title)
  keep_pd <- pd[["cell type:ch1"]] == "Monocytes" & grepl("_unstim_", title, ignore.case = TRUE)
  pd <- pd[keep_pd, , drop = FALSE]; title <- as.character(pd$title)
  ix <- match(title, sample_key)
  if (sum(!is.na(ix)) < 40) stop("Too few unstimulated monocyte samples map to GSE205672 FPKM")
  pd <- pd[!is.na(ix), , drop = FALSE]; ix <- ix[!is.na(ix)]

  symbols <- as.character(x$SymbolID)
  selected_cols <- sample_cols[ix]
  mat <- as.matrix(x[, selected_cols, with = FALSE]); storage.mode(mat) <- "numeric"
  mat <- log2(mat + 1)
  valid <- !is.na(symbols) & nzchar(symbols) & symbols != "NA"
  mat <- rowsum(mat[valid, , drop = FALSE], group = symbols[valid], reorder = FALSE)
  denom <- table(symbols[valid])[rownames(mat)]
  mat <- mat / as.numeric(denom)
  colnames(mat) <- title[!is.na(match(title, sample_key))]

  sig <- read.csv(file.path("07_non_circular_signature", "independent_signature_genes.csv"), check.names = FALSE)
  sig_genes <- unique(sig$gene[sig$signature == "Independent inflammatory, no S100A8/S100A9"])
  hla_genes <- c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA")
  nfkb_genes <- c("NFKBIA", "NFKB1", "RELA", "TNFAIP3", "IL1B", "CXCL8")
  score <- function(gs) {
    gs <- intersect(gs, rownames(mat)); if (length(gs) < 2) return(rep(NA_real_, ncol(mat)))
    colMeans(t(scale(t(mat[gs, , drop = FALSE]))), na.rm = TRUE)
  }
  donor <- data.frame(sample_id = rownames(pd), title = as.character(pd$title),
                      donor_id = sub("_.*$", "", as.character(pd$title)),
                      group = ifelse(pd[["disease state:ch1"]] == "Healthy control", "Control", "Sepsis"),
                      S100A8 = as.numeric(mat[TARGET_GENE, as.character(pd$title)]),
                      independent_state_score = score(sig_genes),
                      antigen_presentation_score = score(hla_genes),
                      NFkB_score = score(nfkb_genes))
  write_csv_utf8(donor, file.path("07_non_circular_signature", "GSE205672_monocyte_scores.csv"))

  measures <- c("S100A8", "independent_state_score", "antigen_presentation_score", "NFkB_score")
  tests <- do.call(rbind, lapply(measures, function(m) {
    wt <- wilcox.test(donor[[m]] ~ donor$group, exact = FALSE)
    data.frame(dataset = "GSE205672", measure = m, n_control = sum(donor$group == "Control"),
               n_sepsis = sum(donor$group == "Sepsis"),
               median_control = median(donor[[m]][donor$group == "Control"], na.rm = TRUE),
               median_sepsis = median(donor[[m]][donor$group == "Sepsis"], na.rm = TRUE),
               p_value = wt$p.value)
  }))
  tests$fdr <- p.adjust(tests$p_value, "BH")
  write_csv_utf8(tests, file.path("07_non_circular_signature", "Table_GSE205672_MonocyteValidation.csv"))

  make_plot <- function(yvar, ylab, stem) {
    p <- ggplot(donor, aes(group, .data[[yvar]], colour = group)) +
      geom_boxplot(width = .45, outlier.shape = NA, colour = "#5F6368", fill = NA, linewidth = .4) +
      geom_point(size = 1.75, alpha = .85, position = position_jitter(width = .065, height = 0)) +
      scale_colour_manual(values = c(Control = palette_nature[["Normal"]], Sepsis = palette_nature[["Sepsis"]])) +
      labs(x = NULL, y = ylab, colour = NULL) + theme_nature() + theme(legend.position = "none")
    save_nature_plot(p, file.path("09_figures", stem), 89, 70,
                     source_data = donor[, c("sample_id", "donor_id", "group", yvar)])
  }
  make_plot("S100A8", expression(log[2](S100A8~FPKM+1)), "Fig_GSE205672_S100A8_Monocytes")
  make_plot("independent_state_score", "Independent emergency-myeloid score", "Fig_GSE205672_IndependentStateScore")
  make_plot("antigen_presentation_score", "Antigen-presentation score", "Fig_GSE205672_AntigenPresentation")
  make_plot("NFkB_score", "NF-kB transcriptional score", "Fig_GSE205672_NFkBScore")

  corr_vars <- c("independent_state_score", "antigen_presentation_score", "NFkB_score")
  cors <- do.call(rbind, lapply(corr_vars, function(v) {
    ct <- cor.test(donor$S100A8, donor[[v]], method = "spearman", exact = FALSE)
    data.frame(dataset = "GSE205672", x = "S100A8", y = v, n = nrow(donor),
               rho = unname(ct$estimate), p_value = ct$p.value)
  }))
  cors$fdr <- p.adjust(cors$p_value, "BH")
  write_csv_utf8(cors, file.path("08_bulk_singlecell_integration", "GSE205672_mechanistic_correlations.csv"))
  p_cor <- ggplot(donor, aes(S100A8, antigen_presentation_score, colour = group)) +
    geom_point(size = 1.9, alpha = .9) + geom_smooth(method = "lm", se = TRUE, colour = "#5F6368", fill = "#D8D8D8", linewidth = .55) +
    scale_colour_manual(values = c(Control = palette_nature[["Normal"]], Sepsis = palette_nature[["Sepsis"]])) +
    labs(x = expression(log[2](S100A8~FPKM+1)), y = "Antigen-presentation score", colour = NULL) + theme_nature()
  save_nature_plot(p_cor, file.path("09_figures", "Fig_GSE205672_S100A8_HLAII_Correlation"), 89, 70, source_data = donor)
})
