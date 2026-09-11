options(stringsAsFactors = FALSE)
invisible(try(Sys.setlocale("LC_CTYPE", ".UTF-8"), silent = TRUE))
suppressPackageStartupMessages({
  library(magick)
  library(grid)
  library(ggplot2)
  library(svglite)
  library(ragg)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Final PLOS ONE submission figure builder. All visual work is performed in R.
# Original submitted figures and raw data are read-only; only the dedicated
# revision output folder is overwritten.
project_root <- "D:/\u684c\u9762/sepsis"
out_root <- file.path(project_root, "\u6587\u7ae0", "\u8fd4\u4fee\u6295\u7a3f\u65b0\u4e3b\u56fe")
work_root <- file.path(tempdir(), "sepsis_plos_figure_panels")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(work_root, recursive = TRUE, showWarnings = FALSE)

fig_root <- file.path(project_root, "review_revision", "09_figures")
sc_root <- file.path(project_root, "sepsis-s100a8-multiomics", "data",
                     "09_scRNA-seq & Bulk Immune Infiltration", "output")
exp_root <- file.path(project_root, "\u5b9e\u9a8c")
exp_ready <- file.path(exp_root, "PLOS_submission_ready_final")
exp_add <- file.path(exp_root, "\u8865\u5b9e\u9a8c", "\u8865\u5b9e\u9a8c", "PDF")
dock_root <- file.path(project_root,
  "S100A8_triptolide_16OH_comparative_docking_topjournal",
  "publication_figure_v2_complete")
md_root <- file.path(project_root, "S100A8_triptolide_MD_topjournal",
  "comparative_analysis_v3", "publication_figures_python_v2")

W <- 2250L
G <- 24L
ink <- "#2B2B2B"
blue <- "#3C6E9E"
red <- "#C5534E"
teal <- "#3A8F86"
gold <- "#D59B35"
group_cols <- c("Control" = "#697386", "TPL" = teal,
                "LPS" = red, "LPS+TPL" = blue)

theme_plos <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = ink),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, colour = ink),
      axis.line = element_line(linewidth = 0.35, colour = ink),
      axis.ticks = element_line(linewidth = 0.35, colour = ink),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(5, 7, 5, 5)
    )
}

save_plot_panel <- function(plot, stem, width = 1000, height = 700) {
  path <- file.path(work_root, paste0(stem, ".png"))
  ragg::agg_png(path, width = width, height = height, units = "px", res = 300,
                background = "white")
  print(plot)
  dev.off()
  path
}

read_any <- function(path) {
  if (!file.exists(path)) stop("Missing panel source: ", path)
  ext <- tolower(tools::file_ext(path))
  img <- if (ext == "pdf") image_read_pdf(path, density = 300)[1] else image_read(path)[1]
  img <- image_background(img, "white", flatten = TRUE)
  tryCatch(image_trim(img, fuzz = 3), error = function(e) img)
}

fit_panel <- function(img, width, height, inset = 22) {
  info <- image_info(img)
  scale <- min((width - 2 * inset) / info$width, (height - 2 * inset) / info$height)
  nw <- max(1, floor(info$width * scale))
  nh <- max(1, floor(info$height * scale))
  img <- image_resize(img, sprintf("%dx%d!", nw, nh), filter = "Lanczos")
  image_extent(img, sprintf("%dx%d", width, height), gravity = "center", color = "white")
}

assemble_rows <- function(rows, labels = LETTERS, width = W, gutter = G) {
  total_h <- sum(vapply(rows, function(x) x$height, numeric(1))) +
    gutter * (length(rows) - 1)
  if (total_h > 2625) stop("Figure exceeds PLOS ONE maximum height: ", total_h)
  canvas <- image_blank(width, total_h, color = "white")
  y <- 0L
  k <- 1L
  for (row in rows) {
    n <- length(row$paths)
    weights <- row$weights %||% rep(1, n)
    usable <- width - gutter * (n - 1)
    widths <- floor(usable * weights / sum(weights))
    widths[n] <- usable - sum(widths[-n])
    x <- 0L
    for (j in seq_len(n)) {
      panel <- fit_panel(read_any(row$paths[j]), widths[j], row$height)
      panel <- image_annotate(panel, labels[k], gravity = "northwest",
                              location = "+8+4", size = 42,
                              font = "Arial", weight = 700, color = ink)
      canvas <- image_composite(canvas, panel, offset = sprintf("+%d+%d", x, y))
      x <- x + widths[j] + gutter
      k <- k + 1L
    }
    y <- y + row$height + gutter
  }
  image_extent(canvas, sprintf("%dx%d", width, total_h), gravity = "center", color = "white")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

write_raster_container <- function(img, path, type = c("pdf", "svg")) {
  type <- match.arg(type)
  info <- image_info(img)
  wi <- info$width / 300
  hi <- info$height / 300
  if (type == "pdf") {
    grDevices::cairo_pdf(path, width = wi, height = hi, family = "Arial")
  } else {
    svglite::svglite(path, width = wi, height = hi, bg = "white")
  }
  grid.newpage()
  grid.raster(as.raster(img), width = unit(1, "npc"), height = unit(1, "npc"),
              interpolate = FALSE)
  dev.off()
}

export_figure <- function(img, stem) {
  img <- image_background(img, "white", flatten = TRUE)
  paths <- file.path(out_root, paste0(stem, c(".png", ".tif", ".pdf", ".svg")))
  image_write(img, paths[1], format = "png", density = "300x300")
  image_write(img, paths[2], format = "tiff", compression = "lzw", density = "300x300")
  write_raster_container(img, paths[3], "pdf")
  write_raster_container(img, paths[4], "svg")
  info <- image_info(img)
  data.frame(figure = stem, format = c("png", "tif", "pdf", "svg"),
             width_px = info$width, height_px = info$height,
             dpi = 300, bytes = file.info(paths)$size,
             path = normalizePath(paths, winslash = "/"))
}

# Remove only old generated submission figures in the dedicated output folder.
old <- list.files(out_root, full.names = TRUE)
old <- old[grepl("^(Fig[1-8]|S[1-8]_fig|FigS[0-9]+|figure_export_manifest|PLOS_figure_QA|interaction_effects|Figure_captions)",
                  basename(old), ignore.case = TRUE)]
if (length(old)) unlink(old, force = TRUE)

# ---------------------------------------------------------------------------
# Regenerated experimental panels from the deposited source-data table.
# ---------------------------------------------------------------------------
exp_dat <- suppressMessages(readr::read_csv(
  file.path(exp_ready, "Fig7_source_data.csv"), show_col_types = FALSE,
  locale = readr::locale(encoding = "UTF-8")))
exp_dat$value <- as.numeric(exp_dat$value)
exp_dat$replicate <- factor(exp_dat$replicate)

qpcr <- exp_dat %>% filter(panel %in% c("B", "C", "D", "E")) %>%
  mutate(group = factor(group, levels = c("Control", "LPS", "LPS+TPL")),
         endpoint = factor(endpoint, levels = c("IL-6", "TNF-\u03b1", "IL-1\u03b2", "S100A8")))
p_qpcr <- ggplot(qpcr, aes(group, value, colour = group)) +
  stat_summary(fun = mean, geom = "col", aes(fill = group), colour = NA,
               width = 0.66, alpha = 0.78) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18,
               linewidth = 0.4, colour = ink) +
  geom_point(position = position_jitter(width = 0.08, height = 0),
             size = 1.2, alpha = 0.85) +
  facet_wrap(~endpoint, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = group_cols) + scale_fill_manual(values = group_cols) +
  labs(x = NULL, y = "Relative mRNA expression") +
  theme_plos() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
qpcr_path <- save_plot_panel(p_qpcr, "Fig7_qPCR", 1800, 620)

wb_old <- exp_dat %>% filter(panel %in% c("G", "H", "I")) %>%
  mutate(endpoint = factor(endpoint,
    levels = c("S100A8 / \u03b2-actin", "NF-\u03baB / \u03b2-actin", "p-NF-\u03baB / NF-\u03baB")),
    group = factor(group, levels = c("Control", "LPS", "LPS+TPL")))
p_wb_old <- ggplot(wb_old, aes(group, value, colour = group)) +
  stat_summary(fun = mean, geom = "col", aes(fill = group), colour = NA,
               width = 0.66, alpha = 0.78) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18,
               linewidth = 0.4, colour = ink) +
  geom_point(position = position_jitter(width = 0.07), size = 1.35) +
  facet_wrap(~endpoint, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = group_cols) + scale_fill_manual(values = group_cols) +
  labs(x = NULL, y = "Normalized band intensity") + theme_plos() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
wb_old_path <- save_plot_panel(p_wb_old, "Fig7_WB_quant", 1500, 650)

kd8 <- exp_dat %>% filter(panel %in% c("N", "O")) %>%
  mutate(siRNA = factor(siRNA, levels = c("siNC", "siS100A8")),
         treatment = factor(treatment, levels = c("Control", "TPL", "LPS", "LPS+TPL")),
         display = factor(paste(siRNA, treatment, sep = "\n"),
           levels = c("siNC\nControl", "siNC\nTPL", "siNC\nLPS", "siNC\nLPS+TPL",
                      "siS100A8\nControl", "siS100A8\nTPL", "siS100A8\nLPS", "siS100A8\nLPS+TPL")))
plot_kd8 <- function(endpoint_value, stem, ylab) {
  d <- kd8 %>% filter(endpoint == endpoint_value)
  p <- ggplot(d, aes(display, value, colour = siRNA, fill = siRNA)) +
    stat_summary(fun = mean, geom = "col", width = 0.66, alpha = 0.65, colour = NA) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18,
                 linewidth = 0.4, colour = ink) +
    geom_point(position = position_jitter(width = 0.07), size = 1.5) +
    scale_colour_manual(values = c("siNC" = blue, "siS100A8" = gold)) +
    scale_fill_manual(values = c("siNC" = blue, "siS100A8" = gold)) +
    labs(x = NULL, y = ylab) + theme_plos() +
    theme(legend.position = "top", axis.text.x = element_text(angle = 42, hjust = 1, size = 6.8))
  save_plot_panel(p, stem, 1000, 740)
}
kd8_pnf <- plot_kd8("p-NF-kB / total NF-kB", "Fig8_pNFkB", "p-NF-\u03baB / total NF-\u03baB")
kd8_s100 <- plot_kd8("S100A8 / beta-actin", "Fig8_S100A8_8group", "S100A8 / \u03b2-actin")

kd4 <- exp_dat %>% filter(panel %in% c("Q", "R", "S")) %>%
  mutate(siRNA = factor(siRNA, levels = c("siNC", "siS100A8")),
         treatment = factor(treatment, levels = c("Control", "LPS")),
         display = factor(paste(siRNA, treatment, sep = "\n"),
           levels = c("siNC\nControl", "siNC\nLPS", "siS100A8\nControl", "siS100A8\nLPS")))
plot_kd4 <- function(endpoint_value, stem, ylab) {
  d <- kd4 %>% filter(endpoint == endpoint_value)
  p <- ggplot(d, aes(display, value, colour = siRNA, fill = siRNA)) +
    stat_summary(fun = mean, geom = "col", width = 0.64, alpha = 0.65, colour = NA) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.18,
                 linewidth = 0.4, colour = ink) +
    geom_point(position = position_jitter(width = 0.065), size = 1.6) +
    scale_colour_manual(values = c("siNC" = blue, "siS100A8" = gold)) +
    scale_fill_manual(values = c("siNC" = blue, "siS100A8" = gold)) +
    labs(x = NULL, y = ylab) + theme_plos() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1, size = 7))
  save_plot_panel(p, stem, 760, 650)
}
kd4_cd74 <- plot_kd4("CD74 / beta-actin", "Fig8_CD74", "CD74 / \u03b2-actin")
kd4_hla <- plot_kd4("HLA-DRalpha / beta-actin", "Fig8_HLADRA", "HLA-DR\u03b1 / \u03b2-actin")
kd4_s100 <- plot_kd4("S100A8 / beta-actin", "Fig8_S100A8_4group", "S100A8 / \u03b2-actin")

# Effect-size summary: TPL rescue within LPS, separately in siNC and siS100A8.
paired_effects <- kd8 %>% filter(treatment %in% c("LPS", "LPS+TPL")) %>%
  select(endpoint, siRNA, replicate, treatment, value) %>%
  pivot_wider(names_from = treatment, values_from = value) %>%
  mutate(effect = `LPS+TPL` - LPS)
within_effects <- paired_effects %>% group_by(endpoint, siRNA) %>%
  summarise(estimate = mean(effect), se = sd(effect) / sqrt(n()), df = n() - 1,
            lower = estimate - qt(0.975, df) * se,
            upper = estimate + qt(0.975, df) * se, .groups = "drop") %>%
  mutate(contrast = as.character(siRNA),
         label = paste(endpoint, contrast, sep = "  |  "))
dependency_effects <- paired_effects %>% select(endpoint, replicate, siRNA, effect) %>%
  pivot_wider(names_from = siRNA, values_from = effect) %>%
  mutate(effect = siS100A8 - siNC) %>% group_by(endpoint) %>%
  summarise(estimate = mean(effect), se = sd(effect) / sqrt(n()), df = n() - 1,
            lower = estimate - qt(0.975, df) * se,
            upper = estimate + qt(0.975, df) * se, .groups = "drop") %>%
  mutate(siRNA = factor("Dependency contrast", levels = c("siNC", "siS100A8", "Dependency contrast")),
         contrast = "Dependency contrast",
         label = paste(endpoint, contrast, sep = "  |  "))
interaction_effects <- bind_rows(within_effects, dependency_effects) %>%
  mutate(siRNA = factor(as.character(siRNA),
                        levels = c("siNC", "siS100A8", "Dependency contrast")),
         endpoint_short = ifelse(grepl("p-NF", endpoint),
                                 "p-NF-\u03baB/NF-\u03baB", "S100A8/\u03b2-actin"),
         label = paste(endpoint_short, contrast, sep = "  |  "))
write.csv(interaction_effects, file.path(out_root, "interaction_effects.csv"), row.names = FALSE)
p_int <- ggplot(interaction_effects, aes(estimate, label, colour = siRNA)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.18, linewidth = 0.5) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c("siNC" = blue, "siS100A8" = gold,
                                 "Dependency contrast" = "#6B6B6B")) +
  labs(x = "TPL effect under LPS (mean difference)\n95% CI", y = NULL) +
  theme_plos() + theme(legend.position = "none", axis.text.y = element_text(size = 7))
int_path <- save_plot_panel(p_int, "Fig8_interaction_effect", 1600, 560)

# ---------------------------------------------------------------------------
# Quantitative computational bridge panels.
# ---------------------------------------------------------------------------
sim <- read.csv(file.path(dock_root, "source_data", "compound_similarity.csv"))
names(sim)[1] <- "metric"
morgan <- as.numeric(sim$value[sim$metric == "morgan_radius2_2048_tanimoto_chiral"])
maccs <- as.numeric(sim$value[sim$metric == "maccs_tanimoto"])
structure_path <- file.path(dock_root, "assets", "compound_structures_horizontal.png")
compound_path <- structure_path

contact <- read.csv(file.path(dock_root, "source_data", "representative_pose_cross_compound_comparison.csv"))
contact$label <- paste(contact$receptor, contact$scoring, sep = " / ")
p_contact <- ggplot(contact, aes(label, contact_jaccard, fill = receptor)) +
  geom_col(width = 0.65, colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", contact_jaccard)), vjust = -0.35, size = 2.7) +
  scale_fill_manual(values = c("AC" = blue, "BD" = teal)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = NULL, y = "Contact-set Jaccard similarity") + theme_plos() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
contact_path <- save_plot_panel(p_contact, "Fig6_contact_similarity", 800, 620)

lig <- read.csv(file.path(dock_root, "source_data", "ligand_pocket.csv"))
names(lig)[1] <- "replicate"
lig$replicate <- factor(lig$replicate, levels = c(1, 2, 3), labels = c("Rep 1", "Rep 2", "Rep 3"))
p_lig <- ggplot(lig, aes(time_ns, value, colour = replicate)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_line(linewidth = 0.42, alpha = 0.9) +
  scale_colour_manual(values = c("Rep 1" = blue, "Rep 2" = teal, "Rep 3" = red)) +
  labs(x = "Time (ns)", y = "Pocket displacement (nm)") + theme_plos() +
  theme(legend.position = "top")
lig_path <- save_plot_panel(p_lig, "Fig6_ligand_displacement", 1100, 620)

# ---------------------------------------------------------------------------
# Original single-panel source map. No panels are cropped from composites.
# ---------------------------------------------------------------------------
P <- function(...) file.path(...)
src <- list(
  f1a = P(project_root, "sepsis-s100a8-multiomics", "data", "02_DEG", "output", "Disease_vs_Normal_volcano_padj_S100A8_Fixed.pdf"),
  f1b = P(project_root, "sepsis-s100a8-multiomics", "data", "03_WGCNA_analysis", "output", "module_group_correlation_heatmap_all.pdf"),
  f1c = P(project_root, "sepsis-s100a8-multiomics", "data", "03_WGCNA_analysis", "output", "module_black_boxplot_half_border.pdf"),
  f1d = P(project_root, "sepsis-s100a8-multiomics", "data", "04_Extract overlapping genes", "output", "Venn_red_blue_horizontal.pdf"),
  f1e = P(project_root, "sepsis-s100a8-multiomics", "data", "05_Enrichment Analysis", "output", "Figure_A_GO.pdf"),
  f1f = P(project_root, "sepsis-s100a8-multiomics", "data", "05_Enrichment Analysis", "output", "Figure_B_KEGG.pdf"),
  f2a = P(fig_root, "Fig_ClassBalance_StudyDesign.svg"),
  f2b = P(fig_root, "Fig_ModelPerformance_BalanceStrategies.svg"),
  f2c = P(fig_root, "Fig_S100A8_SelectionFrequency.svg"),
  f2d = P(fig_root, "Fig_S100A8_RankDistribution.svg"),
  f2e = P(fig_root, "Fig_AUC_Distribution_Downsampling.svg"),
  f2f = P(fig_root, "Fig_Bootstrap_EffectSize.svg"),
  f3a = P(project_root, "sepsis-s100a8-multiomics", "data", "08_External Validation", "output", "Combined_ROC_Validation_Final.pdf"),
  f3b = P(fig_root, "Fig_ClinicalComparator_EffectForest.svg"),
  f3c = P(fig_root, "Fig_ClinicalComparator_Expression.svg"),
  f3d = P(fig_root, "Fig_ClinicalComparator_ROC.svg"),
  f3e = P(fig_root, "Fig_DecisionCurve.svg"),
  f3f = P(fig_root, "Fig_S100A8_Prognosis_Extended.svg"),
  f4a = P(sc_root, "Plot_03_Cleaned_UMAP.pdf"),
  f4b = P(sc_root, "Plot_05_S100A8_FeaturePlot.pdf"),
  f4c = P(sc_root, "Plot_02_Annotation_DotPlot.pdf"),
  f4d = P(fig_root, "Fig_DonorLevel_TFActivity.svg"),
  f4e = P(sc_root, "Figure_4I_KO_Distance_Scatter.pdf"),
  f4f = P(sc_root, "Figure_4J_KO_Distance_Barplot.pdf"),
  f4g = P(sc_root, "Figure_4L_KO_KEGG_Enrichment.pdf"),
  f4h = P(sc_root, "Figure_4K_KO_GO_Enrichment.pdf"),
  f5a = P(fig_root, "Fig_Monocyte_UnsupervisedClusters.svg"),
  f5b = P(fig_root, "Fig_Monocyte_StatePrograms.svg"),
  f5c = P(fig_root, "Fig_SCP548_DonorPseudobulk_S100A8.svg"),
  f5d = P(fig_root, "Fig_SCP548_MS1_Abundance.svg"),
  f5e = P(fig_root, "Fig_SCP548_AntigenPresentation.svg"),
  f5f = P(fig_root, "Fig_GSE205672_S100A8_Monocytes.svg"),
  f5g = P(fig_root, "Fig_GSE205672_AntigenPresentation.svg"),
  f5h = P(fig_root, "Fig_Mechanistic_Triangulation.svg"),
  f6a = P(project_root, "\u836f\u7269\u7b5b\u9009", "Figure9A_S100A8_Heatmap.pdf"),
  f6b = compound_path,
  f6c = P(dock_root, "assets", "docking_overlay_overview_v3.png"),
  f6d = P(dock_root, "assets", "docking_overlay_closeup_v3.png"),
  f6e = P(fig_root, "Fig_TPL_16OH_MatchedDocking.svg"),
  f6f = contact_path,
  f6g = P(md_root, "RMSD_latest.pdf"),
  f6h = lig_path,
  f7a = qpcr_path,
  f7b = P(exp_root, "\u56fe\u516d WB.pdf"),
  f7c = wb_old_path,
  f7d = P(exp_root, "S100A8_CETSA.pdf"),
  f7e = P(exp_root, "CETSA.png"),
  f8a = P(exp_add, "\u8865WB 1.pdf"),
  f8b = kd8_pnf,
  f8c = kd8_s100,
  f8d = P(exp_add, "\u8865WB 2.pdf"),
  f8e = kd4_cd74,
  f8f = kd4_hla,
  f8g = kd4_s100,
  f8h = int_path
)
missing <- names(src)[!file.exists(unlist(src))]
if (length(missing)) stop("Missing mapped panels: ", paste(missing, collapse = ", "))

row <- function(keys, height, weights = NULL) list(paths = unlist(src[keys]), height = height, weights = weights)

figures <- list(
  # Fig1 is assembled title-free from the six original vector panels by
  # scripts/46_reassemble_fig1_titleless.py.
  Fig2 = assemble_rows(list(
    row("f2a", 580), row(c("f2b", "f2c"), 640),
    row(c("f2d", "f2e"), 650), row("f2f", 580))),
  # Fig3 is assembled from six independent Python vector panels by
  # scripts/45_assemble_fig3.py to prevent nested mini-panel composites.
  Fig4 = assemble_rows(list(
    row(c("f4a", "f4b"), 590), row(c("f4c", "f4d"), 560, c(0.9, 1.1)),
    row(c("f4e", "f4f"), 560), row(c("f4g", "f4h"), 560))),
  Fig5 = assemble_rows(list(
    row(c("f5a", "f5b"), 690), row(c("f5c", "f5d", "f5e"), 690),
    row(c("f5f", "f5g"), 690), row("f5h", 420))),
  Fig6 = assemble_rows(list(
    row(c("f6a", "f6b"), 620, c(0.9, 1.1)), row(c("f6c", "f6d"), 650),
    row(c("f6e", "f6f"), 650, c(1.2, 0.8)), row(c("f6g", "f6h"), 620))),
  # Fig7 is assembled by scripts/38_assemble_fig7_fig8.py after its analyte
  # names are placed as axis labels rather than plot titles.
  Fig8 = assemble_rows(list(
    row(c("f8a", "f8b", "f8c"), 650, c(0.9, 1.05, 1.05)),
    row(c("f8d", "f8e", "f8f", "f8g"), 650, c(0.95, 1, 1, 1)),
    row("f8h", 520)))
)

manifest <- do.call(rbind, Map(export_figure, figures, names(figures)))

# Supporting figures: secondary diagnostics, full MD pages and viability.
crop_outer_whitespace <- function(path, stem, top = 0, bottom = 0) {
  img <- read_any(path)
  info <- image_info(img)
  y0 <- floor(info$height * top)
  hh <- floor(info$height * (1 - top - bottom))
  img <- image_crop(img, sprintf("%dx%d+0+%d", info$width, hh, y0), repage = TRUE)
  dst <- file.path(work_root, paste0(stem, ".png"))
  image_write(img, dst, format = "png")
  dst
}
supp_src <- list(
  S1 = c(
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "2_LASSO_CV_Classic_Labeled_Text.pdf"),
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "1_LASSO_Path_Classic_Labeled_Text.pdf"),
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "1_LASSO_Coefficients_Barplot.pdf"),
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "3C_GLMBOOST_Final_Coefficients.pdf"),
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "4A_XGBoost_Importance_Top5.pdf"),
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "3A_RF_Importance_Top5.pdf"),
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "8B_SHAP_BeeSwarm_Plot.pdf"),
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "8D_SHAP_Sample_Waterfall.pdf"),
    P(project_root, "sepsis-s100a8-multiomics", "data", "07_Multi-Model Feature Selection", "output", "8C_SHAP_Feature_Dependence.pdf")),
  S2 = c(P(fig_root, "Fig_ClinicalComparator_PR.svg"), P(fig_root, "Fig_S100A8_APACHEII.svg"),
         P(fig_root, "Fig_S100A8_Mortality.svg"), P(fig_root, "Fig_S100A8_Longitudinal_Trajectory.svg"),
         P(fig_root, "Fig_S100A8_TreatmentResponse.svg"), P(fig_root, "Fig_Prognosis_Forest.svg")),
  S3 = c(P(fig_root, "Fig_Donor_UMAP.svg"), P(fig_root, "Fig_Donor_Composition.svg"),
         P(fig_root, "Fig_Donor_CellNumber.svg"), P(fig_root, "Fig_Donor_StateAbundance.svg"),
         P(fig_root, "Fig_Monocyte_Clustree.svg"), P(fig_root, "Fig_Monocyte_UMAP_ByDonor.svg")),
  S4 = c(P(fig_root, "Fig_Pseudobulk_S100A8_Donor.svg"), P(fig_root, "Fig_IndependentSignature_DonorScores.svg"),
         P(fig_root, "Fig_Bulk_IndependentSignature.svg"), P(fig_root, "Fig_SCP548_S100A8_HLAII_Correlation.svg"),
         P(fig_root, "Fig_GSE205672_S100A8_HLAII_Correlation.svg"), P(fig_root, "Fig_GSE205672_NFkBScore.svg")),
  S5 = P(md_root, "Figure_1_MD_system_metrics.pdf"),
  S6 = P(md_root, "Figure_2_ligand_behavior.pdf"),
  S7 = P(md_root, "Figure_3_PCA_FEL.pdf"),
  S8 = P(exp_root, "\u56fe\u4e00 CCK-8.pdf")
)
if (any(!file.exists(unlist(supp_src)))) stop("One or more supporting-figure sources are missing")

support_figures <- list(
  S1_fig = assemble_rows(list(list(paths = supp_src$S1[1:3], height = 760),
                              list(paths = supp_src$S1[4:6], height = 760),
                              list(paths = supp_src$S1[7:9], height = 760))),
  S2_fig = assemble_rows(list(list(paths = supp_src$S2[1:2], height = 740),
                              list(paths = supp_src$S2[3:4], height = 740),
                              list(paths = supp_src$S2[5:6], height = 740))),
  S3_fig = assemble_rows(list(list(paths = supp_src$S3[1:2], height = 740),
                              list(paths = supp_src$S3[3:4], height = 740),
                              list(paths = supp_src$S3[5:6], height = 740))),
  S4_fig = assemble_rows(list(list(paths = supp_src$S4[1:2], height = 740),
                              list(paths = supp_src$S4[3:4], height = 740),
                              list(paths = supp_src$S4[5:6], height = 740))),
  S5_fig = assemble_rows(list(list(paths = supp_src$S5, height = 2100))),
  S6_fig = assemble_rows(list(list(paths = supp_src$S6, height = 1700))),
  S7_fig = assemble_rows(list(list(paths = supp_src$S7, height = 1700))),
  S8_fig = assemble_rows(list(list(paths = supp_src$S8, height = 1600)))
)
manifest <- rbind(manifest, do.call(rbind, Map(export_figure, support_figures, names(support_figures))))

raw_src <- file.path(exp_ready, "S1_raw_images.pdf")
raw_dst <- file.path(out_root, "S1_raw_images.pdf")
file.copy(raw_src, raw_dst, overwrite = TRUE)

manifest$plos_width_ok <- manifest$width_px >= 789 & manifest$width_px <= 2250
manifest$plos_height_ok <- manifest$height_px <= 2625
manifest$under_10mb <- manifest$bytes <= 10 * 1024^2
write.csv(manifest, file.path(out_root, "figure_export_manifest.csv"), row.names = FALSE)

qa <- manifest %>% filter(format == "tif") %>%
  select(figure, width_px, height_px, dpi, bytes, plos_width_ok, plos_height_ok, under_10mb) %>%
  mutate(size_mb = round(bytes / 1024^2, 2),
         status = ifelse(plos_width_ok & plos_height_ok & under_10mb, "PASS", "CHECK"))
write.csv(qa, file.path(out_root, "PLOS_figure_QA.csv"), row.names = FALSE)

captions <- c(
  "Fig 1. Transcriptomic discovery and network-level prioritization of S100A8 in sepsis. (A) Differential-expression volcano plot. (B) WGCNA module-trait correlations. (C) Disease-associated module eigengene distribution. (D) Overlap between disease-associated module genes and differentially expressed genes. (E,F) GO and KEGG enrichment of overlapping genes.",
  "Fig 2. Class-balance-aware and uncertainty-aware assessment of S100A8 prioritization. (A) Analysis design. (B) Model performance under weighted, unweighted and downsampled strategies. (C) S100A8 selection frequency across resampling strategies. (D) Rank distribution. (E) AUC distribution in balanced resamples. (F) Stratified-bootstrap distributions for AUC and standardized mean difference.",
  "Fig 3. External validation, clinically relevant discrimination boundaries and prognostic context of S100A8. (A) ROC curves in four independent healthy-versus-sepsis validation cohorts. (B) Standardized mean differences with 95% confidence intervals across healthy and clinically similar comparators; total sample sizes are shown at right. (C) Patient-level S100A8 distributions standardized within each comparison. Faint points represent individual samples, thick segments represent interquartile ranges and larger points represent medians. (D) Comparator-specific AUC estimates with 95% confidence intervals; the shaded region indicates weak-to-moderate discrimination and the dashed line denotes chance performance. (E) Net-benefit advantage of cross-validated S100A8 over the better of the treat-all and treat-none strategies across threshold probabilities; grey cells denote thresholds not evaluated in the source data. (F) Odds ratios per 1-SD higher S100A8 for mortality and septic-shock outcomes. Points are estimates, horizontal lines are 95% confidence intervals and exact P values are shown. Intervals crossing the null indicate uncertainty or absence of a clearly supported association.",
  "Fig 4. Single-cell localization of S100A8 and donor-level regulatory context in sepsis. (A) Major PBMC cell classes. (B) S100A8 expression on the UMAP. (C) Horizontal canonical-marker dot plot; dot size denotes the percentage of expressing cells and colour denotes scaled average expression. (D) Donor-level transcription-factor activity differences. Points are effect estimates and horizontal lines are 95% confidence intervals.",
  "Fig 5. Donor-aware and independent single-cell validation of the S100A8-associated myeloid state. (A) Unsupervised monocyte clusters. (B) State-program scores. (C-E) SCP548 donor-level S100A8, MS1-state abundance and antigen-presentation program. (F,G) Independent GSE205672 validation of S100A8 and antigen-presentation changes. (H) Cross-dataset mechanistic triangulation.",
  "Fig 6. Rationale for transition from computational screening to triptolide and structural assessment. (A) Candidate screen. (B) Structural similarity between triptolide and 16-hydroxytriptolide; similarity is used only as a screening bridge. (C,D) Docking-pose overview and binding-site close-up. (E) Matched docking-score differences. (F) Contact-set similarity. (G) Protein backbone RMSD. (H) Ligand pocket displacement across all three simulations, showing replicate heterogeneity.",
  "Fig 7. Cellular effects of triptolide and support for S100A8 engagement. (A) Relative IL-6, TNF-\u03b1, IL-1\u03b2 and S100A8 mRNA expression. (B) Representative immunoblots. (C) Quantification of S100A8, NF-\u03baB and p-NF-\u03baB. (D,E) CETSA temperature-response curve and representative blots. Points show available biological replicates and summary marks show mean \u00b1 SEM. Exact adjusted P values are shown for the two prespecified comparisons, Control versus LPS and LPS versus LPS+TPL; multiplicity was controlled using Tukey's procedure following one-way ANOVA.",
  "Fig 8. Experimental perturbation of S100A8 and boundaries of triptolide dependency. (A) Immunoblots across siNC/siS100A8 and Control/TPL/LPS/LPS+TPL conditions. (B,C) p-NF-\u03baB/total NF-\u03baB and S100A8/\u03b2-actin quantification. (D) Antigen-presentation immunoblots after S100A8 knockdown. (E-G) CD74, HLA-DR\u03b1 and S100A8 quantification. (H) Estimated TPL effects under LPS and the between-siRNA dependency contrast. Points show three independent repeats and summary marks show mean \u00b1 SEM. Exact Holm-adjusted P values are shown for prespecified contrasts; effect estimates in H are presented with 95% confidence intervals, and intervals crossing zero indicate an unresolved dependency effect.",
  "Supporting figures S1-S4 contain feature-selection/SHAP diagnostics, extended clinical analyses, donor/single-cell quality-control analyses and non-circular sensitivity analyses, respectively. Supporting figures S5-S7 contain the complete molecular-dynamics system, ligand-behavior and PCA/free-energy-landscape results. Supporting figure S8 shows the CCK-8 dose-response experiment using the deposited replicate-level values; error bars are shown only where replicate variance is estimable, and the LPS 500 ng/mL series is descriptive because only one value per triptolide dose was available in the deposited source table. S1_raw_images.pdf contains the mapped uncropped blot/gel source images."
)
writeLines(captions, file.path(out_root, "Figure_captions_draft.txt"), useBytes = TRUE)
message("Exported final Fig1-Fig8, S1-S8, source-data effect summary and S1_raw_images.pdf to: ", out_root)
