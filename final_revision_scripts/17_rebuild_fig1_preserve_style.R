.libPaths(c("C:/Users/ZCL/AppData/Local/R/win-library/4.5", .libPaths()))
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(stringr)
  library(grid)
  library(svglite)
  library(ragg)
})

root <- Sys.getenv("SEPSIS_PROJECT_ROOT", unset = normalizePath("..", winslash = "/", mustWork = TRUE))
panel_dir <- file.path(root, "work", "final_revision_20260908", "figure_sources", "Fig1_preserved_style")
wdir <- file.path(root, "work", "final_revision_20260908", "wgcna_allgenes_rerun")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

font <- "Arial"
ink <- "#20252A"; mid <- "#7A858D"; red <- "#D66A5E"; blue <- "#4C78A8"; green <- "#59A14F"

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = font) +
    theme(
      text = element_text(colour = ink),
      axis.text = element_text(size = 8, colour = ink),
      axis.title = element_text(size = 8.5, colour = ink),
      axis.line = element_line(linewidth = 0.35, colour = ink),
      axis.ticks = element_line(linewidth = 0.30, colour = ink),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      plot.margin = margin(3, 3, 3, 4)
    )
}

# A: corrected limma volcano, unchanged content.
deg <- read.csv(file.path(wdir, "corrected_DEG_results_used.csv"), check.names = FALSE)
deg$Regulation <- deg$Regulation_final
pA <- ggplot(deg, aes(log2FoldChange, -log10(pmax(padj, 1e-300)))) +
  geom_point(aes(colour = Regulation), size = 1.05, alpha = 0.65) +
  geom_vline(xintercept = c(-1, 1), colour = mid, linewidth = 0.35, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), colour = mid, linewidth = 0.35, linetype = "dashed") +
  geom_point(data = filter(deg, Gene == "S100A8"), shape = 21, size = 2.7,
             stroke = 0.65, colour = "black", fill = NA) +
  geom_text_repel(data = filter(deg, Gene == "S100A8"), aes(label = Gene),
                  size = 8 / .pt, family = font, fontface = "bold",
                  nudge_x = -0.7, nudge_y = -8, box.padding = 0.5,
                  point.padding = 0.3, min.segment.length = 0,
                  segment.size = 0.35, seed = 123) +
  scale_colour_manual(values = c("Down" = blue, "Not significant" = "#C8CDD0", "Up" = red),
                      breaks = c("Up", "Down", "Not significant"), labels = c("Up", "Down", "NS")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Log2FC", y = expression(-log[10](FDR)), colour = NULL) +
  theme_pub() +
  theme(legend.position = c(0.82, 0.46), legend.background = element_blank(),
        legend.key.width = unit(3, "mm"))

# B: new WGCNA results in the established heatmap geometry and style.
assoc <- read.csv(file.path(wdir, "module_trait_all_results.csv"), check.names = FALSE)
assoc <- transform(assoc, Color = module, Correlation = correlation_sepsis,
                   Significance = ifelse(FDR < 0.001, "***", ifelse(FDR < 0.01, "**", ifelse(FDR < 0.05, "*", "NS"))))
assoc$Color <- factor(assoc$Color, levels = rev(assoc$Color))
bdat <- bind_rows(
  transmute(assoc, Color, Group = "Normal", Value = -Correlation, Significance),
  transmute(assoc, Color, Group = "Disease", Value = Correlation, Significance)
)
bdat$Group <- factor(bdat$Group, levels = c("Normal", "Disease"))
bdat$Label <- sprintf("%+.2f %s", bdat$Value, bdat$Significance)
pB <- ggplot(bdat, aes(Group, Color, fill = Value)) +
  geom_tile() +
  geom_text(aes(label = Label, colour = abs(Value) > 0.42), size = 8 / .pt, family = font) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = ink), guide = "none") +
  scale_fill_gradient2(low = blue, mid = "#F3F3F2", high = red, midpoint = 0,
                       limits = c(-0.65, 0.65), name = NULL) +
  guides(fill = guide_colorbar(barheight = unit(32, "mm"), barwidth = unit(3, "mm"))) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 8, base_family = font) +
  theme(axis.text.x = element_text(size = 8, angle = 35, hjust = 1, colour = ink),
        axis.text.y = element_text(size = 8, colour = ink), panel.grid = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = 8), plot.margin = margin(3, 4, 3, 4))

# C: new black-module eigengene values in the established style.
res <- readRDS(file.path(wdir, "WGCNA_final_results.rds"))
sqc <- read.csv(file.path(wdir, "sample_qc_all_802.csv"), check.names = FALSE)
mi <- match(rownames(res$MEs), sqc$sample_id)
stopifnot(!anyNA(mi))
cdat <- data.frame(Group = ifelse(sqc$group[mi] == "Sepsis", "Disease", "Normal"), ME = res$MEs$MEblack)
cdat$Group <- factor(cdat$Group, levels = c("Normal", "Disease"))
c_ymin <- min(cdat$ME); c_ymax <- max(cdat$ME)
c_bracket <- c_ymax + 0.010
c_top <- c_ymax + 0.021
pC <- ggplot(cdat, aes(Group, ME, fill = Group)) +
  geom_boxplot(width = 0.42, outlier.shape = NA, alpha = 0.80, colour = "black", linewidth = 0.45) +
  geom_point(position = position_jitter(width = 0.16, height = 0, seed = 123),
             size = 0.80, alpha = 0.52, colour = "black") +
  annotate("segment", x = 1, xend = 2, y = c_bracket, yend = c_bracket, linewidth = 0.35) +
  annotate("segment", x = 1, xend = 1, y = c_bracket - 0.006, yend = c_bracket, linewidth = 0.35) +
  annotate("segment", x = 2, xend = 2, y = c_bracket - 0.006, yend = c_bracket, linewidth = 0.35) +
  annotate("text", x = 1.5, y = c_top, label = "****", family = font, size = 8 / .pt) +
  scale_fill_manual(values = c("Normal" = "#377EB8", "Disease" = "#E41A1C")) +
  scale_y_continuous(limits = c(c_ymin - 0.008, c_top + 0.006), expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = "Module Eigengene Value") +
  theme_pub() +
  theme(legend.position = "none", axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8), axis.title.y = element_text(size = 8.5),
        plot.margin = margin(4, 10, 4, 12))

# D: new 265-gene overlap with established geometry, colors, and typography.
overlap <- readLines(file.path(wdir, "primary_sepsis_module_corrected_DEG_overlap.txt"), warn = FALSE)
black <- readLines(file.path(wdir, "primary_sepsis_module_black_genes.txt"), warn = FALSE)
sig <- filter(deg, Regulation != "Not significant")$Gene
only_black <- length(setdiff(black, sig)); only_deg <- length(setdiff(sig, black))
stopifnot(length(overlap) == 265, only_black == 448, only_deg == 1243)
theta <- seq(0, 2 * pi, length.out = 500)
venn_df <- bind_rows(
  data.frame(x = 0.40 + 0.27 * cos(theta), y = 0.50 + 0.31 * sin(theta), set = "WGCNA"),
  data.frame(x = 0.60 + 0.27 * cos(theta), y = 0.50 + 0.31 * sin(theta), set = "DEGs")
)
pD <- ggplot() +
  geom_polygon(data = venn_df, aes(x, y, group = set, fill = set), alpha = 0.80,
               colour = "black", linewidth = 0.65) +
  annotate("text", x = 0.23, y = 0.50, label = only_black, family = font, fontface = "bold", size = 8 / .pt) +
  annotate("text", x = 0.50, y = 0.50, label = length(overlap), family = font, fontface = "bold", size = 8 / .pt) +
  annotate("text", x = 0.79, y = 0.50, label = only_deg, family = font, fontface = "bold", size = 8 / .pt) +
  annotate("text", x = 0.24, y = 0.85, label = "WGCNA", family = font, size = 8 / .pt) +
  annotate("text", x = 0.76, y = 0.85, label = "DEGs", family = font, size = 8 / .pt) +
  scale_fill_manual(values = c("WGCNA" = "#D25E5E", "DEGs" = "#4E8AC8")) +
  coord_fixed(xlim = c(0.03, 0.97), ylim = c(0.08, 0.93), clip = "off") +
  theme_void(base_family = font) +
  theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))

# E/F: corrected enrichment results, unchanged content.
go_all <- read.csv(file.path(wdir, "overlap_GO_all_results.csv"), check.names = FALSE)
go_top <- go_all %>% group_by(ONTOLOGY) %>% arrange(pvalue, .by_group = TRUE) %>%
  slice_head(n = 5) %>% ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, levels = c("MF", "CC", "BP")),
         Description_plot = str_wrap(Description, width = 45)) %>%
  arrange(ONTOLOGY, Count) %>%
  mutate(Description_plot = factor(Description_plot, levels = unique(Description_plot)))
pE <- ggplot(go_top, aes(Count, Description_plot, fill = ONTOLOGY)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = Count), hjust = -0.25, size = 8 / .pt, family = font) +
  scale_fill_manual(values = c(BP = red, CC = blue, MF = green)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.17))) +
  labs(title = "GO Enrichment Analysis", x = "Gene count", y = NULL, fill = NULL) +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8, lineheight = 0.75), axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.25, colour = "#E7EBEE"),
        legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.background = element_blank(), legend.key.height = unit(2.6, "mm"),
        legend.spacing.y = unit(0, "mm"),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        plot.margin = margin(4, 4, 4, 4))

kegg_all <- read.csv(file.path(wdir, "overlap_KEGG_all_results.csv"), check.names = FALSE)
kegg_top <- kegg_all %>% arrange(pvalue) %>% slice_head(n = 5) %>%
  mutate(Category = case_when(
    grepl("Complement|Fc epsilon|TNF|IL-17|immune|phagocytosis", Description, ignore.case = TRUE) ~ "Immune",
    grepl("metabolism|biosynthesis|steroidogenesis", Description, ignore.case = TRUE) ~ "Metabolism",
    TRUE ~ "Related"),
    Description_plot = str_wrap(Description, width = 22)) %>%
  arrange(Category, Count) %>%
  mutate(Description_plot = factor(Description_plot, levels = unique(Description_plot)))
pF <- ggplot(kegg_top, aes(Count, Description_plot, fill = Category)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = Count), hjust = -0.25, size = 8 / .pt, family = font) +
  scale_fill_manual(values = c(Immune = red, Metabolism = green, Related = blue), drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(title = "KEGG Enrichment Analysis", x = "Gene count", y = NULL, fill = NULL,
       caption = "No pathway passed FDR < 0.05") +
  theme_pub() +
  theme(axis.text.y = element_text(size = 8), axis.ticks.y = element_blank(),
        legend.position = "none", panel.grid.major.x = element_line(linewidth = 0.25, colour = "#E7EBEE"),
        plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        plot.caption = element_text(size = 8, hjust = 0, colour = mid),
        plot.margin = margin(4, 4, 3, 4))

# Preserve the established 540 x 552 pt geometry exactly.
page_w <- 540; page_h <- 552
rects <- list(
  c(18, 14, 264, 176), c(280, 14, 526, 176),
  c(18, 192, 264, 350), c(280, 192, 526, 350),
  c(18, 366, 318, 540), c(330, 366, 526, 540)
)
plots <- list(pA, pB, pC, pD, pE, pF)

draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = "white", col = NA))
  for (i in seq_along(plots)) {
    r <- rects[[i]]
    pushViewport(viewport(
      x = unit((r[1] + r[3]) / 2, "pt"),
      y = unit(page_h - (r[2] + r[4]) / 2, "pt"),
      width = unit(r[3] - r[1], "pt"), height = unit(r[4] - r[2], "pt")
    ))
    grid.draw(ggplotGrob(plots[[i]]))
    popViewport()
    grid.text(LETTERS[i], x = unit(r[1] - 11, "pt"), y = unit(page_h - r[2] - 10, "pt"),
              just = c("left", "centre"), gp = gpar(fontfamily = font, fontface = "bold", fontsize = 12, col = ink))
  }
}

save_single <- function(plot, stem, width, height) {
  cairo_pdf(file.path(panel_dir, paste0(stem, ".pdf")), width = width, height = height, family = font)
  print(plot); dev.off()
  svglite(file.path(panel_dir, paste0(stem, ".svg")), width = width, height = height, bg = "white")
  print(plot); dev.off()
}
save_single(pC, "Fig1C_boxplot_font_standardized", 3.42, 2.19)
save_single(pD, "Fig1D_overlap_265_font_standardized", 3.42, 2.19)

candidate_pdf <- file.path(panel_dir, "Fig1_final_preserved_style.pdf")
candidate_svg <- file.path(panel_dir, "Fig1_final_preserved_style.svg")
qa_ascii <- "C:/Users/ZCL/final_revision_outputs"
dir.create(qa_ascii, recursive = TRUE, showWarnings = FALSE)
candidate_tif <- file.path(qa_ascii, "Fig1_preserved_style.tif")
cairo_pdf(candidate_pdf, width = page_w / 72, height = page_h / 72, family = font)
draw_figure(); dev.off()
svglite(candidate_svg, width = page_w / 72, height = page_h / 72, bg = "white")
draw_figure(); dev.off()
agg_tiff(candidate_tif, width = page_w / 72, height = page_h / 72, units = "in",
         res = 300, compression = "lzw", background = "white")
draw_figure(); dev.off()

cat("WILCOX_P\t", format(wilcox.test(ME ~ Group, data = cdat)$p.value, scientific = TRUE), "\n", sep = "")
cat("CANDIDATE_PDF\t", candidate_pdf, "\n", sep = "")
cat("CANDIDATE_SVG\t", candidate_svg, "\n", sep = "")
cat("CANDIDATE_TIF\t", candidate_tif, "\n", sep = "")
