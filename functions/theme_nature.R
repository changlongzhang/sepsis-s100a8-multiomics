# Nature/PLOS兼容的统一R绘图主题与导出函数（R为唯一绘图后端）
palette_nature <- c(
  Normal = "#3C6E8F", Sepsis = "#C65A46", Disease = "#C65A46",
  weighted = "#315F7D", unweighted = "#8B8B8B", downsample = "#4F8A76",
  bootstrap = "#8A6F91", accent = "#D39B44", dark = "#252525", light = "#D8D8D8",
  pale_blue = "#DCE8EE", pale_red = "#F2DED8", pale_teal = "#DDEAE5",
  pale_gold = "#F3E8D2", reference = "#5F6368"
)

theme_nature <- function(base_size = 7.5, base_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.35, colour = "#252525"),
      axis.ticks = ggplot2::element_line(linewidth = 0.35, colour = "#252525"),
      axis.ticks.length = grid::unit(1.4, "mm"),
      axis.title = ggplot2::element_text(size = base_size),
      axis.text = ggplot2::element_text(size = base_size - 0.5, colour = "black"),
      legend.title = ggplot2::element_text(size = base_size - 0.5),
      legend.text = ggplot2::element_text(size = base_size - 1),
      strip.text = ggplot2::element_text(size = base_size - 0.3, face = "bold"),
      strip.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.tag = ggplot2::element_text(size = 8, face = "bold"),
      panel.grid = ggplot2::element_blank(), legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 6, 5, 5)
    )
}

theme_forest <- function(base_size = 7.5, base_family = "Arial") {
  theme_nature(base_size, base_family) +
    ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(colour = "#EFEFEF", linewidth = 0.25))
}

theme_curve <- function(base_size = 7.5, base_family = "Arial") {
  theme_nature(base_size, base_family) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(colour = "#ECECEC", linewidth = 0.25),
                   panel.grid.minor = ggplot2::element_blank(),
                   aspect.ratio = 1)
}

theme_umap <- function(base_size = 7.5, base_family = "Arial") {
  theme_nature(base_size, base_family) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = base_size - 0.2))
}

theme_heatmap <- function(base_size = 7.5, base_family = "Arial") {
  theme_nature(base_size, base_family) +
    ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank())
}

save_nature_plot <- function(plot, stem, width_mm = 89, height_mm = 70, dpi = 600,
                             source_data = NULL) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  w <- width_mm / 25.4; h <- height_mm / 25.4
  svglite::svglite(paste0(stem, ".svg"), width = w, height = h)
  print(plot); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = w, height = h, family = "Arial")
  print(plot); grDevices::dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = w, height = h, units = "in", res = dpi,
                 compression = "lzw")
  print(plot); grDevices::dev.off()
  if (!is.null(source_data)) write.csv(source_data, paste0(stem, "_source_data.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  files <- paste0(stem, c(".svg", ".pdf", ".tiff"))
  if (!all(file.exists(files)) || any(file.info(files)$size == 0)) stop("图形导出失败: ", stem, call. = FALSE)
  invisible(files)
}
