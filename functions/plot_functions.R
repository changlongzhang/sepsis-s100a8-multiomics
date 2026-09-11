# 绘图辅助函数
plot_effect_ci <- function(df, x, estimate, lower, upper, colour = NULL) {
  aes_args <- list(x = rlang::.data[[x]], y = rlang::.data[[estimate]],
                   ymin = rlang::.data[[lower]], ymax = rlang::.data[[upper]])
  if (!is.null(colour)) aes_args$colour <- rlang::.data[[colour]]
  ggplot2::ggplot(df, do.call(ggplot2::aes, aes_args)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey70", linewidth = 0.4) +
    ggplot2::geom_pointrange(linewidth = 0.45, fatten = 1.5) +
    theme_nature()
}

