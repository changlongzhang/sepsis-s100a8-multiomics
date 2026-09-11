# 输入：已注释单核细胞状态及全 PBMC donor 细胞数。
# 输出：每位供者的状态比例、精确置换检验、效应量和95% bootstrap区间。
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")
source(file.path("functions", "scrna_functions.R"), encoding = "UTF-8")

perm_p <- function(y, group) {
  obs <- mean(y[group == "Sepsis"]) - mean(y[group == "Control"])
  n_ctrl <- sum(group == "Control"); cmb <- combn(seq_along(y), n_ctrl)
  null <- apply(cmb, 2, function(ix) mean(y[-ix]) - mean(y[ix]))
  c(effect = obs, p = mean(abs(null) >= abs(obs) - 1e-12))
}
boot_ci <- function(y, group, B = 2000L) {
  vals <- replicate(B, {
    a <- sample(y[group == "Sepsis"], replace = TRUE); b <- sample(y[group == "Control"], replace = TRUE)
    mean(a) - mean(b)
  })
  stats::quantile(vals, c(.025, .975), na.rm = TRUE, names = FALSE)
}

with_script_log("19_donor_level_abundance", {
  suppressPackageStartupMessages({library(Seurat); library(ggplot2)})
  state_file <- file.path("06_monocyte_state", "GSE167363_monocyte_state_annotated.rds")
  donor_file <- file.path("05_single_cell_pseudobulk", "donor_metadata_clean.csv")
  assert_file(state_file); assert_file(donor_file)
  mono <- readRDS(state_file); dm <- read.csv(donor_file, check.names = FALSE)
  mono_tab <- as.data.frame(table(donor_id = mono$donor_id, state_label = mono$state_label), stringsAsFactors = FALSE)
  names(mono_tab)[3] <- "state_cells"
  mono_total <- aggregate(state_cells ~ donor_id, mono_tab, sum); names(mono_total)[2] <- "monocyte_cells"
  all_total <- aggregate(n_cells ~ donor_id, read.csv(file.path("05_single_cell_pseudobulk", "donor_cell_counts.csv")), sum)
  names(all_total)[2] <- "pbmc_cells"
  prop <- Reduce(function(x,y) merge(x,y,by="donor_id",all.x=TRUE), list(mono_tab, mono_total, all_total, dm[, c("donor_id","group")]))
  prop$proportion_monocytes <- prop$state_cells / prop$monocyte_cells
  prop$proportion_pbmc <- prop$state_cells / prop$pbmc_cells
  write_csv_utf8(prop, file.path("06_monocyte_state", "donor_state_proportions.csv"))

  res <- list()
  for (denom in c("proportion_monocytes", "proportion_pbmc")) for (st in unique(prop$state_label)) {
    d <- prop[prop$state_label == st, ]; pp <- perm_p(d[[denom]], d$group); ci <- boot_ci(d[[denom]], d$group)
    res[[paste(st, denom)]] <- data.frame(state_label = st, denominator = denom, n_control = sum(d$group == "Control"), n_sepsis = sum(d$group == "Sepsis"),
      effect_sepsis_minus_control = pp[["effect"]], ci_low = ci[1], ci_high = ci[2], permutation_p = pp[["p"]])
  }
  results <- do.call(rbind, res); results$fdr <- p.adjust(results$permutation_p, "BH")
  write_csv_utf8(results, file.path("06_monocyte_state", "Table_DonorLevel_StateAbundance.csv"))
  pd <- prop[prop$state_label != "Lymphoid-like contamination", ]
  eff<-results[results$denominator=="proportion_monocytes" & results$state_label!="Lymphoid-like contamination",]
  eff$state_label<-factor(eff$state_label,levels=eff$state_label[order(eff$effect_sepsis_minus_control)])
  p <- ggplot(eff,aes(effect_sepsis_minus_control,state_label))+
    geom_vline(xintercept=0,linetype=2,colour="#8C8C8C",linewidth=.35)+
    geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=0,linewidth=.55,colour="#5F6368")+
    geom_point(aes(fill=effect_sepsis_minus_control>0),shape=21,size=2.4,colour="white",stroke=.35)+
    scale_fill_manual(values=c(`TRUE`="#C65A46",`FALSE`="#315F7D"),guide="none")+
    scale_x_continuous(labels=scales::percent_format(accuracy=1))+
    labs(x="Sepsis − control difference\nin monocyte fraction (95% bootstrap CI)",y=NULL)+theme_forest()
  save_nature_plot(p, file.path("09_figures", "Fig_Donor_StateAbundance"), 120, 62, source_data = eff)
})
