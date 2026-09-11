# Purpose: extend the revision with direct SIRS-vs-sepsis discrimination,
# longitudinal severity, mortality/shock and treatment-response analyses.
# Statistical unit is the patient; repeated time points are never treated as
# independent observations. Figures are exported individually (no composites).
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")

read_eset_local <- function(gse) {
  suppressPackageStartupMessages(library(GEOquery))
  f <- list.files(file.path("01_data_inventory", "downloaded_data"),
                  pattern = paste0("^", gse, ".*series_matrix\\.txt\\.gz$"),
                  full.names = TRUE)
  if (!length(f)) stop("Missing series matrix for ", gse)
  GEOquery::getGEO(filename = f[[1]], getGPL = FALSE)
}

extract_array_gene <- function(es, gene, annotation_pkg) {
  suppressPackageStartupMessages({library(Biobase); library(AnnotationDbi)})
  suppressPackageStartupMessages(library(annotation_pkg, character.only = TRUE))
  db <- get(annotation_pkg)
  mp <- AnnotationDbi::select(db, keys = gene, keytype = "SYMBOL",
                              columns = c("PROBEID", "SYMBOL"))
  probes <- intersect(unique(na.omit(mp$PROBEID)), rownames(Biobase::exprs(es)))
  if (!length(probes)) stop("No probe for ", gene, " in ", annotation_pkg)
  x <- Biobase::exprs(es)[probes, , drop = FALSE]
  if (max(x, na.rm = TRUE) > 50) x <- log2(x + 1)
  colMeans(x, na.rm = TRUE)
}

z_or_model <- function(d, outcome, dataset, endpoint) {
  d <- d[complete.cases(d[, c("expression", outcome)]), , drop = FALSE]
  d$expression_z <- as.numeric(scale(d$expression))
  if (nrow(d) < 15 || length(unique(d[[outcome]])) != 2) return(NULL)
  fit <- glm(stats::as.formula(paste(outcome, "~ expression_z")), d,
             family = binomial())
  co <- summary(fit)$coefficients["expression_z", ]
  data.frame(dataset = dataset, endpoint = endpoint, n = nrow(d),
             events = sum(d[[outcome]] == 1), effect_type = "OR per 1 SD",
             estimate = exp(co[["Estimate"]]),
             ci_low = exp(co[["Estimate"]] - 1.96 * co[["Std. Error"]]),
             ci_high = exp(co[["Estimate"]] + 1.96 * co[["Std. Error"]]),
             p_value = co[["Pr(>|z|)"]], log_effect = co[["Estimate"]],
             se = co[["Std. Error"]])
}

with_script_log("27_clinical_extension", {
  suppressPackageStartupMessages({library(Biobase); library(ggplot2); library(openxlsx); library(edgeR)})
  out_samples <- list(); effects <- list()

  # GSE63042: direct sterile SIRS versus sepsis comparison in the same cohort.
  es63042 <- read_eset_local("GSE63042")
  pd63042 <- pData(es63042)
  gz <- file.path("01_data_inventory", "downloaded_data",
                  "GSE63042_capsod_seq_rel_RPM_060314.xlsx.gz")
  xlsx <- sub("\\.gz$", "", gz)
  if (!file.exists(xlsx)) {
    src <- gzfile(gz, "rb"); dst <- file(xlsx, "wb")
    repeat { b <- readBin(src, "raw", 1024^2); if (!length(b)) break; writeBin(b, dst) }
    close(src); close(dst)
  }
  rpm <- openxlsx::read.xlsx(xlsx)
  row <- which(toupper(as.character(rpm$gene_name)) == TARGET_GENE)
  if (length(row) != 1) stop("GSE63042 S100A8 row is not unique")
  ids <- as.character(pd63042$title)
  vals <- as.numeric(rpm[row, match(ids, colnames(rpm))])
  d63042 <- data.frame(dataset = "GSE63042", sample_id = rownames(pd63042),
                       patient_id = pd63042[["patient id:ch1"]],
                       group = pd63042[["sirs vs sepsis:ch1"]],
                       detailed_outcome = pd63042[["sirs outcomes:ch1"]],
                       expression = log2(vals + 1))
  d63042$sepsis <- as.integer(d63042$group == "Sepsis")
  d63042$death <- as.integer(grepl("death", d63042$detailed_outcome, ignore.case = TRUE))
  out_samples[["GSE63042"]] <- d63042
  effects[["GSE63042_sepsis_vs_SIRS"]] <- z_or_model(d63042, "sepsis", "GSE63042", "Sepsis vs sterile SIRS")
  sepsis63042 <- subset(d63042, group == "Sepsis")
  effects[["GSE63042_death"]] <- z_or_model(sepsis63042, "death", "GSE63042", "In-hospital sepsis death")

  p_mimic <- ggplot(d63042, aes(group, expression, colour = group)) +
    geom_boxplot(width = .46, outlier.shape = NA, colour = "#5F6368", fill = NA, linewidth = .4) +
    geom_point(size = 1.65, alpha = .82, position = position_jitter(width = .07, height = 0)) +
    scale_colour_manual(values = c(SIRS = palette_nature[["Normal"]], Sepsis = palette_nature[["Sepsis"]])) +
    labs(x = NULL, y = expression(log[2](S100A8~RPM+1)), colour = NULL) +
    theme_nature() + theme(legend.position = "none")
  save_nature_plot(p_mimic, file.path("09_figures", "Fig_S100A8_Sepsis_vs_SterileSIRS"),
                   89, 70, source_data = d63042)

  # GSE54514: APACHE II and repeated day 1-5 trajectory, patient as unit.
  es54514 <- read_eset_local("GSE54514")
  pd54514 <- pData(es54514); v54514 <- extract_array_gene(es54514, TARGET_GENE, "illuminaHumanv3.db")
  gd <- as.character(pd54514[["group_day:ch1"]])
  d54514 <- data.frame(dataset = "GSE54514", sample_id = rownames(pd54514),
                       patient_id = sub("^[A-Z]+_", "", as.character(pd54514[["group_id:ch1"]])),
                       outcome_group = sub("_D[0-9]+$", "", gd),
                       day = suppressWarnings(as.numeric(sub("^.*_D", "", gd))),
                       apache_ii = suppressWarnings(as.numeric(pd54514[["severity (apacheii):ch1"]])),
                       expression = as.numeric(v54514[rownames(pd54514)]))
  d54514$death <- ifelse(d54514$outcome_group == "NS", 1,
                         ifelse(d54514$outcome_group == "S", 0, NA))
  out_samples[["GSE54514"]] <- d54514
  base54514 <- subset(d54514, day == 1 & outcome_group %in% c("S", "NS"))
  effects[["GSE54514_death"]] <- z_or_model(base54514, "death", "GSE54514", "Mortality (day 1)")
  cor_apache <- cor.test(base54514$expression, base54514$apache_ii,
                         method = "spearman", exact = FALSE)
  apache_result <- data.frame(dataset = "GSE54514", n = sum(complete.cases(base54514[, c("expression", "apache_ii")])),
                              rho = unname(cor_apache$estimate), p_value = cor_apache$p.value)
  write_csv_utf8(apache_result, file.path("04_prognosis_severity", "S100A8_APACHEII_association.csv"))
  p_apache <- ggplot(base54514, aes(apache_ii, expression, colour = factor(death))) +
    geom_point(size = 1.9, alpha = .85) + geom_smooth(method = "lm", se = TRUE, colour = "#5F6368", fill = "#D8D8D8", linewidth = .55) +
    scale_colour_manual(values = c(`0` = palette_nature[["Normal"]], `1` = palette_nature[["Sepsis"]]),
                        labels = c(`0` = "Survivor", `1` = "Non-survivor")) +
    labs(x = "APACHE II score", y = expression(S100A8~expression), colour = NULL) + theme_nature()
  save_nature_plot(p_apache, file.path("09_figures", "Fig_S100A8_APACHEII"), 89, 70, source_data = base54514)

  traj <- subset(d54514, outcome_group %in% c("S", "NS") & day %in% 1:5)
  traj$outcome <- factor(traj$outcome_group, levels = c("S", "NS"), labels = c("Survivor", "Non-survivor"))
  p_traj <- ggplot(traj, aes(day, expression, group = interaction(outcome, patient_id), colour = outcome)) +
    geom_line(alpha = .18, linewidth = .35) + geom_point(alpha = .22, size = .75) +
    stat_summary(aes(group = outcome), fun = mean, geom = "line", linewidth = 1.15) +
    stat_summary(aes(group = outcome), fun = mean, geom = "point", size = 2) +
    scale_colour_manual(values = c(Survivor = palette_nature[["Normal"]], `Non-survivor` = palette_nature[["Sepsis"]])) +
    scale_x_continuous(breaks = 1:5) + labs(x = "Day after enrollment", y = expression(S100A8~expression), colour = NULL) + theme_nature()
  save_nature_plot(p_traj, file.path("09_figures", "Fig_S100A8_Longitudinal_Trajectory"), 120, 76, source_data = traj)

  # GSE272769: independent 30-day mortality and shock endpoints.
  es272769 <- read_eset_local("GSE272769")
  pd272769 <- pData(es272769); v272769 <- extract_array_gene(es272769, TARGET_GENE, "hugene21sttranscriptcluster.db")
  d272769 <- data.frame(dataset = "GSE272769", sample_id = rownames(pd272769),
                        patient_id = sub("_.*$", "", as.character(pd272769$title)),
                        mort30 = as.integer(pd272769[["mort30:ch1"]] == "Yes"),
                        shock = as.integer(pd272769[["shock:ch1"]] == "Yes"),
                        expression = as.numeric(v272769[rownames(pd272769)]))
  out_samples[["GSE272769"]] <- d272769
  effects[["GSE272769_mort30"]] <- z_or_model(d272769, "mort30", "GSE272769", "30-day mortality")
  effects[["GSE272769_shock"]] <- z_or_model(d272769, "shock", "GSE272769", "Septic shock")

  # GSE110487: paired T1/T2 whole-blood response. Paired delta is the analysis unit.
  es110 <- lapply(list.files(file.path("01_data_inventory", "downloaded_data"),
                            pattern = "^GSE110487.*series_matrix", full.names = TRUE),
                  function(f) GEOquery::getGEO(filename = f, getGPL = FALSE))
  pd110 <- do.call(rbind, lapply(es110, pData))
  raw110 <- openxlsx::read.xlsx(file.path("01_data_inventory", "downloaded_data", "GSE110487_rawcounts.xlsx"))
  suppressPackageStartupMessages({library(org.Hs.eg.db); library(AnnotationDbi)})
  ens <- sub("\\..*$", "", as.character(raw110[[1]]))
  mp <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = unique(ens), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  rows <- which(unname(mp[ens]) == TARGET_GENE)
  if (!length(rows)) stop("S100A8 not mapped in GSE110487")
  count_mat <- as.matrix(raw110[, -1, drop = FALSE]); storage.mode(count_mat) <- "numeric"
  lib <- colSums(count_mat, na.rm = TRUE)
  counts <- colSums(count_mat[rows, , drop = FALSE], na.rm = TRUE)
  logcpm <- log2((counts + 2) / (lib + 4) * 1e6)
  d110 <- data.frame(dataset = "GSE110487", sample_id = rownames(pd110), patient_id = pd110[["patient:ch1"]],
                     response = pd110[["clinical classification:ch1"]], timepoint = pd110[["timepoint:ch1"]],
                     expression = as.numeric(logcpm[as.character(pd110$title)]))
  out_samples[["GSE110487"]] <- d110
  wide <- reshape(d110[, c("patient_id", "response", "timepoint", "expression")], idvar = c("patient_id", "response"),
                  timevar = "timepoint", direction = "wide")
  wide$delta_T2_minus_T1 <- wide$expression.T2 - wide$expression.T1
  resp_test <- wilcox.test(delta_T2_minus_T1 ~ response, data = wide, exact = FALSE)
  response_result <- data.frame(dataset = "GSE110487", n_patients = nrow(wide),
                                n_responders = sum(wide$response == "R"), n_nonresponders = sum(wide$response == "NR"),
                                p_value = resp_test$p.value)
  write_csv_utf8(response_result, file.path("04_prognosis_severity", "S100A8_TreatmentResponse_paired.csv"))
  p_resp <- ggplot(wide, aes(response, delta_T2_minus_T1, colour = response)) +
    geom_hline(yintercept = 0, colour = "#A8A8A8", linewidth = .35, linetype = 2) +
    geom_boxplot(width = .45, outlier.shape = NA, colour = "#5F6368", fill = NA, linewidth = .4) +
    geom_point(size = 1.8, position = position_jitter(width = .055, height = 0)) +
    scale_colour_manual(values = c(NR = palette_nature[["Sepsis"]], R = palette_nature[["Normal"]]),
                        labels = c(NR = "Non-responder", R = "Responder")) +
    scale_x_discrete(labels = c(NR = "Non-responder", R = "Responder")) +
    labs(x = NULL, y = expression(Delta*S100A8~"logCPM (T2 - T1)"), colour = NULL) + theme_nature() + theme(legend.position = "none")
  save_nature_plot(p_resp, file.path("09_figures", "Fig_S100A8_TreatmentResponse"), 89, 70, source_data = wide)

  samples <- do.call(rbind, lapply(out_samples, function(x) {
    alln <- unique(unlist(lapply(out_samples, names))); miss <- setdiff(alln, names(x)); for (m in miss) x[[m]] <- NA
    x[, alln, drop = FALSE]
  }))
  write_csv_utf8(samples, file.path("04_prognosis_severity", "prognosis_extension_samples.csv"))
  eff <- do.call(rbind, Filter(Negate(is.null), effects))
  write_csv_utf8(eff, file.path("04_prognosis_severity", "Table_Prognosis_Extended.csv"))
  mimic <- subset(eff, endpoint == "Sepsis vs sterile SIRS")
  write_csv_utf8(mimic, file.path("03_clinical_validation", "Table_ClinicalMimic_Extended.csv"))

  forest <- subset(eff, endpoint != "Sepsis vs sterile SIRS")
  forest$label <- paste0(forest$dataset, "  ", forest$endpoint)
  forest$label <- factor(forest$label, levels = rev(forest$label))
  p_forest <- ggplot(forest, aes(estimate, label)) + geom_vline(xintercept = 1, colour = "#A8A8A8", linewidth = .4) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = .15, linewidth = .45, colour = "#5F6368") +
    geom_point(size = 2.1, colour = palette_nature[["Sepsis"]]) + scale_x_log10() +
    labs(x = "Odds ratio per 1-SD higher S100A8", y = NULL) + theme_forest()
  save_nature_plot(p_forest, file.path("09_figures", "Fig_S100A8_Prognosis_Extended"), 120, 76, source_data = forest)
})
