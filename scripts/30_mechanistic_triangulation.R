# Purpose: triangulate the S100A8/emergency-myeloid/HLA-II relationship across
# independent donor-level datasets and summarize matched TPL-vs-16OH docking.
# These analyses support association/structural plausibility, not causality or
# direct target engagement.
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "theme_nature.R"), encoding = "UTF-8")

fisher_ci <- function(r, n) {
  if (!is.finite(r) || n <= 3 || abs(r) >= 1) return(c(NA_real_, NA_real_))
  z <- atanh(r); se <- 1 / sqrt(n - 3)
  tanh(z + c(-1, 1) * 1.96 * se)
}

with_script_log("30_mechanistic_triangulation", {
  suppressPackageStartupMessages({library(ggplot2); library(data.table)})
  c205 <- read.csv(file.path("08_bulk_singlecell_integration", "GSE205672_mechanistic_correlations.csv"), check.names = FALSE)
  c205$evidence <- ifelse(c205$y == "independent_state_score", "Emergency-myeloid program",
                         ifelse(c205$y == "antigen_presentation_score", "Antigen-presentation program", "NF-kB program"))
  c205 <- c205[, c("dataset", "evidence", "n", "rho", "p_value")]
  c548 <- read.csv(file.path("08_bulk_singlecell_integration", "SCP548_S100A8_HLAII_correlation.csv"), check.names = FALSE)
  c548 <- data.frame(dataset = c548$dataset, evidence = "Antigen-presentation program",
                     n = c548$n_donors, rho = c548$rho, p_value = c548$p_value)

  # Small original GSE167363 dataset retained as supportive, not independent proof.
  a8 <- read.csv(file.path("05_single_cell_pseudobulk", "S100A8_donor_logCPM.csv"), check.names = FALSE)
  st <- read.csv(file.path("06_monocyte_state", "monocyte_donor_state_scores.csv"), check.names = FALSE)
  ind <- read.csv(file.path("07_non_circular_signature", "donor_independent_signature_scores.csv"), check.names = FALSE)
  old <- Reduce(function(x, y) merge(x, y, by = c("donor_id", "group")), list(a8, st, ind))
  old_rows <- lapply(c(antigen_presentation_score = "Antigen-presentation program",
                       Independent_NoA8A9_UCell = "Emergency-myeloid program"), function(label) NULL)
  old_rows <- list()
  for (v in c("antigen_presentation_score", "Independent_NoA8A9_UCell")) {
    ct <- cor.test(old$logCPM, old[[v]], method = "spearman", exact = FALSE)
    old_rows[[v]] <- data.frame(dataset = "GSE167363", evidence = ifelse(v == "antigen_presentation_score", "Antigen-presentation program", "Emergency-myeloid program"),
                                n = nrow(old), rho = unname(ct$estimate), p_value = ct$p.value)
  }
  tri <- rbind(c205, c548, do.call(rbind, old_rows))
  ci <- t(mapply(fisher_ci, tri$rho, tri$n)); tri$ci_low <- ci[, 1]; tri$ci_high <- ci[, 2]
  tri$interpretation <- ifelse(tri$evidence == "Antigen-presentation program" & tri$rho < 0,
                               "Direction supports S100A8-high/HLA-II-low association",
                               ifelse(tri$evidence == "Emergency-myeloid program" & tri$rho > 0,
                                      "Direction supports S100A8-linked emergency myeloid state",
                                      "Direction does not support the proposed mechanism"))
  tri$claim_boundary <- "Association only; no perturbation or causal inference"
  write_csv_utf8(tri, file.path("08_bulk_singlecell_integration", "Table_Mechanistic_Triangulation.csv"))

  plot_tri <- subset(tri, evidence != "NF-kB program")
  plot_tri$label <- paste(plot_tri$dataset, plot_tri$evidence, sep = " | ")
  plot_tri$label <- factor(plot_tri$label, levels = rev(plot_tri$label))
  p <- ggplot(plot_tri, aes(rho, label, colour = evidence)) +
    geom_vline(xintercept = 0, colour = "#A8A8A8", linewidth = .4) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = .12, linewidth = .45) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(`Emergency-myeloid program` = palette_nature[["Sepsis"]],
                                   `Antigen-presentation program` = palette_nature[["Normal"]])) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, .5)) +
    labs(x = "Spearman correlation with donor-level S100A8", y = NULL, colour = NULL) + theme_forest() +
    theme(legend.position = "bottom")
  save_nature_plot(p, file.path("09_figures", "Fig_Mechanistic_Triangulation"), 140, 82, source_data = plot_tri)

  # Matched receptor/scoring/seed comparison; positive delta means 16OH is less favorable.
  dock_root <- file.path(PROJECT_ROOT, "S100A8_triptolide_16OH_comparative_docking_topjournal")
  reps <- read.csv(file.path(dock_root, "results", "comparative_analysis", "matched_local_replicates_with_rmsd.csv"), check.names = FALSE)
  wide <- reshape(reps[, c("receptor", "scoring", "seed", "compound", "best_score_kcal_mol")],
                  idvar = c("receptor", "scoring", "seed"), timevar = "compound", direction = "wide")
  names(wide) <- gsub("best_score_kcal_mol\\.", "", names(wide))
  wide$delta_16OH_minus_TPL <- wide[["16-hydroxytriptolide"]] - wide$triptolide
  ds <- aggregate(delta_16OH_minus_TPL ~ receptor + scoring, wide, function(x) c(mean = mean(x), sd = sd(x), n = length(x)))
  dock <- data.frame(receptor = ds$receptor, scoring = ds$scoring,
                     mean_delta = ds$delta_16OH_minus_TPL[, "mean"], sd = ds$delta_16OH_minus_TPL[, "sd"],
                     n_pairs = ds$delta_16OH_minus_TPL[, "n"])
  dock$se <- dock$sd / sqrt(dock$n_pairs); dock$ci_low <- dock$mean_delta - qt(.975, dock$n_pairs - 1) * dock$se
  dock$ci_high <- dock$mean_delta + qt(.975, dock$n_pairs - 1) * dock$se
  sim <- read.csv(file.path(dock_root, "results", "comparative_analysis", "compound_similarity.csv"), check.names = FALSE)
  contacts <- read.csv(file.path(dock_root, "results", "comparative_analysis", "representative_pose_cross_compound_comparison.csv"), check.names = FALSE)
  dock$morgan_tanimoto <- sim$value[sim$metric == "morgan_radius2_2048_tanimoto_chiral"]
  dock$maccs_tanimoto <- sim$value[sim$metric == "maccs_tanimoto"]
  dock$contact_jaccard <- contacts$contact_jaccard[match(paste(dock$receptor, dock$scoring), paste(contacts$receptor, contacts$scoring))]
  dock$conclusion <- "Closely matched pose/contact pattern; not proof of equal pharmacology or direct binding"
  write_csv_utf8(dock, file.path("08_bulk_singlecell_integration", "Table_TPL_16OH_MatchedDocking.csv"))
  dock$condition <- factor(paste(dock$receptor, toupper(dock$scoring), sep = " | "),
                           levels = rev(paste(dock$receptor, toupper(dock$scoring), sep = " | ")))
  pd <- ggplot(dock, aes(mean_delta, condition)) +
    geom_vline(xintercept = 0, colour = "#A8A8A8", linewidth = .4) +
    geom_errorbar(aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = .13, linewidth = .45, colour = "#5F6368") +
    geom_point(size = 2.2, colour = palette_nature[["Sepsis"]]) +
    labs(x = expression(Delta*" docking score (16OH - triptolide), kcal/mol"), y = NULL) + theme_forest()
  save_nature_plot(pd, file.path("09_figures", "Fig_TPL_16OH_MatchedDocking"), 110, 70, source_data = dock)

  # Update the existing decision file in place.
  dec_file <- file.path("12_manuscript_revision", "final_decision_answers.csv")
  dec <- read.csv(dec_file, check.names = FALSE)
  repl <- c(
    `4` = "进一步支持但仍非完全特异：GSE63042同队列中S100A8区分脓毒症与无菌性SIRS（OR 2.27/每1 SD，P=0.0017）；SCP548供者梯度中感染无器官功能障碍到尿源性脓毒症方向升高，但ICU-SEP与ICU-NoSEP差异不显著。因此可回应泛炎症质疑，但不能声称普适临床特异性。",
    `5` = "预后与严重程度已扩展：GSE272769中30天死亡和休克仅呈边缘方向性；GSE54514中基线S100A8与APACHE II及死亡不显著；GSE63042院内死亡不显著；GSE110487配对治疗反应已分析。结论仍不支持把S100A8写成稳定预后标志物。",
    `6` = "独立供者级证据已显著加强：SCP548按65名供者聚合单核细胞，GSE205672以15名健康对照和36名脓毒症患者的未刺激纯化单核细胞独立验证，避免细胞伪重复。",
    `7` = "跨独立队列支持S100A8-high/HLA-II-low关联：GSE205672抗原呈递分数在脓毒症下降且与S100A8负相关（rho=-0.462，P=0.00065）；SCP548供者级相关同方向（rho=-0.299，P=0.0156）。该证据仍属关联，不能替代蛋白或因果实验。",
    `8` = "外部验证支持应急髓系状态：GSE205672中不含S100A8/S100A9的30基因程序显著升高；SCP548中MS1比例及独立程序随部分临床严重度梯度升高。",
    `9` = "是：非循环签名已在独立GSE205672纯化单核细胞中验证，并在SCP548供者级单核细胞中获得方向一致证据，不再局限于原小样本内部稳健性。",
    `10` = "S100A8与应急髓系程序及HLA-II低表达的关联在GSE205672、SCP548和原GSE167363中进行三角验证；NF-kB汇总转录分数在独立GSE205672中未获支持，故NF-kB机制必须保持探索性。",
    `12` = "临床特异性、预后/严重程度、独立单细胞、非循环状态和TPL/16OH比较对接均已明显补强；但预后结果多为阴性或边缘，不能转化为阳性主张。",
    `14` = "计算层面已基本穷尽：临床相似对照、纵向/结局、独立供者单细胞、纯化单核细胞、非循环签名、HLA-II关联和TPL/16OH匹配对接均已完成。仍无法由生物信息学解决直接结合、S100A8扰动/救援及NF-kB功能因果性，需湿实验或明确降级结论。")
  for (q in names(repl)) dec$answer[dec$question == as.integer(q)] <- repl[[q]]
  write_csv_utf8(dec, dec_file)

  ev_file <- file.path("13_response_to_reviewers", "reviewer_evidence_matrix.csv")
  ev <- read.csv(ev_file, check.names = FALSE)
  upd <- function(id, analysis, scripts, figure, table, result, status) {
    i <- grep(paste0("^R1\\.", id, ":"), ev$reviewer_comment)
    if (!length(i)) return()
    ev$analysis_performed[i] <<- analysis; ev$script[i] <<- scripts; ev$figure[i] <<- figure; ev$table[i] <<- table
    ev$main_result[i] <<- result; ev$status[i] <<- status
  }
  upd(3, "GSE63042 sterile-SIRS comparison plus SCP548 clinical cohort gradients", "27,28", "Fig_S100A8_Sepsis_vs_SterileSIRS; Fig_SCP548_ClinicalGradient", "Table_ClinicalMimic_Extended; Table_SCP548_DonorComparisons", "Direct same-cohort SIRS comparison is positive, while ICU mimic discrimination remains heterogeneous.", "Fully addressed")
  upd(4, "Independent mortality, shock, APACHE II, longitudinal and paired treatment-response analyses", "27", "Fig_S100A8_Prognosis_Extended; Fig_S100A8_APACHEII; Fig_S100A8_Longitudinal_Trajectory; Fig_S100A8_TreatmentResponse", "Table_Prognosis_Extended", "Most endpoints are null or borderline; this answers the request but does not support a prognostic claim.", "Fully addressed")
  upd(9, "Independent SCP548 monocyte pseudobulk and donor-level state abundance across 65 donors", "28", "Fig_SCP548_DonorPseudobulk_S100A8; Fig_SCP548_MS1_Abundance", "Table_SCP548_DonorComparisons", "Independent donor-level analysis removes cell-level pseudoreplication.", "Fully addressed")
  upd(10, "Independent emergency-myeloid/MS1 and antigen-presentation programs in SCP548 and GSE205672", "28,29", "Fig_SCP548_ClinicalGradient; Fig_GSE205672_IndependentStateScore", "Table_GSE205672_MonocyteValidation", "The state is independently supported as an emergency-myeloid/MS1-like program.", "Fully addressed")
  upd(11, "Cross-dataset donor/sample-level HLA-II transcriptional validation", "28-30", "Fig_GSE205672_AntigenPresentation; Fig_Mechanistic_Triangulation", "Table_Mechanistic_Triangulation", "Independent transcriptomic evidence supports HLA-II suppression but does not replace protein validation.", "Partially addressed")
  upd(12, "External validation of the S100A8/S100A9-excluded program", "28,29", "Fig_GSE205672_IndependentStateScore; Fig_SCP548_ClinicalGradient", "Table_GSE205672_MonocyteValidation", "Non-circular program reproduces in purified monocytes and an independent single-cell cohort.", "Fully addressed")
  upd(14, "Matched TPL-versus-16OH chemical similarity and 5HLO docking across receptor dimers, scoring functions and ten seeds", "30", "Fig_TPL_16OH_MatchedDocking", "Table_TPL_16OH_MatchedDocking", "High structural similarity and matched contact geometry support a shared docking rationale, not pharmacological equivalence.", "Fully addressed")
  upd(19, "Independent donor-level NF-kB transcriptional score test", "29,30", "Fig_GSE205672_NFkBScore", "Table_GSE205672_MonocyteValidation; Table_Mechanistic_Triangulation", "NF-kB score was not supported in GSE205672; causal NF-kB depth remains unresolved.", "Partially addressed")
  write_csv_utf8(ev, ev_file)
})
