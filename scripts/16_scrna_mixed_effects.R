# 输入：GSE167363 单核细胞 log-normalized 表达。
# 输出：重点基因的 donor 随机截距混合模型，仅作为方向敏感性分析。
rm(list = ls()); options(stringsAsFactors = FALSE, warn = 1)
source(file.path("config", "project_config.R"), encoding = "UTF-8")
source(file.path("functions", "io_functions.R"), encoding = "UTF-8")
source(file.path("functions", "scrna_functions.R"), encoding = "UTF-8")

with_script_log("16_scrna_mixed_effects", {
  suppressPackageStartupMessages({library(Seurat); library(lme4)})
  in_file <- file.path("06_monocyte_state", "GSE167363_monocyte_state_annotated.rds")
  assert_file(in_file); mono <- readRDS(in_file); DefaultAssay(mono) <- "RNA"
  focus <- c("S100A8","S100A9","S100A12","IL1B","FCN1","CTSS","NFKBIA","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","CD74","CIITA","IL1R2","RETN","LILRB1","LILRB3","SOCS3","ARG1","MPO","ELANE")
  focus <- intersect(focus, rownames(mono))
  mat <- SeuratObject::LayerData(mono[["RNA"]], layer = "data")[focus, , drop = FALSE]
  base <- data.frame(donor_id = factor(mono$donor_id), group = relevel(factor(mono$group), "Control"))
  rows <- lapply(focus, function(g) {
    d <- base; d$expression <- as.numeric(mat[g, ])
    fit <- tryCatch(lme4::lmer(expression ~ group + (1 | donor_id), data = d, REML = FALSE,
                              control = lmerControl(check.nobs.vs.nRE = "ignore")), error = identity)
    if (inherits(fit, "error")) return(data.frame(gene=g,n_cells=nrow(d),n_donors=nlevels(d$donor_id),estimate=NA,se=NA,ci_low=NA,ci_high=NA,singular=NA,status=conditionMessage(fit)))
    cf <- summary(fit)$coefficients["groupSepsis", ]
    data.frame(gene=g,n_cells=nrow(d),n_donors=nlevels(d$donor_id),estimate=cf[["Estimate"]],se=cf[["Std. Error"]],ci_low=cf[["Estimate"]]-1.96*cf[["Std. Error"]],ci_high=cf[["Estimate"]]+1.96*cf[["Std. Error"]],singular=lme4::isSingular(fit),status="Sensitivity only; donor-level pseudo-bulk is primary")
  })
  out <- do.call(rbind, rows)
  write_csv_utf8(out, file.path("05_single_cell_pseudobulk", "Table_MixedEffects_FocusGenes.csv"))
})
