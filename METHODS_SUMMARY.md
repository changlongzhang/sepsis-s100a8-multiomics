# Methods summary aligned with the final revision

## Whole-blood discovery and feature prioritization

The primary GSE65682 whole-blood analysis used a log2 RMA expression matrix. Probe sets were collapsed to gene level with the maximum-probe rule for the primary analysis; a mean-probe analysis was retained as a sensitivity analysis. Differential expression used limma. The corrected differential-expression set contained 1,508 genes, and 265 genes overlapped the primary WGCNA module used for downstream interpretation.

The primary signed WGCNA network included 802 samples and 11,719 finite, nonzero-variance genes. It used Pearson correlation, soft-threshold power 19, minimum module size 200, `deepSplit = 2`, `pamRespectsDendro = FALSE`, and `mergeCutHeight = 0.25`. The sepsis-associated black module had module–trait correlation 0.409 (FDR 3.72 × 10^-33); S100A8 module membership was 0.700. The 50% MAD and mean-probe networks were separate prespecified sensitivity analyses; S100A8 was not manually restored when it fell below a variability filter.

Feature prioritization was leakage controlled: feature screening and model fitting were performed within training folds. Comparator, class-balance, temporal, endotype, prognosis, and treatment-response analyses are retained as context-specific validations. Higher S100A8 expression was oriented toward sepsis in ROC analyses; no automatic direction reversal was used.

## Donor-aware single-cell analysis

The primary single-cell pseudobulk analysis used donor-level aggregation in a prespecified GSE167363 subset of two controls and four sepsis donors. Cell-level visualizations were not treated as independent donor-level inference. GSE205672 and SCP548 were used as independent validation resources. The analyses evaluate S100A8 alongside emergency-myeloid, HLA-DR-low, antigen-presentation, and NF-κB-associated programs.

## In vitro and computational structural analyses

Initial three-group in vitro comparisons used one-way ANOVA. Loss-of-function factorial models included biological repeat as a blocking factor; planned contrasts used multiplicity control as reported in the final S7 Table. RT-qPCR statistics use archived relative-expression values because sample-level Ct/ΔCt values were unavailable.

For molecular dynamics, the protein and water were modeled using AMBER99SB-ILDN and TIP3P; triptolide was modeled using OpenFF Sage 2.2.1 with AM1-BCC ELF10 charges. Docking and molecular-dynamics outputs are computational summaries and do not establish direct binding. The triptolide-dependency experiment was negative and is reported explicitly.

## Software versions recorded for the final analysis

R 4.5.1; limma 3.64.3; edgeR 4.6.3; caret 7.0-1; glmnet 4.1-10; mboost 2.9-11; xgboost 1.7.11.1; ranger 0.17.0; pROC 1.19.0.1; PRROC 1.4; dcurves 0.5.0; metafor 4.8-0; Seurat 5.3.1.9999; decoupleR 2.14.0; UCell 2.12.0; GSVA 2.2.0; GROMACS 2023.2; AutoDock Vina 1.2.7; Meeko 0.7.1; OpenFF Sage 2.2.1. The archived materials did not record an ImageJ version.

Detailed parameters, model specifications, exclusions, and limitations remain in the manuscript, Supporting Tables, scripts, and machine-readable outputs; this summary does not replace those sources.
