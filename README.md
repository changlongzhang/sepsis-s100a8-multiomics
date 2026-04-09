# Sepsis Multi-omics Diagnostic Biomarker Discovery and Immune Microenvironment Analysis: Full Pipeline Reproducibility Guide

This open-source repository contains all R and Bash scripts utilized in our study. To meet the most rigorous reproducibility standards, we provide highly detailed execution instructions, input/output specifications, and core methodological highlights for every essential module within our data analysis pipeline.

## 🌐 Data Availability

Adhering to the principles of Open Science, all raw transcriptomic and single-cell omics data involved in the full-pipeline analysis are hosted in the **NCBI Gene Expression Omnibus (GEO)** database. Researchers can directly access the relevant raw files and metadata via the following links:

### 1. Bulk RNA-seq Discovery Cohort
* **[GSE65682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65682)**: The core training set of this study. It contains whole-blood transcriptomic data from sepsis patients based on the Affymetrix GPL13667 platform, utilized for differential expression analysis, WGCNA module construction, and machine learning model training.

### 2. scRNA-seq Mechanistic Cohort
* **[GSE167363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167363)**: Contains single-cell transcriptomic data of peripheral blood mononuclear cells (PBMCs) from sepsis patients utilizing the 10X Genomics platform. Used to identify subpopulation-specific marker localization, infer transcription factor activities, and execute `scTenifoldKnk` in silico knockout simulations.

### 3. Independent External Validation Cohorts
To ensure the robustness of the identified biomarkers, multiple independent datasets were introduced for cross-platform validation:
* **[GSE9692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9692)** (Affymetrix GPL201)
* **[GSE26440](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26440)** (Illumina GPL10558)
* **[GSE28750](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28750)** (Illumina GPL10558)
* **[GSE69528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69528)** (Illumina GPL10558)

### 4. Auxiliary Annotation Files
* Platform probe mapping files (e.g., `GPL13667-15572.txt`) can be downloaded from the corresponding GEO platform pages or directly cited from the pre-organized version in the `Annotations/` folder within this repository.

### 5. Large Processed Datasets & Objects
Due to GitHub repository file size constraints (single files < 25MB), large pre-processed expression matrices, fully clustered single-cell Seurat objects (`.rds`), and continuous trajectory files from molecular dynamics simulations (`.xtc/.trr`) required for reproduction have been packaged and hosted on the international open-science archiving platform **Zenodo**. Readers must first download this data package, extract it, and place the contents into the corresponding `data/` directory of this repository for seamless script execution.
* **Zenodo Exclusive Download Link**: [https://doi.org/10.5281/zenodo.19480952](https://doi.org/10.5281/zenodo.19480952)
* **Global Unique DOI**: `10.5281/zenodo.19480952`

## 💻 Software Environment & Reproducibility Statement
To ensure absolute reproducibility of all computational results in this repository (particularly modules involving machine learning feature selection, non-linear dimensionality reduction, and permutation testing), we detail the underlying execution environment, core dependency package versions, and globally hardcoded random seeds. By configuring an identical local environment and running the scripts, operators will obtain statistical outputs and high-fidelity figures exactly consistent with the original manuscript.

### 1. Global Random Seeds
Throughout our analytical pipeline, all algorithms involving random sampling or heuristic optimization have strictly locked random seeds at the code level:
* **Machine Learning & Statistical Pipelines** (Including 10-fold cross-validation, LASSO/RF/XGBoost training, CIBERSORT permutation testing, and Bootstrap calibration resampling): Globally set to `set.seed(123)`.
* **High-Dimensional Single-Cell Dimensionality Reduction & Clustering** (Including UMAP projection, t-SNE, and K-means algorithms): Globally set to `set.seed(1234)`.

### 2. Core R Environment & Package Dependencies
* **Underlying R Environment**: R version 4.5.2

The following are the core R packages invoked in this study along with their precise version numbers:

| Domain | Package Name | Version | Functional Description |
| :--- | :--- | :--- | :--- |
| **scRNA-seq Analysis & QC** | `Seurat` | 5.3.1.9999 | Core framework for single-cell data processing |
| | `DoubletFinder` | 2.0.6 | Prediction and removal of doublets |
| | `harmony` | 1.2.4 | Multi-sample batch effect correction |
| | `clustree` | 0.5.1 | Clustering resolution dendrogram evaluation |
| **Transcriptomics & Network Construction** | `WGCNA` | 1.73 | Weighted Gene Co-expression Network Analysis |
| | `DESeq2` | 1.48.2 | Differential expression analysis |
| | `limma` | 3.64.3 | Microarray data processing and differential analysis |
| | `flashClust` | 1.1.2 | Highly efficient hierarchical clustering |
| | `dynamicTreeCut` | 1.63.1 | Dynamic branch cutting and module identification |
| **Functional Mechanisms & Virtual Knockout** | `clusterProfiler` | 4.16.0 | High-throughput GO/KEGG enrichment analysis |
| | `org.Hs.eg.db` | 3.21.0 | Underlying database for Gene ID mapping |
| | `scTenifoldKnk` | 1.0.2 | Single-cell in silico gene knockout simulation |
| | `decoupleR` | 2.14.0 | Transcription factor activity network inference |
| | `OmnipathR` | 3.16.2 | Prior biological knowledge base retrieval |
| **Machine Learning & Interpretability** | `caret` | 7.0.1 | Machine learning ensemble pipeline and resampling |
| | `glmnet` | 4.1.10 | LASSO penalized linear regression |
| | `randomForest` | 4.7.1.2 | Random Forest algorithm foundation |
| | `xgboost` | 1.7.11.1 | eXtreme Gradient Boosting tree algorithm foundation |
| | `mboost` | 2.9.11 | Generalized linear boosting regression |
| | `kernelshap` | 0.9.1 | Robust model-agnostic SHAP value estimation |
| | `shapviz` | 0.10.2 | Advanced SHAP visualization analysis |
| **Clinical Efficacy & Multidimensional Validation** | `pROC` | 1.19.0.1 | Traditional ROC curve and AUC calculation |
| | `precrec` | 0.14.5 | PR curve and advanced AUPRC calculation |
| | `rms` | 8.0.0 | Logistic regression and Bootstrap calibration analysis |
| | `tigeR` | 1.0.0 | Tumor/Immune microenvironment assessment foundation |
| **Data Retrieval & Engineering Acceleration** | `GEOquery` | 2.76.0 | Automated public database parsing |
| | `linkET` | 0.0.7.4 | Mantel test matrix correlation plots |
| | `doParallel` | 1.0.17 | Multi-core parallel computing acceleration support |

---

## 📦 Module 01: Microarray Data Preprocessing & Probe Annotation (`01_GSE65682_preprocessing.R`)

### 🎯 Module Overview
Acting as the starting point of the entire analysis pipeline, this module automates the downloading of the GSE65682 raw microarray dataset and executes high-quality probe filtering, ID conversion, and expression matrix reconstruction.

### 🛠 Environment & Dependencies
* `GEOquery`, `data.table`, `dplyr`, `tibble`
* **Note: For exact R and package versions, refer to the Supplementary Materials of this paper.**

### 📥 Required Input Files
1. `GPL13667-15572.txt`: Customized microarray platform probe annotation file.

### ⚙️ Core Analytical Workflow
1. **Automated Retrieval**: Uses `getGEO` to fetch GSE65682 clinical and expression objects.
2. **Scale Reversion**: Subjects the default Log2 expression matrices generated by the RMA algorithm to exponential anti-log transformation ($2^x$) to revert them to a linear scale.
3. **Probe Mapping**: Reads GPL annotations, filtering invalid Gene Symbols (e.g., rows containing null values or "---").
4. **Probe Deduplication**: Associates matrices via `left_join` and resolves redundancy ("multiple probes mapping to one gene") using a strict `max` aggregation strategy to retain maximum transcript abundance.

### 📤 Complete Output Details
**📊 Data & Table Files:**
* `GSE65682_gene_matrix.csv`: Cleaned, deduplicated pure gene expression matrix (row names as standard Gene Symbols).
* `GSE65682_clinical.csv`: Complete clinical phenotypes and metadata extracted from raw objects.

### 👣 Detailed Reproduction Guide
1. Place `GPL13667-15572.txt` into the working directory.
2. Modify line 8 in the script: `setwd("Your Local Path")`.
3. Execute the script to wait for network download and automatic generation of the two CSV files above.
4. If a network timeout occurs during execution, researchers are advised to manually download the raw matrix via the GEO link provided and place it in the directory.

### 🔬 Methodological & Statistical Highlights
> **Scale Reversion & Negative Value Evasion Strategy:**
> In this module, we applied an exponential transformation ($2^x$) to the Log2-scaled expression matrices outputted by the RMA algorithm, returning them to a linear scale. The underlying rationale is that **negative expression values** can trigger mathematical exceptions in input-sensitive machine learning algorithms (e.g., non-negative matrix factorization). Furthermore, to handle multi-probe redundancies, we strictly employed a `max` aggregation strategy to represent the final gene expression value, ensuring maximum transcript sensitivity.

---

## 📦 Module 02: DESeq2-Based Differential Expression Analysis (`02_DEG.R`)

### 🎯 Module Overview
Executes rigorous Disease vs. Normal differentially expressed gene (DEG) analysis and utilizes highly customized graphical parameters to render publication-quality volcano plots.

### 🛠 Environment & Dependencies
* `DESeq2`, `dplyr`, `tibble`, `stringr`, `ggplot2`, `ggrepel`

### 📥 Required Input Files
1. `GSE65682_gene.csv`: (Renamed from the expression matrix generated in Module 01).
2. `GSE65682_Groups.csv`: Pre-organized sample grouping table (containing Sample and Group columns).

### ⚙️ Core Analytical Workflow
1. **Data Cleaning & Integer Conversion**: Intercepts negative values in the matrix (replacing with 0) and forces the continuous signals into the non-negative integer matrix required by DESeq2 using the `round()` function.
2. **Sample Alignment**: Strictly matches expression matrix column names with clinical grouping files via intersection operations to ensure correct factor level mapping (Normal as Control, Disease as Experimental).
3. **Model Fitting**: Constructs the `DESeqDataSet` object and performs generalized linear model fitting based on a negative binomial distribution.
4. **Result Extraction & Repulsion Rendering**: Identifies up/down-regulated genes using thresholds of padj < 0.05 and |Log2FC| > 1, extracting coordinates for target genes (e.g., S100A8) for graphical rendering.

### 📤 Complete Output Details
**📊 Data & Table Files:**
* `Disease_vs_Normal_all_results_padj.csv`: Complete table of differential analysis statistics across the genome.
* `Disease_vs_Normal_sig_genes_padj.csv`: List of significantly differentially expressed genes after strict threshold filtering.
**🖼️ Figures:**
* `Disease_vs_Normal_volcano_padj_S100A8_Fixed.pdf`: Highly customized vector volcano plot (S100A8 highlighted with a bold red box and anti-collision layout).

### 👣 Detailed Reproduction Guide
1. Ensure the two CSV input files are present in the directory.
2. Modify the working directory in line 9 of the script.
3. Run the script with one click. The code will automatically perform model calculations and output statistical tables and the PDF volcano plot.

### 🔬 Methodological & Statistical Highlights
> **DESeq2 Integer Adaptation & Customized Label Avoidance:**
> 1. **DESeq2 `round()` Adaptation**: Because we extracted linearized continuous signals from microarrays, we innovatively introduced non-negative truncation and `round()` integer approximation strategies. This perfectly integrates the signal matrix with `DESeq2`'s discrete dispersion estimation engine, circumventing underlying model crashes.
> 2. **`ggrepel` Customized Label Repulsion**: To prevent core genes from being obscured in the volcano plot, we deeply engaged the `ggrepel` algorithmic engine. By hardcoding physical collision parameters and precise leader-line micro-adjustments, we established a customized label spatial avoidance mechanism.

---

## 📦 Module 03: Weighted Gene Co-expression Network Analysis (WGCNA) (`03_WGCNA_Analysis.R`)

### 🎯 Module Overview
Mines synergistic gene expression patterns, clusters high-dimensional genes into co-expression modules, and pinpoints core modules significantly associated with the sepsis phenotype.

### 🛠 Environment & Dependencies
* `WGCNA`, `flashClust`, `dynamicTreeCut`, `limma`, `factoextra`, `ggpubr`

### 📥 Required Input Files
1. `GSE65682_gene.csv`: Gene expression matrix.
2. `GSE65682_Groups.csv`: Clinical phenotype grouping file.

### ⚙️ Core Analytical Workflow
1. **Preprocessing & Dimensionality Cleaning**: Filters low-expression genes and performs Z-score sample-level normalization on the full matrix.
2. **Rigorous Outlier Sample Removal**: Executes PCA dimensionality reduction, calculates the Mahalanobis distance of each sample in the principal component space, and automatically strips outlier points deviating >3 SD from the centroid.
3. **Adaptive Network Optimization**: Iterates through powers 1–20 to automatically identify the minimum Soft Power meeting scale-free topology fit ($R^2 > 0.9$).
4. **Module Division & Merging**: Constructs a Signed Network, generates a Topological Overlap Matrix (TOM), performs dynamic tree cutting, and merges sub-modules with similarity >0.75.
5. **Phenotype Association Calculation**: Calculates correlation P-values and FDR between module Eigengenes and binary clinical labels to screen for significant modules and generate charts.

### 📤 Complete Output Details
**📊 Data & Table Files:**
* `module_group_association.csv`: Complete statistical correlation table between modules and clinical phenotypes.
* `module_[color]_genes.txt`: The full internal gene set for each clustered module.
* `DIFF_MODULE_[color]_GENES.txt`: Candidate gene sets outputted **only for significantly associated modules**.
* `preprocessed_data.RData` & `wgcna_analysis_complete.RData`: Complete environment and network topology snapshots.
**🖼️ Figures:**
* QC Plots: `sample_clustering_before_filtering.pdf`, `pca_before_outlier_removal.pdf`, `pca_after_outlier_removal.pdf`, `sample_clustering_final.pdf`
* Topology & Clustering Plots: `soft_threshold_selection.pdf`, `gene_dendrogram_modules.pdf`, `merged_modules.pdf`
* Correlation Heatmaps: `module_group_correlation_heatmap_all.pdf`, `module_group_correlation_heatmap_significant.pdf`
* Boxplots: `module_[color]_boxplot_half_border.pdf` (Automatically generates distribution differences for significant modules).

### 👣 Detailed Reproduction Guide
1. Prepare the two input files and modify the working path.
2. Run the script. WGCNA is computationally intensive (multi-threading is enabled), so please be patient. All dendrograms and heatmaps will be sequentially saved to the directory.

### 🔬 Methodological & Statistical Highlights
> **Multi-dimensional Outlier Removal & Signed Adaptive Network Optimization:**
> 1. **Mahalanobis Distance-Based Strict Outlier Removal (>3 SD)**: Innovatively combines PCA dimensionality reduction with Mahalanobis distance metrics to accurately calculate deviations in low-dimensional space, automatically stripping outliers exceeding 3 standard deviations.
> 2. **Signed Network & Adaptive Soft-Thresholding**: Explicitly sets `networkType = "signed"` to construct a signed network, strictly distinguishing between positive synergistic regulation and negative antagonistic effects among genes. Built-in adaptive loops automatically lock onto the optimal power.

---

## 📦 Module 04: Core Overlapping Gene Extraction (`04_Extract overlapping genes.R`)

### 🎯 Module Overview
Extracts overlapping genes derived from multi-dimensional analysis strategies (WGCNA modules vs. DEGs) and renders high-precision Venn diagrams using a physical coordinate engine.

### 🛠 Environment & Dependencies
* `VennDiagram`, `grid`

### 📥 Required Input Files
Extract the significant genes outputted from Module 02 into a single column and rename to `DEGs.txt`; rename the target color module from Module 03 to `WGCNA.txt`.
1. `WGCNA.txt`: Gene sets from specified significant modules (see Materials and Methods).
2. `DEGs.txt`: Significant differentially expressed gene sets from Module 02.

### ⚙️ Core Analytical Workflow
1. **List Reading**: Reads the two TXT files without headers into 1D vectors.
2. **Set Operations**: Executes `intersect()` and `setdiff()` to calculate overlaps and unique subsets.
3. **Coordinate Rewriting**: Disables the native scaling and labeling of `VennDiagram`, activates the `grid` engine, and uses hardcoded absolute coordinates to write labels directly into the vector circles.

### 📤 Complete Output Details
**🖼️ Figures:**
* `Venn_red_blue_horizontal.pdf`: High-contrast custom Venn diagram with perfectly centered labels and symmetrical, equally-sized circles.

### 👣 Detailed Reproduction Guide
1. Ensure the two gene set TXT files are prepared. Modify the path and run the code to generate the PDF. Manually save the overlapping genes generated in the R environment as `Overlap_Genes.txt` for downstream use.

### 🔬 Methodological & Statistical Highlights
> **Underlying Rendering Override & Absolute Coordinate Alignment:**
> Disabled the package's built-in text rendering engine (hardcoded `cex = 0` and `scaled = FALSE`) to invoke R's foundational plotting system, the `grid` package. By utilizing `grid.text()`, we redrew all set cardinalities and intersection quantities directly onto the absolute coordinate system of the physical canvas, achieving perfect physical centering and ensuring symmetrical representation.

---

## 📦 Module 05: Functional Enrichment Analysis & Custom Visualization (`05_Enrichment Analysis.R`)

### 🎯 Module Overview
Relies on the `clusterProfiler` engine to perform high-throughput GO and KEGG semantic classification enrichment analysis on overlapping genes.

### 🛠 Environment & Dependencies
* `clusterProfiler`, `org.Hs.eg.db`, `ggplot2`, `stringr`

### 📥 Required Input Files
1. `Overlap_Genes.txt`: The intersection result from Module 04 (one Gene Symbol per line).

### ⚙️ Core Analytical Workflow
1. **ID Conversion**: Uses `org.Hs.eg.db` to losslessly map Gene Symbols to ENTREZ IDs required for underlying computations.
2. **GO Analysis & Typography**: Executes `enrichGO`, extracts the Top 5 items for BP/CC/MF, and implements RGB color mapping with left-aligned typography.
3. **Regex Semantic Matching (KEGG)**: Extracts KEGG results, uses `grepl` with a pre-defined core keyword matrix to forcefully dimensionally-reduce and classify pathway descriptions (into Infection, Metabolism, Immunity, etc.), and re-levels factors to construct a visual color gradient.

### 📤 Complete Output Details
**🖼️ Figures:**
* `Figure_A_GO.pdf`: Minimalist bar plot for GO enrichment.
* `Figure_B_KEGG.pdf`: Customized KEGG bar plot featuring regex-driven secondary semantic classification and visual gradient mapping.

### 👣 Detailed Reproduction Guide
1. Verify `Overlap_Genes.txt` exists in the working directory. Run the script; the engine will auto-download online databases, complete alignments, and render PDFs.

### 🔬 Methodological & Statistical Highlights
> **Regex-Driven KEGG Semantic Reclassification & Visual Gradient Mapping:**
> Innovatively introduces a secondary semantic classification algorithm based on regular expressions. The algorithm hard-codes scans for core keywords within pathway descriptions, intelligently dimensionally-reducing flat pathways into three macro-pathological modules. Building on this, factor re-leveling is utilized to achieve a flawless bottom-up visual gradient mapping (e.g., from green to blue to red).

---

## 📦 Module 06: Machine Learning Training Set Reconstruction (`06_Prepare ML Data.R`)

### 🎯 Module Overview
Executes precise data cleaning and physical transposition on high-dimensional expression matrices, transforming text-based clinical phenotypes into tensor structures and binary numerical labels required by downstream algorithms (e.g., Scikit-Learn/caret).

### 🛠 Environment & Dependencies
* `openxlsx`, R Base Functions

### 📥 Required Input Files
1. `Overlap_Genes.txt`: Feature gene name list.
2. `GSE65682_gene.csv`: Cleaned full microarray matrix.
3. `GSE65682_Groups.csv`: Metadata containing text-based group labels.

### ⚙️ Core Analytical Workflow
1. **Space Truncation**: Filters the expression matrix using overlapping genes, drastically reducing the feature space.
2. **Matrix Transposition**: Executes the `t()` function on the filtered matrix, flipping "Rows: Genes, Cols: Samples" to "Rows: Samples, Cols: Features".
3. **Label Encoding**: Extracts clinical groupings and forces "Control/Disease" translation into `1/0` numerical response variables using `ifelse`.
4. **Attribute Merging**: Performs precise inner joins based on SampleID and reconstructs dataframe column order (ID first, Labels last).

### 📤 Complete Output Details
**📊 Data & Table Files:**
* `ML_training_data_prepared.csv`: A highly standardized 2D machine learning training set (SampleID + Feature Genes + Binary Group) ready for seamless algorithm integration.

### 👣 Detailed Reproduction Guide
1. Prepare the 3 dependency files and set paths. Run the script. The console will print a data preview and locally generate a standardized CSV file.

### 🔬 Methodological & Statistical Highlights
> **Structured Feature Transposition & Target Variable Label Encoding:**
> Strictly implements target variable numerical encoding, safely translating clinical phenotypes into binary, machine-readable numbers. This preprocessing step completely prevents type-crashing during downstream LASSO or XGBoost loss function calculations and provides a highly standardized binary dependent variable (Y) for subsequent continuous-probability metric evaluations (ROC, PR curves, and Brier scores).

---

## 📦 Module 07: Multi-Model Ensemble Selection, Single-Gene Validation & Deep Interpretability (`07_Multi-Model Feature Selection.R`)

### 🎯 Module Overview
This module acts as the core engine for feature reduction and biomarker discovery. Addressing high-dimensional omics data, this script constructs a multi-dimensional machine learning ensemble pipeline containing LASSO, Random Forest, glmBOOST, and XGBoost. After extracting consensus biomarkers via 10-fold cross-validation (with 3 repeats), it performs rigorous single-gene independent ROC validation on core targets (MAFG, S100A8). Furthermore, it deeply integrates the `kernelshap` explanation algorithm, shattering the "black box" nature of complex tree models by strictly locking the Top 6 key genes for local and global feature attribution.

### 🛠 Environment & Dependencies
* **Core Algorithm Libraries**: `glmnet`, `caret`, `randomForest`, `mboost`, `xgboost`
* **Model Explanation & Validation**: `shapviz`, `kernelshap`, `fastshap`, `pROC`
* **Advanced Visualization**: `ggplot2`, `cowplot`, `ggbeeswarm`, `patchwork`, `UpSetR`, `reshape2`, `ggsci`

### 📥 Required Input Files
1. `Train_Data.csv`: Highly standardized expression matrix combined with `1/0` label data generated by upstream feature engineering.

### ⚙️ Core Analytical Workflow
1. **Multicollinearity Processing & Ensemble Selection**: Removes Near-Zero Variance (NZV) predictors, imputes missing values via medians, and performs Z-score normalization. Aligns the evaluation baselines of four heterogeneous models.
2. **Consensus Feature Extraction & Benchmarking**: Extracts feature coefficients/Gini importance/Information gain across models, calculates cross-validated ROC curves, and locks intersection genes using the UpSet algorithm.
3. **Targeted Biomarker Validation**: Detaches from the multi-feature system to isolate core target genes (MAFG, S100A8), independently calculating and restoring customized dual-gene ROC curves matching journal visual standards.
4. **Hardcoded Top-6 Targeted SHAP Interpretability**: Utilizes the robust model-agnostic estimator `kernelshap` to estimate global marginal contributions based on background sample databases. Locks the Top 6 highest-ranked variables, forcefully constructing targeted explanation objects to render mean-based importance barplots, non-linear dependence plots, and patient-level waterfall attribution plots.

### 📤 Complete Output Details
**🖼️ Figures (All outputted to `FinalPlots/` directory):**
* **Independent Model QC Plots**:
  * `1_LASSO_Path_Classic_Labeled_Text.pdf`
  * `2_LASSO_CV_Classic_Labeled_Text.pdf`
  * `3B_RF_OOB_Error_Curve.pdf`
  * `3A_GLMBOOST_CV_Risk_Curve.pdf`
  * `3B_GLMBOOST_Coefficient_Path.pdf`
* **Uniform Gradient Style Importance Barplots (Top 5)**:
  * `1_LASSO_Coefficients_Barplot.pdf` 
  * `3A_RF_Importance_Top5.pdf` 
  * `3C_GLMBOOST_Final_Coefficients.pdf`
  * `4A_XGBoost_Importance_Top5.pdf` 
* **Ensemble & Independent Diagnostic Efficacy**:
  * `5A_MultiModel_CV_ROC_Comparison.pdf` (4-model CV ROC Comparison)
  * `5B_UpSet_Feature_Overlap.pdf` (Ensemble feature intersection)
  * `7_Single_Gene_ROC_MAFG_S100A8.pdf` (Targeted single-gene ROCs)
* **Deep Interpretability SHAP Atlas**:
  * `8A_SHAP_Feature_Importance_Barplot.pdf` (Absolute importance ranking of Top 6)
  * `8B_SHAP_BeeSwarm_Plot.pdf` (High/Low expression push distribution beeswarm plot)
  * `8C_SHAP_Feature_Dependence.pdf` (Core feature non-linear dependence plot)
  * `8D_SHAP_Sample_Waterfall.pdf` (Micro-individual predictive deconstruction waterfall plot)

### 👣 Detailed Reproduction Guide
1. Ensure `Train_Data.csv` is prepared and modify the absolute working directory (`setwd()`) on line 11.
2. Run the script. Due to the deep `kernelshap` background bootstrapping calculations (dynamically locked to `min(520, nrow(X_train))` for balance between computation and precision), execution may take several minutes depending on parallel CPU power.
3. Retrieve all high-resolution vector figures in the `FinalPlots/` directory upon completion.

### 🔬 Methodological & Statistical Highlights
> **Targeted Interpretability Reconstruction & Granular Validation:**
>
> 1. **Locally Stripped Single-Gene Efficacy Validation**: Complex machine learning networks often evaluate the predictive limits of large Feature Panels. To meet clinical laboratory medicine's urgent need for single or dual-biomarker translation, this module hardcodes an independent non-parametric ROC extraction pipeline specifically for `MAFG` and `S100A8` at the end of the macro-ensemble model. This ensures practical diagnostic value is uncoupled and proven while maintaining group omics depth.
> 2. **KernelSHAP-Based Locked-Domain Attribution Analysis**: To entirely penetrate XGBoost's complex ensemble-tree black box, we bypassed crude global explainers, employing the computationally expensive but highly precise `KernelSHAP` algorithm to extract Shapley values. After rigorous average pooling of marginal contributions, the model "locks in" the Top 6 variables. We then deployed non-linear Dependence Plots and micro-push Waterfall Plots, granting pure computer algorithms profound molecular and pathophysiological interpretability from macro-cohorts down to micro-individuals.

## 📦 Module 08: Core Biomarker Independent External Cohort Validation & Advanced Performance Evaluation (`08_External Validation.R`)

### 🎯 Module Overview
This module represents the final hurdle in testing the clinical translational potential of machine-learning-discovered core diagnostic biomarkers (e.g., S100A8). The script builds a highly automated external validation pipeline containing an exceptionally robust microarray probe adaptive-matching engine and a regex-based clinical phenotype smart-parser. With zero selection bias explicitly guaranteed, the code outputs traditional ROC curves and forcefully introduces Precision-Recall (PR) curves, Bootstrap internal calibration curves, and Brier scores to comprehensively evaluate model prediction efficacy and calibration under the strictest SCI standards for imbalanced clinical scenarios.

### 🛠 Environment & Dependencies
Ensure the following core packages are installed:
* `GEOquery`, `limma` (For online GEO retrieval and matrix processing)
* `pROC`, `precrec`, `rms` (For advanced performance curves, AUPRC, and Logistic Bootstrap calibration)
* `ggplot2`, `ggpubr`, `RColorBrewer` (For cohort merging and publication-level rendering)

### 📥 Required Input Files
This module possesses extremely high automation and **requires NO local input files**.
1. The built-in `GEOquery` engine will automatically download the hardcoded external validation cohorts (`GSE9692`, `GSE26440`, `GSE28750`, `GSE69528`) from line 130. Please ensure a stable internet connection.

### ⚙️ Core Analytical Workflow
1. **Forced Hardcoded Loading**: Rejecting manual cohort picking, the script bulk-fetches the four pre-defined GEO external validation cohorts.
2. **Robust Probe Finder**: Targeting different sequencing platforms (GPL), the code first matches against hardcoded common probe libraries. If missed, it automatically retrieves feature annotations (fData) to perform regex fuzzy matching, precisely locating the target gene's (S100A8) expression signal.
3. **Smart Clinical Phenotype Parsing**: Reads metadata (pData) and utilizes a built-in regex keyword library (e.g., `sepsis`, `healthy`, `sirs`) to automatically strip irrelevant or confounding samples (e.g., uninfected SIRS patients), strictly establishing pure `Normal` vs `Sepsis` groups.
4. **Efficacy Calculation & Multi-Dimensional Evaluation**:
   * *Baseline:* Calculates non-parametric P-values (Wilcoxon) and traditional ROC curves (AUC) per dataset.
   * *Advanced:* Computes AUPRC using `precrec`; constructs binomial logistic regression with `rms` and resamples (Bootstrap = 400) to plot Calibration Curves; calculates Brier Scores for prediction probability mean squared errors.

### 📤 Complete Output Details
**🖼️ Figures (Generated in the current working directory):**
* `Boxplot_[GSE]_S100A8.pdf`: Differential distribution boxplots of target gene expression across independent cohorts (with P-values).
* `Combined_ROC_Validation_Final.pdf`: Traditional Receiver Operating Characteristic (ROC) curve set across four cohorts.
* `Combined_PR_Validation_Final.pdf`: Precision-Recall (PR) curve set with marked AUPRC values.
* `Combined_Calibration_Final.pdf`: Bootstrap Calibration curves comparing model predictive probabilities with actual incidence rates.
* `Combined_Brier_Final.pdf`: Bar plot of Brier scores quantifying predictive precision across cohorts.

### 👣 Detailed Reproduction Guide
1. Ensure a stable network environment for real-time NCBI Series Matrix downloads.
2. Modify `work_dir` on line 10 to your desired local absolute path for saving figures.
3. Run the script. The code automatically catches and bypasses single-cohort download crashes (`tryCatch`), ensuring the macro-loop reaches completion and auto-renders all evaluation charts.

### 🔬 Methodological & Statistical Highlights
> **Zero Selection Bias & Imbalance-Aware Advanced Metrics System:**
> 
> 1. **Zero Selection Bias with Hardcoded Cohorts**: During multi-omics feature validation, studies often quietly discard poorly performing validation cohorts. To demonstrate ultimate academic transparency and model robustness, this script strictly "hardcodes" the target cohorts (`GSE9692`, `GSE26440`, `GSE28750`, `GSE69528`) at the very frontend. Coupled with underlying `tryCatch` mechanisms, the model faithfully parses and reports evaluation results entirely, completely eradicating "data snooping" and selection bias.
> 2. **Advanced Clinical Metrics (PR, Calibration, Brier) for Imbalanced Samples**: Traditional ROC curves (AUC) often present overly optimistic illusions when facing extremely imbalanced positive/negative ratios in real-world clinical settings. This module forcefully introduces three-dimensional advanced metrics:
>    * **PR Curves (AUPRC)**: Applies stricter penalties for minority class (Sepsis) prediction failures in the tradeoff between recall and precision, reflecting true positive-capture capabilities.
>    * **Bootstrap Calibration Curves**: Through 400 resamples, evaluates whether the biomarker's "absolute risk probability" closely aligns with the true "clinical incidence rate" (diagonal fit).
>    * **Brier Scores**: Calculates the mean squared error between predictive probabilities and true binary labels, flawlessly quantifying not just discrimination, but the clinical predictive calibration of the biomarker.

## 📦 Module 09: Single-Cell Multi-Omics Integration & Virtual Knockout Loop (`09_scRNA-seq & Bulk Immune Infiltration.R`)

### 🎯 Module Overview
This module represents the pinnacle of mechanistic investigation and multi-omics validation in this study. The script fully processes scRNA-seq sequencing data (GSE167363) from ground-up QC and batch correction to fine cellular annotation. Crucially, it utilizes the `scTenifoldKnk` algorithm to execute an in silico perturbation (virtual knockout) of S100A8 at the single-cell level, unveiling its downstream cascade targets. Furthermore, it shatters traditional Bulk immune infiltration feature-set barriers by innovatively purifying monocyte feature matrices directly from scRNA-seq results, driving CIBERSORT to complete a flawless multi-omics cross-validation loop.

### 🛠 Environment & Dependencies
* **Core Framework**: `Seurat` (V5 compatible), `tidyverse`, `Matrix`, `data.table`
* **QC & Dim. Reduction**: `DoubletFinder` (doublet removal), `harmony` (batch correction), `clustree`
* **Advanced Analysis**: `scTenifoldKnk` (Virtual KO), `decoupleR` & `OmnipathR` (TF activity), `clusterProfiler` (Enrichment), `patchwork`, `ggpubr`
* **Validation & Correlation**: `CIBERSORT.R` (underlying script), `pheatmap`, `ComplexHeatmap`, `linkET`, `corrplot`

### 📥 Required Input Files
1. **Raw scRNA Matrix**: The `GSE167363_RAW/` directory must contain `matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz` files for all samples.
2. **Bulk Validation Arrays**: The `GSE65682_gene.csv` and `GSE65682_Groups.csv` generated from earlier modules.
3. **Algorithm Dependency**: Standard `CIBERSORT.R` script must be present in the working directory. (Note: Due to Stanford University copyright, researchers must independently request academic authorization, download the source code, rename it to `CIBERSORT.R`, and place it in the directory).

### ⚙️ Core Analytical Workflow
1. **Deep QC & Harmony Integration**: Merges 10X Genomics sparse matrices, executing strict MT/RBC filtering. Applies `DoubletFinder` to strip doublets, regresses cell-cycle variation via `SCTransform`, and invokes `Harmony` to eliminate patient-level batch effects.
2. **Targeted Dim. Reduction & Annotation**: Clusters cells based on high-dimensional spaces, localizing Monocytes, T/B/NK lineages via multi-gene marker joint-validation, stripping background noise (DCs/Platelets), and outputting pristine UMAP panoramas.
3. **Subpopulation Contrast & K-means TF Inference**: Subdivides disease monocytes into S100A8 High/Low sub-groups (Top 33% vs Bottom 33%) for differential GO enrichment. Invokes `decoupleR` to calculate core TF activities, introducing rigid K-means clustering during UMAP re-projection to forcefully trim outliers for extremely clear activity projections.
4. **scTenifoldKnk Virtual KO (In Silico Perturbation)**: Extracts disease monocytes, cleans ribosome/MT noise, and constructs non-linear manifold networks. Executes a single-thread virtual knockout of S100A8, calculating manifold distance perturbations (Z-scores) for all genes pre- and post-knockout to render clean GO/KEGG mechanistic enrichment charts.
5. **Multi-Omics CIBERSORT Loop**: Discarding public immune matrices (e.g., LM22), the algorithm isolates the disease monocyte cluster based on current single-cell clustering and S100A8 expression (Top 40%), forcefully purifying an absolutely specific "Custom Single-Cell CIBERSORT Feature Matrix". This is fed back into the underlying Bulk matrix executing 1,000 permutation tests, rendering publication-grade horizontal/vertical reference heatmaps, custom stacked bar plots, and multi-dimensional correlation networks.

### 📤 Complete Output Details
**📊 Core Data Snapshots (RDS & CSV):**
* `Sepsis_GSE167363_Prepared.rds` & `Sepsis_GSE167363_Annotated_Full.rds`.
* `Sepsis_GSE167363_Final_Clean.rds`: Final high-order Seurat object stripped of interference clusters (DCs/Platelets).
* `Result_Analysis4_S100A8_KO.csv`: Genome-wide manifold perturbation distance list post-virtual KO.
* `Sepsis_Custom_Reference_Matrix.csv`: Custom immune feature reference matrix extracted from scRNA-seq.
* `Final_Deconvolution_Results.csv` & `Figure_5D_Correlation_Stats.csv`.

**🖼️ Core Mechanistic & Infiltration Atlas (Plot 01-09 & Figures 4-5):**
* **scRNA Panorama**: UMAPs, Multi-gene DotPlots, FeaturePlots, Violin Plots.
* **Subpopulation & TF Activity**: Grouping Boxplots, Symmetric Bi-directional GO Barplots, TF Activity FeaturePlots (with/without axes).
* **Virtual KO Results (Fig 4)**: Manifold Distance Scatter Plots (`Figure_4I_KO_Distance_Scatter.pdf`), Targeted Distance Barplots (excluding S100A9 heterodimer interference), GO/KEGG Enrichment.
* **Custom Immune Infiltration (Fig 5)**: Minimalist Horizontal/Vertical Reference Heatmaps (`Figure_5A`), Re-scaled Patchwork StackPlots (`Figure_5B`), Pan-cellular Proportion Boxplots (`Figure_5C`), Sepsis-Specific Immune Interaction Corrplots (`Figure_5E`).

### 👣 Detailed Reproduction Guide
1. Ensure the massive raw sequencing data is correctly placed in `GSE167363_RAW` and Bulk matrices are ready.
2. Due to massive tensor calculations and 1,000 deconvolution permutations, a single workstation may require **2~4 hours** to complete execution.
3. Memory peaks may exceed 32GB (especially during manifold tensor construction); ensure sufficient computational buffering. All figures and core data dictionaries will auto-aggregate upon completion.

### 🔬 Methodological & Statistical Highlights
> **In Silico Perturbation Technology & Breaking LM22 Multi-Omics Barriers:**
> 
> 1. **`scTenifoldKnk` In Silico Knockout**: Traditional mechanistic inference is often hindered by the ultra-long cycles of wet-lab animal models. This module introduces cutting-edge Manifold Tensor Analysis. By constructing high-dimensional non-linear single-cell transcriptional regulatory networks, the algorithm forcefully zeros out S100A8 expression in a silicon-based environment, simulating and calculating the collapse and response of the entire transcriptome topology. This strategy of locking downstream targets via Distance Perturbation provides an entirely new *in silico* reductionist approach for biological mechanism exploration.
> 2. **Breaking Traditional Feature Matrix Limits**: Most immune infiltration deconvolutions heavily rely on pre-defined databases (e.g., LM22). However, these public matrices were not constructed under severe systemic inflammation (sepsis) environments, suffering from massive baseline drift. This study completely **breaks this limitation**, innovatively extracting true environment-specific cell markers directly from our self-purified single-cell population via strict `logfc.threshold` metrics. Dynamically assembling this highly adaptive `Sepsis_Custom_Reference_Matrix` and reverse-projecting it into the Bulk transcriptome cohorts for proportional deduction flawlessly bridges single-cell high-resolution with Bulk high-statistical-power, achieving an irrefutable validation loop.

## 📦 Module 10: Core Target Validation & Molecular Dynamics Simulation Loop (`10_Molecular_Docking_&_MD_Simulation.sh`)

### 🎯 Module Overview
This module bridges the macroscopic multi-omics screening with microscopic atomic-level validation. After pinpointing S100A8 as the core driver in the sepsis monocyte immune microenvironment, this script utilizes Computer-Aided Drug Design (CADD) to evaluate the feasibility of targeted intervention by the natural small molecule 16-hydroxytriptolide. The pipeline deeply integrates molecular docking preprocessing with up to 100 ns of all-atom Molecular Dynamics (MD) simulations, providing formidable physicochemical evidence for our prognostic model and potential clinical translation through dimensions of dynamic conformational evolution, thermodynamic stability, and free energy landscapes.

### 🛠 Environment & Dependencies
* **Core Engine**: `GROMACS` (V2023.2, GPU CUDA acceleration highly recommended)
* **OS**: Bash script. Windows users must configure the runtime environment via **WSL 2 (Ubuntu distribution recommended)** or submit jobs to a Linux High-Performance Computing (HPC) cluster.
* **Docking & Electronic Analysis**: `AutoDock Vina` (V1.2.2), `Multiwfn`, `VMD`
* **Force Field Systems**: `AMBER14SB` (Receptor All-Atom Force Field) / `GAFF` (General Amber Force Field for Ligands)
* **Advanced Vis. & Dim. Reduction Plotting**: `QtGrace` (V2.6), Custom Python 3.1.3 script (`xpm2png.py`) 
* **Note: Massive high-precision continuous trajectory files (.xtc/.trr, >10GB total) and core processed data packages have been hosted on Zenodo (DOI: 10.5281/zenodo.19480952). Researchers can access [https://doi.org/10.5281/zenodo.19480952](https://doi.org/10.5281/zenodo.19480952) for download or pull directly via wget.**

### 📥 Required Input Files
1. **Docking Target & Ligand**: S100A8 high-resolution crystal structure (PDB: `4GGF`) and 3D conformation of 16-hydroxytriptolide (PubChem CID: `107985`). Grid Box coordinate configurations are detailed in the accompanying `config.txt`.
2. **MD Initial System**: The assembled initial coordinate file `merge_gmx.gro` and full system topology files incorporating AMBER/GAFF parameters (`merge_topol.top` and `lig.itp`).
3. **Algorithm Dependency Parameters**: 5 core control scripts located in the `MDP_Files/` directory (`ref.em_steep.mdp`, `ref.em_cg.mdp`, `ref.nvt.mdp`, `ref.npt.mdp`, `ref.md.mdp`) and the custom grouping index `index.ndx`.

### ⚙️ Core Analytical Workflow
1. **Multi-Algorithm Energy Minimization**: Solvates the system in a 1.0 nm TIP3P cubic water box and neutralizes charges. Subsequently, sequentially invokes Steepest Descent and Conjugate Gradient methods for rigorous dual energy minimization, eliminating steric hindrances and unreasonable atomic contacts from initial docking conformations.
2. **Rigid & Flexible Equilibration**: Applies position restraints to protein backbones and ligands, executing NVT (constant volume/temperature) and NPT (constant pressure/temperature) pre-equilibration at 298 K and 1 bar standard physiological conditions, ensuring smooth transitions of system density and thermodynamic properties.
3. **100 ns Unrestrained Production MD**: Removes all position restraints, applies the LINCS algorithm to fix high-frequency bond length vibrations, and executes 100 ns of all-atom production dynamics simulation.
4. **Advanced Trajectory Cleaning**: To address Periodic Boundary Condition (PBC) artifacts, implements a 5-step trajectory correction method: cross-boundary molecule splicing -> forced clustering to prevent flying -> removing coordinate jumps -> center translation -> backbone rotational fitting, generating a research-grade pristine trajectory (`FINAL_CLEAN.xtc`).
5. **Conformational Dynamics & Thermodynamics Extraction**: Concurrently extracts RMSD, RMSF, Radius of Gyration (Rg), Intermolecular Hydrogen Bond Count, and Buried SASA based on the cleaned trajectory. Furthermore, utilizes PCA dimensionality reduction to extract core motion modes and construct a Gibbs Free Energy Landscape.

### 📤 Complete Output Details
**📊 Core Data Snapshots (.xtc, .tpr, .xvg):**
* `FINAL_CLEAN.xtc`: Final 100 ns continuous analysis trajectory completely devoid of PBC artifacts.
* `pic.rmsd.xvg` / `pic.rmsf.xvg`: Raw matrices containing per-picosecond (ps) backbone deviations and residue flexibility fluctuations.
* `pic.hbnum.xvg` / `pic.area.xvg`: Dynamic hydrogen bond occupancy lists and Solvent Accessible Surface Area statistics.
* `eigenvalues.xvg` / `2dproj.xvg`: PCA-reduced eigenvector extractions and 2D spatial projection matrices.

**🖼️ Core Mechanisms & Binding Atlas (Figures 6 & 7):**
* **System Stability Panorama**: `RMSD_Backbone_Fluctuation.pdf`, `RMSF_Residue_Flexibility.pdf` (Precisely localizing rigidity enhancement regions induced by small molecule binding).
* **Binding Compactness Dissection**: `Radius_of_Gyration_Compactness.pdf`, `Intermolecular_H-bonds_Count.pdf`.
* **Free Energy Landscape (FEL) (Fig 7)**: `gibbs.png` (2D free energy basin rendering generated by Python script, displaying the stable binding dominant conformation with the lowest global energy).

### 👣 Detailed Reproduction Guide
1. Create corresponding folder directories. Windows users should operate within a WSL Linux terminal and ensure GROMACS environment variables and CUDA acceleration support are correctly configured.
2. 100 ns all-atom calculations demand extreme computational power. A single consumer-grade GPU (e.g., RTX 3090/4090) may require **24~48 hours** of continuous operation.
3. Strictly follow the execution sequence in the `10_Molecular_Docking_&_MD_Simulation.sh` script. Skipping the advanced trajectory correction steps in phase 6 is strictly prohibited, as it will inevitably cause subsequent PCA and FEL analyses to crash or distort.

### 🔬 Methodological & Statistical Highlights
> **Breaking Static Docking Limitations via Atomic-Level Validation & Non-Linear Free Energy Landscapes:**
> 
> 1. **Cross-Dimensional Validation from Macro-Omics to Micro-Atomic Dynamics**: If the previous scRNA-seq and Bulk joint infiltration analyses inferred mechanisms at the "tissue and cellular network" level, this module provides profound validation from a microscopic perspective. We abandoned static molecular docking limitations, directly dropping the screened target (S100A8) complex into a 100 ns dynamic sandbox simulating realistic physiological solvent environments. Capturing true physical atomic interactions provides the most convincing pharmacodynamic computational evidence for the massive upstream machine-learning screening.
> 2. **PCA-Based Free Energy Landscape (FEL) Revealing Binding States**: This module introduces Gibbs Free Energy Landscape calculations based on Principal Component Analysis (PCA). The algorithm extracts extremely complex high-dimensional protein conformational changes, dimensionally-reducing and projecting them into a low-dimensional phase space formed by PC1 and PC2, calculating probability densities based on the Boltzmann distribution. The resulting FEL precisely pinpoints the "Global Minimum Basin" the system falls into during ligand binding, rigorously proving the high stability and specificity of the binding between 16-hydroxytriptolide and S100A8.
