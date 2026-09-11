# Execution guide

## Scope

The numbered filenames record the revision workflow. Scripts `02`-`24` contain
the main data QC, model, external-validation, single-cell/pseudobulk, integration,
and table-generation analyses. Scripts `27`-`30` contain added external and
mechanistic triangulation analyses. Scripts `31` onward mainly assemble or audit
submission figures. Unnumbered scripts are document, workbook, or Illustrator
utilities and are not part of the primary statistical workflow.

## Software families

- R/Bioconductor: R, GEOquery, Biobase, limma, edgeR, Seurat,
  SingleCellExperiment, SummarizedExperiment, Matrix, lme4, GSVA, UCell,
  decoupleR, dorothea, AnnotationDbi, org.Hs.eg.db, clusterProfiler, pROC,
  caret, glmnet, randomForest, xgboost, mboost, metafor, data.table, dplyr,
  tidyr, readr, stringr, ggplot2, ggrepel, patchwork, svglite, and ragg as
  called by the individual scripts.
- Python: numpy, pandas, scipy, statsmodels, matplotlib, Pillow, PyMuPDF, and
  python-docx as called by the individual scripts.
- Optional production tools: Node.js with `@oai/artifact-tool` for workbook
  production and Adobe Illustrator for `.jsx` export scripts.

Exact package versions from the author environment are reported in the software
and environment table supplied with the manuscript. Install only the packages
required by the script being evaluated.

## Input and run order

1. Download the public datasets cited in the manuscript and preserve their
   accession-level sample metadata.
2. Set a local project root and update input/output roots at the top of the
   script you intend to run. The portable pattern is:

   - Python: `PROJECT_ROOT = Path(os.environ["SEPSIS_PROJECT_ROOT"])`
   - R: `project_root <- Sys.getenv("SEPSIS_PROJECT_ROOT")`

3. Run the principal scripts in numeric order only for the analysis branch of
   interest. Later scripts consume saved objects or tables produced by earlier
   scripts; they are not all required for every figure.
4. Compare regenerated numerical outputs with `tables/`,
   `figure_source_data/`, and the analysis-specific folders in the archive.

## Legacy absolute paths

Some figure assembly, Illustrator, audit, and document-production scripts retain
the original Windows path `D:\\桌面\\sepsis` (or its forward-slash equivalent).
This is deliberate provenance, not a hidden data dependency. Before running one
of those scripts on another computer, replace only that root with the local
unpacked project root and keep the remaining relative path unchanged. Do not
change data values, group labels, statistics, or panel selections during path
migration.

The primary numerical tables are included so that a reader can inspect the
reported values without Adobe Illustrator or the author's workstation layout.

