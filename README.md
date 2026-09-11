# sepsis-s100a8-multiomics

Code, numerical source data, and reproducibility documentation for the revised manuscript:

**Multi-cohort transcriptomics and donor-aware single-cell analyses prioritize S100A8 within an emergency-myeloid, HLA-DR-low state in sepsis**

Manuscript ID: **PONE-D-26-24904**

Authors: Changlong Zhang, Qianqian Cui, Qiyun Wu, Wei Tang, Liang Chu, and Lei Zhang.

## Scope

This repository is the release candidate aligned with the final revised manuscript. It contains analysis scripts, numerical source data, intermediate summary tables, and execution notes. It is an auditable research record; it is not a single-command workflow and it does not redistribute the large raw public datasets.

The study integrates whole-blood transcriptomic discovery and external validation, leakage-controlled feature prioritization, donor-aware single-cell analyses, bulk–single-cell integration, in vitro perturbation experiments, docking, and molecular-dynamics summaries. The results prioritize S100A8 within an emergency-myeloid, HLA-DR-low state; they do not establish a clinical diagnostic, prognostic biomarker, direct drug target, or direct triptolide–S100A8 binding.

## Public datasets

Raw human transcriptomic data are available under GEO accessions [GSE65682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65682), [GSE9692](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9692), [GSE26440](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE26440), [GSE28750](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28750), [GSE69528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69528), [GSE134347](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134347), [GSE131411](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131411), [GSE272769](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272769), [GSE54514](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54514), [GSE63042](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63042), [GSE110487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110487), [GSE167363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167363), [GSE205672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205672), and [SCP548](https://singlecell.broadinstitute.org/single_cell/study/SCP548).

See [DATASET_MANIFEST.csv](DATASET_MANIFEST.csv) for study roles and cohort definitions. The manuscript's complete cohort definitions are also reproduced in `tables/published_supporting/S1_Table_Cohorts.xlsx`.

## Repository map

- `scripts/`: primary analysis, validation, statistics, figure assembly, and audit scripts.
- `corrected_DEG_Fig1_S1/`: the final 265-gene intersection and Fig 1/S1 Fig source data.
- `wgcna_final_rerun/`: primary maximum-probe WGCNA results.
- `wgcna_50pct_MAD_sensitivity/` and `wgcna_meanprobe_sensitivity/`: prespecified sensitivity analyses.
- `clinical_validation/`: external diagnostic, comparator, prognosis, longitudinal, and treatment-response summaries.
- `single_cell_pseudobulk/` and `bulk_singlecell_integration/`: donor-aware single-cell and integrated evidence outputs.
- `class_imbalance_summaries/`: class-weighting, downsampling, and bootstrap analyses.
- `in_vitro_source_data/`: archived quantitative in vitro measurements, CETSA values, and provenance notes.
- `figure_source_data/`: machine-readable panel-level data and corrected repeat-blocked Fig 7/8 statistical outputs.
- `docking_summaries/` and `MD_summaries/`: docking and molecular-dynamics summaries.
- `tables/published_supporting/`: the seven final Supporting Tables.

Read [METHODS_SUMMARY.md](METHODS_SUMMARY.md), [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md), and `scripts/README_EXECUTION.md` before reuse.

## Reproducing the revised analysis

The accompanying Zenodo v2 deposit is split into ten checksummed ZIP files. Extract archives 01-04 into one common parent directory. Place this GitHub repository snapshot in that same parent directory, so the following sibling paths exist:

```text
release_root/
├── GitHub_Repository_Snapshot/
├── GSE65682/
├── sepsis-s100a8-multiomics/data/
├── review_revision/01_data_inventory/downloaded_data/
├── 单细胞测序/Sepsis_GSE167363_Final_Clean.rds
├── work/final_revision_20260908/
└── S100A8_triptolide_MD_topjournal/
```

Set `SEPSIS_PROJECT_ROOT` to `release_root`, make `GitHub_Repository_Snapshot` the working directory, and run only the analysis branch needed. The repository includes frozen numerical outputs so readers can verify the reported values without rerunning every high-cost model. Delete or move a cached output only when deliberately requesting a full recalculation.

- `REPRODUCIBILITY_LEVELS.md` distinguishes raw-public-data reruns, frozen-object reruns, numerical-output verification, and visual assembly provenance.
- `FIGURE_TABLE_REPRODUCIBILITY_MAP.csv` maps each main/supplementary figure and table to its input, script, and deposited output.
- `ZENODO_ARCHIVE_LAYOUT.md` identifies which Zenodo upload archive supplies each required input.
- `run_analysis.R` is the portable cached runner for numerical scripts 01-24. Structural simulations and final WGCNA reruns have separate commands because of their computational cost.

The six MD trajectory archives are not required to inspect the manuscript's deposited summary tables, but they permit independent re-analysis of the reported trajectories with the matching TPR/NDX files.

## Analysis module index

| Module | Final-revision purpose | Main entry points |
|---|---|---|
| 01 | Environment, inventory, and discovery-scale QC | `scripts/01_package_environment.R`, `02_discovery_data_qc.R` |
| 02 | Leakage-controlled modeling and class-imbalance sensitivity | `scripts/03_class_weighted_models.R` through `07_class_imbalance_summary.py` |
| 03 | Clinical comparators, calibration, decision curves, prognosis, and response | `scripts/08_external_metadata_audit.R` through `13_clinical_meta_analysis.R`; `27_clinical_extension.R` |
| 04 | Donor-aware single-cell pseudobulk and monocyte states | `scripts/14_scrna_metadata_qc.R` through `23_integrated_evidence_summary.R` |
| 05 | Independent GSE205672/SCP548 and mechanistic triangulation | `scripts/28_scp548_independent_validation.R` through `30_mechanistic_triangulation.R` |
| 06 | Formal WGCNA reruns and sensitivities | `final_revision_scripts/04_wgcna_allgenes_rerun.R`, `01_wgcna_final_rerun.R`, `20_wgcna_meanprobe_sensitivity.R` |
| 07 | Matched docking and six 100-ns MD trajectories | `structural_reproducibility/` |
| 08 | In vitro statistics and figure source data | `scripts/36_fig7_fig8_statistics.py`; `in_vitro_source_data/` |

Scripts numbered 31 onward mainly preserve figure assembly and quality-control provenance. They are not the authoritative input for scientific statistics when a machine-readable table is supplied.

## Reproducibility boundaries

Public datasets must be downloaded from their repositories. Large public raw data, licensed software, editable Adobe Illustrator files, and proprietary instrument-native blot acquisition files are not redistributed. Some legacy figure-assembly scripts retain author-workstation paths solely as provenance; replace those roots before execution.

RT-qPCR inference is based on archived relative-expression values. Sample-level Ct/ΔCt tables were not present in the project archive, so no Ct values were reconstructed and independent re-analysis on the preferred ΔCt scale is not possible from the archived material.

Author-provided full-field blot exports and crop mappings are supplied separately as PLOS `S1 Raw images`; no brightness, contrast, gamma, or background adjustment was applied in preparing that file. The instrument-native acquisition files were not present in the archive.

## Versioning and citation

This is the **v2.0.0 final-revision release candidate**. The existing Zenodo record [10.5281/zenodo.19480952](https://doi.org/10.5281/zenodo.19480952) is retained as the earlier archived version. The version-independent Zenodo concept DOI is [10.5281/zenodo.19480951](https://doi.org/10.5281/zenodo.19480951). A new version DOI must be inserted into `CITATION.cff`, `.zenodo.json`, and the GitHub release notes only after Zenodo mints it.

Code authored for this project is released under the repository's MIT license. The Zenodo license applies to the authors' original documentation, derived outputs, and packaging; redistributed GEO/SCP snapshots remain subject to the source repositories' terms and must be cited by accession. No new rights are asserted over third-party source data.
