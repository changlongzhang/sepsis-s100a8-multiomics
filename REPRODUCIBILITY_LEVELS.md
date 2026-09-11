# Reproducibility levels and boundaries

## Level A — raw public transcriptomic data

GEO/SCP accession identifiers and exact cohort rules are listed in `DATASET_MANIFEST.csv` and S1 Table. Raw third-party data remain available from their source repositories. The Zenodo v2 transcriptomic archive also preserves the exact downloaded matrices, annotations, metadata, and analysis caches used for the final revision, so later repository updates do not silently change the analysis input.

## Level B — frozen analysis objects

The Zenodo v2 single-cell archive supplies `Sepsis_GSE167363_Final_Clean.rds`, the exact object consumed by the final donor-aware scripts, plus donor pseudobulk and state-analysis outputs. The raw GSE167363 files are public at GEO; regeneration of the frozen Seurat object requires the original QC/integration workflow and compatible package versions. The frozen object therefore provides the shortest faithful route to the reported donor-level results.

The maximum-probe and mean-probe WGCNA RDS objects are supplied losslessly as PLOS S3-S10 files and in the Zenodo core archive. Use `reassemble_large_rds.py` and verify the SHA-256 values in `large_rds_manifest.json` before loading them.

## Level C — numerical result verification

All main and supplementary figure-level numerical source data, final Supporting Tables, clinical-validation summaries, class-balance outputs, donor-level single-cell results, in vitro data, docking summaries, and MD summaries are deposited. These files permit independent checking of manuscript values without proprietary figure software.

## Level D — computational reruns

The repository includes R, Python, PowerShell, Node, and Illustrator-provenance scripts. Primary numerical scripts use the deposited directory layout and the `SEPSIS_PROJECT_ROOT` environment variable. Figure-assembly scripts that still name the author's workstation are retained as provenance and must have only their root path changed before use. Adobe Illustrator scripts are optional and are not needed to verify numerical claims.

The final MD route uses PDB 5HLO, structural calcium, three independently seeded apo trajectories, and three independently seeded S100A8-triptolide trajectories. Inputs, MDP files, ligand parameters, scripts, TPR/NDX files, canonical trajectories, and summary tables are deposited. Rerunning six 100-ns trajectories requires GROMACS 2023.2 and substantial compute; exact floating-point identity can depend on hardware and parallel execution, while summary-level scientific comparison should remain reproducible.

## Non-reconstructable source boundary

The available qPCR archive lacks sample-level Ct/ΔCt values. No Ct values were reconstructed from relative expression. The available blot archive contains full-field exports and crop mappings, but proprietary instrument-native acquisition files were not present. These limitations are disclosed rather than filled with inferred data.
