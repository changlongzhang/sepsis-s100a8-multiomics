# Zenodo v2 archive layout

| Upload file | Purpose |
|---|---|
| `01_PONE-D-26-24904_v2_core_code_source_data_and_final_SI.zip` | GitHub snapshot, final S1 raw images, S1-S7 Tables, S1-S8 Figures, S2-S10 code/source-data files |
| `02_PONE-D-26-24904_v2_transcriptomic_input_snapshots_and_caches.zip` | Exact GSE65682 inputs, four original external matrices, added clinical datasets, numerical caches, and three WGCNA branches |
| `03_PONE-D-26-24904_v2_GSE167363_frozen_analysis_object.zip` | Frozen GSE167363 Seurat object and donor-aware derived objects/results; includes the downloaded Reyes object used in the revision workspace |
| `04_PONE-D-26-24904_v2_docking_MD_inputs_parameters_and_summaries.zip` | PDB 5HLO route, docking inputs/results, MD topology/parameters/MDP/scripts, TPR/NDX files, and numerical MD summaries |
| `MD_apo_rep1_100ns_canonical_trajectory.zip` to `MD_apo_rep3_...` | Three independent apo canonical trajectories |
| `MD_complex_rep1_100ns_canonical_trajectory.zip` to `MD_complex_rep3_...` | Three independent S100A8-triptolide canonical trajectories |

Every upload file has a SHA-256 entry in `SHA256SUMS.txt` and a row in `ZENODO_UPLOAD_MANIFEST.csv`. The earlier Zenodo `data.zip` is retained as an immutable historical version but is not an input to the final revised analysis.
