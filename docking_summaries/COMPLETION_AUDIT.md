# Completion audit

- Completed: 2026-08-21 (Asia/Shanghai).
- Existing triptolide docking and molecular-dynamics project directories were not overwritten.
- 16-hydroxytriptolide source: PubChem CID 126556, downloaded through the PubChem PUG REST SDF endpoint.
- Downloaded SDF SHA256: `596B5D689CAA759B5DC47F3340AE814404FD94E18B91B409106F2F67EAE4A0B6`.
- Ligand preparation: RDKit 2025.09.5 and Meeko 0.7.1, deterministic ETKDGv3/MMFF94s conformer workflow, Gasteiger PDBQT charges.
- Docking engine: AutoDock Vina 1.2.7, executable SHA256 `E0C4B2715E0C1A74F6E92D0F3BE0328AC97542EAFBC111E6B1EFAD897A73CCE5`.
- Receptors, protonation, grids and local-pocket centres were reused from the accepted 5HLO triptolide project.
- Blind docking: 24/24 jobs returned code 0; no missing best score.
- Local docking: 40/40 jobs returned code 0; every job returned at least 19 modes and a valid best score.
- Expected local outputs: 40 PDBQT pose files and 40 logs; all present.
- Error scan: no `FATAL`, `ERROR`, `failed`, `NaN` or `exception` patterns in successful docking logs.
- Analysis completion marker: `results/comparative_analysis/MATCHED_DOCKING_ANALYSIS_COMPLETE.txt`.
- No molecular-dynamics or wet-lab workflow was created or run for 16-hydroxytriptolide.

## Recorded recoverable preparation issue

The first Meeko invocation attempted to pass an absolute Windows path containing Chinese characters to RDKit and failed with a bad-input-file error. The failure is preserved in `meeko_ligand_preparation.log`. The preparation was safely rerun from the ligand directory using ASCII relative filenames; the molecular input, stereochemistry and preparation parameters were unchanged.
