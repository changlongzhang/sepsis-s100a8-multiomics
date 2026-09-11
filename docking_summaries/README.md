# S100A8 matched comparative docking

This directory contains a non-destructive matched comparison of triptolide (PubChem CID 107985) and 16-hydroxytriptolide (PubChem CID 126556). The existing triptolide docking project was not overwritten.

## Start here

1. `RESULTS_INTERPRETATION_CN.md` — Chinese interpretation and replacement decision.
2. `MANUSCRIPT_REPLACEMENT_TEXT_EN.md` — Methods, Results and reviewer-response wording.
3. `results/comparative_analysis/matched_replicate_statistics.csv` — matched score and pose-convergence summary.
4. `results/comparative_analysis/representative_pose_cross_compound_comparison.csv` — common-scaffold RMSD and contact overlap.
5. `results/comparative_analysis/matched_blind_tiling_comparison.csv` — complete blind-search comparison.
6. `input/prepared/compound_comparison_metadata.json` — compound identity, physicochemical properties and similarity audit.

## Scope boundary

No molecular dynamics or wet-lab experiment was performed for 16-hydroxytriptolide. This package supports predicted shared-pocket compatibility only; it does not establish functional equivalence or permit triptolide experiments to be described as validation of 16-hydroxytriptolide.
