# In vitro source-data provenance

`Fig7_Fig8_quantitative_source_data_original.csv` is the unchanged author-provided table. In its CCK-8 panel, the `TPL_nM` cells for biological repeats 2 and 3 contain normalized measurement-like values rather than dose labels. `Fig7_Fig8_quantitative_source_data.csv` preserves every original row and value, retains the original field as `source_TPL_nM`, and reconstructs the intended ordered dose labels (0, 10, 25, 50, and 100 nM) within each LPS group and repeat. Missing measurements remain missing; none were imputed.

`Western_blot_raw_densitometry.xlsx` is the archived densitometry workbook with a portable filename. The spelling-only label correction `Tritolide` → `Triptolide` was applied; numerical cells and workbook structure were not changed.

`CETSA_normalized_densitometry.csv` contains the available normalized temperature-gradient values. `CETSA_plot.R` reproduces the archived plot from those values.

RT-qPCR material contains relative-expression values, but the project archive did not contain sample-level Ct or ΔCt tables. No Ct/ΔCt values were inferred or reconstructed.
