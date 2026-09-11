# Complete figure quality-control report

- Revision: the original 18-candidate S100A8 computational-screening matrix was added as panel A. It was compacted from 18 x 1 to two 9-row blocks solely for legibility; all 18 identities, their original order and all binary values were preserved in `source_data/S100A8_compound_screening_binary.csv`.
- Target highlight: 16-hydroxytriptolide is labelled as 16-OH-TPL in panel A and distinguished by bold text plus a purple keyline; its source value remains 1.
- Panel map: A-T. Former panels A-S were advanced to B-T without changing their underlying docking or molecular-dynamics data.
- Visual review: passed after iterative high-resolution PDF rendering. No clipped labels, overlapping text, broken axes, missing panels or displaced legends were observed in the final render.
- Layout revision: the existing five-row/four-column arrangement was retained. Panels A and B were extended vertically. Panel A uses wider separation between its two drug blocks, narrower binary tiles and the display label `PKCβ inhibitor`; panel B uses a vertical two-structure arrangement to increase molecular-structure size and separate its labels. The bottom row of A, the bottom similarity label of B and the legend baseline of C were visually aligned.
- Panel-label alignment: A-D share one y baseline and use the exact same four column-start coordinates as E-H, I-L, M-P and Q-T.
- Typography: Arial throughout; panel labels are bold uppercase to match the manuscript's established figure convention. PDF fonts are embedded as Arial TrueType-compatible Type 0/Type 42 output.
- SVG editability: 233 live text nodes are retained; text is not converted to outlines.
- PDF: one page; MediaBox 576 x 680.4 pt. The vector PDF was directly rendered at high resolution for final visual inspection.
- TIFF: 4800 x 5670 pixels, RGB, 600 x 600 dpi, LZW compression, no alpha channel.
- AI: PDF-compatible Illustrator working copy generated from the final vector PDF.
- Screening source audit: 18 rows, 14 positive values, and 16-HYDROXYTRIPTOLIDE = 1.
- Titles: no overall title and no narrative plot-panel titles.
- MD style: all quantitative MD axes retain complete four-sided boxes with ticks on all sides.
- Replicate display: G-L show the median and full range across three independent replicas; M-Q and R show complex replicas separately.
- Scientific scope: docking includes triptolide and 16-hydroxytriptolide; MD includes triptolide complex and apo controls only.
- Data integrity: only accepted v3 PBC-corrected MD source tables were used; no failed, superseded or extrapolated trajectories were included.
- Interpretation guardrail: the figure does not imply MD or experimental validation of 16-hydroxytriptolide.
