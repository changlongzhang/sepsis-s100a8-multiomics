import fs from "node:fs/promises";
import path from "node:path";
import { Workbook, SpreadsheetFile } from "file:///C:/Users/ZCL/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules/@oai/artifact-tool/dist/artifact_tool.mjs";

const root = process.cwd();
const previewDir = path.join(root, "11_logs", "spreadsheet_previews");
await fs.mkdir(previewDir, { recursive: true });

const specs = [
  ["10_tables/Table_ClassBalance_Performance.xlsx", [["Performance","02_class_imbalance/Table_ClassBalance_Performance.csv"],["S100A8 stability","02_class_imbalance/Table_S100A8_FeatureStability.csv"],["Four-model common","02_class_imbalance/Table_S100A8_FourModelCommonFrequency.csv"]]],
  ["10_tables/Table_S100A8_FeatureStability.xlsx", [["Feature stability","02_class_imbalance/Table_S100A8_FeatureStability.csv"],["Empirical intervals","02_class_imbalance/S100A8_empirical_intervals.csv"]]],
  ["10_tables/sample_group_mapping_all_cohorts.xlsx", [["Sample mapping","03_clinical_validation/sample_group_mapping_all_cohorts.csv"]]],
  ["10_tables/clinical_dataset_screening.xlsx", [["Dataset screening","03_clinical_validation/clinical_dataset_screening.csv"]]],
  ["10_tables/Table_ClinicalValidation.xlsx", [["Validation metrics","03_clinical_validation/Table_ClinicalValidation_long.csv"],["DCA calibration","03_clinical_validation/Table_DCA_Calibration.csv"],["Meta analysis","03_clinical_validation/clinical_effect_random_effects_meta.csv"]]],
  ["10_tables/Table_Prognosis_Associations.xlsx", [["Associations","04_prognosis_severity/Table_Prognosis_Associations.csv"],["Mortality meta","04_prognosis_severity/prognosis_random_effects_meta.csv"],["Variable availability","04_prognosis_severity/severity_variable_availability.csv"]]],
  ["10_tables/Table_Biomarker_Comparison.xlsx", [["Availability","03_clinical_validation/Table_DCA_Calibration.csv"]]],
  ["10_tables/Table_DonorLevel_TFActivity.xlsx", [["TF activity","05_single_cell_pseudobulk/Table_DonorLevel_TFActivity.csv"],["Donor values","05_single_cell_pseudobulk/donor_level_TF_activity.csv"]]],
  ["13_response_to_reviewers/reviewer_evidence_matrix.xlsx", [["Evidence matrix","13_response_to_reviewers/reviewer_evidence_matrix.csv"]]],
  ["11_logs/run_status.xlsx", [["Run status","11_logs/run_status.csv"]]],
];

for (const [outRel, sheets] of specs) {
  const wb = Workbook.create();
  let first = true;
  for (const [sheetName, csvRel] of sheets) {
    const csvPath = path.join(root, csvRel);
    const csv = await fs.readFile(csvPath, "utf8");
    await wb.fromCSV(csv, { sheetName, renameFirstIfOnlyNewSpreadsheet: first });
    first = false;
    const ws = wb.worksheets.getItem(sheetName);
    ws.showGridLines = false;
    ws.freezePanes.freezeRows(1);
    const used = ws.getUsedRange();
    if (used) {
      used.format.font = { name: "Arial", size: 9, color: "#222222" };
      used.format.autofitColumns();
      used.format.columnWidth = 18;
      used.format.wrapText = true;
      used.format.autofitRows();
      const header = used.getRow(0);
      header.format = { fill: "#3D7EA6", font: { name: "Arial", size: 9, bold: true, color: "#FFFFFF" }, verticalAlignment: "center", wrapText: true, borders: { bottom: { style: "medium", color: "#2F617D" } } };
      header.format.rowHeight = 32;
    }
  }
  const outPath = path.join(root, outRel);
  await fs.mkdir(path.dirname(outPath), { recursive: true });
  const xlsx = await SpreadsheetFile.exportXlsx(wb);
  await xlsx.save(outPath);
  const firstSheet = sheets[0][0];
  // 仅预览首个工作表的前20行，避免样本映射表（数百行）生成超高位图。
  const preview = await wb.render({ sheetName: firstSheet, range: "A1:H20", scale: 0.8, format: "png" });
  const stem = path.basename(outRel, ".xlsx").replace(/[^A-Za-z0-9_-]/g, "_");
  await fs.writeFile(path.join(previewDir, `${stem}.png`), new Uint8Array(await preview.arrayBuffer()));
  const inspect = await wb.inspect({ kind: "sheet", include: "id,name", range: "A1:F8" });
  await fs.writeFile(path.join(previewDir, `${stem}_inspect.json`), JSON.stringify(inspect, null, 2), "utf8");
  const autoInspect = `${outPath}.inspect.ndjson`;
  try { await fs.rename(autoInspect, path.join(previewDir, `${stem}_artifact_inspect.ndjson`)); } catch {}
}
