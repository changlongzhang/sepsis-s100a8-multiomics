import { FileBlob, SpreadsheetFile } from "@oai/artifact-tool";

const inputs = process.argv.slice(2);
if (!inputs.length) throw new Error("Provide one or more workbook paths.");

for (const inputPath of inputs) {
  const workbook = await SpreadsheetFile.importXlsx(await FileBlob.load(inputPath));
  const sheets = await workbook.inspect({
    kind: "workbook,sheet,table,region",
    maxChars: 18000,
    tableMaxRows: 40,
    tableMaxCols: 20,
    tableMaxCellChars: 80,
  });
  process.stdout.write(`\n=== ${inputPath} ===\n${sheets.ndjson}\n`);
}
