import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const root = path.resolve(process.env.FANLAB_RESEARCH2_ROOT || process.cwd());
const pkg = path.join(root, "final_manuscript_planning_20260605");
const outDir = path.join(pkg, "manuscript_tables");
const xlsxPath = path.join(outDir, "manuscript_figure_table_package.xlsx");

function csvEscape(value) {
  const s = value === null || value === undefined ? "" : String(value);
  return /[",\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

function toCsv(rows) {
  return rows.map((row) => row.map(csvEscape).join(",")).join("\r\n") + "\r\n";
}

function toMarkdown(headers, rows) {
  const line = (cells) => `| ${cells.map((x) => String(x ?? "").replace(/\|/g, "\\|")).join(" | ")} |`;
  return [
    line(headers),
    line(headers.map(() => "---")),
    ...rows.map((row) => line(row)),
  ].join("\n") + "\n";
}

function splitCsvLine(line) {
  const out = [];
  let cur = "";
  let inQuotes = false;
  for (let i = 0; i < line.length; i += 1) {
    const ch = line[i];
    if (ch === '"') {
      if (inQuotes && line[i + 1] === '"') {
        cur += '"';
        i += 1;
      } else {
        inQuotes = !inQuotes;
      }
    } else if (ch === "," && !inQuotes) {
      out.push(cur);
      cur = "";
    } else {
      cur += ch;
    }
  }
  out.push(cur);
  return out;
}

async function readCsv(file) {
  const txt = await fs.readFile(file, "utf8");
  const lines = txt.replace(/^\uFEFF/, "").trim().split(/\r?\n/);
  if (lines.length === 0 || !lines[0]) return [];
  const headers = splitCsvLine(lines[0]);
  return lines.slice(1).filter(Boolean).map((line) => {
    const cells = splitCsvLine(line);
    return Object.fromEntries(headers.map((h, i) => [h, cells[i] ?? ""]));
  });
}

async function writeTable(name, headers, rows) {
  await fs.writeFile(path.join(outDir, `${name}.csv`), toCsv([headers, ...rows]), "utf8");
  await fs.writeFile(path.join(outDir, `${name}.md`), toMarkdown(headers, rows), "utf8");
}

async function listFilesRecursively(dir) {
  const entries = await fs.readdir(dir, { withFileTypes: true });
  const files = [];
  for (const entry of entries) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) files.push(...await listFilesRecursively(full));
    else files.push(full);
  }
  return files;
}

function rel(file) {
  return path.relative(root, file).replace(/\\/g, "/");
}

function addSheet(workbook, name, headers, rows, widths = []) {
  const sheet = workbook.worksheets.add(name);
  const values = [headers, ...rows];
  sheet.getRangeByIndexes(0, 0, values.length, headers.length).values = values;
  const endCol = columnName(headers.length);
  const table = sheet.tables.add(`A1:${endCol}${values.length}`, true, name.replace(/[^A-Za-z0-9_]/g, "_").slice(0, 240));
  sheet.freezePanes.freezeRows(1);
  headers.forEach((_, i) => {
    const width = widths[i] ?? Math.min(260, Math.max(80, Math.max(...values.map((r) => String(r[i] ?? "").length)) * 6 + 24));
    sheet.getRangeByIndexes(0, i, values.length, 1).format.columnWidthPx = width;
  });
  sheet.getRangeByIndexes(0, 0, 1, headers.length).format.font.bold = true;
  sheet.getRangeByIndexes(0, 0, values.length, headers.length).format.wrapText = true;
  sheet.getRangeByIndexes(0, 0, 1, headers.length).format.rowHeightPx = 28;
  if (values.length > 1) {
    sheet.getRangeByIndexes(1, 0, values.length - 1, headers.length).format.rowHeightPx = 72;
  }
  return sheet;
}

function columnName(n) {
  let s = "";
  while (n > 0) {
    const rem = (n - 1) % 26;
    s = String.fromCharCode(65 + rem) + s;
    n = Math.floor((n - 1) / 26);
  }
  return s;
}

await fs.mkdir(outDir, { recursive: true });

const figureOrderHeaders = [
  "DisplayOrder",
  "ManuscriptPlacement",
  "FigureLabel",
  "FigureTitle",
  "PrimaryConclusion",
  "EvidenceRole",
  "AssetStatus",
  "RecommendedAction",
];

const figureOrderRows = [
  ["GA", "Graphical abstract", "Graphical abstract", "Association-based milk-gut-SCFA maturation model in early infancy", "Milk bioactives, infant gut succession, SCFA maturation, maternal gut context and exploratory clinical modifiers form an association-based model.", "Conceptual summary, not a result figure", "Conclusion-integrated PNG/SVG available", "Use as graphical abstract or first-page visual summary"],
  ["1", "Main text", "Fig. 1", "Qinghai-Tibetan Plateau mother-infant cohort and modular multi-omics design", "The study design supports modular longitudinal integration rather than strict complete-case all-omics modeling.", "Cohort/design/coverage anchor", "Stable PNG copied into figures_redrawn/fig1_cohort_design_overview", "Use as Fig. 1 after any final journal-size resizing"],
  ["2", "Main text", "Fig. 2", "Human milk bioactive component remodeling", "Lactoferrin and selected HMO/LCFA features show distinct early postnatal dynamics.", "Exposure-side longitudinal dynamics", "Redrawn PNG/PDF/SVG/TIFF and source CSV available", "Keep as main Fig. 2"],
  ["3", "Main text", "Fig. 3", "Infant gut microbiome maturation", "Infant fecal microbiome maturation is Bifidobacterium-centered, with B. longum as the most stable high-abundance species label.", "Infant microbiome maturation", "Redrawn PNG/PDF/SVG/TIFF and source CSV available", "Keep as main Fig. 3; retain species-label caution"],
  ["4", "Main text", "Fig. 4", "Fecal SCFA metabolic maturation", "Infant fecal SCFA profiles form an age-associated maturation axis from D05 to D90.", "Infant metabolome anchor", "Redrawn PNG/PDF/SVG/TIFF and source CSV available", "Keep as main Fig. 4"],
  ["5", "Main text", "Fig. 5", "Candidate lagged milk bioactive links to infant fecal SCFAs", "Early milk LTF and selected HMO features show candidate lagged associations with later infant fecal SCFA outcomes.", "Core integration figure", "Redrawn PNG/PDF/SVG/TIFF and source CSV available", "Keep as final main result Fig. 5; avoid causal wording"],
  ["S1", "Supplementary", "Fig. S1", "Maternal fecal metagenomic context", "Mother-infant fecal paired-similarity signals are sparse and contextual rather than proof of vertical transmission.", "Contextual dyad ecology", "Redrawn Fig. 6 package available", "Move current Fig. 6 to Supplementary Fig. S1 or Extended Data"],
  ["S2", "Supplementary", "Fig. S2", "Clinical modifier evidence audit", "Infant age remains the dominant correlate of SCFA maturation; maternal BMI-by-age is nominal; feeding/delivery/birthweight are unclear.", "Exploratory clinical modifier audit", "Redrawn Fig. 7 package available", "Move current Fig. 7 to Supplementary Fig. S2 or Extended Data"],
];

const suppFigureHeaders = ["SupplementaryLabel", "FormerLabel", "Title", "ReasonForSupplementaryPlacement", "MainTextReference"];
const suppFigureRows = [
  ["Fig. S1", "Fig. 6", "Maternal fecal metagenomic context", "Useful context but sparse own-dyad similarity and no strain-level transmission evidence.", "Mention briefly after Fig. 5 as maternal gut ecological background."],
  ["Fig. S2", "Fig. 7", "Clinical modifier evidence audit", "Clinical covariate coverage is incomplete; BMI-by-age is nominal and other modifiers are not clear.", "Mention in a short exploratory clinical modifiers paragraph."],
];

const mainTableHeaders = ["Module", "RowsOrSamples", "Dyads", "Timepoints", "ManuscriptUse", "SourcePath", "CoverageInterpretation"];
const mainTableRows = [
  ["Baseline clinical data", "52 dyads", "52", "Baseline/prenatal", "Cohort description and modifier availability", "Derived from final planning report and clinical modifier matrix", "Baseline clinical information anchors the cohort but is not complete for every downstream model."],
  ["Human milk bioactive module", "324 records", "60", "D05, D14, D30, D90", "Exposure-side milk dynamics and lagged integration", "p2_bridge_analyses/tables/p2_milk_module_clinical_matrix.csv", "Milk coverage is broadest for lactoferrin; HMO and LCFA coverage are window-specific."],
  ["Infant fecal SCFA module", "91 records", "32", "D05, D14, D30, D90", "SCFA maturation score and infant metabolome trends", "p0_multiomics_extension/tables/p0_scfa_maturation_score.csv", "Sufficient for age-associated SCFA maturation, but the score remains exploratory."],
  ["Infant fecal metagenome module", "73 mother-infant comparison records; metagenome coverage varies by window", "29", "D05, D14, D30, D90", "Infant gut microbiome maturation and maternal fecal context", "p0_multiomics_extension/tables/p0_mf_bf_own_vs_other_similarity_summary.csv", "D90 metagenomic coverage is limited; avoid over-reading late species patterns."],
  ["Bifidobacterium-SCFA matched module", "19 matched records", "14", "D05, D14, D30, D90", "Supplementary species-metabolite bridge", "p2_bridge_analyses/tables/p2_bifidobacterium_scfa_same_time_matrix.csv", "Small matched N; no FDR-supported bridge claims."],
  ["Clinical modifier matrix", "91 records", "32", "D05, D14, D30, D90", "Exploratory modifier audit", "p1_clinical_modifiers/tables/p1_scfa_maturation_clinical_analysis_matrix.csv", "Use as exploratory support because clinical variables have incomplete coverage."],
];

const tableListHeaders = ["TableLabel", "Placement", "Title", "SourceOrLinkedCSV", "Purpose", "SubmissionStatus"];
const tableListRows = [
  ["Table 1", "Main text", "Cohort and modular multi-omics coverage", "Table1_cohort_module_coverage.csv", "Defines the available analysis modules, dyads, timepoints and interpretation boundaries.", "Ready as formal manuscript table; numeric coverage should be checked against final locked cohort metadata before submission."],
  ["Table S1", "Supplementary", "Figure order and source-data crosswalk", "TableS1_figure_order_crosswalk.csv", "Maps final figure order, manuscript placement, asset status and action.", "Ready for internal manuscript assembly and source-data tracking."],
  ["Table S2", "Supplementary", "Evidence strength and claim boundary table", "final_evidence_strength_table.csv", "Keeps each result block aligned with claim strength and caution language.", "Ready as supplementary planning table; can be shortened for journal supplement."],
  ["Table S3", "Supplementary", "Analysis atlas and module source paths", "final_analysis_atlas.csv", "Documents row counts, dyads, timepoints and source paths for reproducibility.", "Ready as internal/supplementary methods table."],
  ["Table S4", "Supplementary", "Lagged milk-to-SCFA candidate model highlights", "figures_redrawn/fig5_milk_scfa_lagged_candidates/Fig5_source_formal_model_candidates.csv", "Provides the statistical backbone for Fig. 5.", "Use as source data; keep exploratory candidate wording."],
  ["Table S5", "Supplementary", "Maternal fecal context source summary", "figures_redrawn/fig6_maternal_fecal_context/Fig6_source_own_vs_other_summary.csv", "Provides source data for supplementary maternal fecal context figure.", "Supplementary/source-data ready."],
  ["Table S6", "Supplementary", "Clinical modifier model and contrast source data", "figures_redrawn/fig7_clinical_modifier_layer/Fig7_source_model_forest_terms.csv", "Provides source data for exploratory clinical modifier audit.", "Supplementary/source-data ready."],
];

const sourceFiles = await listFilesRecursively(path.join(pkg, "figures_redrawn"));
const sourceDataRows = sourceFiles
  .filter((file) => /\.(csv|png|pdf|svg|tiff|R|md)$/i.test(file))
  .map((file) => {
    const base = path.basename(file);
    const dir = path.basename(path.dirname(file));
    let linkedFigure = "";
    if (dir.startsWith("fig1_")) linkedFigure = "Fig. 1";
    else if (dir.startsWith("fig2_")) linkedFigure = "Fig. 2";
    else if (dir.startsWith("fig3_")) linkedFigure = "Fig. 3";
    else if (dir.startsWith("fig4_")) linkedFigure = "Fig. 4";
    else if (dir.startsWith("fig5_")) linkedFigure = "Fig. 5";
    else if (dir.startsWith("fig6_")) linkedFigure = "Fig. S1 (former Fig. 6)";
    else if (dir.startsWith("fig7_")) linkedFigure = "Fig. S2 (former Fig. 7)";
    const ext = path.extname(base).replace(".", "").toUpperCase();
    const role = /\.csv$/i.test(base) ? "Source data" : /\.(png|pdf|svg|tiff)$/i.test(base) ? "Figure export" : /\.R$/i.test(base) ? "Reproducible redraw script" : "Notes";
    return [linkedFigure, dir, role, ext, base, rel(file)];
  });

const sourceDataHeaders = ["LinkedFigure", "FigureFolder", "FileRole", "Format", "FileName", "RelativePath"];

await writeTable("Table1_cohort_module_coverage", mainTableHeaders, mainTableRows);
await writeTable("TableS1_figure_order_crosswalk", figureOrderHeaders, figureOrderRows);
await writeTable("TableS2_supplementary_figure_plan", suppFigureHeaders, suppFigureRows);
await writeTable("TableS3_table_package_index", tableListHeaders, tableListRows);
await writeTable("TableS4_source_data_index", sourceDataHeaders, sourceDataRows);

const oldEvidence = await readCsv(path.join(pkg, "tables", "final_evidence_strength_table.csv"));
const evidenceHeaders = ["ResultBlock", "EvidenceLayer", "MainFinding", "KeyStatistic", "ClaimStrength", "ManuscriptUse", "Caution"];
const evidenceRows = oldEvidence.map((r) => evidenceHeaders.map((h) => r[h] ?? ""));
await writeTable("TableS5_evidence_strength_claim_boundaries", evidenceHeaders, evidenceRows);

const analysisAtlas = await readCsv(path.join(pkg, "tables", "final_analysis_atlas.csv"));
const atlasHeaders = ["AnalysisSet", "Rows", "Dyads", "Timepoints", "SourcePath"];
const atlasRows = analysisAtlas.map((r) => atlasHeaders.map((h) => r[h] ?? ""));
await writeTable("TableS6_analysis_atlas", atlasHeaders, atlasRows);

const workbook = Workbook.create();
addSheet(workbook, "README", ["Item", "Value"], [
  ["Package", "Manuscript figure and table ordering package"],
  ["Generated", "2026-06-07"],
  ["Main figure order", "Graphical abstract; Fig. 1 design overview; Fig. 2 milk bioactive dynamics; Fig. 3 infant microbiome maturation; Fig. 4 SCFA maturation; Fig. 5 lagged milk-to-SCFA candidate links"],
  ["Supplementary figure order", "Fig. S1 maternal fecal context; Fig. S2 clinical modifier evidence audit"],
  ["Claim boundary", "Longitudinal associations and candidate temporal links only; avoid causal, mediation, prediction, regulation, and strain-transmission wording."],
], [140, 520]);
addSheet(workbook, "MainFigures", figureOrderHeaders, figureOrderRows, [72, 120, 90, 220, 320, 180, 210, 220]);
addSheet(workbook, "Table1", mainTableHeaders, mainTableRows, [170, 100, 60, 120, 210, 260, 300]);
addSheet(workbook, "SuppFigures", suppFigureHeaders, suppFigureRows, [110, 90, 230, 320, 300]);
addSheet(workbook, "TableIndex", tableListHeaders, tableListRows, [85, 105, 230, 260, 320, 280]);
addSheet(workbook, "SourceDataIndex", sourceDataHeaders, sourceDataRows, [105, 190, 135, 70, 230, 420]);
addSheet(workbook, "EvidenceStrength", evidenceHeaders, evidenceRows, [200, 80, 340, 320, 150, 220, 280]);
addSheet(workbook, "AnalysisAtlas", atlasHeaders, atlasRows, [220, 70, 70, 130, 340]);

const inspect = await workbook.inspect({
  kind: "workbook,sheet,table",
  maxChars: 6000,
  tableMaxRows: 4,
  tableMaxCols: 8,
});
console.log(inspect.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "formula error scan",
});
console.log(errors.ndjson);

for (const sheetName of ["README", "MainFigures", "Table1", "SourceDataIndex"]) {
  const preview = await workbook.render({ sheetName, autoCrop: "all", scale: 1, format: "png" });
  await fs.writeFile(path.join(outDir, `${sheetName}_preview.png`), new Uint8Array(await preview.arrayBuffer()));
}

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(xlsxPath);
console.log(`WROTE ${xlsxPath}`);
