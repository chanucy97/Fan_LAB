import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const WORKSPACE_ROOT = path.resolve(process.env.FANLAB_RESEARCH2_ROOT || process.cwd());
const ROOT = path.join(WORKSPACE_ROOT, "final_manuscript_planning_20260605");
const OUT_DIR = path.join(ROOT, "gut_microbes_revision", "supplementary_tables");
const OUT_XLSX = path.join(OUT_DIR, "gut_microbes_supplementary_tables_20260613.xlsx");

const CODE_REPOSITORY = "https://github.com/chanucy97/Fan_LAB";

function csvParse(text) {
  const rows = [];
  let row = [];
  let cell = "";
  let inQuotes = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    const next = text[i + 1];
    if (inQuotes) {
      if (ch === '"' && next === '"') {
        cell += '"';
        i += 1;
      } else if (ch === '"') {
        inQuotes = false;
      } else {
        cell += ch;
      }
      continue;
    }
    if (ch === '"') {
      inQuotes = true;
    } else if (ch === ",") {
      row.push(cell);
      cell = "";
    } else if (ch === "\n") {
      row.push(cell);
      rows.push(row);
      row = [];
      cell = "";
    } else if (ch !== "\r") {
      cell += ch;
    }
  }
  if (cell.length || row.length) {
    row.push(cell);
    rows.push(row);
  }
  return rows;
}

function coerceCell(value) {
  if (value === "") return null;
  if (value === "NA" || value === "NaN" || value === "Inf" || value === "-Inf") return value;
  if (/^-?\d+(\.\d+)?([eE][+-]?\d+)?$/.test(value)) return Number(value);
  return value;
}

async function readCsv(relPath) {
  const full = path.join(ROOT, relPath);
  const text = await fs.readFile(full, "utf8");
  const rows = csvParse(text).map((r) => r.map(coerceCell));
  return { rows, full, relPath };
}

function colLetter(n) {
  let s = "";
  while (n > 0) {
    const m = (n - 1) % 26;
    s = String.fromCharCode(65 + m) + s;
    n = Math.floor((n - 1) / 26);
  }
  return s;
}

function safeTableName(sheetName) {
  return `${sheetName.replace(/[^A-Za-z0-9_]/g, "_").slice(0, 24)}_Tbl`;
}

function writeMatrixSheet(workbook, sheetName, rows, options = {}) {
  const sheet = workbook.worksheets.add(sheetName);
  sheet.showGridLines = false;
  const maxCols = Math.max(...rows.map((r) => r.length));
  const normalized = rows.map((r) => {
    const copy = [...r];
    while (copy.length < maxCols) copy.push(null);
    return copy;
  });
  const range = sheet.getRangeByIndexes(0, 0, normalized.length, maxCols);
  range.values = normalized;
  const header = sheet.getRangeByIndexes(0, 0, 1, maxCols);
  header.format = {
    fill: "#1F4E79",
    font: { bold: true, color: "#FFFFFF" },
    wrapText: true,
  };
  range.format.borders = { preset: "all", style: "thin", color: "#D9E2F3" };
  range.format.wrapText = true;
  try {
    sheet.freezePanes.freezeRows(1);
  } catch {}
  try {
    range.format.autofitColumns();
    range.format.autofitRows();
  } catch {}
  const end = `${colLetter(maxCols)}${normalized.length}`;
  if (normalized.length > 1 && maxCols > 1 && options.table !== false) {
    try {
      const table = sheet.tables.add(`A1:${end}`, true, safeTableName(sheetName));
      table.showFilterButton = true;
      table.showBandedRows = true;
    } catch {}
  }
  return sheet;
}

function makeRows(headers, objects) {
  return [headers, ...objects.map((obj) => headers.map((h) => obj[h] ?? null))];
}

const tableSpecs = [
  {
    sheet: "S1_Module_Coverage",
    relPath: "manuscript_tables/Table1_cohort_module_coverage.csv",
    use: "Cohort module overview for Figure 1 and Methods coverage description",
  },
  {
    sheet: "S2_Coverage_By_Window",
    relPath: "figures_redrawn/fig1_cohort_design_overview/Fig1_source_coverage_module_counts.csv",
    use: "Dyad counts by assay module and time window for Figure 1B",
  },
  {
    sheet: "S3_Dyad_Coverage",
    relPath: "figures_redrawn/fig1_cohort_design_overview/Fig1_source_coverage_dyad_counts.csv",
    use: "Per-dyad module completeness supporting the coverage heatmap",
  },
  {
    sheet: "S4_Milk_Atlas",
    relPath: "figures_redrawn/fig2_milk_bioactive_remodeling/Fig2_source_atlas.csv",
    use: "Milk LTF/HMO/LCFA-FFA longitudinal summaries for Figure 2A",
  },
  {
    sheet: "S5_Milk_Temporal_Rank",
    relPath: "figures_redrawn/fig2_milk_bioactive_remodeling/Fig2_source_ranked_dynamic_features.csv",
    use: "Ranked milk temporal effects and candidate exposure flags for Figure 2B/E",
  },
  {
    sheet: "S6_Microbiome_Counts",
    relPath: "figures_redrawn/fig3_infant_microbiome_maturation/Fig3_source_sample_counts.csv",
    use: "Infant fecal metagenome sample counts by timepoint for Figure 3",
  },
  {
    sheet: "S7_Microbiome_PERMANOVA",
    relPath: "figures_redrawn/fig3_infant_microbiome_maturation/Fig3_source_permanova.csv",
    use: "Community-level age-associated separation statistic for Figure 3B",
  },
  {
    sheet: "S8_Bifido_GTDB_Labels",
    relPath: "figures_redrawn/fig3_infant_microbiome_maturation/Fig3_source_bifidobacterium_species_summary.csv",
    use: "Taxonomy-reviewed sylph/GTDB Bifidobacterium labels for Figure 3E/F",
  },
  {
    sheet: "S9_SCFA_Time_Trends",
    relPath: "figures_redrawn/fig4_scfa_maturation/Fig4_source_scfa_time_trends.csv",
    use: "Feature-level SCFA age trends and score directions for Figure 4C",
  },
  {
    sheet: "S10_SCFA_Score_Model",
    relPath: "figures_redrawn/fig4_scfa_maturation/Fig4_source_score_model.csv",
    use: "Composite SCFA maturation score age model for Figure 4B",
  },
  {
    sheet: "S11_SCFA_Score_Values",
    relPath: "figures_redrawn/fig4_scfa_maturation/Fig4_source_score_values.csv",
    use: "Sample-level SCFA maturation score values supporting Figure 4",
  },
  {
    sheet: "S12_Lagged_Screen",
    relPath: "figures_redrawn/fig5_milk_scfa_lagged_candidates/Fig5_source_lagged_screen_hmo_ltf_scfa.csv",
    use: "Full lagged HMO/LTF-to-SCFA screening table for Figure 5A",
  },
  {
    sheet: "S13_Formal_Models",
    relPath: "figures_redrawn/fig5_milk_scfa_lagged_candidates/Fig5_source_formal_model_candidates.csv",
    use: "Formal lagged candidate model table for Figure 5B/C",
  },
  {
    sheet: "S14_Evidence_Counts",
    relPath: "figures_redrawn/fig5_milk_scfa_lagged_candidates/Fig5_source_evidence_counts.csv",
    use: "Screen-to-model-to-adjusted signal attrition counts for Figure 5C",
  },
  {
    sheet: "S15_LODO_Summary",
    relPath: "gut_microbes_revision/gut_microbes_lagged_leave_one_dyad_out_summary.csv",
    use: "Leave-one-dyad-out stability summary for prioritized lagged candidates",
  },
  {
    sheet: "S16_LODO_Details",
    relPath: "gut_microbes_revision/gut_microbes_lagged_leave_one_dyad_out_details.csv",
    use: "Leave-one-dyad-out refit-level details",
  },
  {
    sheet: "S17_Maternal_Context",
    relPath: "figures_redrawn/fig6_maternal_fecal_context/Fig6_source_own_vs_other_summary.csv",
    use: "Maternal-infant fecal similarity context for Supplementary Figure S1",
  },
  {
    sheet: "S18_Clinical_Coverage",
    relPath: "figures_redrawn/fig7_clinical_modifier_layer/Fig7_source_coverage_by_time.csv",
    use: "Clinical and SCFA covariate coverage by timepoint for Supplementary Figure S2D",
  },
  {
    sheet: "S19_Clinical_Models",
    relPath: "figures_redrawn/fig7_clinical_modifier_layer/Fig7_source_model_forest_terms.csv",
    use: "Exploratory clinical modifier model terms for Supplementary Figure S2C",
  },
  {
    sheet: "S20_Clinical_Matrix",
    relPath: "figures_redrawn/fig7_clinical_modifier_layer/Fig7_source_clinical_analysis_matrix.csv",
    use: "Sample-level clinical covariate matrix; clinical covariates are not fully available",
  },
];

async function main() {
  await fs.mkdir(OUT_DIR, { recursive: true });
  const workbook = Workbook.create();

  const indexRows = [
    {
      Item: "Manuscript",
      Detail: "Human milk bioactives and infant gut metabolic maturation in a Qinghai-Tibetan Plateau mother-infant cohort",
      Notes: "Supplementary tables assembled from figure/model source data generated for the expert-review draft.",
    },
    {
      Item: "Code availability",
      Detail: CODE_REPOSITORY,
      Notes: "Analysis and figure-generation code will be maintained in the Fan_LAB repository.",
    },
    {
      Item: "Taxonomy review",
      Detail: "Completed",
      Notes: "sylph/GTDB Bifidobacterium labels have been reviewed and are treated as taxonomy-supported labels, not strain-level transmission evidence.",
    },
    {
      Item: "Clinical covariates",
      Detail: "Not fully available",
      Notes: "Clinical modifier analyses are exploratory because feeding, delivery, BMI, birth-weight and related covariates are not complete across all SCFA/metagenome samples.",
    },
    {
      Item: "Data access boundary",
      Detail: "Controlled access for participant-level data",
      Notes: "Derived tables are assembled here; participant-level clinical/sequencing data should be shared only under ethics/institutional approval.",
    },
  ];
  writeMatrixSheet(workbook, "README_Index", makeRows(["Item", "Detail", "Notes"], indexRows), { table: false });

  const sourceIndex = [];
  for (const spec of tableSpecs) {
    const { rows } = await readCsv(spec.relPath);
    writeMatrixSheet(workbook, spec.sheet, rows);
    sourceIndex.push({
      Sheet: spec.sheet,
      SourcePath: spec.relPath,
      Rows: Math.max(rows.length - 1, 0),
      Columns: rows[0]?.length ?? 0,
      ManuscriptUse: spec.use,
    });
  }
  writeMatrixSheet(workbook, "Source_Index", makeRows(["Sheet", "SourcePath", "Rows", "Columns", "ManuscriptUse"], sourceIndex));

  const output = await SpreadsheetFile.exportXlsx(workbook);
  await output.save(OUT_XLSX);
  console.log(OUT_XLSX);
}

await main();
