// Table 3. Prevalence ratios for reported prior diagnosis: osteoporosis vs other measured
// chronic conditions, adults 50+ with health care contact (primary sample, measured definitions).
// Reads the CSV outputs of nhanes_crosscondition_analysis.R and writes Table3.docx (JAMA style:
// title above, horizontal rules only, footnotes with superscript letters, editable Word table).
// Usage: node make_table3.js [results_dir] [out.docx]
const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
  BorderStyle, AlignmentType, VerticalAlign, PageOrientation,
} = require("docx");

const DIR = process.argv[2] || ".";
const OUT = process.argv[3] || "Table3.docx";

// ---- CSV reader (handles quoted fields) ------------------------------------
function readCsv(file) {
  const text = fs.readFileSync(path.join(DIR, file), "utf8").replace(/\r/g, "");
  const lines = text.split("\n").filter((l) => l.length);
  const parse = (line) => {
    const out = []; let cur = ""; let q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (q) { if (c === '"' && line[i + 1] === '"') { cur += '"'; i++; } else if (c === '"') q = false; else cur += c; }
      else if (c === '"') q = true; else if (c === ",") { out.push(cur); cur = ""; } else cur += c;
    }
    out.push(cur); return out;
  };
  const head = parse(lines[0]);
  return lines.slice(1).map((l) => Object.fromEntries(parse(l).map((v, i) => [head[i], v])));
}
const SAMPLE = "Health care contact (primary)", DEF = "Primary (measured)";
const a1 = readCsv("crosscond_A1_prevalence_ratios.csv").filter((r) => r.sample === SAMPLE && r.definition === DEF);
const a1x = readCsv("crosscond_A1_sex_specificity.csv").filter((r) => r.sample === SAMPLE && r.definition === DEF);
const a2 = readCsv("crosscond_A2_prevalence_ratios.csv")
  .filter((r) => r.sample === SAMPLE && r.definition === DEF && r.pairs === "All three conditions");

const one = (arr, pred, what) => {
  const m = arr.filter(pred);
  if (m.length !== 1) throw new Error(`Expected 1 row for ${what}, found ${m.length}`);
  return m[0];
};
const f2 = (x) => Number(x).toFixed(2);
const ci = (r) => `${f2(r.PR)} (${f2(r.lcl)}–${f2(r.ucl)})`;
const n0 = (x) => Number(x).toLocaleString("en-US");

// ---- content ---------------------------------------------------------------
const COMPS = [["Diabetes", "diabetes"], ["Hypertension", "hypertension"], ["High cholesterol", "high cholesterol"]];
const sexes = ["Women", "Men", "Both"];

const a1n = sexes.map((s) => one(a1, (r) => r.sex === s && r.contrast === "Osteoporosis vs Diabetes", `A1 n ${s}`));
const a1rows = COMPS.map(([key, label]) => ({
  label: `Osteoporosis vs ${label}`,
  cells: [
    ...sexes.map((s) => ci(one(a1, (r) => r.sex === s && r.contrast === `Osteoporosis vs ${key}`, `A1 ${s} ${key}`))),
    ci(one(a1x, (r) => r.contrast === `Sex PR (men/women): osteoporosis relative to ${key}`, `A1x ${key}`)),
  ],
}));
const a2by = Object.fromEntries(["Women", "Men", "Both", "Men vs women"].map((s) => [s, one(a2, (r) => r.sex === s, `A2 ${s}`)]));
const pairs = (r) => Number(r.n) / 2;   // each pair contributes 2 rows

// ---- styling helpers -------------------------------------------------------
const FONT = "Arial", SZ = 18;                 // 9 pt
const W = [3060, 1575, 1575, 1575, 1575];      // sums to 9360 DXA = 6.5 in
const TOTAL = W.reduce((a, b) => a + b, 0);
const RULE = { style: BorderStyle.SINGLE, size: 6, color: "000000" };
const NONE = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };

function run(text, opts = {}) { return new TextRun({ text, font: FONT, size: opts.size || SZ, bold: opts.bold, italics: opts.italics, superScript: opts.sup }); }
function para(runs, align = AlignmentType.LEFT, indent = 0) {
  return new Paragraph({ children: runs, alignment: align, indent: indent ? { left: indent } : undefined, spacing: { before: 30, after: 30 } });
}
function cell(runs, i, { align, top, bottom, span, indent } = {}) {
  return new TableCell({
    children: [para(runs, align ?? (i === 0 ? AlignmentType.LEFT : AlignmentType.CENTER), indent)],
    width: { size: span ? TOTAL : W[i], type: WidthType.DXA },
    columnSpan: span ? 5 : undefined,
    verticalAlign: VerticalAlign.BOTTOM,
    margins: { top: 20, bottom: 20, left: 80, right: 80 },
    borders: { top: top ? RULE : NONE, bottom: bottom ? RULE : NONE, left: NONE, right: NONE },
  });
}
const dataRow = (label, values, opts = {}) => new TableRow({
  children: [cell([run(label)], 0, { indent: 180, bottom: opts.last }),
             ...values.map((v, i) => cell([run(v)], i + 1, { bottom: opts.last }))],
});
const sectionRow = (title, sup) => new TableRow({
  children: [cell([run(title, { bold: true }), run(sup, { sup: true, bold: true })], 0, { span: true })],
});

const header = new TableRow({
  tableHeader: true,
  children: [
    cell([run("Comparison", { bold: true })], 0, { top: true }),
    ...["Women", "Men", "Both sexes"].map((h, i) =>
      cell([run(h, { bold: true })], i + 1, { top: true })),
    cell([run("Ratio of PRs, men vs women", { bold: true }), run("c", { sup: true, bold: true })], 4, { top: true }),
  ],
});
const subhead = new TableRow({
  tableHeader: true,
  children: [cell([run("")], 0, { bottom: true }),
             ...[0, 1, 2].map((i) => cell([run("PR (95% CI)")], i + 1, { bottom: true })),
             cell([run("Ratio (95% CI)")], 4, { bottom: true })],
});

const rows = [
  header,
  subhead,
  sectionRow("Population-level comparison", "a"),
  dataRow("No. of observations (reporting a diagnosis)",
          [...a1n.map((r) => `${n0(r.n)} (${n0(r.events)})`), "NA"]),
  ...a1rows.map((r) => dataRow(r.label, r.cells)),
  sectionRow("Within-person comparison", "b"),
  dataRow("No. of participants (pairs)",
          [...["Women", "Men", "Both"].map((s) => `${n0(a2by[s].persons)} (${n0(pairs(a2by[s]))})`), "NA"]),
  dataRow("Osteoporosis vs comparison condition, same participant",
          [ci(a2by.Women), ci(a2by.Men), ci(a2by.Both), ci(a2by["Men vs women"])], { last: true }),
];

const table = new Table({ rows, width: { size: TOTAL, type: WidthType.DXA }, columnWidths: W,
  borders: { top: NONE, bottom: NONE, left: NONE, right: NONE, insideHorizontal: NONE, insideVertical: NONE } });

const note = (runs) => new Paragraph({ children: runs, spacing: { before: 40, after: 40 } });
const NS = 16; // 8 pt footnotes
const footnotes = [
  note([run("Abbreviations: NA, not applicable; PR, prevalence ratio.", { size: NS })]),
  note([run("Data are from the National Health and Nutrition Examination Survey, 2005-2010, 2013-2014, and 2017-March 2020, among adults 50 years or older with health care contact (health insurance, a usual source of care, and at least 1 health care visit in the past 12 months). Each condition was defined by measurement at the survey examination (definitions in Figure 1). PRs compare the survey-weighted prevalence of a reported prior diagnosis of osteoporosis with that of the comparison condition; values below 1 indicate that osteoporosis was reported less often. Variance estimation accounts for the complex survey design, including multiple observations from the same participant.", { size: NS })]),
  note([run("a", { size: NS, sup: true }), run(" Adults with each condition; participants contributed 1 observation per condition they had. Modified Poisson regression adjusted for age group and race and ethnicity (and sex in the both-sexes model).", { size: NS })]),
  note([run("b", { size: NS, sup: true }), run(" Adults with osteoporosis who also had diabetes, hypertension, or high cholesterol. Each pair contributed the participant’s osteoporosis report and comparison condition report; pairs were weighted by the survey weight divided by the participant’s number of pairs, so each participant counted once. Unadjusted, because both reports come from the same participant.", { size: NS })]),
  note([run("c", { size: NS, sup: true }), run(" The PR for men divided by the PR for women, estimated from a condition (or report type) × sex interaction term. Equivalently, the sex PR for osteoporosis divided by the sex PR for the comparison condition; values below 1 indicate a larger sex difference for osteoporosis.", { size: NS })]),
];

const title = new Paragraph({
  spacing: { after: 120 },
  children: [run("Table 3. Prevalence Ratios for Reported Prior Diagnosis of Osteoporosis Compared With Diabetes, Hypertension, and High Cholesterol Among US Adults 50 Years or Older With Health Care Contact", { bold: true, size: 20 })],
});

const doc = new Document({
  sections: [{
    properties: { page: { size: { width: 12240, height: 15840, orientation: PageOrientation.PORTRAIT },
                          margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    children: [title, table, ...footnotes],
  }],
});
Packer.toBuffer(doc).then((buf) => { fs.writeFileSync(OUT, buf); console.log("Wrote", OUT); });
