#!/usr/bin/env python3
"""
Figure 1. Reported prior diagnosis of osteoporosis compared with diabetes,
hypertension and high cholesterol among US adults 50+ with health care contact,
by sex. NHANES 2005-2010, 2013-2014 and 2017-March 2020.

Input : crosscond_A1_percentages.csv (from nhanes_crosscondition_analysis.R)
Output: Figure1.pdf (vector, fonts embedded), Figure1.tiff (600 dpi, LZW), Figure1.png (preview)

Usage : python3 make_figure1.py [path/to/crosscond_A1_percentages.csv] [out_dir]

Design notes (JAMA Network Open artwork):
  - Sans-serif type (Arial; Liberation Sans is the metric-identical fallback), 8-9 pt at final size
  - No title inside the artwork (title and footnotes go in the figure legend)
  - Full-page width 6.8 in; colors checked for color-vision deficiency, and sex is also
    encoded by marker shape so the figure reads in grayscale
  - Unweighted No.; weighted % with Korn-Graubard 95% CI
"""
import sys, os, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
from matplotlib.lines import Line2D

CSV = sys.argv[1] if len(sys.argv) > 1 else "crosscond_A1_percentages.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."
os.makedirs(OUT, exist_ok=True)

SAMPLE = "Health care contact (primary)"
DEF = "Primary (measured)"
CONDITIONS = ["Osteoporosis", "Diabetes", "Hypertension", "High cholesterol"]
SEXES = ["Women", "Men"]

# ---- read the analysis output -------------------------------------------------
rows = {}
with open(CSV, newline="") as f:
    for r in csv.DictReader(f):
        if r["sample"] == SAMPLE and r["definition"] == DEF and r["sex"] in SEXES:
            rows[(r["condition"], r["sex"])] = r
missing = [(c, s) for c in CONDITIONS for s in SEXES if (c, s) not in rows]
if missing:
    sys.exit(f"Missing rows in {CSV}: {missing}")
for k, r in rows.items():
    if r["nchs_flag"] != "Reliable" or r["nchs_flag_complement"] != "Reliable":
        sys.exit(f"Estimate not reliable under NCHS standards, check before plotting: {k} {r['nchs_flag']}")

# ---- style --------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8.5,
    "pdf.fonttype": 42, "ps.fonttype": 42,     # embed TrueType so text stays editable
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "xtick.major.size": 3,
})
INK, INK2, GRID, BAND = "#1a1a1a", "#555555", "#d9d9d9", "#f1f1f1"
COLOR = {"Women": "#2a78d6", "Men": "#eb6834"}     # validated pair (CVD dE 24.7, contrast >= 3:1)
MARKER = {"Women": "o", "Men": "s"}
MSIZE = {"Women": 5.6, "Men": 5.0}                 # square drawn slightly smaller to match visual area

# ---- layout: a heading line per condition, then one row per sex ---------------
ROW, GAP = 1.0, 0.45
ys, head_y, group_span, y = {}, {}, {}, 0.0
for c in CONDITIONS:
    top = y + 0.5
    head_y[c] = y
    y -= ROW
    for s in SEXES:
        ys[(c, s)] = y
        y -= ROW
    group_span[c] = (top, y + 0.5)
    y -= GAP
y_min, y_max = y + GAP + 0.35, 1.35

fig = plt.figure(figsize=(6.8, 4.1))
ax = fig.add_axes([0.50, 0.135, 0.475, 0.79])
ax.set_xlim(0, 100)
ax.set_ylim(y_min, y_max)
tr = blended_transform_factory(ax.transAxes, ax.transData)   # x in axes fraction, y in data

# column positions (axes-fraction; negative = left of the plot)
X_LAB, X_IND, X_N, X_EST, X_END = -1.03, -0.97, -0.60, -0.05, 1.0

# shaded band behind the osteoporosis group
lo, hi = group_span["Osteoporosis"][1], group_span["Osteoporosis"][0]
ax.add_patch(plt.Rectangle((X_LAB - 0.015, lo), X_END - X_LAB + 0.015, hi - lo, transform=tr,
                           color=BAND, zorder=0, clip_on=False, lw=0))

# vertical grid
for gx in (25, 50, 75, 100):
    ax.axvline(gx, color=GRID, lw=0.5, zorder=1)

for c in CONDITIONS:
    ax.text(X_LAB, head_y[c], c, transform=tr, ha="left", va="center", fontweight="bold", color=INK)
    for s in SEXES:
        r = rows[(c, s)]
        yy = ys[(c, s)]
        p, l, u, n = float(r["pct"]), float(r["lcl"]), float(r["ucl"]), int(r["n"])
        ax.plot([l, u], [yy, yy], color=COLOR[s], lw=1.4, solid_capstyle="butt", zorder=3)
        for cap in (l, u):
            ax.plot([cap, cap], [yy - 0.16, yy + 0.16], color=COLOR[s], lw=1.0, zorder=3)
        ax.plot(p, yy, marker=MARKER[s], ms=MSIZE[s], color=COLOR[s],
                mec="white", mew=0.7, zorder=4, linestyle="none")
        ax.text(X_IND, yy, s, transform=tr, ha="left", va="center", color=INK)
        ax.text(X_N, yy, f"{n:,}", transform=tr, ha="right", va="center", color=INK)
        ax.text(X_EST, yy, f"{p:.1f} ({l:.1f}\u2013{u:.1f})", transform=tr, ha="right",
                va="center", color=INK)

# column headers
HEAD_Y = y_max - 0.15
for x, txt, ha in [(X_LAB, "Condition and sex", "left"), (X_N, "No.", "right"),
                   (X_EST, "Reported prior diagnosis,\n% (95% CI)", "right")]:
    ax.text(x, HEAD_Y, txt, transform=tr, ha=ha, va="bottom", fontweight="bold", color=INK,
            linespacing=1.1)
ax.plot([X_LAB, X_END], [HEAD_Y - 0.15] * 2, transform=tr, color=INK, lw=0.6, clip_on=False)

# axes
ax.set_yticks([])
for side in ("left", "right", "top"):
    ax.spines[side].set_visible(False)
ax.spines["bottom"].set_color(INK)
ax.set_xticks([0, 25, 50, 75, 100])
ax.tick_params(axis="x", colors=INK, labelsize=8.5)
ax.set_xlabel("Adults reporting a prior diagnosis, %", color=INK, labelpad=5)

# legend (identity is also given by the Sex column, so color is never the only cue)
handles = [Line2D([0], [0], marker=MARKER[s], color=COLOR[s], lw=1.4, ms=MSIZE[s],
                  mec="white", mew=0.7, label=s) for s in SEXES]
ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.995), ncol=2, frameon=False,
          handlelength=2.0, columnspacing=1.4, borderaxespad=0.3, fontsize=8.5)

for ext, kw in [("pdf", {}), ("png", {"dpi": 300}),
                ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})]:
    fig.savefig(os.path.join(OUT, f"Figure1.{ext}"), facecolor="white", **kw)
print("Wrote Figure1.pdf, Figure1.png, Figure1.tiff to", os.path.abspath(OUT))
