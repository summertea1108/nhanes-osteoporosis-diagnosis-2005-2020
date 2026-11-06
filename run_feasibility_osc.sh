#!/bin/bash
# =============================================================================
# Run the feasibility and reliability checks on OSC
#
# Usage (from an OSC login node):
#   bash run_feasibility_osc.sh
#
# What it does:
#   1. Checks that both R scripts and the existing NHANES data are in place
#   2. Downloads the 25 extra NHANES files (login nodes have internet;
#      compute nodes usually do not)
#   3. Submits a short batch job that runs the R script
# Results: results_feasibility/run_log.txt (summary at the end) plus CSV tables
# =============================================================================
set -euo pipefail

# ---- settings (edit if your layout differs) ---------------------------------
PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # folder containing this script
DATA_DIR="${PROJ_DIR}/nhanes_data"          # folder holding the .xpt files used by the main script
OUT_DIR="${PROJ_DIR}/results_feasibility"
MAIN_R="${PROJ_DIR}/nhanes_cascade_analysis.R"
FEAS_R="${PROJ_DIR}/nhanes_feasibility_crosscondition.R"
ACCOUNT="YOUR_OSC_PROJECT"   # replace with your OSC project code
R_MODULE=""                                 # leave empty to use OSC's default R; or e.g. "R/4.4.0-gnu11.2"
# ------------------------------------------------------------------------------

cd "$PROJ_DIR"

echo "== 1. Checking files =="
for f in "$MAIN_R" "$FEAS_R"; do
  [[ -f "$f" ]] || { echo "Missing: $f  (upload it to $PROJ_DIR)"; exit 1; }
done
[[ -d "$DATA_DIR" ]] || { echo "Data folder not found: $DATA_DIR  (edit DATA_DIR at the top)"; exit 1; }
[[ -f "$DATA_DIR/DEMO_D.xpt" || -f "$DATA_DIR/P_DEMO.xpt" ]] || \
  echo "Warning: no DEMO file found in $DATA_DIR. Check DATA_DIR points to the main analysis data."
echo "OK"

echo "== 2. Downloading extra NHANES files (skips files already present) =="
declare -A YEAR=( [_D]=2005 [_E]=2007 [_F]=2009 [_H]=2013 )
fetch() {   # $1 = file name without .xpt, $2 = cycle start year
  local f="$1" y="$2" dest="$DATA_DIR/$1.xpt"
  if [[ -s "$dest" ]]; then echo "  have $f"; return; fi
  echo "  downloading $f"
  curl -fsSL --retry 3 -o "$dest.part" \
    "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/${y}/DataFiles/${f}.xpt" \
    && mv "$dest.part" "$dest" \
    || { rm -f "$dest.part"; echo "  FAILED: $f"; FAILED=1; }
}
FAILED=0
for suf in _D _E _F _H; do
  for base in GHB DIQ BPQ BPX TCHOL; do fetch "${base}${suf}" "${YEAR[$suf]}"; done
done
for base in P_GHB P_DIQ P_BPQ P_BPXO P_TCHOL; do fetch "$base" 2017; done

# a real XPT file starts with a SAS header; catch HTML error pages saved by mistake
for f in "$DATA_DIR"/{GHB,DIQ,BPQ,BPX,TCHOL}_{D,E,F,H}.xpt "$DATA_DIR"/P_{GHB,DIQ,BPQ,BPXO,TCHOL}.xpt; do
  [[ -f "$f" ]] || continue
  if ! head -c 80 "$f" | grep -q "HEADER RECORD"; then
    echo "  Not a valid XPT file (probably an error page): $f"; rm -f "$f"; FAILED=1
  fi
done
if [[ "$FAILED" -ne 0 ]]; then
  echo "Some downloads failed. Re-run this script, or download those files manually from"
  echo "https://wwwn.cdc.gov/nchs/nhanes/ into $DATA_DIR"
  exit 1
fi
echo "All 25 files present."

echo "== 3. Submitting the batch job =="
mkdir -p "$OUT_DIR"
JOB=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=op_feasibility
#SBATCH --account=${ACCOUNT}
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=${OUT_DIR}/slurm_%j.out

set -e
module purge
if [[ -n "${R_MODULE}" ]]; then module load ${R_MODULE}; else module load R; fi
cd "${PROJ_DIR}"
echo "R: \$(Rscript --version 2>&1)"
Rscript "${FEAS_R}" "${DATA_DIR}" "${OUT_DIR}" "${MAIN_R}"
EOF
)
echo "Submitted job $JOB"
echo
echo "Check progress:   squeue -u \$USER"
echo "Watch the output: tail -f ${OUT_DIR}/slurm_${JOB}.out"
echo "When finished, the summary is at the end of: ${OUT_DIR}/run_log.txt"
echo "   tail -n 45 ${OUT_DIR}/run_log.txt"
