#!/bin/bash
# =============================================================================
# Run the cross-condition analysis on OSC (Ascend)
# Usage (from a login node):   bash run_crosscondition_osc.sh
# Results: results_crosscondition/run_log.txt (summary at the end) plus CSV tables
# The extra NHANES files were already downloaded by run_feasibility_osc.sh.
# =============================================================================
set -euo pipefail

PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # folder containing this script
DATA_DIR="${PROJ_DIR}/nhanes_data"
OUT_DIR="${PROJ_DIR}/results_crosscondition"
MAIN_R="${PROJ_DIR}/nhanes_cascade_analysis.R"
ANALYSIS_R="${PROJ_DIR}/nhanes_crosscondition_analysis.R"
ACCOUNT="YOUR_OSC_PROJECT"   # replace with your OSC project code
R_MODULE="gcc/12.3.0 R/4.4.0"

cd "$PROJ_DIR"
for f in "$MAIN_R" "$ANALYSIS_R"; do
  [[ -f "$f" ]] || { echo "Missing: $f  (upload it to $PROJ_DIR)"; exit 1; }
done
for f in P_GHB P_DIQ P_BPQ P_BPXO P_TCHOL GHB_D BPX_H; do
  [[ -s "$DATA_DIR/$f.xpt" ]] || { echo "Missing $f.xpt in $DATA_DIR: run bash run_feasibility_osc.sh first"; exit 1; }
done
mkdir -p "$OUT_DIR"

JOB=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=op_crosscond
#SBATCH --account=${ACCOUNT}
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=${OUT_DIR}/slurm_%j.out

set -e
module reset
module load ${R_MODULE}
cd "${PROJ_DIR}"
Rscript "${ANALYSIS_R}" "${DATA_DIR}" "${OUT_DIR}" "${MAIN_R}"
EOF
)
echo "Submitted job $JOB"
echo "Check progress:  squeue -u \$USER"
echo "When it is gone from the queue, read the results:"
echo "   tail -n 40 ${OUT_DIR}/run_log.txt"
echo "If run_log.txt is missing or ends early:"
echo "   tail -n 25 ${OUT_DIR}/slurm_${JOB}.out"
