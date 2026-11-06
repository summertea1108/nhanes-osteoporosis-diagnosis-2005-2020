#!/bin/bash
# Pre-submission checks for the osteoporosis sex-difference paper.  Usage:  bash run_submission_checks_osc.sh
set -euo pipefail
PROJ_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # folder containing this script
OUT_DIR="${PROJ_DIR}/results_submission_checks"
cd "$PROJ_DIR"
for f in nhanes_cascade_analysis.R nhanes_crosscondition_analysis.R nhanes_submission_checks.R; do
  [[ -f "$f" ]] || { echo "Missing: $f  (upload it to $PROJ_DIR)"; exit 1; }
done
mkdir -p "$OUT_DIR"
JOB=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=sub_checks
#SBATCH --account=YOUR_OSC_PROJECT
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --output=${OUT_DIR}/slurm_%j.out
set -e
module reset
module load gcc/12.3.0 R/4.4.0
cd "${PROJ_DIR}"
Rscript nhanes_submission_checks.R "${PROJ_DIR}/nhanes_data" "${OUT_DIR}" nhanes_cascade_analysis.R nhanes_crosscondition_analysis.R
EOF
)
echo "Submitted job $JOB"
echo "Check progress:  squeue -u \$USER"
echo "When done:   sed -n '/SUMMARY/,\$p' ${OUT_DIR}/run_log.txt"
echo "If it fails: tail -n 25 ${OUT_DIR}/slurm_${JOB}.out"
