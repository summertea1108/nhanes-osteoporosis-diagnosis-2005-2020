#!/bin/bash
#SBATCH --job-name=nhanes_cascade
#SBATCH --account=YOUR_OSC_PROJECT
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --output=nhanes_cascade_%j.log

cd "$SLURM_SUBMIT_DIR"
module load gcc/12.3.0 R/4.6.1
Rscript nhanes_cascade_analysis.R nhanes_data results
