# Sex differences in reported diagnosis of osteoporosis vs other chronic conditions (NHANES 2005–March 2020)

Analysis code for:

> Liu J, Li Q, Wu Q. Sex Differences in Reported Diagnosis of Osteoporosis vs Other Chronic Conditions Among US Adults With Health Care Contact. *Submitted to JAMA Network Open.*

The study compares reported prior diagnosis of DXA-defined osteoporosis with reported diagnosis of measured diabetes, hypertension, and high cholesterol among US adults 50 years or older, overall and among adults with health care contact. It uses the National Health and Nutrition Examination Survey (NHANES) cycles with hip DXA: 2005-2006, 2007-2008, 2009-2010, 2013-2014, and the 2017-March 2020 prepandemic file.

## Data

All data are public NHANES files from the National Center for Health Statistics (https://wwwn.cdc.gov/nchs/nhanes/). No data are included in this repository. `download_nhanes.sh` downloads every file used into `nhanes_data/`:

DEMO, DXXFEM, OSQ, RXQ_RX, HIQ, HUQ, BMX, GHB, DIQ, BPQ, BPX (BPXO in 2017-March 2020), and TCHOL.

## Requirements

- R 4.6.1 with `haven`, `survey`, `dplyr`, `ggplot2` (`Rscript install_packages.R`)
- Python 3 with `matplotlib` (Figure)
- Node.js with the `docx` package (optional; Word formatting of Table 3)

The analyses were run on the Ohio Supercomputer Center (Slurm). The `run_*_osc.sh` scripts submit the jobs; before using them, replace `YOUR_OSC_PROJECT` with your own project code. Every R script can also be run directly with `Rscript` on any machine.

## Run order

| Step | Command | Output folder | Used for |
|---|---|---|---|
| 0 | `bash download_nhanes.sh` | `nhanes_data/` | Input data |
| 1 | `Rscript nhanes_cascade_analysis.R nhanes_data results` (or `sbatch run_analysis_osc.sh`) | `results/` | Analytic sample, osteoporosis definition, survey design; Table 2 and adjusted prevalence ratio models (eTable 2), subgroup estimates (eTable 3), participant flow |
| 2 | `Rscript nhanes_feasibility_crosscondition.R` (or `bash run_feasibility_osc.sh`) | `results_feasibility/` | Prespecified feasibility counts and NCHS reliability checks |
| 3 | `Rscript nhanes_crosscondition_analysis.R` (or `bash run_crosscondition_osc.sh`) | `results_crosscondition/` | Primary cross-condition analyses: Figure, Table 3, population-level and within-person comparisons |
| 4 | `Rscript nhanes_revision_analyses.R` (or `bash run_revision_osc.sh`) | `results_revision/` | Table 1; alternative definitions, visit gradient, discordant pairs, missing-scan weighting, medication (eTables 5-10). These analyses were added after the 16 October 2026 analysis plan and are reported as post hoc. |
| 5 | `Rscript nhanes_submission_checks.R` (or `bash run_submission_checks_osc.sh`) | `results_submission_checks/` | Prescription-inclusive definitions for all 4 conditions, paired outcome counts, person-level estimator, survey-period checks |
| 6 | `python3 make_figure1.py results_crosscondition/crosscond_A1_percentages.csv .` | — | Figure (PDF and 600-dpi TIFF) |
| 7 | `node make_table3.js results_crosscondition Table3.docx` (optional) | — | Table 3 in Word |

Steps 2-5 reuse the data-building code in `nhanes_cascade_analysis.R` (sections 1-6, plus the `fit_pr()` function), so keep all scripts in the same folder. Each run writes a `run_log.txt` summarizing its results.

## Methods notes

- Osteoporosis: femoral neck or total hip T-score of −2.5 or less (NHANES III reference for non-Hispanic White women aged 20-29 years); a male reference is used in a sensitivity analysis.
- Comparison conditions were defined by measurement at the examination (HbA1c, blood pressure, total cholesterol). Medication-inclusive definitions use the 30-day prescription inventory, which NHANES asks every participant regardless of reported diagnosis.
- All estimates account for strata, primary sampling units, and combined examination weights. Percentages use Korn-Graubard CIs and are flagged against NCHS data presentation standards.

## License

MIT (see `LICENSE`).
