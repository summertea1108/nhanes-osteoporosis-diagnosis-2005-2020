#!/usr/bin/env Rscript
# =============================================================================
# Pre-analysis checks for the revised osteoporosis paper
#
#   Step 1  Feasibility counts for the two new analyses
#           A. Cross-condition benchmark: reported diagnosis of osteoporosis vs
#              diabetes, hypertension and high cholesterol (population level and
#              within the same person)
#           B. Sex x health care contact model for reported osteoporosis diagnosis
#   Step 2  NCHS reliability of the headline reported-diagnosis estimates
#           (20.7% men / 49.4% women with health care contact; men 70+; etc.)
#
# Run:   Rscript nhanes_feasibility_crosscondition.R [data_dir] [out_dir] [main_script]
#        defaults: nhanes_data  results_feasibility  nhanes_cascade_analysis.R
#
# This script runs sections 1-6 of the main script (data build, definitions,
# survey design and helpers est_cell / nchs_flag / fit_pr) so every definition
# is identical to the main analysis. It then adds the comparison conditions.
#
# IMPORTANT (analysis plan): Step 1 reports denominators, degrees of freedom and
# adequacy flags only. It does NOT print how often the comparison conditions are
# reported, so the analysis plan for A and B can be written and dated before
# those results are seen. Osteoporosis counts are shown because they are
# already known from the main analysis.
#
# New files needed per cycle (downloaded from CDC if missing and the node has
# internet; otherwise download them on a login node):
#   GHB (HbA1c), DIQ (diabetes questionnaire), BPQ (BP/cholesterol questionnaire),
#   BPX (BP exam; BPXO = oscillometric in 2017-2020), TCHOL (total cholesterol)
#
# Comparison-condition definitions (primary = measurement only, mirroring the
# osteoporosis definition):
#   Diabetes        HbA1c >= 6.5%                         told: DIQ010 = 1 (3 "borderline" = not told)
#   Hypertension    mean SBP >= 140 or mean DBP >= 90     told: BPQ020 = 1
#   High cholesterol total cholesterol >= 240 mg/dL       told: BPQ080 = 1 (never checked = not told)
# Medication-inclusive sensitivity definitions add: insulin or diabetes pills
# (DIQ050 / DIQ070-DID070); now taking BP medicine (BPQ050A); now taking
# cholesterol medicine (BPQ100D).
# BP: mean of available readings (2005-2016 BPXSY1-4/BPXDI1-4 auscultatory;
# 2017-2020 BPXOSY1-3/BPXODI1-3 oscillometric); diastolic 0 treated as missing.
# =============================================================================

args_f     <- commandArgs(trailingOnly = TRUE)
DATA_DIR_F <- if (length(args_f) >= 1) args_f[1] else "nhanes_data"
OUT_DIR_F  <- if (length(args_f) >= 2) args_f[2] else "results_feasibility"
MAIN_F     <- if (length(args_f) >= 3) args_f[3] else "nhanes_cascade_analysis.R"
if (!file.exists(MAIN_F)) stop("Cannot find the main script: ", MAIN_F,
                               " (pass its path as the third argument)", call. = FALSE)
if (normalizePath(OUT_DIR_F, mustWork = FALSE) == normalizePath("results", mustWork = FALSE))
  stop("Use a separate output folder so the main run log is not overwritten.", call. = FALSE)

# -----------------------------------------------------------------------------
# 0. Make sure the extra NHANES files are present
# -----------------------------------------------------------------------------
CYC_F <- data.frame(
  cycle = c("2005-2006", "2007-2008", "2009-2010", "2013-2014", "2017-2020"),
  pre   = c("", "", "", "", "P_"),
  suf   = c("_D", "_E", "_F", "_H", ""),
  year  = c(2005, 2007, 2009, 2013, 2017),
  stringsAsFactors = FALSE
)
extra_name <- function(base, i) {
  if (base == "BPX" && CYC_F$cycle[i] == "2017-2020") base <- "BPXO"
  paste0(CYC_F$pre[i], base, CYC_F$suf[i])
}
EXTRA <- c("GHB", "DIQ", "BPQ", "BPX", "TCHOL")
need <- unlist(lapply(seq_len(nrow(CYC_F)), function(i) {
  f <- vapply(EXTRA, extra_name, character(1), i = i)
  setNames(rep(CYC_F$year[i], length(f)), f)
}))
dir.create(DATA_DIR_F, showWarnings = FALSE, recursive = TRUE)
missing_f <- names(need)[!file.exists(file.path(DATA_DIR_F, paste0(names(need), ".xpt")))]
for (f in missing_f) {
  url <- sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s.xpt", need[[f]], f)
  cat("Downloading ", f, " ...\n", sep = "")
  tryCatch(utils::download.file(url, file.path(DATA_DIR_F, paste0(f, ".xpt")), mode = "wb", quiet = TRUE),
           error = function(e) cat("  could not download ", f, ": ", conditionMessage(e), "\n", sep = ""))
}
still <- names(need)[!file.exists(file.path(DATA_DIR_F, paste0(names(need), ".xpt")))]
if (length(still))
  stop("Missing data files (download them from the CDC NHANES site into ", DATA_DIR_F, "): ",
       paste0(still, ".xpt", collapse = ", "), call. = FALSE)

# -----------------------------------------------------------------------------
# 1. Run sections 1-6 of the main script (definitions, data build, design, helpers)
# -----------------------------------------------------------------------------
main_lines <- readLines(MAIN_F, warn = FALSE)
cut_at <- grep("^# 7\\. Tables", main_lines)
if (length(cut_at) != 1) stop("Could not find the '# 7. Tables' header in ", MAIN_F, call. = FALSE)
commandArgs <- function(trailingOnly = FALSE) c(DATA_DIR_F, OUT_DIR_F)   # feeds the main script's arguments
eval(parse(text = main_lines[seq_len(cut_at - 1)]), envir = globalenv())
rm(commandArgs)
# fit_pr() is defined later in the main script (Table 4 block); load just that function
fp_start <- grep("^fit_pr <- function", main_lines)
fp_end   <- grep("^M <- list\\(", main_lines)
if (length(fp_start) != 1 || length(fp_end) != 1 || fp_end <= fp_start)
  stop("Could not locate fit_pr() in ", MAIN_F, call. = FALSE)
eval(parse(text = main_lines[fp_start:(fp_end - 1)]), envir = globalenv())
# objects now available: all, exam, D, A, CYCLES, read_nh, getcol, num, yn, say,
# write_out, sub_des, domain_df, est_cell, nchs_flag, fit_pr, SEXES, sex_cond

say("")
say("=== Feasibility and reliability checks ===")
say("Main script: ", normalizePath(MAIN_F), "  (sections 1-6)")

# -----------------------------------------------------------------------------
# 2. Add the comparison conditions
# -----------------------------------------------------------------------------
yn_dm <- function(x) { x <- num(x); ifelse(x == 1, 1, ifelse(x %in% c(2, 3), 0, NA_real_)) }
cols_num <- function(d, pattern) {
  cc <- grep(pattern, names(d), value = TRUE)
  if (!length(cc)) return(matrix(NA_real_, nrow(d), 1))
  do.call(cbind, lapply(cc, function(k) num(d[[k]])))
}
row_mean <- function(m) { r <- rowMeans(m, na.rm = TRUE); r[is.nan(r)] <- NA; r }

build_extra <- function(i) {
  cy  <- CYC_F[i, ]
  ghb <- read_nh(extra_name("GHB", i)); diq <- read_nh(extra_name("DIQ", i))
  bpq <- read_nh(extra_name("BPQ", i)); bpx <- read_nh(extra_name("BPX", i))
  tch <- read_nh(extra_name("TCHOL", i))

  g  <- data.frame(SEQN = ghb$SEQN, hba1c = num(getcol(ghb, "LBXGH")))
  dq <- data.frame(SEQN = diq$SEQN,
                   told_dm    = yn_dm(getcol(diq, "DIQ010")),
                   dm_insulin = yn(getcol(diq, "DIQ050")),
                   dm_pills   = yn(getcol(diq, "DIQ070", "DID070")))
  checked   <- yn(getcol(bpq, "BPQ060"))
  told_chol <- yn(getcol(bpq, "BPQ080"))
  told_chol <- ifelse(is.na(told_chol) & checked %in% 0, 0, told_chol)   # never checked -> never told
  bq <- data.frame(SEQN = bpq$SEQN,
                   told_htn  = yn(getcol(bpq, "BPQ020")),
                   htn_med   = yn(getcol(bpq, "BPQ050A")),
                   told_chol = told_chol,
                   chol_med  = yn(getcol(bpq, "BPQ100D")))
  if (cy$cycle == "2017-2020") {
    sy <- cols_num(bpx, "^BPXOSY[1-3]$"); di <- cols_num(bpx, "^BPXODI[1-3]$")
  } else {
    sy <- cols_num(bpx, "^BPXSY[1-4]$");  di <- cols_num(bpx, "^BPXDI[1-4]$")
  }
  di[!is.na(di) & di == 0] <- NA
  bp <- data.frame(SEQN = bpx$SEQN, sbp = row_mean(sy), dbp = row_mean(di))
  tc <- data.frame(SEQN = tch$SEQN, tchol = num(getcol(tch, "LBXTC")))

  out <- Reduce(function(a, b) merge(a, b, by = "SEQN", all = TRUE), list(g, dq, bq, bp, tc))
  out$cycle <- cy$cycle
  say("  ", cy$cycle, ": HbA1c ", sum(!is.na(out$hba1c)), ", BP ", sum(!is.na(out$sbp)),
      ", cholesterol ", sum(!is.na(out$tchol)), " records")
  out
}
say("Reading comparison-condition files")
extra <- do.call(rbind, lapply(seq_len(nrow(CYC_F)), build_extra))

all <- merge(all, extra, by = c("SEQN", "cycle"), all.x = TRUE)
all <- all %>% mutate(
  dm_meas   = ifelse(is.na(hba1c), NA, as.numeric(hba1c >= 6.5)),
  htn_meas  = ifelse(is.na(sbp), NA, as.numeric(sbp >= 140 | (!is.na(dbp) & dbp >= 90))),
  chol_meas = ifelse(is.na(tchol), NA, as.numeric(tchol >= 240)),
  dm_medinc   = ifelse(is.na(dm_meas) & !(dm_insulin %in% 1 | dm_pills %in% 1), NA,
                       as.numeric(dm_meas %in% 1 | dm_insulin %in% 1 | dm_pills %in% 1)),
  htn_medinc  = ifelse(is.na(htn_meas) & !(htn_med %in% 1), NA,
                       as.numeric(htn_meas %in% 1 | htn_med %in% 1)),
  chol_medinc = ifelse(is.na(chol_meas) & !(chol_med %in% 1), NA,
                       as.numeric(chol_meas %in% 1 | chol_med %in% 1))
)

# rebuild the designs on the enlarged file (same weights, strata, PSUs)
exam <- all[all$examined, ]
D  <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_pooled, nest = TRUE, data = exam)
A  <- subset(D, analytic)
OP <- subset(A, op == 1)
ENG <- sub_des(OP, OP$variables$engaged %in% "Engaged in care")

COND <- data.frame(
  condition = c("Osteoporosis", "Diabetes", "Hypertension", "High cholesterol"),
  measured  = c("op", "dm_meas", "htn_meas", "chol_meas"),
  medinc    = c(NA, "dm_medinc", "htn_medinc", "chol_medinc"),
  told      = c("told_op", "told_dm", "told_htn", "told_chol"),
  stringsAsFactors = FALSE
)
MIN_EVENTS <- 10   # each outcome (told / not told) should have at least this many people
MIN_N      <- 30   # NCHS minimum sample size for a proportion

# counts in one cell; outcome split shown only when show_split = TRUE
count_cell <- function(v, has, told, show_split) {
  den   <- v[[has]] %in% 1
  known <- den & !is.na(v[[told]])
  yes   <- sum(known & v[[told]] == 1); no <- sum(known) - yes
  data.frame(
    n_with_condition = sum(den),
    n_report_known   = sum(known),
    told_yes = if (show_split) yes else NA,
    told_no  = if (show_split) no else NA,
    domain_df = if (sum(known) > 0) domain_df(v[known, ]) else NA,
    n_ge_30 = sum(known) >= MIN_N,
    both_outcomes_ge_10 = yes >= MIN_EVENTS & no >= MIN_EVENTS
  )
}

# -----------------------------------------------------------------------------
# STEP 1A(i). Population-level benchmark: who has each condition, by sex and contact
# -----------------------------------------------------------------------------
say("")
say("--- Step 1A(i): population-level benchmark denominators (age 50+, valid hip DXA) ---")
va <- A$variables
f1 <- list()
for (smp in c("All analytic adults 50+", "Health care contact")) {
  in_smp <- if (smp == "Health care contact") va$engaged %in% "Engaged in care" else rep(TRUE, nrow(va))
  for (s in SEXES) {
    v <- va[which(in_smp & sex_cond(va, s)), ]
    for (k in seq_len(nrow(COND))) for (def in c("Primary (measured)", "Medication-inclusive")) {
      has <- if (def == "Primary (measured)") COND$measured[k] else COND$medinc[k]
      if (is.na(has)) next
      f1[[length(f1) + 1]] <- cbind(
        data.frame(sample = smp, sex = s, condition = COND$condition[k], definition = def),
        count_cell(v, has, COND$told[k], show_split = COND$condition[k] == "Osteoporosis"))
    }
  }
}
f1 <- do.call(rbind, f1)
write_out(f1, "feasibility_1a_population_benchmark.csv")

# -----------------------------------------------------------------------------
# STEP 1A(ii). Within-person comparison: adults with osteoporosis who ALSO have
# a measured comparison condition, with both reports known
# -----------------------------------------------------------------------------
say("")
say("--- Step 1A(ii): within-person comparison (osteoporosis + comparison condition) ---")
vo <- OP$variables
f2 <- list()
for (smp in c("All with osteoporosis", "Osteoporosis + health care contact")) {
  in_smp <- if (smp == "All with osteoporosis") rep(TRUE, nrow(vo)) else vo$engaged %in% "Engaged in care"
  for (s in SEXES) {
    v <- vo[which(in_smp & sex_cond(vo, s)), ]
    rows_total <- 0; persons_any <- rep(FALSE, nrow(v))
    for (k in 2:nrow(COND)) {
      has_both  <- v[[COND$measured[k]]] %in% 1
      known     <- has_both & !is.na(v$told_op) & !is.na(v[[COND$told[k]]])
      discord   <- known & v$told_op != v[[COND$told[k]]]
      rows_total  <- rows_total + sum(known)
      persons_any <- persons_any | known
      f2[[length(f2) + 1]] <- data.frame(
        sample = smp, sex = s, comparison = COND$condition[k],
        n_osteoporosis = nrow(v),
        n_also_has_condition = sum(has_both),
        n_both_reports_known = sum(known),
        n_discordant_pairs = sum(discord),   # pairs that carry the within-person information
        domain_df = if (sum(known) > 0) domain_df(v[known, ]) else NA,
        n_ge_30 = sum(known) >= MIN_N,
        discordant_ge_10 = sum(discord) >= MIN_EVENTS)
    }
    f2[[length(f2) + 1]] <- data.frame(
      sample = smp, sex = s, comparison = "Any of the three (stacked rows)",
      n_osteoporosis = nrow(v), n_also_has_condition = NA,
      n_both_reports_known = rows_total, n_discordant_pairs = NA,
      domain_df = if (any(persons_any)) domain_df(v[persons_any, ]) else NA,
      n_ge_30 = sum(persons_any) >= MIN_N, discordant_ge_10 = NA)
    f2[[length(f2)]]$n_also_has_condition <- sum(persons_any)   # here: persons contributing >= 1 row
  }
}
f2 <- do.call(rbind, f2)
write_out(f2, "feasibility_1a_within_person.csv")

# -----------------------------------------------------------------------------
# STEP 1B. Sex x contact model: cells for reported osteoporosis diagnosis
# -----------------------------------------------------------------------------
say("")
say("--- Step 1B: sex x health care contact cells (adults with osteoporosis) ---")
f3 <- list()
for (s in c("Women", "Men")) for (cg in c("Engaged in care", "Not engaged", "Contact missing")) {
  in_cg <- if (cg == "Contact missing") is.na(vo$engaged) else vo$engaged %in% cg
  v <- vo[which(vo$sex == s & in_cg), ]
  for (pp in c("All periods", "Early", "Late")) {
    vp <- if (pp == "All periods") v else v[which(v$period == pp), ]
    vp$one <- rep(1, nrow(vp))
    f3[[length(f3) + 1]] <- cbind(data.frame(sex = s, contact = cg, period = pp),
                                  count_cell(vp, "one", "told_op", show_split = TRUE)[, -1])
  }
}
f3 <- do.call(rbind, f3)
write_out(f3, "feasibility_1b_sex_by_contact_cells.csv")

# -----------------------------------------------------------------------------
# STEP 2. NCHS reliability of the headline reported-diagnosis estimates
# -----------------------------------------------------------------------------
say("")
say("--- Step 2: NCHS reliability of reported-diagnosis estimates ---")
rel <- list()
add_rel <- function(label, cond) {
  des <- sub_des(OP, cond)
  r <- est_cell(des, "told_op")
  w <- r$ucl - r$lcl
  rel[[length(rel) + 1]] <<- cbind(data.frame(estimate = label), r[, c("n", "events", "pct", "lcl", "ucl")],
    data.frame(ci_width_pts = w,
               rel_ci_width_pct = 100 * w / r$pct,
               rel_ci_width_complement_pct = 100 * w / (100 - r$pct)),
    r[, c("n_eff", "df", "nchs_flag", "nchs_flag_complement")])
}
eng  <- vo$engaged %in% "Engaged in care"
men  <- vo$sex == "Men"; wom <- vo$sex == "Women"
add_rel("All with osteoporosis, both sexes", rep(TRUE, nrow(vo)))
add_rel("All with osteoporosis, women", wom)
add_rel("All with osteoporosis, men", men)
add_rel("HEADLINE: health care contact, women", eng & wom)
add_rel("HEADLINE: health care contact, men", eng & men)
add_rel("Health care contact, both sexes", eng)
add_rel("No health care contact, women", vo$engaged %in% "Not engaged" & wom)
add_rel("No health care contact, men", vo$engaged %in% "Not engaged" & men)
add_rel("Men 50-69", men & vo$age < 70)
add_rel("Men 70+", men & vo$age >= 70)
add_rel("Men 70+, health care contact", men & vo$age >= 70 & eng)
add_rel("Guideline-age (women 65+), health care contact", wom & vo$age >= 65 & eng)
add_rel("Guideline-age (men 70+), health care contact", men & vo$age >= 70 & eng)
for (pp in c("Early", "Late")) for (s in c("Women", "Men"))
  add_rel(paste0("Health care contact, ", s, ", ", pp, " era"), eng & vo$sex == s & vo$period == pp)
for (rc in c("NH White", "NH Black", "Hispanic")) {
  add_rel(paste0(rc, ", all with osteoporosis"), vo$race %in% rc)
  add_rel(paste0(rc, ", health care contact"), vo$race %in% rc & eng)
}
rel <- do.call(rbind, rel)
write_out(rel, "reliability_step2_reported_diagnosis.csv")

# sex prevalence ratio in the contact subgroup (the fallback if the men's % is not reportable)
pr_eng <- fit_pr(ENG, "told_op", "sex + age_grp + race_m", "Sex gap, age/race adjusted", "Osteoporosis, health care contact")
pr_eng <- pr_eng[!is.na(pr_eng$term) & pr_eng$term == "sexMen", ]
write_out(pr_eng, "reliability_step2_sex_pr_contact.csv")

# -----------------------------------------------------------------------------
# Plain-language summary
# -----------------------------------------------------------------------------
say("")
say("=================== SUMMARY ===================")
verdict <- function(flag) if (grepl("^Reliable$", flag)) "REPORTABLE" else
  if (grepl("^Reliable", flag)) "REPORTABLE WITH FOOTNOTE (df < 8)" else "NOT REPORTABLE as a percentage"
say("Step 2 - reliability (NCHS 2017 standards):")
for (k in seq_len(nrow(rel))) {
  r <- rel[k, ]
  if (!grepl("HEADLINE|Men 70\\+|Guideline", r$estimate)) next
  say(sprintf("  %-50s %5.1f%% (%.1f-%.1f), n=%d, events=%d, width %.1f pts, eff n %.0f, df %d -> %s",
              r$estimate, r$pct, r$lcl, r$ucl, r$n, r$events, r$ci_width_pts, r$n_eff, r$df,
              verdict(paste(r$nchs_flag))))
}
hm <- rel[rel$estimate == "HEADLINE: health care contact, men", ]
hw <- rel[rel$estimate == "HEADLINE: health care contact, women", ]
say(sprintf("  Matches manuscript (men 20.7%%, women 49.4%%)? men %s, women %s",
            ifelse(abs(hm$pct - 20.7) < 0.05, "YES", sprintf("NO (%.1f%%)", hm$pct)),
            ifelse(abs(hw$pct - 49.4) < 0.05, "YES", sprintf("NO (%.1f%%)", hw$pct))))
if (nrow(pr_eng))
  say(sprintf("  Sex PR (men vs women) with health care contact: %.2f (%.2f-%.2f), n=%d",
              pr_eng$PR, pr_eng$lcl, pr_eng$ucl, pr_eng$n))

say("")
say("Step 1A - cross-condition benchmark, men with health care contact:")
for (cn in COND$condition[-1]) {
  r <- f1[f1$sample == "Health care contact" & f1$sex == "Men" & f1$condition == cn &
            f1$definition == "Primary (measured)", ]
  say(sprintf("  Population level, %-17s n=%5d  df=%s  -> %s", cn, r$n_report_known, r$domain_df,
              ifelse(r$n_ge_30 & r$both_outcomes_ge_10, "FEASIBLE", "THIN")))
}
for (s in c("Men", "Women")) for (cn in c(COND$condition[-1], "Any of the three (stacked rows)")) {
  r <- f2[f2$sample == "Osteoporosis + health care contact" & f2$sex == s & f2$comparison == cn, ]
  ok <- if (is.na(r$discordant_ge_10)) r$n_ge_30 else r$n_ge_30 & r$discordant_ge_10
  say(sprintf("  Within person, %-5s %-32s persons=%3d discordant=%s  -> %s", s, cn,
              if (grepl("^Any", cn)) r$n_also_has_condition else r$n_both_reports_known,
              ifelse(is.na(r$n_discordant_pairs), "-", r$n_discordant_pairs),
              ifelse(ok, "FEASIBLE", "THIN")))
}

say("")
say("Step 1B - sex x contact cells (all periods):")
for (k in which(f3$period == "All periods"))
  say(sprintf("  %-5s %-16s n=%3d  told yes=%3s  no=%3s  -> %s", f3$sex[k], f3$contact[k],
              f3$n_report_known[k], f3$told_yes[k], f3$told_no[k],
              ifelse(f3$both_outcomes_ge_10[k], "OK", "THIN")))

say("")
say("How to read this:")
say("  REPORTABLE         -> lead with the percentage.")
say("  NOT REPORTABLE     -> lead with the sex prevalence ratio; give the % only with a caution.")
say("  FEASIBLE / THIN    -> THIN within-person cells for men mean the within-person analysis")
say("                        rests on hypertension (or both sexes combined) only; plan accordingly.")
say("  Comparison-condition reporting rates are intentionally NOT shown here (write the plan first).")
say("Done. Files in ", normalizePath(OUT_DIR))
