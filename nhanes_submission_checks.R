#!/usr/bin/env Rscript
# =============================================================================
# Pre-submission checks requested in review (osteoporosis sex-difference paper)
#
#   C1  Parallel medication-inclusive definitions for ALL FOUR conditions, from the
#       30-day prescription inventory (RXQ_RX), which every participant is asked
#       regardless of diagnosis. Replaces the questionnaire-based medication items,
#       which NHANES asks only of people who already report the diagnosis.
#       Reruns the population-level (A1) and within-person (A2) comparisons.
#   C2  Skip-pattern audit: who answered each medication item, by reported diagnosis
#   C3  Missing-data audit: missing measurements or answers are excluded, never coded no
#   C4  Unique participants vs condition records (A1) and persons vs pairs (A2), by sex;
#       paired outcome counts (both / neither / osteoporosis only / comparison only)
#   C5  Variance check: within-person PR from the stacked pair records vs a
#       person-level ratio estimator (one record per person)
#   C6  Survey-period adjustment and era-specific estimates for the cross-condition models
#   C7  Denominators: screening-age groups and headline counts (total, known status, events)
#
# Run:   Rscript nhanes_submission_checks.R [data_dir] [out_dir] [main_script] [crosscondition_script]
#        defaults: nhanes_data  results_submission_checks  nhanes_cascade_analysis.R
#                  nhanes_crosscondition_analysis.R
# =============================================================================

a_c <- commandArgs(trailingOnly = TRUE)
DATA_C  <- if (length(a_c) >= 1) a_c[1] else "nhanes_data"
OUT_C   <- if (length(a_c) >= 2) a_c[2] else "results_submission_checks"
MAIN_C  <- if (length(a_c) >= 3) a_c[3] else "nhanes_cascade_analysis.R"
CROSS_C <- if (length(a_c) >= 4) a_c[4] else "nhanes_crosscondition_analysis.R"
for (f in c(MAIN_C, CROSS_C)) if (!file.exists(f)) stop("Cannot find ", f, call. = FALSE)

# -----------------------------------------------------------------------------
# 0. Sections 0-3 (setup and helpers) of the cross-condition script
# -----------------------------------------------------------------------------
cl <- readLines(CROSS_C, warn = FALSE)
a1_at  <- grep("^# A1\\. Population-level comparison", cl)
sum_at <- grep("^# 4\\. Results summary", cl)
if (length(a1_at) != 1 || length(sum_at) != 1) stop("Unexpected layout of ", CROSS_C, call. = FALSE)
commandArgs <- function(trailingOnly = FALSE) c(DATA_C, OUT_C, MAIN_C)
eval(parse(text = cl[seq_len(a1_at - 2)]), envir = globalenv())
if (exists("commandArgs", envir = globalenv(), inherits = FALSE)) rm(commandArgs, envir = globalenv())

say("")
say("=== Pre-submission checks ===")

# -----------------------------------------------------------------------------
# 1. Prescription-inventory medication flags (RXQ_RX; asked of all participants)
# -----------------------------------------------------------------------------
RX_DM <- paste0("\\b(METFORMIN|GLIPIZIDE|GLYBURIDE|GLIMEPIRIDE|CHLORPROPAMIDE|TOLBUTAMIDE|TOLAZAMIDE|",
  "PIOGLITAZONE|ROSIGLITAZONE|SITAGLIPTIN|SAXAGLIPTIN|LINAGLIPTIN|ALOGLIPTIN|INSULIN|EXENATIDE|",
  "LIRAGLUTIDE|DULAGLUTIDE|SEMAGLUTIDE|LIXISENATIDE|ALBIGLUTIDE|CANAGLIFLOZIN|DAPAGLIFLOZIN|",
  "EMPAGLIFLOZIN|ERTUGLIFLOZIN|ACARBOSE|MIGLITOL|REPAGLINIDE|NATEGLINIDE|PRAMLINTIDE)\\b")
RX_HTN <- paste0("\\b(BENAZEPRIL|CAPTOPRIL|ENALAPRIL|FOSINOPRIL|LISINOPRIL|MOEXIPRIL|PERINDOPRIL|QUINAPRIL|",
  "RAMIPRIL|TRANDOLAPRIL|CANDESARTAN|EPROSARTAN|IRBESARTAN|LOSARTAN|OLMESARTAN|TELMISARTAN|VALSARTAN|",
  "AZILSARTAN|ALISKIREN|AMLODIPINE|FELODIPINE|ISRADIPINE|NICARDIPINE|NIFEDIPINE|NISOLDIPINE|DILTIAZEM|",
  "VERAPAMIL|HYDROCHLOROTHIAZIDE|CHLORTHALIDONE|CHLOROTHIAZIDE|INDAPAMIDE|METOLAZONE|FUROSEMIDE|",
  "BUMETANIDE|TORSEMIDE|SPIRONOLACTONE|EPLERENONE|TRIAMTERENE|AMILORIDE|ACEBUTOLOL|ATENOLOL|BETAXOLOL|",
  "BISOPROLOL|CARVEDILOL|LABETALOL|METOPROLOL|NADOLOL|NEBIVOLOL|PINDOLOL|PROPRANOLOL|DOXAZOSIN|",
  "PRAZOSIN|TERAZOSIN|CLONIDINE|GUANFACINE|METHYLDOPA|HYDRALAZINE|MINOXIDIL|RESERPINE)\\b")
RX_LIP <- paste0("\\b(ATORVASTATIN|FLUVASTATIN|LOVASTATIN|PITAVASTATIN|PRAVASTATIN|ROSUVASTATIN|SIMVASTATIN|",
  "CERIVASTATIN|EZETIMIBE|FENOFIBRATE|FENOFIBRIC|GEMFIBROZIL|NIACIN|COLESEVELAM|COLESTIPOL|",
  "CHOLESTYRAMINE|EVOLOCUMAB|ALIROCUMAB|ICOSAPENT|OMEGA-3-ACID|BEMPEDOIC)\\b")

rx_flags <- do.call(rbind, lapply(seq_len(nrow(CYCLES)), function(i) {
  cy <- CYCLES[i, ]
  rx <- read_nh(paste0(cy$pre, "RXQ_RX", cy$suf))
  rx$drug <- toupper(ifelse(is.na(rx$RXDDRUG), "", as.character(rx$RXDDRUG)))
  rx$use  <- yn(rx$RXDUSE)
  p <- rx %>% group_by(SEQN) %>%
    summarise(use = first(use), dm = any(grepl(RX_DM, drug, perl = TRUE)),
              htn = any(grepl(RX_HTN, drug, perl = TRUE)), lip = any(grepl(RX_LIP, drug, perl = TRUE))) %>%
    as.data.frame()
  f <- function(x) ifelse(is.na(p$use), NA, ifelse(p$use == 0, 0, as.numeric(x)))
  data.frame(SEQN = p$SEQN, cycle = cy$cycle, rx_dm = f(p$dm), rx_htn = f(p$htn), rx_lip = f(p$lip))
}))
all <- merge(all, rx_flags, by = c("SEQN", "cycle"), all.x = TRUE)
inc <- function(meas, rx) ifelse(is.na(meas) & !(rx %in% 1), NA, as.numeric(meas %in% 1 | rx %in% 1))
all <- all %>% mutate(
  op_rxinc   = ifelse(is.na(op), NA, as.numeric(op %in% 1 | rx_tx %in% 1)),   # DXA still required
  dm_rxinc   = inc(dm_meas, rx_dm),
  htn_rxinc  = inc(htn_meas, rx_htn),
  chol_rxinc = inc(chol_meas, rx_lip)
)
exam <- all[all$examined, ]
D  <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_pooled, nest = TRUE, data = exam)
A  <- subset(D, analytic)
OP <- subset(A, op == 1)

COND <- data.frame(
  condition = c("Osteoporosis", "Diabetes", "Hypertension", "High cholesterol"),
  measured  = c("op", "dm_meas", "htn_meas", "chol_meas"),
  medinc    = c("op_rxinc", "dm_rxinc", "htn_rxinc", "chol_rxinc"),
  told      = c("told_op", "told_dm", "told_htn", "told_chol"),
  stringsAsFactors = FALSE
)
DEFS <- c("Primary (measured)", "Prescription-inclusive, all 4 conditions")
PRIMARY_DEF <- DEFS[1]
KEEP <- c("SEQN", "psu", "strata", "w_pooled", "sex", "age_grp", "race_m", "analytic", "op", "engaged", "period", "cycle")
COND_VARS <- unique(na.omit(c(COND$measured, COND$medinc, COND$told)))
base <- exam[, c(KEEP, COND_VARS)]
COND_LEV <- COND$condition; COMPS <- COND$condition[-1]

# -----------------------------------------------------------------------------
# C1. Rerun A1 and A2 with the parallel prescription-inclusive definitions.
#     In A2 the osteoporosis side of each pair also uses the definition in force.
# -----------------------------------------------------------------------------
say("")
say("--- C1: parallel prescription-inclusive definitions (A1 and A2 rerun) ---")
blk <- cl[(a1_at - 1):(sum_at - 2)]
blk <- gsub("base$op %in% 1", "base[[has_var(1, df_)]] %in% 1", blk, fixed = TRUE)
blk <- gsub("\"crosscond_", "\"check_C1_crosscond_", blk, fixed = TRUE)
eval(parse(text = blk), envir = globalenv())
say("  Added by the prescription inventory (analytic sample, unweighted): osteoporosis ",
    sum(A$variables$op_rxinc %in% 1) - sum(A$variables$op %in% 1),
    "; diabetes ", sum(A$variables$dm_rxinc %in% 1) - sum(A$variables$dm_meas %in% 1),
    "; hypertension ", sum(A$variables$htn_rxinc %in% 1) - sum(A$variables$htn_meas %in% 1),
    "; high cholesterol ", sum(A$variables$chol_rxinc %in% 1) - sum(A$variables$chol_meas %in% 1))

# -----------------------------------------------------------------------------
# C2. Skip-pattern audit (analytic sample, unweighted)
# -----------------------------------------------------------------------------
va <- A$variables
skip_row <- function(med, told, label) data.frame(
  item = label,
  answered_if_told_yes = sum(!is.na(va[[med]]) & va[[told]] %in% 1), told_yes = sum(va[[told]] %in% 1),
  answered_if_told_no = sum(!is.na(va[[med]]) & va[[told]] %in% 0), told_no = sum(va[[told]] %in% 0),
  med_yes_if_told_no = sum(va[[med]] %in% 1 & va[[told]] %in% 0))
c2 <- rbind(
  skip_row("htn_med", "told_htn", "BPQ050A now taking BP medicine (questionnaire)"),
  skip_row("chol_med", "told_chol", "BPQ100D now taking cholesterol medicine (questionnaire)"),
  skip_row("dm_pills", "told_dm", "DIQ070/DID070 diabetes pills (questionnaire)"),
  skip_row("dm_insulin", "told_dm", "DIQ050 insulin (questionnaire)"),
  skip_row("rx_htn", "told_htn", "Antihypertensive in prescription inventory"),
  skip_row("rx_lip", "told_chol", "Lipid-lowering drug in prescription inventory"),
  skip_row("rx_dm", "told_dm", "Antidiabetic drug in prescription inventory"),
  skip_row("rx_tx", "told_op", "Osteoporosis drug in prescription inventory"))
write_out(c2, "check_C2_skip_patterns.csv")

# -----------------------------------------------------------------------------
# C3. Missing-data audit: primary sample (hip DXA + health care contact)
# -----------------------------------------------------------------------------
prim <- va$engaged %in% "Engaged in care"
c3 <- list()
for (k in seq_len(nrow(COND))) for (sx in SEXES) {
  s <- prim & sex_cond(va, sx); m <- va[[COND$measured[k]]]; t <- va[[COND$told[k]]]
  c3[[length(c3) + 1]] <- data.frame(condition = COND$condition[k], sex = sx, in_sample = sum(s),
    measurement_missing = sum(s & is.na(m)), measured_positive = sum(s & m %in% 1),
    positive_answer_missing_excluded = sum(s & m %in% 1 & is.na(t)),
    analyzed = sum(s & m %in% 1 & !is.na(t)), analyzed_yes = sum(s & m %in% 1 & t %in% 1),
    analyzed_no = sum(s & m %in% 1 & t %in% 0))
}
c3 <- do.call(rbind, c3)
write_out(c3, "check_C3_missing_audit.csv")

# -----------------------------------------------------------------------------
# C4. Participants vs records (A1), persons vs pairs and paired outcome counts (A2)
# -----------------------------------------------------------------------------
S1 <- mk_design(stack_A1(PRIMARY_DEF), "w_pooled"); vs <- S1$variables
dom1 <- in_sample(vs, SAMPLES[1]) & vs$has %in% 1 & !is.na(vs$told)
c4a <- do.call(rbind, lapply(SEXES, function(sx) {
  d <- dom1 & sex_cond(vs, sx)
  data.frame(sex = sx, condition_records = sum(d), unique_participants = length(unique(vs$SEQN[d])),
             records_reporting_diagnosis = sum(d & vs$told == 1),
             participants_with_1_condition = sum(table(vs$SEQN[d]) == 1),
             participants_with_2plus = sum(table(vs$SEQN[d]) >= 2))
}))
write_out(c4a, "check_C4a_A1_records_vs_participants.csv")

b <- base; bc <- in_sample(b, SAMPLES[1])
c4b <- list()
for (cn in c(COMPS, "Any")) for (sx in SEXES) {
  ks <- if (cn == "Any") 2:4 else which(COND$condition == cn)
  rows <- do.call(rbind, lapply(ks, function(k) {
    e <- bc & b$op %in% 1 & !is.na(b$told_op) & b[[COND$measured[k]]] %in% 1 & !is.na(b[[COND$told[k]]]) & sex_cond(b, sx)
    data.frame(SEQN = b$SEQN[e], o = b$told_op[e], c = b[[COND$told[k]]][e])
  }))
  c4b[[length(c4b) + 1]] <- data.frame(comparison = cn, sex = sx, persons = length(unique(rows$SEQN)), pairs = nrow(rows),
    both_reported = sum(rows$o == 1 & rows$c == 1), neither_reported = sum(rows$o == 0 & rows$c == 0),
    osteoporosis_only = sum(rows$o == 1 & rows$c == 0), comparison_only = sum(rows$o == 0 & rows$c == 1))
}
c4b <- do.call(rbind, c4b)
write_out(c4b, "check_C4b_paired_outcome_counts.csv")

# -----------------------------------------------------------------------------
# C5. Variance check for the within-person comparison: person-level ratio estimator
#     (osteoporosis report / mean of the person's comparison reports, one record per
#     person) vs the stacked pair-record model. Point estimates must be identical.
# -----------------------------------------------------------------------------
elig_k <- lapply(2:4, function(k) b$analytic %in% TRUE & b$op %in% 1 & !is.na(b$told_op) &
                   b[[COND$measured[k]]] %in% 1 & !is.na(b[[COND$told[k]]]))
kc <- Reduce(`+`, elig_k)
cs <- Reduce(`+`, lapply(1:3, function(j) ifelse(elig_k[[j]], b[[COND$told[j + 1]]], 0)))
exam$comp_mean <- ifelse(kc > 0, cs / pmax(kc, 1), NA)
exam$in_pairs  <- kc > 0 & in_sample(exam, SAMPLES[1])
Dp <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_pooled, nest = TRUE, data = exam)
c5 <- list()
P2 <- mk_design(stack_A2(PRIMARY_DEF, COMPS), "w_pair"); vp <- P2$variables
for (sx in SEXES) {
  dd <- sub_des(Dp, exam$in_pairs & sex_cond(exam, sx))
  rr <- svyratio(~told_op, ~comp_mean, dd)
  est <- as.numeric(coef(rr)); se <- as.numeric(SE(rr)); df <- degf(Dp)
  st <- fit_model(sub_des(P2, in_sample(vp, SAMPLES[1]) & vp$in_pair & sex_cond(vp, sx)), told ~ type, "stacked", SAMPLES[1], PRIMARY_DEF)
  st <- st[!is.na(st$term) & st$term == "typeOsteoporosis", ]
  c5[[length(c5) + 1]] <- data.frame(sex = sx, persons = sum(exam$in_pairs & sex_cond(exam, sx)),
    person_level_PR = est, person_level_lcl = est * exp(-qt(0.975, df) * se / est),
    person_level_ucl = est * exp(qt(0.975, df) * se / est),
    stacked_PR = st$PR, stacked_lcl = st$lcl, stacked_ucl = st$ucl)
}
c5 <- do.call(rbind, c5)
write_out(c5, "check_C5_variance_person_vs_stacked.csv")

# -----------------------------------------------------------------------------
# C6. Survey period: adjusted models, condition x period test, era-specific PRs
# -----------------------------------------------------------------------------
c6 <- list()
grab <- function(r, lab) { r <- r[!is.na(r$term) & grepl("^cond", r$term) & !grepl(":", r$term), , drop = FALSE]
  if (!nrow(r)) return(NULL); r <- invert(r); r$contrast <- paste0("Osteoporosis vs ", sub("^cond", "", r$term)); r$analysis <- lab; r }
for (sx in SEXES) {
  f <- if (sx == "Both") told ~ cond + sex + age_grp + race_m + period else told ~ cond + age_grp + race_m + period
  c6[[length(c6) + 1]] <- grab(cbind(fit_model(sub_des(S1, dom1 & sex_cond(vs, sx)), f, sx, SAMPLES[1], PRIMARY_DEF)),
                               paste("Adjusted for survey period,", sx))
}
rx6 <- fit_model(sub_des(S1, dom1), told ~ cond * sex + age_grp + race_m + period, "Condition x sex", SAMPLES[1], PRIMARY_DEF)
rx6 <- rx6[!is.na(rx6$term) & grepl("^cond.*:sexMen$", rx6$term), , drop = FALSE]
if (nrow(rx6)) { rx6 <- invert(rx6); rx6$contrast <- paste0("Ratio of PRs (men/women), osteoporosis vs ",
  sub(":sexMen$", "", sub("^cond", "", rx6$term))); rx6$analysis <- "Adjusted for survey period, condition x sex"; c6[[length(c6) + 1]] <- rx6 }
for (pp in c("Early", "Late")) for (sx in SEXES) {
  f <- if (sx == "Both") told ~ cond + sex + age_grp + race_m else told ~ cond + age_grp + race_m
  c6[[length(c6) + 1]] <- grab(fit_model(sub_des(S1, dom1 & sex_cond(vs, sx) & vs$period == pp), f, sx, SAMPLES[1], PRIMARY_DEF),
                               paste0("Within era ", pp, ", ", sx))
}
c6 <- do.call(rbind, c6)
dd6 <- S1[dom1, ]
for (tv in c("age_grp", "race_m")) dd6$variables[[tv]] <- droplevels(dd6$variables[[tv]])
f6  <- told ~ cond * period + sex + age_grp + race_m
st6 <- c(log(mean(dd6$variables$told)), rep(0, ncol(model.matrix(f6, dd6$variables)) - 1))
m6  <- tryCatch(do.call(svyglm, list(formula = f6, design = dd6, family = quasipoisson(link = "log"), start = st6)),
                error = function(e) NULL)
p_cp <- if (is.null(m6)) NA_real_ else as.numeric(regTermTest(m6, ~cond:period)$p)
c6$condition_x_period_P <- p_cp
write_out(c6, "check_C6_survey_period.csv")
# paired comparison within each era
c6b <- list()
for (pp in c("Early", "Late")) for (sx in SEXES) {
  r <- fit_model(sub_des(P2, in_sample(vp, SAMPLES[1]) & vp$in_pair & sex_cond(vp, sx) & vp$period == pp), told ~ type,
                 "Paired", SAMPLES[1], PRIMARY_DEF)
  r <- r[!is.na(r$term) & r$term == "typeOsteoporosis", , drop = FALSE]
  if (nrow(r)) { r$era <- pp; r$sex <- sx; c6b[[length(c6b) + 1]] <- r }
}
write_out(do.call(rbind, c6b), "check_C6b_paired_by_era.csv")

# -----------------------------------------------------------------------------
# C7. Denominators: total, known diagnosis status, events, weighted %
# -----------------------------------------------------------------------------
vo <- OP$variables
GROUPS <- list(
  "All with osteoporosis" = rep(TRUE, nrow(vo)),
  "Health care contact" = vo$engaged %in% "Engaged in care",
  "No health care contact" = vo$engaged %in% "Not engaged",
  "Contact status missing" = is.na(vo$engaged),
  "Screening age (women >=65, men >=70)" = (vo$sex == "Women" & vo$age >= 65) | (vo$sex == "Men" & vo$age >= 70),
  "Screening age, health care contact" = ((vo$sex == "Women" & vo$age >= 65) | (vo$sex == "Men" & vo$age >= 70)) &
                                          vo$engaged %in% "Engaged in care")
c7 <- list()
for (g in names(GROUPS)) for (sx in SEXES) {
  s <- GROUPS[[g]] & sex_cond(vo, sx)
  e <- est_cell(sub_des(OP, s), "told_op")
  c7[[length(c7) + 1]] <- data.frame(group = g, sex = sx, total = sum(s), known_status = sum(s & !is.na(vo$told_op)),
    missing_status = sum(s & is.na(vo$told_op)), reported_diagnosis = sum(s & vo$told_op %in% 1),
    pct = e$pct, lcl = e$lcl, ucl = e$ucl, nchs_flag = e$nchs_flag)
}
c7 <- do.call(rbind, c7)
write_out(c7, "check_C7_denominators.csv")

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
say("")
say("=================== SUMMARY ===================")
say("C1 - prescription-inclusive (all 4 conditions), health care contact: % reporting a diagnosis")
for (sx in c("Women", "Men")) for (cn in COND$condition) {
  r <- a1p[a1p$sample == SAMPLES[1] & a1p$definition == DEFS[2] & a1p$sex == sx & a1p$condition == cn, ][1, ]
  say(sprintf("  %-5s %-17s %s", sx, cn, fmt_pct(r)))
}
for (sx in SEXES) for (cn in COMPS) {
  r <- a1m[a1m$sample == SAMPLES[1] & a1m$definition == DEFS[2] & a1m$sex == sx & a1m$contrast == paste0("Osteoporosis vs ", cn), ]
  say(sprintf("  PR %-5s osteoporosis vs %-17s %s", sx, cn, fmt_pr(r)))
}
for (k in seq_len(nrow(a1x))) if (a1x$definition[k] == DEFS[2] && a1x$sample[k] == SAMPLES[1])
  say(sprintf("  %-70s %s", a1x$contrast[k], fmt_pr(a1x[k, ])))
r2 <- a2m[a2m$pairs == "All three conditions" & a2m$sample == SAMPLES[1] & a2m$definition == DEFS[2], ]
for (k in seq_len(nrow(r2))) say(sprintf("  Paired %-12s persons=%s PR %s", r2$sex[k], r2$persons[k], fmt_pr(r2[k, ])))
say("C2 - medication items answered among adults NOT reporting the diagnosis (0 = asked only after a diagnosis):")
for (k in seq_len(nrow(c2))) say(sprintf("  %-58s answered %5d of %5d told no | %5d of %5d told yes",
  c2$item[k], c2$answered_if_told_no[k], c2$told_no[k], c2$answered_if_told_yes[k], c2$told_yes[k]))
say("C3 - excluded because the diagnosis answer was missing (primary sample, both sexes):")
for (k in which(c3$sex == "Both")) say(sprintf("  %-17s measured %5d, answer missing %3d, analyzed %5d (yes %5d, no %5d)",
  c3$condition[k], c3$measured_positive[k], c3$positive_answer_missing_excluded[k], c3$analyzed[k], c3$analyzed_yes[k], c3$analyzed_no[k]))
say("C4a - A1 records vs participants:")
for (k in seq_len(nrow(c4a))) say(sprintf("  %-5s %5d records from %5d participants", c4a$sex[k], c4a$condition_records[k], c4a$unique_participants[k]))
say("C4b - paired outcome counts (pairs):")
for (k in seq_len(nrow(c4b))) say(sprintf("  %-17s %-5s persons %4d pairs %4d | both %4d neither %4d osteoporosis only %3d comparison only %4d",
  c4b$comparison[k], c4b$sex[k], c4b$persons[k], c4b$pairs[k], c4b$both_reported[k], c4b$neither_reported[k],
  c4b$osteoporosis_only[k], c4b$comparison_only[k]))
say("C5 - within-person PR, person-level ratio vs stacked records (should match):")
for (k in seq_len(nrow(c5))) say(sprintf("  %-5s person-level %.3f (%.2f-%.2f) | stacked %.3f (%.2f-%.2f)", c5$sex[k],
  c5$person_level_PR[k], c5$person_level_lcl[k], c5$person_level_ucl[k], c5$stacked_PR[k], c5$stacked_lcl[k], c5$stacked_ucl[k]))
say("C6 - survey period (condition x period P = ", signif(p_cp, 3), "):")
for (k in seq_len(nrow(c6))) say(sprintf("  %-50s %-40s %s", c6$analysis[k], c6$contrast[k], fmt_pr(c6[k, ])))
say("C7 - denominators:")
for (k in seq_len(nrow(c7))) say(sprintf("  %-38s %-5s total %4d known %4d missing %2d reported %3d  %s", c7$group[k], c7$sex[k],
  c7$total[k], c7$known_status[k], c7$missing_status[k], c7$reported_diagnosis[k],
  ifelse(is.na(c7$pct[k]), "NA", sprintf("%.1f%% (%.1f-%.1f) [%s]", c7$pct[k], c7$lcl[k], c7$ucl[k], c7$nchs_flag[k]))))
say("Done. Files in ", normalizePath(OUT_DIR))
