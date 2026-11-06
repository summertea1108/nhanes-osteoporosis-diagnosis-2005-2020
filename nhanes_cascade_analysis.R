#!/usr/bin/env Rscript
# =============================================================================
# Osteoporosis care cascade in US adults 50+, NHANES 2005-2006 to 2017-March 2020
# Main analysis: survey-weighted estimates, Korn-Graubard CIs, NCHS reliability
# flags, prevalence-ratio models, early-vs-late trends, sensitivity analyses.
#
# Run:   Rscript nhanes_cascade_analysis.R [data_dir] [out_dir]
#        defaults: nhanes_data  results
# Needs: haven, survey, dplyr  (ggplot2 optional, for the figure)
#
# Definitions
#   Analytic sample : age >= 50, examined (MEC weight > 0), valid femoral neck or total hip BMD
#   Osteoporosis    : lower of femoral-neck / total-hip T-score <= -2.5, both sexes scored
#                     against NHANES III white women aged 20-29 (ISCD convention)
#   Diagnosed       : OSQ060 = 1 ("doctor ever told you had osteoporosis")
#   Treated         : any osteoporosis drug in the 30-day prescription file (RXQ_RX)
#   Fragility fx    : hip or spine fracture at age >= 50 from a fall from standing height or less
#   Indicated       : osteoporosis OR fragility hip/spine fracture
#   Engaged in care : insured AND usual place of care AND >= 1 clinician visit in past 12 months
#
# Reliability: NCHS Data Presentation Standards for Proportions (Parker et al. 2017).
#   Suppress if n < 30 or effective n < 30; reliable if CI width <= 5 points; suppress if
#   width >= 30 points or relative CI width > 130%; flag if degrees of freedom < 8.
#   Each proportion AND its complement are judged separately (e.g. "treated" vs "untreated").
#
# Reference values: VERIFY against Looker et al., Osteoporos Int 1998;8:468-89 before submission.
# =============================================================================

user_lib <- Sys.getenv("R_LIBS_USER")
if (nzchar(user_lib) && dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  library(haven)
  library(survey)
  library(dplyr)
})
options(survey.lonely.psu = "adjust", dplyr.summarise.inform = FALSE, warn = 1)

args     <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(args) >= 1) args[1] else "nhanes_data"
OUT_DIR  <- if (length(args) >= 2) args[2] else "results"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
LOG <- file.path(OUT_DIR, "run_log.txt")
cat("", file = LOG)
say <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n")
  cat(msg, "\n", file = LOG, append = TRUE)
}
write_out <- function(df, name) {
  f <- file.path(OUT_DIR, name)
  write.csv(df, f, row.names = FALSE, na = "")
  say("  wrote ", f, " (", nrow(df), " rows)")
}

# -----------------------------------------------------------------------------
# 1. Configuration
# -----------------------------------------------------------------------------
CYCLES <- data.frame(
  cycle  = c("2005-2006", "2007-2008", "2009-2010", "2013-2014", "2017-2020"),
  pre    = c("", "", "", "", "P_"),
  suf    = c("_D", "_E", "_F", "_H", ""),
  wt     = c("WTMEC2YR", "WTMEC2YR", "WTMEC2YR", "WTMEC2YR", "WTMECPRP"),
  years  = c(2, 2, 2, 2, 3.2),
  period = c("Early", "Early", "Early", "Late", "Late"),
  stringsAsFactors = FALSE
)
TOTAL_YEARS <- sum(CYCLES$years)   # 11.2

REF_F <- list(FN = c(0.858, 0.120), TH = c(0.942, 0.122))   # NHANES III white women 20-29
REF_M <- list(FN = c(0.930, 0.136), TH = c(1.041, 0.144))   # NHANES III white men 20-29

OP_DRUGS <- c("ALENDRONATE", "RISEDRONATE", "IBANDRONATE", "ZOLEDRONIC", "ETIDRONATE",
              "RALOXIFENE", "BAZEDOXIFENE", "TERIPARATIDE", "ABALOPARATIDE",
              "DENOSUMAB", "ROMOSOZUMAB", "CALCITONIN")
OP_PATTERN        <- paste(OP_DRUGS, collapse = "|")
OP_PATTERN_NOCALC <- paste(setdiff(OP_DRUGS, "CALCITONIN"), collapse = "|")

# -----------------------------------------------------------------------------
# 2. Helpers for reading and recoding
# -----------------------------------------------------------------------------
read_nh <- function(name) {
  f <- file.path(DATA_DIR, paste0(name, ".xpt"))
  if (!file.exists(f)) stop("Missing data file: ", f, call. = FALSE)
  d <- as.data.frame(haven::zap_labels(haven::read_xpt(f)))
  names(d) <- toupper(names(d))
  d
}
getcol <- function(d, ...) {
  for (n in c(...)) if (n %in% names(d)) return(d[[n]])
  rep(NA_real_, nrow(d))
}
num <- function(x) suppressWarnings(as.numeric(x))
# yes/no items: 1 = yes, 2 = no, 7/9 = refused/don't know
yn <- function(x) { x <- num(x); ifelse(x == 1, 1, ifelse(x == 2, 0, NA_real_)) }
# multi-digit items: drop only 77/99 style codes (7 and 9 can be real answers)
drop_77_99 <- function(x) { x <- num(x); x[x %in% c(77, 99, 777, 999, 7777, 9999, 77777, 99999)] <- NA; x }

# Fragility fracture at a site ("A" = hip, "C" = spine): any fracture episode
# at age >= 50 (OSD030 age, or OSQ040 = "50 or older" if age unknown) caused by
# a fall from standing height or less (OSD050 = 1).
fragility_fx <- function(o, site) {
  out <- rep(FALSE, nrow(o))
  age_cols <- grep(paste0("^OSD030", site, "[A-Z]$"), names(o), value = TRUE)
  for (a in age_cols) {
    k     <- substring(a, nchar(a))
    age   <- num(o[[a]]); age[!is.na(age) & age >= 7777] <- NA
    u50   <- num(getcol(o, paste0("OSQ040", site, k)))
    cause <- num(getcol(o, paste0("OSD050", site, k)))
    ge50  <- (!is.na(age) & age >= 50) | (is.na(age) & !is.na(u50) & u50 == 2)
    out   <- out | (ge50 & !is.na(cause) & cause == 1)
  }
  as.numeric(out)
}

# -----------------------------------------------------------------------------
# 3. Build one person-level file per cycle
# -----------------------------------------------------------------------------
build_cycle <- function(i) {
  cy <- CYCLES[i, ]
  nm <- function(f) paste0(cy$pre, f, cy$suf)
  say("Reading ", cy$cycle)
  demo <- read_nh(nm("DEMO")); fem <- read_nh(nm("DXXFEM")); osq <- read_nh(nm("OSQ"))
  rx   <- read_nh(nm("RXQ_RX")); hiq <- read_nh(nm("HIQ")); huq <- read_nh(nm("HUQ"))

  d <- data.frame(
    SEQN     = demo$SEQN,
    cycle    = cy$cycle,
    period   = cy$period,
    sex      = ifelse(demo$RIAGENDR == 1, "Men", ifelse(demo$RIAGENDR == 2, "Women", NA)),
    age      = pmin(num(demo$RIDAGEYR), 80),          # 2005-06 top-codes at 85
    race_code = num(demo$RIDRETH1),
    pir      = num(demo$INDFMPIR),
    educ_code = num(getcol(demo, "DMDEDUC2")),
    psu      = num(demo$SDMVPSU),
    strata   = paste0(cy$cycle, "_", demo$SDMVSTRA),
    w_cycle  = num(demo[[cy$wt]]),
    stringsAsFactors = FALSE
  )
  d$w_pooled <- d$w_cycle * cy$years / TOTAL_YEARS

  d <- merge(d, fem[, c("SEQN", "DXXNKBMD", "DXXOFBMD")], by = "SEQN", all.x = TRUE)

  # body measures (BMI and weight predict a missing scan: the scanner has a weight limit)
  bmx <- read_nh(nm("BMX"))
  d <- merge(d, data.frame(SEQN = bmx$SEQN, bmi = num(getcol(bmx, "BMXBMI")),
                           wt_kg = num(getcol(bmx, "BMXWT"))), by = "SEQN", all.x = TRUE)

  o <- data.frame(
    SEQN        = osq$SEQN,
    told_op     = yn(getcol(osq, "OSQ060")),
    q_tx        = yn(getcol(osq, "OSQ072", "OSQ070")),   # wording differs: 2005-10 ever treated; 2013+ prescribed
    fx_hip_any  = yn(getcol(osq, "OSQ010A")),
    fx_spine_any = yn(getcol(osq, "OSQ010C")),
    frag_hip    = fragility_fx(osq, "A"),
    frag_spine  = fragility_fx(osq, "C"),
    steroid     = yn(getcol(osq, "OSQ130"))
  )
  d <- merge(d, o, by = "SEQN", all.x = TRUE)

  rx$drug   <- toupper(ifelse(is.na(rx$RXDDRUG), "", as.character(rx$RXDDRUG)))
  rx$use    <- yn(rx$RXDUSE)
  rx$op     <- grepl(OP_PATTERN, rx$drug)
  rx$op_nc  <- grepl(OP_PATTERN_NOCALC, rx$drug)
  p <- rx %>% group_by(SEQN) %>%
    summarise(use = first(use), op = any(op), op_nc = any(op_nc)) %>% as.data.frame()
  p$rx_tx        <- ifelse(is.na(p$use), NA, ifelse(p$use == 0, 0, as.numeric(p$op)))
  p$rx_tx_nocalc <- ifelse(is.na(p$use), NA, ifelse(p$use == 0, 0, as.numeric(p$op_nc)))
  d <- merge(d, p[, c("SEQN", "rx_tx", "rx_tx_nocalc")], by = "SEQN", all.x = TRUE)

  ins <- yn(getcol(hiq, "HIQ011"))
  # insurance-type items: "yes" is coded as the item's own number (14, 15, ... in 2005-14;
  # 1, 2, ... in 2017-20), so test for a non-missing answer rather than == 1
  covered <- function(...) {
    x <- drop_77_99(getcol(hiq, ...))
    ifelse(is.na(ins), NA, ifelse(ins == 0, 0, as.numeric(!is.na(x))))
  }
  h <- data.frame(
    SEQN         = hiq$SEQN,
    insured      = ins,
    ins_private  = covered("HIQ032A", "HIQ031A"),
    ins_medicare = covered("HIQ032B", "HIQ031B"),
    ins_medicaid = covered("HIQ032D", "HIQ031D")
  )
  d <- merge(d, h, by = "SEQN", all.x = TRUE)

  usual  <- num(getcol(huq, "HUQ030"))
  visits <- drop_77_99(getcol(huq, "HUQ051", "HUQ050"))   # 0 = none in both versions
  u <- data.frame(
    SEQN        = huq$SEQN,
    gen_health  = { g <- num(getcol(huq, "HUQ010")); g[g %in% c(7, 9)] <- NA; g },
    usual_place = ifelse(usual %in% c(1, 3), 1, ifelse(usual == 2, 0, NA)),
    any_visit   = ifelse(is.na(visits), NA, as.numeric(visits >= 1))
  )
  d <- merge(d, u, by = "SEQN", all.x = TRUE)
  d
}

all <- do.call(rbind, lapply(seq_len(nrow(CYCLES)), build_cycle))

# -----------------------------------------------------------------------------
# 4. Derived variables
# -----------------------------------------------------------------------------
tsc <- function(bmd, ref) (bmd - ref[1]) / ref[2]
min2 <- function(a, b) ifelse(is.na(a), b, ifelse(is.na(b), a, pmin(a, b)))

all <- all %>% mutate(
  t_fn_f = tsc(DXXNKBMD, REF_F$FN), t_th_f = tsc(DXXOFBMD, REF_F$TH),
  t_fn_m = tsc(DXXNKBMD, REF_M$FN), t_th_m = tsc(DXXOFBMD, REF_M$TH),
  t_min_f = min2(t_fn_f, t_th_f),
  t_min_m = min2(t_fn_m, t_th_m),
  has_dxa = !is.na(t_min_f),
  examined = !is.na(w_cycle) & w_cycle > 0,
  analytic = examined & !is.na(age) & age >= 50 & has_dxa,

  op         = ifelse(has_dxa, as.numeric(t_min_f <= -2.5), NA),
  op_maleref = ifelse(has_dxa, as.numeric(ifelse(sex == "Men", t_min_m, t_min_f) <= -2.5), NA),
  op_hiponly = ifelse(!is.na(t_th_f), as.numeric(t_th_f <= -2.5), NA),
  frag_any   = as.numeric(frag_hip == 1 | frag_spine == 1),
  frag_any   = ifelse(is.na(frag_any), 0, frag_any),
  indicated  = ifelse(has_dxa, as.numeric(op == 1 | frag_any == 1), NA),
  indicated_maleref = ifelse(has_dxa, as.numeric(op_maleref == 1 | frag_any == 1), NA),

  undiagnosed = 1 - told_op,
  untreated   = 1 - rx_tx,
  q_or_rx_tx  = ifelse(told_op == 1, as.numeric((!is.na(q_tx) & q_tx == 1) | (!is.na(rx_tx) & rx_tx == 1)), NA),

  sex = factor(sex, levels = c("Women", "Men")),
  period = factor(period, levels = c("Early", "Late")),
  age_grp = cut(age, c(50, 60, 65, 70, 80, Inf), right = FALSE,
                labels = c("50-59", "60-64", "65-69", "70-79", "80+")),
  age_grp3 = cut(age, c(50, 65, 80, Inf), right = FALSE, labels = c("50-64", "65-79", "80+")),
  race = factor(case_when(race_code %in% c(1, 2) ~ "Hispanic", race_code == 3 ~ "NH White",
                          race_code == 4 ~ "NH Black", race_code == 5 ~ "Other/multiracial"),
                levels = c("NH White", "NH Black", "Hispanic", "Other/multiracial")),
  pir3 = factor(case_when(is.na(pir) ~ "Missing", pir < 1.3 ~ "<1.3", pir < 3.5 ~ "1.3-<3.5",
                          TRUE ~ ">=3.5"), levels = c(">=3.5", "1.3-<3.5", "<1.3", "Missing")),
  pir2 = factor(case_when(is.na(pir) ~ NA_character_, pir < 2 ~ "<2.0", TRUE ~ ">=2.0"),
                levels = c(">=2.0", "<2.0")),
  educ = factor(case_when(educ_code %in% c(1, 2) ~ "<High school", educ_code == 3 ~ "High school/GED",
                          educ_code %in% c(4, 5) ~ "Some college+", TRUE ~ "Missing"),
                levels = c("Some college+", "High school/GED", "<High school", "Missing")),
  ins_type = factor(case_when(insured == 0 ~ "Uninsured", ins_medicaid == 1 ~ "Medicaid (any)",
                              ins_private == 1 ~ "Private (any)", insured == 1 ~ "Medicare/other public only",
                              TRUE ~ "Missing"),
                    levels = c("Private (any)", "Medicare/other public only", "Medicaid (any)",
                               "Uninsured", "Missing")),
  insured_f = factor(case_when(insured == 1 ~ "Insured", insured == 0 ~ "Uninsured"),
                     levels = c("Insured", "Uninsured")),
  usual_f = factor(case_when(usual_place == 1 ~ "Has usual place", usual_place == 0 ~ "No usual place"),
                   levels = c("Has usual place", "No usual place")),
  visit_f = factor(case_when(any_visit == 1 ~ "Visit in past 12 mo", any_visit == 0 ~ "No visit"),
                   levels = c("Visit in past 12 mo", "No visit")),
  engaged = factor(case_when(insured == 1 & usual_place == 1 & any_visit == 1 ~ "Engaged in care",
                             !is.na(insured) & !is.na(usual_place) & !is.na(any_visit) ~ "Not engaged"),
                   levels = c("Engaged in care", "Not engaged")),
  access_barrier = factor(case_when(insured == 0 | usual_place == 0 ~ "Uninsured or no usual place",
                                    insured == 1 & usual_place == 1 ~ "Insured with usual place"),
                          levels = c("Insured with usual place", "Uninsured or no usual place")),
  frag_f = factor(ifelse(frag_any == 1, "Fragility hip/spine fx", "No fragility fx"),
                  levels = c("No fragility fx", "Fragility hip/spine fx")),
  steroid_f = factor(case_when(steroid == 1 ~ "Ever daily steroid", steroid == 0 ~ "No",
                               TRUE ~ "Missing"), levels = c("No", "Ever daily steroid", "Missing")),
  all_f = factor("All"),
  health_f = factor(case_when(gen_health %in% c(1, 2) ~ "Excellent/very good", gen_health == 3 ~ "Good",
                              gen_health %in% c(4, 5) ~ "Fair/poor", TRUE ~ "Missing"),
                    levels = c("Excellent/very good", "Good", "Fair/poor", "Missing")),
  bmi_f = factor(case_when(is.na(bmi) ~ "Missing", bmi < 25 ~ "<25", bmi < 30 ~ "25-<30",
                           bmi < 35 ~ "30-<35", bmi < 40 ~ "35-<40", TRUE ~ ">=40"),
                 levels = c("<25", "25-<30", "30-<35", "35-<40", ">=40", "Missing")),
  heavy_f = factor(case_when(is.na(wt_kg) ~ "Missing", wt_kg >= 136 ~ ">=136 kg", TRUE ~ "<136 kg"),
                   levels = c("<136 kg", ">=136 kg", "Missing"))
)
# model versions of factors with explicit Missing level (so models keep everyone)
add_missing <- function(f) { f <- addNA(f); levels(f)[is.na(levels(f))] <- "Missing"; f }
all <- all %>% mutate(
  race_m = add_missing(race), usual_m = add_missing(usual_f), visit_m = add_missing(visit_f)
)

# -----------------------------------------------------------------------------
# 5. Survey design (full examined sample; analyses use subsets = domains)
# -----------------------------------------------------------------------------
exam <- all[all$examined, ]
D <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_pooled, nest = TRUE, data = exam)
A <- subset(D, analytic)                      # age 50+, valid hip DXA

va <- A$variables
say("")
say("=== Check against feasibility run (unweighted) ===")
say("Analytic sample: ", nrow(va))
say("Osteoporosis: ", sum(va$op == 1, na.rm = TRUE), "  women ", sum(va$op == 1 & va$sex == "Women", na.rm = TRUE),
    "  men ", sum(va$op == 1 & va$sex == "Men", na.rm = TRUE))
say("Diagnosed among osteoporosis: women ", sum(va$op == 1 & va$told_op == 1 & va$sex == "Women", na.rm = TRUE),
    "  men ", sum(va$op == 1 & va$told_op == 1 & va$sex == "Men", na.rm = TRUE))
say("Treated among osteoporosis:   women ", sum(va$op == 1 & va$rx_tx == 1 & va$sex == "Women", na.rm = TRUE),
    "  men ", sum(va$op == 1 & va$rx_tx == 1 & va$sex == "Men", na.rm = TRUE))
say("Indicated (osteoporosis or fragility fx): women ", sum(va$indicated == 1 & va$sex == "Women", na.rm = TRUE),
    "  men ", sum(va$indicated == 1 & va$sex == "Men", na.rm = TRUE))
say("Fragility hip/spine fx (analytic sample): ", sum(va$frag_any == 1, na.rm = TRUE))
say("Expected from feasibility run: 12,240 / 779 / 627 / 152; diagnosed 251 / 22; treated 98 / 3")
say("")

# -----------------------------------------------------------------------------
# 6. Estimation helpers
# -----------------------------------------------------------------------------
sub_des <- function(des, cond) {
  cond[is.na(cond)] <- FALSE
  if (sum(cond) == 0) return(NULL)
  des[cond, ]
}
domain_df <- function(v) length(unique(paste(v$strata, v$psu))) - length(unique(v$strata))

nchs_flag <- function(p, l, u, n, n_eff, df) {
  if (is.na(p)) return("No data")
  if (n < 30 || is.na(n_eff) || n_eff < 30) return("Suppress: n or effective n < 30")
  if (is.na(l) || is.na(u)) return("Review: CI not computable")
  w <- u - l
  out <- if (w <= 0.05) "Reliable"
         else if (w >= 0.30) "Suppress: CI width >= 30 points"
         else if (p <= 0 || w / p > 1.30) "Suppress: relative CI width > 130%"
         else "Reliable"
  if (df < 8 && out == "Reliable") out <- "Reliable, but df < 8 (flag)"
  out
}

empty_cell <- function(n = 0) data.frame(
  n = n, events = NA, pct = NA, lcl = NA, ucl = NA, se = NA, n_eff = NA, df = NA,
  pop_events_millions = NA, nchs_flag = "No data", nchs_flag_complement = "No data")

# One proportion in one domain: % with y == 1 among non-missing y
est_cell <- function(des, y) {
  if (is.null(des)) return(empty_cell())
  v  <- des$variables
  ok <- !is.na(v[[y]]) & weights(des) > 0
  n  <- sum(ok)
  if (n == 0) return(empty_cell())
  d2 <- des[ok, ]
  v2 <- v[ok, ]
  events <- sum(v2[[y]] == 1)
  df <- max(domain_df(v2), 1)
  f  <- as.formula(paste0("~", y))
  m  <- svymean(f, d2)
  p  <- as.numeric(coef(m)); se <- as.numeric(SE(m))
  ci <- if (events == 0 || events == n) c(NA, NA) else
        tryCatch(as.numeric(attr(svyciprop(f, d2, method = "beta", df = df), "ci")),
                 error = function(e) c(NA, NA))
  n_eff <- if (se > 0) min(n, p * (1 - p) / se^2) else n
  tot <- as.numeric(coef(svytotal(f, d2))) / 1e6
  data.frame(
    n = n, events = events, pct = 100 * p, lcl = 100 * ci[1], ucl = 100 * ci[2], se = 100 * se,
    n_eff = round(n_eff, 1), df = df, pop_events_millions = round(tot, 3),
    nchs_flag = nchs_flag(p, ci[1], ci[2], n, n_eff, df),
    nchs_flag_complement = nchs_flag(1 - p, 1 - ci[2], 1 - ci[1], n, n_eff, df)
  )
}

# The cascade for one denominator-domain
cascade <- function(des, labels) {
  v <- if (is.null(des)) NULL else des$variables
  dx_des <- if (is.null(des)) NULL else sub_des(des, v$told_op == 1)
  rows <- list(
    cbind(labels, step = "1 Denominator size (n with outcome data not required)",
          empty_cell(if (is.null(des)) 0 else nrow(v))),
    cbind(labels, step = "2 Diagnosed (told by doctor)",      est_cell(des, "told_op")),
    cbind(labels, step = "2b Undiagnosed",                     est_cell(des, "undiagnosed")),
    cbind(labels, step = "3 Currently treated (Rx, 30 days)",  est_cell(des, "rx_tx")),
    cbind(labels, step = "3b Untreated",                       est_cell(des, "untreated")),
    cbind(labels, step = "4 Treated, among diagnosed",         est_cell(dx_des, "rx_tx"))
  )
  if (!is.null(des)) rows[[1]]$pop_events_millions <- round(sum(weights(des)) / 1e6, 3)
  do.call(rbind, rows)
}

SEXES <- c("Both", "Women", "Men")
sex_cond <- function(v, s) if (s == "Both") rep(TRUE, nrow(v)) else v$sex == s

# -----------------------------------------------------------------------------
# 7. Tables
# -----------------------------------------------------------------------------
say("=== Estimation ===")

# --- Table 0: who is missing a DXA scan (age 50+, examined) --------------------
E50 <- subset(D, !is.na(age) & age >= 50)
ve <- E50$variables
t0 <- list()
for (g in c(TRUE, FALSE)) {
  dg <- sub_des(E50, ve$has_dxa == g)
  for (var in c("sex", "age_grp3", "race", "pir3", "insured_f", "usual_f")) {
    vv <- dg$variables
    for (lv in levels(factor(vv[[var]]))) {
      dg$variables$ind_tmp <- as.numeric(as.character(vv[[var]]) == lv)
      t0[[length(t0) + 1]] <- cbind(data.frame(has_valid_dxa = g, variable = var, level = lv),
                                    est_cell(dg, "ind_tmp")[, c("n", "pct", "lcl", "ucl")])
    }
  }
}
write_out(do.call(rbind, t0), "table0_missing_dxa_comparison.csv")

# --- Table 1: characteristics of adults with osteoporosis ---------------------
OP <- subset(A, op == 1)
t1 <- list()
for (s in SEXES) {
  ds <- sub_des(OP, sex_cond(OP$variables, s))
  vv <- ds$variables
  for (var in c("age_grp", "race", "pir3", "educ", "ins_type", "usual_f", "visit_f", "engaged",
                "frag_f", "steroid_f", "period")) {
    for (lv in levels(factor(vv[[var]]))) {
      ds$variables$ind_tmp <- as.numeric(as.character(vv[[var]]) == lv)
      t1[[length(t1) + 1]] <- cbind(data.frame(sex = s, variable = var, level = lv),
                                    est_cell(ds, "ind_tmp")[, c("n", "events", "pct", "lcl", "ucl", "nchs_flag")])
    }
  }
}
write_out(do.call(rbind, t1), "table1_characteristics_osteoporosis.csv")

# --- Table 2: the cascade, by sex, for each denominator definition ------------
DENOMS <- c(op = "Osteoporosis (primary: female reference, neck/total hip)",
            op_maleref = "Osteoporosis (male reference for men)",
            op_hiponly = "Osteoporosis (total hip only)",
            indicated = "Treatment indicated (osteoporosis or fragility hip/spine fx)",
            indicated_maleref = "Treatment indicated (male reference for men)")
t2 <- list()
for (dn in names(DENOMS)) {
  dd <- sub_des(A, A$variables[[dn]] == 1)
  for (s in SEXES) {
    t2[[length(t2) + 1]] <- cascade(sub_des(dd, sex_cond(dd$variables, s)),
                                    data.frame(denominator = DENOMS[[dn]], sex = s))
  }
}
t2 <- do.call(rbind, t2)
write_out(t2, "table2_cascade_by_sex.csv")

# sensitivity: treatment definitions among diagnosed (primary osteoporosis domain)
t2b <- list()
for (s in SEXES) {
  dx <- sub_des(OP, OP$variables$told_op == 1 & sex_cond(OP$variables, s))
  for (pp in c("All", "Early", "Late")) {
    dp <- if (pp == "All") dx else sub_des(dx, dx$variables$period == pp)
    t2b[[length(t2b) + 1]] <- rbind(
      cbind(data.frame(sex = s, period = pp, definition = "Current Rx (30 days), primary"), est_cell(dp, "rx_tx")),
      cbind(data.frame(sex = s, period = pp, definition = "Current Rx excluding calcitonin"), est_cell(dp, "rx_tx_nocalc")),
      cbind(data.frame(sex = s, period = pp, definition = "Questionnaire: ever treated (2005-10) / prescribed (2013-20)"), est_cell(dp, "q_tx")),
      cbind(data.frame(sex = s, period = pp, definition = "Either current Rx or questionnaire"), est_cell(dp, "q_or_rx_tx")))
  }
}
write_out(do.call(rbind, t2b), "table2b_treatment_definitions_among_diagnosed.csv")

# --- Table 3: who is missed -- cascade by subgroup ----------------------------
GROUPS <- c("age_grp", "race", "pir3", "pir2", "educ", "ins_type", "insured_f", "usual_f",
            "visit_f", "engaged", "access_barrier", "frag_f", "period")
t3 <- list()
for (s in SEXES) {
  ds <- sub_des(OP, sex_cond(OP$variables, s))
  for (g in GROUPS) {
    vv <- ds$variables
    for (lv in levels(factor(vv[[g]]))) {
      t3[[length(t3) + 1]] <- cascade(sub_des(ds, as.character(vv[[g]]) == lv),
                                      data.frame(sex = s, subgroup = g, level = lv))
    }
  }
}
t3 <- do.call(rbind, t3)
t3 <- t3[!grepl("^1 ", t3$step), ]
write_out(t3, "table3_cascade_by_subgroup.csv")

# --- Table 4: adjusted prevalence ratios (modified Poisson, log link) ---------
fit_pr <- function(des, y, rhs, label, sample) {
  v  <- des$variables
  terms_in <- trimws(strsplit(rhs, "\\+")[[1]])
  ok <- complete.cases(v[, c(y, terms_in), drop = FALSE])
  d2 <- des[ok, ]
  # drop empty factor levels, and drop covariates left with a single level in this sample
  for (tv in terms_in) if (is.factor(d2$variables[[tv]])) d2$variables[[tv]] <- droplevels(d2$variables[[tv]])
  keep <- terms_in[vapply(terms_in, function(tv) length(unique(d2$variables[[tv]])) > 1, logical(1))]
  if (length(keep) == 0) keep <- "1"
  f  <- as.formula(paste(y, "~", paste(keep, collapse = " + ")))
  mm <- model.matrix(f, d2$variables)
  st <- c(log(max(mean(d2$variables[[y]]), 1e-4)), rep(0, ncol(mm) - 1))
  # do.call passes the starting values by value (svyglm evaluates glm() elsewhere,
  # so a plain `start = st` is not found)
  m  <- tryCatch(do.call(svyglm, list(formula = f, design = d2,
                                      family = quasipoisson(link = "log"), start = st)),
                 error = function(e) e)
  if (inherits(m, "error"))
    return(data.frame(sample = sample, outcome = y, model = label, term = NA, PR = NA, lcl = NA,
                      ucl = NA, p = NA, n = sum(ok), events = sum(v[[y]][ok] == 1),
                      note = conditionMessage(m)))
  ct <- summary(m)$coefficients
  ci <- suppressMessages(confint(m))[rownames(ct), , drop = FALSE]
  data.frame(sample = sample, outcome = y, model = label, term = rownames(ct),
             PR = exp(ct[, 1]), lcl = exp(ci[, 1]), ucl = exp(ci[, 2]), p = ct[, 4],
             n = sum(ok), events = sum(v[[y]][ok] == 1), note = "", row.names = NULL)
}
M <- list(
  M1 = "sex + age_grp + race_m",
  M2 = "sex + age_grp + race_m + frag_f + steroid_f",
  M3 = "sex + age_grp + race_m + frag_f + steroid_f + pir3 + educ",
  M4 = "sex + age_grp + race_m + frag_f + steroid_f + pir3 + educ + ins_type + usual_m + visit_m"
)
t4 <- list()
for (y in c("told_op", "rx_tx")) {
  for (mn in names(M)) {
    t4[[length(t4) + 1]] <- fit_pr(OP, y, M[[mn]], mn, "Osteoporosis, both sexes")
    t4[[length(t4) + 1]] <- fit_pr(sub_des(OP, OP$variables$sex == "Women"), y,
                                   sub("sex \\+ ", "", M[[mn]]), mn, "Osteoporosis, women")
  }
}
# treated among diagnosed (both sexes): is there a second barrier after diagnosis?
DXD <- sub_des(OP, OP$variables$told_op == 1)
for (mn in c("M1", "M3", "M4"))
  t4[[length(t4) + 1]] <- fit_pr(DXD, "rx_tx", M[[mn]], mn, "Osteoporosis and diagnosed, both sexes")
t4 <- do.call(rbind, t4)
t4 <- t4[is.na(t4$term) | t4$term != "(Intercept)", ]
write_out(t4, "table4_prevalence_ratios.csv")

# --- Table 5: missed opportunities -- people already engaged in care ----------
ENG <- sub_des(OP, OP$variables$engaged == "Engaged in care")
t5 <- list()
for (s in SEXES)
  t5[[length(t5) + 1]] <- cascade(sub_des(ENG, sex_cond(ENG$variables, s)),
                                  data.frame(sample = "Osteoporosis, engaged in care", sex = s))
IND_ENG <- sub_des(A, A$variables$indicated == 1 & A$variables$engaged == "Engaged in care")
for (s in SEXES)
  t5[[length(t5) + 1]] <- cascade(sub_des(IND_ENG, sex_cond(IND_ENG$variables, s)),
                                  data.frame(sample = "Treatment indicated, engaged in care", sex = s))
write_out(do.call(rbind, t5), "table5_cascade_engaged_in_care.csv")
t5m <- rbind(
  fit_pr(ENG, "told_op", "sex + age_grp + race_m", "Sex gap, age/race adjusted", "Osteoporosis, engaged"),
  fit_pr(ENG, "told_op", "sex + age_grp + race_m + frag_f + steroid_f + pir3 + educ", "Sex gap, fully adjusted", "Osteoporosis, engaged"),
  fit_pr(ENG, "rx_tx", "sex + age_grp + race_m", "Sex gap, age/race adjusted", "Osteoporosis, engaged"),
  fit_pr(ENG, "rx_tx", "sex + age_grp + race_m + frag_f + steroid_f + pir3 + educ", "Sex gap, fully adjusted", "Osteoporosis, engaged"))
write_out(t5m[t5m$term != "(Intercept)" | is.na(t5m$term), ], "table5b_sex_gap_engaged_models.csv")

# --- Table 6: trends, early (2005-10) vs late (2013-20) ------------------------
t6 <- list(); t6d <- list(); t6m <- list()
for (dn in c("op", "op_hiponly", "indicated")) {
  dd <- sub_des(A, A$variables[[dn]] == 1)
  for (s in SEXES) {
    ds <- sub_des(dd, sex_cond(dd$variables, s))
    for (pp in c("Early", "Late"))
      t6[[length(t6) + 1]] <- cascade(sub_des(ds, ds$variables$period == pp),
                                      data.frame(denominator = DENOMS[[dn]], sex = s, period = pp))
    for (y in c("told_op", "rx_tx")) {
      vv <- ds$variables
      dy <- sub_des(ds, !is.na(vv[[y]]))
      b  <- tryCatch(svyby(as.formula(paste0("~", y)), ~period, dy, svymean, covmat = TRUE), error = function(e) NULL)
      if (!is.null(b) && nrow(b) == 2) {
        cc <- svycontrast(b, c(-1, 1))
        est <- as.numeric(coef(cc)); se <- as.numeric(SE(cc))
        dfd <- max(domain_df(dy$variables), 1)
        t6d[[length(t6d) + 1]] <- data.frame(denominator = DENOMS[[dn]], sex = s, outcome = y,
          early_pct = 100 * coef(b)[1], late_pct = 100 * coef(b)[2],
          diff_pts = 100 * est, lcl = 100 * (est - qt(0.975, dfd) * se),
          ucl = 100 * (est + qt(0.975, dfd) * se), p = 2 * pt(-abs(est / se), dfd), row.names = NULL)
      }
      rhs <- if (s == "Both") "period + sex + age_grp + race_m" else "period + age_grp + race_m"
      t6m[[length(t6m) + 1]] <- fit_pr(ds, y, rhs, "Late vs early, adjusted",
                                       paste0(DENOMS[[dn]], " | ", s))
    }
  }
}
write_out(do.call(rbind, t6), "table6_cascade_by_period.csv")
write_out(do.call(rbind, t6d), "table6b_period_differences.csv")
t6m <- do.call(rbind, t6m)
write_out(t6m[!is.na(t6m$term) & grepl("^period", t6m$term), ], "table6c_period_prevalence_ratios.csv")

# descriptive five-cycle series (women and both sexes; primary definition)
t6e <- list()
for (s in c("Both", "Women", "Men")) {
  ds <- sub_des(OP, sex_cond(OP$variables, s))
  for (cy in CYCLES$cycle)
    t6e[[length(t6e) + 1]] <- cascade(sub_des(ds, ds$variables$cycle == cy),
                                      data.frame(sex = s, cycle = cy))
}
write_out(do.call(rbind, t6e), "table6d_cascade_by_cycle_descriptive.csv")

# --- Table 7: sensitivity -- inverse-probability weighting for missing DXA ----
# Model the probability of a valid hip scan among examined adults 50+, then
# divide the survey weight of each scanned participant by that probability, so the
# scanned sample also stands in for similar people whose scan was missing or invalid.
# (Variance does not account for estimating the weights; this is a sensitivity analysis.)
say("")
say("=== Sensitivity: inverse-probability weighting for missing DXA ===")
E50$variables$has_dxa_n <- as.numeric(E50$variables$has_dxa)
dxa_rhs <- "sex + age_grp + race_m + pir3 + educ + ins_type + usual_m + visit_m + health_f + bmi_f + heavy_f"
m_dxa <- svyglm(as.formula(paste("has_dxa_n ~", dxa_rhs)), design = E50, family = quasibinomial())
ve2 <- E50$variables
p_hat <- as.numeric(predict(m_dxa, newdata = ve2, type = "response", se.fit = FALSE))
ipw_tab <- data.frame(SEQN = ve2$SEQN, p_scan = p_hat)
exam2 <- merge(exam, ipw_tab, by = "SEQN", all.x = TRUE)
exam2$w_ipw <- ifelse(exam2$analytic & !is.na(exam2$p_scan), exam2$w_pooled / exam2$p_scan, exam2$w_pooled)
say("  Predicted probability of a valid scan (age 50+): min ", round(min(p_hat), 3),
    ", median ", round(median(p_hat), 3), ", max ", round(max(p_hat), 3))
ct_dxa <- summary(m_dxa)$coefficients
write_out(data.frame(term = rownames(ct_dxa), OR = exp(ct_dxa[, 1]), p = ct_dxa[, 4], row.names = NULL),
          "table7a_missing_dxa_model_odds_ratios.csv")

D_ipw  <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_ipw, nest = TRUE, data = exam2)
A_ipw  <- subset(D_ipw, analytic)
OP_ipw <- subset(A_ipw, op == 1)
t7 <- list()
for (s in SEXES) {
  t7[[length(t7) + 1]] <- cascade(sub_des(OP_ipw, sex_cond(OP_ipw$variables, s)),
                                  data.frame(weighting = "IPW for missing DXA", sample = "Osteoporosis", sex = s))
  eng <- sub_des(OP_ipw, OP_ipw$variables$engaged == "Engaged in care" & sex_cond(OP_ipw$variables, s))
  t7[[length(t7) + 1]] <- cascade(eng, data.frame(weighting = "IPW for missing DXA",
                                                  sample = "Osteoporosis, engaged in care", sex = s))
}
write_out(do.call(rbind, t7), "table7b_cascade_ipw_missing_dxa.csv")
t7m <- list()
for (y in c("told_op", "rx_tx")) for (mn in c("M1", "M4"))
  t7m[[length(t7m) + 1]] <- fit_pr(OP_ipw, y, M[[mn]], mn, "Osteoporosis, both sexes, IPW for missing DXA")
t7m <- do.call(rbind, t7m)
write_out(t7m[is.na(t7m$term) | t7m$term != "(Intercept)", ], "table7c_prevalence_ratios_ipw_missing_dxa.csv")
t7t <- do.call(rbind, t7)
for (s in SEXES) for (st in c("2 Diagnosed (told by doctor)", "3b Untreated")) {
  r <- t7t[t7t$sample == "Osteoporosis" & t7t$sex == s & t7t$step == st, ][1, ]
  say("  IPW ", s, " ", sub("^[0-9a-z]+ ", "", st), ": ", sprintf("%.1f%% (%.1f-%.1f) [%s]", r$pct, r$lcl, r$ucl, r$nchs_flag))
}

# --- Table 8: revision analyses (reviewer comments) -----------------------------
# 8a  diagnosis x current treatment cross-tab (are the cascade stages nested?)
# 8b  drugs recorded among treated adults who did not report a diagnosis
# 8c  nested cascade: diagnosed AND currently treated, as % of all with osteoporosis
# 8d  screening-age-eligible adults (women >= 65, men >= 70)
# 8e  medication-inclusive denominator (T-score <= -2.5 OR currently on an osteoporosis drug)
# 8f  age summary for the abstract;  8g  race codes behind "Other/multiracial"
say("")
say("=== Revision analyses (Table 8) ===")
fmt8 <- function(r) if (is.na(r$pct)) paste0("NA  n=", r$n) else
  sprintf("%5.1f%% (%.1f-%.1f)  n=%d  [%s | complement: %s]", r$pct, r$lcl, r$ucl, r$n, r$nchs_flag, r$nchs_flag_complement)
mk_xtab <- function(v) {
  both <- !is.na(v$told_op) & !is.na(v$rx_tx)
  data.frame(
    dx1_tx1 = ifelse(both, as.numeric(v$told_op == 1 & v$rx_tx == 1), NA),
    dx1_tx0 = ifelse(both, as.numeric(v$told_op == 1 & v$rx_tx == 0), NA),
    dx0_tx1 = ifelse(both, as.numeric(v$told_op == 0 & v$rx_tx == 1), NA),
    dx0_tx0 = ifelse(both, as.numeric(v$told_op == 0 & v$rx_tx == 0), NA))
}
A$variables <- cbind(A$variables, mk_xtab(A$variables))
OP <- subset(A, op == 1)

# 8a cross-tab, among adults with osteoporosis (complete on both items)
t8a <- list()
for (s in SEXES) for (smp in c("Osteoporosis", "Osteoporosis, engaged in care")) {
  ds <- sub_des(OP, sex_cond(OP$variables, s) &
                  (smp == "Osteoporosis" | OP$variables$engaged %in% "Engaged in care"))
  for (cl in c("dx1_tx1", "dx1_tx0", "dx0_tx1", "dx0_tx0"))
    t8a[[length(t8a) + 1]] <- cbind(data.frame(sample = smp, sex = s, cell = cl), est_cell(ds, cl))
}
t8a <- do.call(rbind, t8a)
write_out(t8a, "table8a_diagnosis_by_treatment_crosstab.csv")
vo <- OP$variables
say("  Unweighted cross-tab (osteoporosis, both sexes): dx+tx+ ", sum(vo$dx1_tx1 == 1, na.rm = TRUE),
    "; dx+tx- ", sum(vo$dx1_tx0 == 1, na.rm = TRUE), "; dx-tx+ ", sum(vo$dx0_tx1 == 1, na.rm = TRUE),
    "; dx-tx- ", sum(vo$dx0_tx0 == 1, na.rm = TRUE),
    "; missing diagnosis ", sum(is.na(vo$told_op)), "; missing Rx ", sum(is.na(vo$rx_tx)))

# 8b which drugs the treated-but-not-diagnosed were taking (unweighted counts)
tnd <- vo[vo$told_op %in% 0 & vo$rx_tx %in% 1, c("SEQN", "cycle", "sex")]
t8b <- list()
for (i in seq_len(nrow(CYCLES))) {
  cy <- CYCLES[i, ]; ids <- tnd$SEQN[tnd$cycle == cy$cycle]
  if (!length(ids)) next
  rx <- read_nh(paste0(cy$pre, "RXQ_RX", cy$suf))
  rx$drug <- toupper(ifelse(is.na(rx$RXDDRUG), "", as.character(rx$RXDDRUG)))
  rx <- rx[rx$SEQN %in% ids & grepl(OP_PATTERN, rx$drug), ]
  for (dg in OP_DRUGS) {
    k <- length(unique(rx$SEQN[grepl(dg, rx$drug)]))
    if (k > 0) t8b[[length(t8b) + 1]] <- data.frame(cycle = cy$cycle, drug = dg, participants = k)
  }
}
if (length(t8b)) {
  t8b <- do.call(rbind, t8b)
  tot <- aggregate(participants ~ drug, data = t8b, FUN = sum); tot$cycle <- "All"
  t8b <- rbind(t8b, tot[, c("cycle", "drug", "participants")])
} else t8b <- data.frame(cycle = "All", drug = "none", participants = 0)
write_out(t8b, "table8b_drugs_treated_not_diagnosed.csv")
write_out(as.data.frame(table(sex = tnd$sex, cycle = tnd$cycle)), "table8b2_treated_not_diagnosed_by_sex_cycle.csv")

# 8c-8e cascades with the nested step added
cascade_nested <- function(des, labels) {
  out <- cascade(des, labels)
  rbind(out, cbind(labels, step = "5 Diagnosed AND currently treated (nested)",
                   if (is.null(des)) empty_cell() else est_cell(des, "dx1_tx1")))
}
t8c <- list()
add <- function(des, ...) t8c[[length(t8c) + 1]] <<- cascade_nested(des, data.frame(...))
for (s in SEXES) {
  ds <- sub_des(OP, sex_cond(OP$variables, s))
  add(ds, sample = "Osteoporosis", sex = s)
  add(sub_des(ds, ds$variables$engaged %in% "Engaged in care"), sample = "Osteoporosis, engaged in care", sex = s)
}
# 8d screening-age eligible: women >= 65, men >= 70
elig <- with(OP$variables, (sex == "Women" & age >= 65) | (sex == "Men" & age >= 70))
SEL <- sub_des(OP, elig)
for (s in SEXES) {
  ds <- sub_des(SEL, sex_cond(SEL$variables, s))
  add(ds, sample = "Screening-age eligible (women 65+, men 70+)", sex = s)
  add(sub_des(ds, ds$variables$engaged %in% "Engaged in care"),
      sample = "Screening-age eligible, engaged in care", sex = s)
}
# men 70+ vs men 50-69 (descriptive contrast)
for (ag in c("50-69", "70+")) {
  cond <- OP$variables$sex == "Men" & (if (ag == "70+") OP$variables$age >= 70 else OP$variables$age < 70)
  add(sub_des(OP, cond), sample = paste0("Men aged ", ag), sex = "Men")
}
# 8e medication-inclusive denominator
A$variables$op_or_rx <- ifelse(is.na(A$variables$op), NA,
                               as.numeric(A$variables$op == 1 | A$variables$rx_tx %in% 1))
OPR <- subset(A, op_or_rx == 1)
for (s in SEXES) add(sub_des(OPR, sex_cond(OPR$variables, s)),
                     sample = "Osteoporosis OR current osteoporosis drug", sex = s)
say("  Medication-inclusive denominator adds ",
    sum(A$variables$op_or_rx == 1, na.rm = TRUE) - sum(A$variables$op == 1, na.rm = TRUE),
    " treated participants with T-score > -2.5")
t8c <- do.call(rbind, t8c)
write_out(t8c, "table8c_nested_and_eligible_cascades.csv")

# sex gap among screening-age-eligible adults (age group and race adjusted)
t8d <- rbind(fit_pr(SEL, "told_op", "sex + age_grp + race_m", "Sex gap, screening-age eligible", "Eligible"),
             fit_pr(SEL, "rx_tx", "sex + age_grp + race_m", "Sex gap, screening-age eligible", "Eligible"))
write_out(t8d[is.na(t8d$term) | t8d$term != "(Intercept)", ], "table8d_sex_gap_screening_eligible.csv")

# 8f age (NHANES top-codes age at 80 in 2005-2006 and later; 80 = "80 or older")
t8f <- list()
for (s in SEXES) {
  ds <- sub_des(OP, sex_cond(OP$variables, s))
  qv <- tryCatch(as.numeric(coef(svyquantile(~age, ds, c(0.25, 0.5, 0.75), ci = FALSE))),
                 error = function(e) rep(NA_real_, 3))
  mn <- svymean(~age, ds)
  sdv <- sqrt(as.numeric(coef(svyvar(~age, ds))))
  ds$variables$age80 <- as.numeric(ds$variables$age >= 80)
  t8f[[length(t8f) + 1]] <- data.frame(sex = s, n = nrow(ds$variables), mean = as.numeric(coef(mn)),
    sd = sdv, p25 = qv[1], median = qv[2], p75 = qv[3],
    pct_age80plus = 100 * as.numeric(coef(svymean(~age80, ds))))
}
t8f <- do.call(rbind, t8f)
write_out(t8f, "table8f_age_summary.csv")
say(sprintf("  Age, both sexes: weighted median %.0f (IQR %.0f-%.0f); mean %.1f (SD %.1f)",
            t8f$median[1], t8f$p25[1], t8f$p75[1], t8f$mean[1], t8f$sd[1]))

# 8g race codes (RIDRETH1: 1 Mexican American, 2 Other Hispanic, 3 NH White, 4 NH Black,
#    5 Other race including multiracial -- this is the "Other/multiracial" group)
write_out(as.data.frame(table(RIDRETH1 = vo$race_code, cycle = vo$cycle)), "table8g_race_codes_osteoporosis.csv")

# headline lines for the log
t8s <- t8c[t8c$step %in% c("2 Diagnosed (told by doctor)", "3b Untreated",
                            "5 Diagnosed AND currently treated (nested)"), ]
for (k in seq_len(nrow(t8s))) say(sprintf("  %-46s %-6s %-40s %s", t8s$sample[k], t8s$sex[k],
                                          sub("^[0-9a-z]+ ", "", t8s$step[k]), fmt8(t8s[k, ])))


# --- Table 9: complete-case sensitivity (missing covariates) ---------------------
# Primary models keep everyone by giving missing covariates their own category.
# Here the fully adjusted (model 4) diagnosis and treatment models are refit among
# participants with no missing covariate, to check the sex and race estimates.
say("")
say("=== Complete-case sensitivity (Table 9) ===")
vcc <- OP$variables
cc <- !is.na(vcc$race) & vcc$pir3 != "Missing" & vcc$educ != "Missing" & vcc$ins_type != "Missing" &
      vcc$steroid_f != "Missing" & !is.na(vcc$usual_f) & !is.na(vcc$visit_f)
say("  Complete cases: ", sum(cc), " of ", nrow(vcc), " adults with osteoporosis")
OPCC <- sub_des(OP, cc)
t9 <- rbind(fit_pr(OPCC, "told_op", M[["M1"]], "M1", "Complete cases"),
            fit_pr(OPCC, "told_op", M[["M4"]], "M4", "Complete cases"),
            fit_pr(OPCC, "rx_tx",   M[["M4"]], "M4", "Complete cases"))
t9 <- t9[is.na(t9$term) | t9$term != "(Intercept)", ]
write_out(t9, "table9_complete_case_prevalence_ratios.csv")
for (k in which(t9$term %in% c("sexMen", "race_mNH Black")))
  say(sprintf("  %-8s %-3s %-15s PR %.2f (%.2f-%.2f)", t9$outcome[k], t9$model[k], t9$term[k], t9$PR[k], t9$lcl[k], t9$ucl[k]))


# --- Table 10: third-review checks ---------------------------------------------
# 10a participant flow by cycle; 10b engaged-in-care results within questionnaire era
# (HUQ050 in 2005-2010 counts ER and home visits; HUQ051 from 2013 excludes them);
# 10c missing-scan weights with survey cycle and cycle-specific scanner weight limits
# (300 lb ~136 kg in 2005-2010; 450 lb ~204 kg from 2013); 10d treatment-model counts.
say("")
say("=== Third-review checks (Table 10) ===")
flow <- all %>% group_by(cycle) %>% summarise(
  demo_records = n(),
  examined = sum(examined),
  examined_50plus = sum(examined & !is.na(age) & age >= 50),
  valid_hip_dxa = sum(examined & !is.na(age) & age >= 50 & has_dxa),
  missing_or_invalid_dxa = sum(examined & !is.na(age) & age >= 50 & !has_dxa),
  osteoporosis = sum(examined & !is.na(age) & age >= 50 & has_dxa & op %in% 1),
  osteoporosis_missing_diagnosis = sum(examined & !is.na(age) & age >= 50 & has_dxa & op %in% 1 & is.na(told_op))) %>%
  as.data.frame()
flow <- rbind(flow, cbind(cycle = "All", as.data.frame(t(colSums(flow[, -1])))))
write_out(flow, "table10a_participant_flow_by_cycle.csv")

# 10b engaged in care, within questionnaire era
t10b <- list(); t10bm <- list()
for (pp in c("Early", "Late")) {
  eng <- sub_des(OP, OP$variables$engaged %in% "Engaged in care" & OP$variables$period == pp)
  for (s in SEXES)
    t10b[[length(t10b) + 1]] <- cbind(data.frame(period = pp, sex = s),
      est_cell(sub_des(eng, sex_cond(eng$variables, s)), "told_op")[, c("n", "events", "pct", "lcl", "ucl", "nchs_flag")])
  t10bm[[length(t10bm) + 1]] <- fit_pr(eng, "told_op", "sex + age_grp + race_m",
                                       paste("Engaged in care,", pp), "Osteoporosis, engaged")
}
write_out(do.call(rbind, t10b), "table10b_engaged_diagnosis_by_era.csv")
t10bm <- do.call(rbind, t10bm)
t10bm <- t10bm[is.na(t10bm$term) | t10bm$term %in% c("sexMen", "race_mNH Black"), ]
write_out(t10bm, "table10b2_engaged_diagnosis_prs_by_era.csv")
for (k in seq_len(nrow(t10bm))) say(sprintf("  %-28s %-15s PR %.2f (%.2f-%.2f) n=%s events=%s",
  t10bm$model[k], t10bm$term[k], t10bm$PR[k], t10bm$lcl[k], t10bm$ucl[k], t10bm$n[k], t10bm$events[k]))

# 10c missing-scan weights: add cycle and the cycle-specific weight limit
E50$variables$cycle_f <- factor(E50$variables$cycle)
lim <- ifelse(E50$variables$period == "Early", 136, 204)
E50$variables$over_limit_f <- factor(case_when(is.na(E50$variables$wt_kg) ~ "Missing",
  E50$variables$wt_kg >= lim ~ "Over limit", TRUE ~ "Under limit"), levels = c("Under limit", "Over limit", "Missing"))
say("  Over cycle-specific scanner limit (examined 50+): ", sum(E50$variables$over_limit_f == "Over limit"),
    "; >=136 kg: ", sum(E50$variables$heavy_f == ">=136 kg"))
rhs2 <- "sex + age_grp + race_m + pir3 + educ + ins_type + usual_m + visit_m + health_f + bmi_f + over_limit_f + cycle_f"
m_dxa2 <- svyglm(as.formula(paste("has_dxa_n ~", rhs2)), design = E50, family = quasibinomial())
p2 <- as.numeric(predict(m_dxa2, newdata = E50$variables, type = "response", se.fit = FALSE))
exam3 <- merge(exam, data.frame(SEQN = E50$variables$SEQN, p_scan2 = p2), by = "SEQN", all.x = TRUE)
exam3$w_ipw2 <- ifelse(exam3$analytic & !is.na(exam3$p_scan2), exam3$w_pooled / exam3$p_scan2, exam3$w_pooled)
an <- exam3$analytic & !is.na(exam3$p_scan2)
say(sprintf("  Revised scan probability among analytic: min %.3f, 1st pct %.3f, median %.3f; max weight ratio %.2f",
  min(p2[E50$variables$has_dxa]), quantile(p2[E50$variables$has_dxa], 0.01), median(p2[E50$variables$has_dxa]),
  max(1 / exam3$p_scan2[an])))
D3 <- svydesign(ids = ~psu, strata = ~strata, weights = ~w_ipw2, nest = TRUE, data = exam3)
OP3 <- subset(subset(D3, analytic), op == 1)
t10c <- list()
for (s in SEXES) {
  ds <- sub_des(OP3, sex_cond(OP3$variables, s))
  t10c[[length(t10c) + 1]] <- cbind(data.frame(weighting = "Revised IPW (cycle, cycle-specific limit)", sex = s,
    pop_millions = round(sum(weights(ds)) / 1e6, 2)), est_cell(ds, "told_op")[, c("n", "events", "pct", "lcl", "ucl")])
}
t10c <- do.call(rbind, t10c)
t10c$undiag_millions <- NA; t10c$untreated_millions <- NA
for (k in seq_len(nrow(t10c))) {
  ds <- sub_des(OP3, sex_cond(OP3$variables, t10c$sex[k]))
  t10c$undiag_millions[k] <- round(sum(weights(ds)[ds$variables$told_op %in% 0]) / 1e6, 2)
  t10c$untreated_millions[k] <- round(sum(weights(ds)[ds$variables$rx_tx %in% 0]) / 1e6, 2)
}
write_out(t10c, "table10c_revised_missing_scan_weights.csv")
m3 <- fit_pr(OP3, "told_op", M[["M4"]], "M4", "Osteoporosis, revised IPW")
m3 <- m3[m3$term %in% c("sexMen", "race_mNH Black"), ]
write_out(m3, "table10c2_revised_ipw_prs.csv")
for (k in seq_len(nrow(t10c))) say(sprintf("  Revised IPW %-5s total %.2f M, undiagnosed %.2f M, no Rx %.2f M; diagnosed %.1f%% (%.1f-%.1f)",
  t10c$sex[k], t10c$pop_millions[k], t10c$undiag_millions[k], t10c$untreated_millions[k], t10c$pct[k], t10c$lcl[k], t10c$ucl[k]))
for (k in seq_len(nrow(m3))) say(sprintf("  Revised IPW M4 %-15s PR %.2f (%.2f-%.2f)", m3$term[k], m3$PR[k], m3$lcl[k], m3$ucl[k]))

# 10d outcome counts behind the treatment models
vo2 <- OP$variables
t10d <- data.frame(
  group = c("All", "Women", "Men", "Diagnosis missing"),
  n = c(nrow(vo2), sum(vo2$sex == "Women"), sum(vo2$sex == "Men"), sum(is.na(vo2$told_op))),
  rx_recorded = c(sum(vo2$rx_tx %in% 1), sum(vo2$rx_tx %in% 1 & vo2$sex == "Women"),
                  sum(vo2$rx_tx %in% 1 & vo2$sex == "Men"), sum(vo2$rx_tx %in% 1 & is.na(vo2$told_op))),
  rx_recorded_men_diag_missing = c(NA, NA, NA, sum(vo2$rx_tx %in% 1 & is.na(vo2$told_op) & vo2$sex == "Men")))
write_out(t10d, "table10d_treatment_counts.csv")
say("  Treated men (any diagnosis status): ", t10d$rx_recorded[3], "; treated with missing diagnosis: ", t10d$rx_recorded[4],
    " (men: ", t10d$rx_recorded_men_diag_missing[4], ")")

# -----------------------------------------------------------------------------
# 8. Figure 1: the cascade by sex (primary definition)
# -----------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  fd <- t2[t2$denominator == DENOMS[["op"]] & t2$sex != "Both" &
             t2$step %in% c("2 Diagnosed (told by doctor)", "3 Currently treated (Rx, 30 days)"), ]
  fd <- rbind(data.frame(sex = c("Women", "Men"), step = "1 Osteoporosis on DXA", pct = 100, lcl = NA, ucl = NA),
              fd[, c("sex", "step", "pct", "lcl", "ucl")])
  fd$step <- factor(sub("^[0-9]+ ", "", fd$step),
                    levels = c("Osteoporosis on DXA", "Diagnosed (told by doctor)", "Currently treated (Rx, 30 days)"))
  fd$sex <- factor(fd$sex, levels = c("Women", "Men"))
  p <- ggplot(fd, aes(x = step, y = pct, fill = sex)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = lcl, ymax = ucl), position = position_dodge(0.8), width = 0.2, na.rm = TRUE) +
    geom_text(aes(label = sprintf("%.0f%%", pct), y = pmax(pct, ifelse(is.na(ucl), 0, ucl)) + 3),
              position = position_dodge(0.8), size = 3.5) +
    scale_fill_manual(values = c(Women = "#2F6690", Men = "#D9822B")) +
    scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25)) +
    labs(x = NULL, y = "% of adults with osteoporosis", fill = NULL,
         title = "Osteoporosis care cascade, US adults 50+, NHANES 2005-2020") +
    theme_minimal(base_size = 12) + theme(legend.position = "top", panel.grid.major.x = element_blank())
  ggsave(file.path(OUT_DIR, "figure1_cascade_by_sex.png"), p, width = 7.5, height = 5, dpi = 300)
  say("  wrote figure1_cascade_by_sex.png")
} else {
  say("  ggplot2 not installed: figure skipped (tables are complete)")
}

# -----------------------------------------------------------------------------
# 9. Short readable summary
# -----------------------------------------------------------------------------
fmt <- function(r) if (is.na(r$pct)) "NA" else
  sprintf("%5.1f%% (%.1f-%.1f)  n=%d  [%s]", r$pct, r$lcl, r$ucl, r$n, r$nchs_flag)
pick <- function(tab, ...) { k <- list(...); r <- tab; for (nm in names(k)) r <- r[r[[nm]] == k[[nm]], ]; r[1, ] }
say("")
say("=== Headline numbers (primary definition, pooled 2005-2020) ===")
for (s in SEXES) {
  say(s, ":")
  for (st in c("2 Diagnosed (told by doctor)", "3 Currently treated (Rx, 30 days)", "3b Untreated",
               "4 Treated, among diagnosed"))
    say("   ", sprintf("%-36s", sub("^[0-9a-z]+ ", "", st)), fmt(pick(t2, denominator = DENOMS[["op"]], sex = s, step = st)))
}
t5t <- do.call(rbind, t5)
say("Engaged in care (insured + usual place + visit):")
for (s in c("Women", "Men")) for (st in c("2 Diagnosed (told by doctor)", "3b Untreated"))
  say("   ", s, " ", sprintf("%-30s", sub("^[0-9a-z]+ ", "", st)),
      fmt(pick(t5t, sample = "Osteoporosis, engaged in care", sex = s, step = st)))
say("")
say("Done. All tables are in ", normalizePath(OUT_DIR))
say("Report only estimates whose nchs_flag (or nchs_flag_complement for the complement) is 'Reliable'.")
