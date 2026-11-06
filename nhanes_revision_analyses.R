#!/usr/bin/env Rscript
# =============================================================================
# Revision analyses requested in review (POST HOC: added after the analysis plan
# dated 16 October 2026; report them as such)
#
#   R1  Comparison-condition definitions: medication-inclusive by condition and sex;
#       stricter thresholds (HbA1c >= 7.0%, BP >= 160/100 mm Hg); cholesterol with
#       never-tested participants excluded
#   R2  Reported diagnosis by number of health care visits in the past year (gradient)
#   R3  Within-person: one comparison condition per person (hypertension > diabetes >
#       cholesterol) vs all pairs
#   R4  Within-person: design-based analysis of discordant pairs
#   R5  Missing-scan weighting: final model, weight diagnostics, 99th-percentile truncation
#   R6  Numbers for manuscript placeholders, Table 1, main-analysis sensitivity summary
#
# Run:   Rscript nhanes_revision_analyses.R [data_dir] [out_dir] [main_script]
#        defaults: nhanes_data  results_revision  nhanes_cascade_analysis.R
# Setup (data build and definitions) is identical to nhanes_crosscondition_analysis.R.
# =============================================================================


args_f     <- commandArgs(trailingOnly = TRUE)
DATA_DIR_F <- if (length(args_f) >= 1) args_f[1] else "nhanes_data"
OUT_DIR_F  <- if (length(args_f) >= 2) args_f[2] else "results_revision"
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

# =============================================================================
# 3. Revision analyses (added after the 16 October 2026 plan, in response to
#    review; these are post hoc and are reported as such)
# =============================================================================
say("")
say("=== Revision analyses (post hoc, added after the 2026-10-16 plan) ===")

# ---- helpers -----------------------------------------------------------------
mk_design <- function(d, wcol)
  svydesign(ids = ~psu, strata = ~strata, weights = as.formula(paste0("~", wcol)), nest = TRUE, data = d)
fit_model <- function(des, f, label, extra = list()) {
  lab <- as.data.frame(c(list(model = label), extra), stringsAsFactors = FALSE)
  empty <- function(n, ev, note) cbind(lab, data.frame(term = NA, PR = NA, lcl = NA, ucl = NA, p = NA,
                                                       n = n, events = ev, note = note))
  if (is.null(des)) return(empty(0, 0, "no data"))
  v <- des$variables; vars <- all.vars(f); y <- vars[1]
  ok <- complete.cases(v[, vars, drop = FALSE]) & weights(des) > 0
  n <- sum(ok); ev <- sum(v[[y]][ok] == 1)
  if (n == 0) return(empty(0, 0, "no data"))
  d2 <- des[ok, ]
  for (tv in vars[-1]) if (is.factor(d2$variables[[tv]])) d2$variables[[tv]] <- droplevels(d2$variables[[tv]])
  mm <- model.matrix(f, d2$variables)
  st <- c(log(max(mean(d2$variables[[y]]), 1e-4)), rep(0, ncol(mm) - 1))
  m  <- tryCatch(do.call(svyglm, list(formula = f, design = d2, family = quasipoisson(link = "log"), start = st)),
                 error = function(e) e)
  if (inherits(m, "error")) return(empty(n, ev, conditionMessage(m)))
  ct <- summary(m)$coefficients
  ci <- suppressMessages(confint(m))[rownames(ct), , drop = FALSE]
  cbind(lab, data.frame(term = rownames(ct), PR = exp(ct[, 1]), lcl = exp(ci[, 1]), ucl = exp(ci[, 2]),
                        p = ct[, 4], n = n, events = ev, note = "", row.names = NULL))
}
invert <- function(r) { l <- 1 / r$ucl; u <- 1 / r$lcl; r$PR <- 1 / r$PR; r$lcl <- l; r$ucl <- u; r }
fmt_pr <- function(r) if (!nrow(r) || is.na(r$PR[1])) "not estimable" else
  sprintf("%.2f (%.2f-%.2f)", r$PR[1], r$lcl[1], r$ucl[1])
fmt_pct <- function(r) if (is.null(r) || !nrow(r) || is.na(r$pct[1])) "NA" else
  sprintf("%.1f%% (%.1f-%.1f) n=%d [%s]", r$pct[1], r$lcl[1], r$ucl[1], r$n[1], r$nchs_flag[1])
pct_row <- function(des, y, ...) cbind(data.frame(..., stringsAsFactors = FALSE),
                                       est_cell(des, y)[, c("n", "events", "pct", "lcl", "ucl", "n_eff", "df",
                                                            "nchs_flag", "nchs_flag_complement")])
cut_codes <- function(v, map) {
  out <- rep(NA_character_, length(v))
  for (nm in names(map)) out[v %in% map[[nm]]] <- nm
  factor(out, levels = names(map))
}

# ---- extra variables: visit counts (harmonized) and cholesterol testing ---------
# HUQ050 (2005-2010): 0 none, 1 one, 2 2-3, 3 4-9, 4 10-12, 5 13+
# HUQ051 (2013-2020): 0 none, 1 one, 2 2-3, 3 4-5, 4 6-7, 5 8-9, 6 10-12, 7 13-15, 8 16+
VISIT_LEV <- c("0", "1", "2-3", "4-9", ">=10")
read_extra2 <- function(i) {
  cy  <- CYC_F[i, ]
  huq <- read_nh(paste0(cy$pre, "HUQ", cy$suf))
  bpq <- read_nh(extra_name("BPQ", i))
  v   <- drop_77_99(getcol(huq, "HUQ051", "HUQ050"))
  vg  <- if (cy$cycle %in% c("2013-2014", "2017-2020"))
           cut_codes(v, list("0" = 0, "1" = 1, "2-3" = 2, "4-9" = 3:5, ">=10" = 6:8)) else
           cut_codes(v, list("0" = 0, "1" = 1, "2-3" = 2, "4-9" = 3, ">=10" = 4:5))
  h <- data.frame(SEQN = huq$SEQN, visit_grp = vg)
  t80 <- yn(getcol(bpq, "BPQ080")); chk <- yn(getcol(bpq, "BPQ060"))
  # never tested in any cycle: 2005-2014 skip BPQ080 when BPQ060 = 2; 2017-2020 ask BPQ080 first,
  # then BPQ060 of those not answering yes, so BPQ060 = 2 identifies never tested in every cycle
  b <- data.frame(SEQN = bpq$SEQN, chol_never_tested = chk %in% 0)
  out <- merge(h, b, by = "SEQN", all = TRUE); out$cycle <- cy$cycle; out
}
ex2 <- do.call(rbind, lapply(seq_len(nrow(CYC_F)), read_extra2))
exam <- merge(exam, ex2, by = c("SEQN", "cycle"), all.x = TRUE)
exam <- exam %>% mutate(
  visit_grp   = factor(as.character(visit_grp), levels = VISIT_LEV),
  visit_num   = as.numeric(visit_grp) - 1,                     # 0..4, for the trend test
  dm_strict   = ifelse(is.na(hba1c), NA, as.numeric(hba1c >= 7.0)),
  htn_strict  = ifelse(is.na(sbp), NA, as.numeric(sbp >= 160 | (!is.na(dbp) & dbp >= 100))),
  told_chol_tested = ifelse(chol_never_tested %in% TRUE, NA, told_chol),
  contact     = engaged %in% "Engaged in care"
)
D  <- mk_design(exam, "w_pooled")
A  <- subset(D, analytic)
OP <- subset(A, op == 1)
say("Visit groups among adults with osteoporosis: ",
    paste(names(table(OP$variables$visit_grp)), table(OP$variables$visit_grp), sep = "=", collapse = ", "))
say("Never tested for cholesterol among adults with measured high cholesterol and contact: ",
    sum(A$variables$chol_meas %in% 1 & A$variables$contact & A$variables$chol_never_tested %in% TRUE))

# =============================================================================
# R1. Condition definitions: medication-inclusive (by condition), stricter
#     thresholds, and cholesterol excluding never tested  (reviewer comment 1, 6)
# =============================================================================
say("")
say("--- R1: comparison-condition definitions ---")
CONDS <- c("Osteoporosis", "Diabetes", "Hypertension", "High cholesterol")
SPECS <- list(
  "Primary (measured)" = list(has = c("op", "dm_meas", "htn_meas", "chol_meas"),
                              told = c("told_op", "told_dm", "told_htn", "told_chol")),
  "Medication-inclusive" = list(has = c("op", "dm_medinc", "htn_medinc", "chol_medinc"),
                                told = c("told_op", "told_dm", "told_htn", "told_chol")),
  "Stricter thresholds (HbA1c >=7.0%; BP >=160/100 mm Hg)" =
    list(has = c("op", "dm_strict", "htn_strict", "chol_meas"),
         told = c("told_op", "told_dm", "told_htn", "told_chol")),
  "Cholesterol: never tested excluded" = list(has = c("op", "dm_meas", "htn_meas", "chol_meas"),
                                              told = c("told_op", "told_dm", "told_htn", "told_chol_tested"))
)
base_cols <- c("SEQN", "psu", "strata", "w_pooled", "sex", "age_grp", "race_m", "analytic", "contact")
r1p <- list(); r1m <- list()
for (sp in names(SPECS)) {
  s <- SPECS[[sp]]
  va <- A$variables
  for (k in seq_along(CONDS)) for (sx in SEXES)
    r1p[[length(r1p) + 1]] <- pct_row(sub_des(A, va$contact & sex_cond(va, sx) & va[[s$has[k]]] %in% 1), s$told[k],
                                      definition = sp, sex = sx, condition = CONDS[k])
  stk <- do.call(rbind, lapply(seq_along(CONDS), function(k) {
    d <- exam[, base_cols]; d$cond <- CONDS[k]; d$has <- exam[[s$has[k]]]; d$told <- exam[[s$told[k]]]; d }))
  stk$cond <- factor(stk$cond, levels = CONDS)
  S  <- mk_design(stk, "w_pooled"); vs <- S$variables
  dom <- vs$analytic %in% TRUE & vs$contact & vs$has %in% 1 & !is.na(vs$told)
  for (sx in SEXES) {
    f <- if (sx == "Both") told ~ cond + sex + age_grp + race_m else told ~ cond + age_grp + race_m
    r <- fit_model(sub_des(S, dom & sex_cond(vs, sx)), f, "Adjusted", list(definition = sp, sex = sx))
    r <- r[!is.na(r$term) & grepl("^cond", r$term), , drop = FALSE]
    if (nrow(r)) { r <- invert(r); r$contrast <- paste0("Osteoporosis vs ", sub("^cond", "", r$term)); r1m[[length(r1m) + 1]] <- r }
  }
  rx <- fit_model(sub_des(S, dom), told ~ cond * sex + age_grp + race_m, "Condition x sex", list(definition = sp, sex = "Men vs women"))
  rx <- rx[!is.na(rx$term) & grepl(":sexMen$", rx$term), , drop = FALSE]
  if (nrow(rx)) { rx <- invert(rx); rx$contrast <- paste0("Ratio of PRs (men/women), osteoporosis vs ",
                                                          sub(":sexMen$", "", sub("^cond", "", rx$term))); r1m[[length(r1m) + 1]] <- rx }
}
r1p <- do.call(rbind, r1p); r1m <- do.call(rbind, r1m)
write_out(r1p, "rev_R1_definitions_percentages.csv")
write_out(r1m, "rev_R1_definitions_prevalence_ratios.csv")

# =============================================================================
# R2. Visit-frequency gradient (reviewer comment 2)
# =============================================================================
say("")
say("--- R2: reported diagnosis by number of health care visits in the past year ---")
r2 <- list(); r2m <- list()
va <- A$variables
for (k in seq_along(CONDS)) {
  hv <- SPECS[[1]]$has[k]; tv <- SPECS[[1]]$told[k]
  for (sx in SEXES) for (vg in VISIT_LEV)
    r2[[length(r2) + 1]] <- pct_row(sub_des(A, va[[hv]] %in% 1 & sex_cond(va, sx) & va$visit_grp %in% vg), tv,
                                    condition = CONDS[k], sex = sx, visits = vg)
  for (sx in SEXES) {
    f <- as.formula(paste(tv, if (sx == "Both") "~ visit_num + sex + age_grp + race_m" else "~ visit_num + age_grp + race_m"))
    r <- fit_model(sub_des(A, va[[hv]] %in% 1 & sex_cond(va, sx)), f, "Trend per visit category",
                   list(condition = CONDS[k], sex = sx))
    r2m[[length(r2m) + 1]] <- r[!is.na(r$term) & r$term == "visit_num", , drop = FALSE]
  }
}
r2 <- do.call(rbind, r2); r2m <- do.call(rbind, r2m)
write_out(r2, "rev_R2_visit_gradient_percentages.csv")
write_out(r2m, "rev_R2_visit_gradient_trend.csv")

# =============================================================================
# R3/R4. Within-person: pair-level data, one condition per person, and
#        design-based analysis of discordant pairs (reviewer comment 3)
# =============================================================================
say("")
say("--- R3/R4: within-person sensitivity analyses ---")
COMP_ORDER <- c("Hypertension", "Diabetes", "High cholesterol")   # prespecified order for one pair per person
HAS  <- c(Hypertension = "htn_meas", Diabetes = "dm_meas", `High cholesterol` = "chol_meas")
TOLD <- c(Hypertension = "told_htn", Diabetes = "told_dm", `High cholesterol` = "told_chol")
elig_k <- lapply(COMP_ORDER, function(cn) exam$analytic %in% TRUE & exam$op %in% 1 & !is.na(exam$told_op) &
                   exam[[HAS[cn]]] %in% 1 & !is.na(exam[[TOLD[cn]]]))
names(elig_k) <- COMP_ORDER
kcount <- Reduce(`+`, elig_k)
first_k <- rep(NA_character_, nrow(exam))
for (cn in rev(COMP_ORDER)) first_k[elig_k[[cn]]] <- cn          # earliest in COMP_ORDER wins

# pair-level rows: one row per person-condition pair (all pairs) or one per person (first pair only)
pair_rows <- function(one_per_person) {
  do.call(rbind, lapply(COMP_ORDER, function(cn) {
    d <- exam[, c(base_cols, "told_op")]
    d$pair_cond <- cn
    d$told_comp <- exam[[TOLD[cn]]]
    d$in_pair   <- if (one_per_person) first_k %in% cn else elig_k[[cn]]
    d$w_pair    <- if (one_per_person) d$w_pooled else d$w_pooled / pmax(kcount, 1)
    d
  }))
}
long_rows <- function(pr) {   # 2 rows per pair for the marginal PR model
  a <- pr; a$type <- "Osteoporosis"; a$told <- pr$told_op
  b <- pr; b$type <- "Comparison";   b$told <- pr$told_comp
  out <- rbind(a, b); out$type <- factor(out$type, levels = c("Comparison", "Osteoporosis")); out
}
r3 <- list(); r4 <- list()
for (opp in c(FALSE, TRUE)) {
  lab <- if (opp) "One pair per person (hypertension > diabetes > cholesterol)" else "All pairs (primary)"
  for (smp in c("Health care contact", "All adults 50+ with hip DXA")) {
    pr <- pair_rows(opp)
    pr$dom <- pr$in_pair & (if (smp == "Health care contact") pr$contact else TRUE)
    # marginal PR (as in the primary analysis)
    L <- mk_design(long_rows(pr), "w_pair"); vl <- L$variables
    for (sx in SEXES) {
      r <- fit_model(sub_des(L, vl$dom & sex_cond(vl, sx)), told ~ type, "Paired marginal PR",
                     list(pairs = lab, sample = smp, sex = sx))
      r <- r[!is.na(r$term) & r$term == "typeOsteoporosis", , drop = FALSE]
      if (nrow(r)) {
        r$persons <- length(unique(vl$SEQN[vl$dom & sex_cond(vl, sx)]))
        r3[[length(r3) + 1]] <- r
      }
    }
    rx <- fit_model(sub_des(L, vl$dom), told ~ type * sex, "Report type x sex", list(pairs = lab, sample = smp, sex = "Men vs women"))
    rx <- rx[!is.na(rx$term) & rx$term == "typeOsteoporosis:sexMen", , drop = FALSE]
    if (nrow(rx)) { rx$persons <- NA; r3[[length(r3) + 1]] <- rx }
    # discordant pairs: among pairs whose two reports disagree, share in which osteoporosis is the
    # unreported one. Ratio = share / (1 - share) = weighted b/c (McNemar-type ratio).
    pr$disc <- !is.na(pr$told_op) & !is.na(pr$told_comp) & pr$told_op != pr$told_comp
    pr$osteo_unreported <- ifelse(pr$disc, as.numeric(pr$told_op == 0 & pr$told_comp == 1), NA)
    P <- mk_design(pr, "w_pair"); vp <- P$variables
    for (sx in SEXES) {
      e <- pct_row(sub_des(P, vp$dom & vp$disc & sex_cond(vp, sx)), "osteo_unreported",
                   pairs = lab, sample = smp, sex = sx)
      e$discordant_pairs <- e$n
      e$ratio <- e$pct / (100 - e$pct); e$ratio_lcl <- e$lcl / (100 - e$lcl); e$ratio_ucl <- e$ucl / (100 - e$ucl)
      e$concordant_pairs <- sum(vp$dom & !vp$disc & sex_cond(vp, sx) & !is.na(vp$told_op) & !is.na(vp$told_comp))
      r4[[length(r4) + 1]] <- e
    }
  }
}
r3 <- do.call(rbind, r3); r4 <- do.call(rbind, r4)
write_out(r3, "rev_R3_within_person_prevalence_ratios.csv")
write_out(r4, "rev_R4_discordant_pairs.csv")

# =============================================================================
# R5. Missing-scan weighting: final model, weight diagnostics, truncation (comment 5)
# =============================================================================
say("")
say("--- R5: missing-scan weighting diagnostics ---")
E50 <- subset(D, !is.na(age) & age >= 50)
ve  <- E50$variables
lim <- ifelse(ve$period == "Early", 136, 204)          # scanner weight limit, kg (300 lb; 450 lb from 2013)
E50$variables$has_dxa_n    <- as.numeric(ve$has_dxa)
E50$variables$cycle_f      <- factor(ve$cycle)
E50$variables$over_limit_f <- factor(case_when(is.na(ve$wt_kg) ~ "Missing", ve$wt_kg >= lim ~ "Over limit",
                                               TRUE ~ "Under limit"), levels = c("Under limit", "Over limit", "Missing"))
RHS_IPW <- "sex + age_grp + race_m + pir3 + educ + ins_type + usual_m + visit_m + health_f + bmi_f + over_limit_f + cycle_f"
m_scan <- svyglm(as.formula(paste("has_dxa_n ~", RHS_IPW)), design = E50, family = quasibinomial())
p_scan <- as.numeric(predict(m_scan, newdata = E50$variables, type = "response"))
exam$p_scan <- NA_real_
exam$p_scan[match(E50$variables$SEQN, exam$SEQN)] <- p_scan
exam$w_ipw  <- ifelse(exam$analytic %in% TRUE & !is.na(exam$p_scan), exam$w_pooled / exam$p_scan, exam$w_pooled)
opw  <- exam$w_ipw[exam$analytic %in% TRUE & exam$op %in% 1]
cut99 <- quantile(opw, 0.99)
exam$w_ipw_trunc <- ifelse(exam$analytic %in% TRUE & exam$op %in% 1, pmin(exam$w_ipw, cut99), exam$w_ipw)
kish <- function(w) sum(w)^2 / sum(w^2)
srt  <- sort(opw, decreasing = TRUE)
r5d <- data.frame(
  scan_probability_min = min(p_scan[E50$variables$has_dxa]), scan_probability_p1 = quantile(p_scan[E50$variables$has_dxa], 0.01),
  scan_probability_median = median(p_scan[E50$variables$has_dxa]),
  op_n = length(opw), op_weight_min = min(opw), op_weight_median = median(opw), op_weight_max = max(opw),
  op_weight_ratio_max_to_median = max(opw) / median(opw),
  op_share_of_weight_top1pct = sum(srt[seq_len(max(1, ceiling(0.01 * length(srt))))]) / sum(srt),
  op_kish_ess_survey_weights = kish(exam$w_pooled[exam$analytic %in% TRUE & exam$op %in% 1]),
  op_kish_ess_ipw = kish(opw), truncation_cut_99th = cut99, row.names = NULL)
write_out(r5d, "rev_R5_ipw_weight_diagnostics.csv")
ct <- summary(m_scan)$coefficients
write_out(data.frame(term = rownames(ct), OR = exp(ct[, 1]), p = ct[, 4], row.names = NULL), "rev_R5_ipw_scan_model.csv")
r5 <- list()
for (wc in c("w_pooled", "w_ipw", "w_ipw_trunc")) {
  Dw <- mk_design(exam, wc); OPw <- subset(subset(Dw, analytic), op == 1); vw <- OPw$variables
  for (smp in c("All with osteoporosis", "Health care contact")) {
    ind <- if (smp == "All with osteoporosis") rep(TRUE, nrow(vw)) else vw$contact
    for (sx in SEXES)
      r5[[length(r5) + 1]] <- pct_row(sub_des(OPw, ind & sex_cond(vw, sx)), "told_op", weights = wc, sample = smp, sex = sx)
    r <- fit_model(sub_des(OPw, ind), told_op ~ sex + age_grp + race_m, "Sex PR", list(weights = wc, sample = smp))
    r <- r[!is.na(r$term) & r$term == "sexMen", , drop = FALSE]
    if (nrow(r)) r5[[length(r5) + 1]] <- data.frame(weights = wc, sample = smp, sex = "PR men vs women", n = r$n, events = r$events,
                                                    pct = r$PR, lcl = r$lcl, ucl = r$ucl, n_eff = NA, df = NA,
                                                    nchs_flag = "PR (not a %)", nchs_flag_complement = "")
  }
}
r5 <- do.call(rbind, r5)
write_out(r5, "rev_R5_ipw_estimates.csv")

# =============================================================================
# R6. Numbers for the manuscript placeholders, Table 1, and the main-analysis
#     sensitivity summary
# =============================================================================
say("")
say("--- R6: manuscript placeholders ---")
va <- A$variables; vo <- OP$variables
m1 <- fit_pr(OP, "told_op", "sex + age_grp + race_m", "M1", "Osteoporosis, all")
m1 <- m1[!is.na(m1$term) & m1$term != "(Intercept)", ]
write_out(m1, "rev_R6_model1_all_osteoporosis.csv")
qv <- as.numeric(coef(svyquantile(~age, OP, c(0.25, 0.5, 0.75), ci = FALSE)))
E50f <- exam$examined %in% TRUE & !is.na(exam$age) & exam$age >= 50
flow <- data.frame(examined_50plus = sum(E50f), missing_or_invalid_dxa = sum(E50f & !exam$has_dxa),
                   pct_missing = 100 * mean(!exam$has_dxa[E50f]), analytic = sum(exam$analytic %in% TRUE),
                   osteoporosis = sum(exam$analytic %in% TRUE & exam$op %in% 1),
                   median_age = qv[2], age_p25 = qv[1], age_p75 = qv[3])
write_out(flow, "rev_R6_flow_and_age.csv")
med <- do.call(rbind, lapply(SEXES, function(sx) pct_row(sub_des(OP, sex_cond(vo, sx)), "rx_tx", sex = sx)))
write_out(med, "rev_R6_medication_recorded.csv")

# Table 1: characteristics by sex and by health care contact (weighted % among adults with osteoporosis)
T1VARS <- c("age_grp", "race", "pir3", "educ", "ins_type", "usual_f", "visit_f", "frag_f", "steroid_f", "period")
t1 <- list()
groups <- list(All = rep(TRUE, nrow(vo)), Women = vo$sex == "Women", Men = vo$sex == "Men",
               `Health care contact` = vo$contact, `No health care contact` = vo$engaged %in% "Not engaged")
for (g in names(groups)) {
  dg <- sub_des(OP, groups[[g]])
  for (var in T1VARS) {
    vv <- dg$variables
    for (lv in levels(factor(vv[[var]]))) {
      dg$variables$ind_tmp <- as.numeric(as.character(vv[[var]]) %in% lv)
      t1[[length(t1) + 1]] <- pct_row(dg, "ind_tmp", group = g, variable = var, level = lv)
    }
  }
}
write_out(do.call(rbind, t1), "rev_R6_table1.csv")

# main-analysis sensitivity summary: reported diagnosis by sex and sex PR
sens <- list()
add_sens <- function(lab, cond_vec, rhs = "sex + age_grp + race_m", des = A) {
  vv <- des$variables
  for (sx in SEXES) sens[[length(sens) + 1]] <<- pct_row(sub_des(des, cond_vec & sex_cond(vv, sx)), "told_op",
                                                         analysis = lab, sex = sx)
  r <- fit_model(sub_des(des, cond_vec), as.formula(paste("told_op ~", rhs)), "Sex PR", list(analysis = lab))
  r <- r[!is.na(r$term) & r$term == "sexMen", , drop = FALSE]
  if (nrow(r)) sens[[length(sens) + 1]] <<- data.frame(analysis = lab, sex = "PR men vs women", n = r$n, events = r$events,
                                                       pct = r$PR, lcl = r$lcl, ucl = r$ucl, n_eff = NA, df = NA,
                                                       nchs_flag = "PR (not a %)", nchs_flag_complement = "")
}
add_sens("Primary: female reference, neck or total hip", va$op %in% 1)
add_sens("Male reference for men", va$op_maleref %in% 1)
add_sens("Total hip only", va$op_hiponly %in% 1)
add_sens("Primary, adjusted for survey period", va$op %in% 1, "sex + age_grp + race_m + period")
add_sens("Health care contact, adjusted for survey period", va$op %in% 1 & va$contact, "sex + age_grp + race_m + period")
cc <- va$op %in% 1 & !is.na(va$race) & va$pir3 != "Missing" & va$educ != "Missing" & va$ins_type != "Missing" &
      va$steroid_f != "Missing" & !is.na(va$usual_f) & !is.na(va$visit_f)
add_sens("Complete cases, fully adjusted", cc,
         "sex + age_grp + race_m + frag_f + steroid_f + pir3 + educ + ins_type + usual_m + visit_m")
for (pp in c("Early", "Late")) add_sens(paste("Health care contact,", pp, "questionnaire era"),
                                        va$op %in% 1 & va$contact & va$period == pp)
write_out(do.call(rbind, sens), "rev_R6_main_sensitivity_summary.csv")

# =============================================================================
# Summary
# =============================================================================
say("")
say("=================== SUMMARY (health care contact unless stated) ===================")
say("R1 - PR, osteoporosis vs each condition (both sexes), and ratio of PRs men/women:")
for (sp in names(SPECS)) for (cn in CONDS[-1]) {
  b <- r1m[r1m$definition == sp & r1m$sex == "Both" & r1m$contrast == paste0("Osteoporosis vs ", cn), ]
  x <- r1m[r1m$definition == sp & r1m$sex == "Men vs women" & r1m$contrast == paste0("Ratio of PRs (men/women), osteoporosis vs ", cn), ]
  say(sprintf("  %-55s %-17s PR %s | ratio %s", sp, cn, fmt_pr(b), fmt_pr(x)))
}
say("R2 - osteoporosis reported, by visits in past year (all adults with osteoporosis):")
for (sx in c("Women", "Men")) for (vg in VISIT_LEV)
  say(sprintf("  %-5s visits %-4s %s", sx, vg, fmt_pct(r2[r2$condition == "Osteoporosis" & r2$sex == sx & r2$visits == vg, ])))
for (k in seq_len(nrow(r2m))) say(sprintf("  Trend PR per visit category, %-17s %-5s %s", r2m$condition[k], r2m$sex[k], fmt_pr(r2m[k, ])))
say("R3 - within-person marginal PR:")
for (k in seq_len(nrow(r3))) say(sprintf("  %-58s %-27s %-12s persons=%s PR %s", r3$pairs[k], r3$sample[k], r3$sex[k],
                                         r3$persons[k], fmt_pr(r3[k, ])))
say("R4 - discordant pairs: % in which osteoporosis was the unreported diagnosis (ratio = b/c):")
for (k in seq_len(nrow(r4))) say(sprintf("  %-58s %-27s %-5s discordant=%3d  %5.1f%% (%.1f-%.1f) ratio %.1f (%.1f-%.1f) [%s]",
  r4$pairs[k], r4$sample[k], r4$sex[k], r4$discordant_pairs[k], r4$pct[k], r4$lcl[k], r4$ucl[k],
  r4$ratio[k], r4$ratio_lcl[k], r4$ratio_ucl[k], r4$nchs_flag[k]))
say("R5 - weights among adults with osteoporosis: max/median ", round(r5d$op_weight_ratio_max_to_median, 1),
    "; top 1% hold ", round(100 * r5d$op_share_of_weight_top1pct, 1), "% of weight; Kish ESS ",
    round(r5d$op_kish_ess_survey_weights), " (survey) vs ", round(r5d$op_kish_ess_ipw), " (IPW)")
for (k in which(r5$sex %in% c("Women", "Men", "PR men vs women")))
  say(sprintf("  %-12s %-22s %-16s %.2f (%.2f-%.2f)", r5$weights[k], r5$sample[k], r5$sex[k], r5$pct[k], r5$lcl[k], r5$ucl[k]))
say("R6 - placeholders: examined 50+ ", flow$examined_50plus, "; missing/invalid DXA ", flow$missing_or_invalid_dxa,
    sprintf(" (%.1f%%)", flow$pct_missing), "; weighted median age ", round(flow$median_age), " (IQR ",
    round(flow$age_p25), "-", round(flow$age_p75), ")")
for (k in seq_len(nrow(m1))) say(sprintf("  M1 %-26s PR %.2f (%.2f-%.2f)", m1$term[k], m1$PR[k], m1$lcl[k], m1$ucl[k]))
for (k in seq_len(nrow(med))) say(sprintf("  Medication recorded, %-5s %s", med$sex[k], fmt_pct(med[k, ])))
say("Done. Files in ", normalizePath(OUT_DIR))
