#!/usr/bin/env Rscript
# =============================================================================
# Cross-condition analysis: reported diagnosis of osteoporosis compared with
# diabetes, hypertension and high cholesterol, US adults 50+, NHANES 2005-2020
#
# Implements the analysis plan dated 16 October 2026
# ("Analysis Plan: Osteoporosis Recognition Compared With Other Chronic Conditions").
#
#   A1  Population-level comparison (primary sample: health care contact)
#       - % reporting a prior diagnosis for each measured condition, by sex
#       - stacked person-condition model: PR osteoporosis vs each condition,
#         by sex and both sexes combined (adjusted for age group, race/ethnicity)
#       - condition x sex: ratio of sex PRs (supportive only)
#   A2  Within-person paired comparison (adults with osteoporosis who also have
#       >= 1 comparison condition): PR osteoporosis report vs comparison report,
#       by sex; report type x sex = Gate 3
#   Sensitivity: medication-inclusive definitions; all adults 50+ with hip DXA;
#       hypertension-only pairs; diabetes/cholesterol pairs in women
#
# Run:   Rscript nhanes_crosscondition_analysis.R [data_dir] [out_dir] [main_script]
#        defaults: nhanes_data  results_crosscondition  nhanes_cascade_analysis.R
# Definitions are shared with nhanes_cascade_analysis.R (sections 1-6 and fit_pr)
# and with nhanes_feasibility_crosscondition.R (comparison conditions).
# =============================================================================

args_f     <- commandArgs(trailingOnly = TRUE)
DATA_DIR_F <- if (length(args_f) >= 1) args_f[1] else "nhanes_data"
OUT_DIR_F  <- if (length(args_f) >= 2) args_f[2] else "results_crosscondition"
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
# 3. Analysis (plan dated 16 October 2026)
# =============================================================================
say("")
say("=== Cross-condition analysis (plan dated 2026-10-16) ===")

KEEP      <- c("SEQN", "psu", "strata", "w_pooled", "sex", "age_grp", "race_m", "analytic", "op", "engaged")
COND_VARS <- unique(na.omit(c(COND$measured, COND$medinc, COND$told)))
base      <- exam[, c(KEEP, COND_VARS)]
COND_LEV  <- COND$condition                         # Osteoporosis first = reference
COMPS     <- COND$condition[-1]

SAMPLES <- c("Health care contact (primary)", "All adults 50+ with hip DXA")
DEFS    <- c("Primary (measured)", "Medication-inclusive")
PRIMARY_SAMPLE <- SAMPLES[1]; PRIMARY_DEF <- DEFS[1]

in_sample <- function(d, s) {
  a <- d$analytic %in% TRUE
  if (s == SAMPLES[1]) a & d$engaged %in% "Engaged in care" else a
}
has_var <- function(k, def) {
  if (def == DEFS[2] && !is.na(COND$medinc[k])) COND$medinc[k] else COND$measured[k]
}
mk_design <- function(d, wcol)
  svydesign(ids = ~psu, strata = ~strata, weights = as.formula(paste0("~", wcol)), nest = TRUE, data = d)

# survey-weighted modified Poisson (log link); returns one row per coefficient
fit_model <- function(des, f, label, sample, def) {
  lab <- data.frame(sample = sample, definition = def, model = label)
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
# turn "comparison vs osteoporosis" into "osteoporosis vs comparison"
invert <- function(r) { l <- 1 / r$ucl; u <- 1 / r$lcl; r$PR <- 1 / r$PR; r$lcl <- l; r$ucl <- u; r }
fmt_pr <- function(r) if (!nrow(r) || is.na(r$PR[1])) "not estimable" else
  sprintf("%.2f (%.2f-%.2f)", r$PR[1], r$lcl[1], r$ucl[1])
fmt_pct <- function(r) if (is.na(r$pct)) paste0("NA (n=", r$n, ")") else
  sprintf("%5.1f%% (%.1f-%.1f) n=%d [%s]", r$pct, r$lcl, r$ucl, r$n, r$nchs_flag)

# -----------------------------------------------------------------------------
# A1. Population-level comparison
# -----------------------------------------------------------------------------
say("")
say("--- A1: population-level comparison ---")
va <- A$variables
a1p <- list()
for (s in SAMPLES) for (df_ in DEFS) for (sx in SEXES) for (k in seq_len(nrow(COND))) {
  if (df_ == DEFS[2] && is.na(COND$medinc[k]) && k != 1) next
  hv <- has_var(k, df_)
  des <- sub_des(A, in_sample(va, s) & sex_cond(va, sx) & va[[hv]] %in% 1)
  a1p[[length(a1p) + 1]] <- cbind(data.frame(sample = s, definition = df_, sex = sx, condition = COND$condition[k]),
                                  est_cell(des, COND$told[k]))
}
a1p <- do.call(rbind, a1p)
a1p <- a1p[, setdiff(names(a1p), "pop_events_millions")]
write_out(a1p, "crosscond_A1_percentages.csv")

stack_A1 <- function(df_) {
  st <- do.call(rbind, lapply(seq_len(nrow(COND)), function(k) {
    d <- base[, KEEP]
    d$cond <- COND$condition[k]
    d$has  <- base[[has_var(k, df_)]]
    d$told <- base[[COND$told[k]]]
    d
  }))
  st$cond <- factor(st$cond, levels = COND_LEV)
  st
}
a1m <- list(); a1x <- list()
for (df_ in DEFS) {
  S1 <- mk_design(stack_A1(df_), "w_pooled")
  vs <- S1$variables
  for (s in SAMPLES) {
    dom <- in_sample(vs, s) & vs$has %in% 1 & !is.na(vs$told)
    for (sx in SEXES) {
      des <- sub_des(S1, dom & sex_cond(vs, sx))
      f <- if (sx == "Both") told ~ cond + sex + age_grp + race_m else told ~ cond + age_grp + race_m
      r <- fit_model(des, f, paste0("Adjusted, ", sx), s, df_)
      r <- r[!is.na(r$term) & grepl("^cond", r$term), , drop = FALSE]
      if (nrow(r)) {
        r <- invert(r)
        r$contrast <- paste0("Osteoporosis vs ", sub("^cond", "", r$term)); r$sex <- sx
        a1m[[length(a1m) + 1]] <- r
      }
    }
    rx <- fit_model(sub_des(S1, dom), told ~ cond * sex + age_grp + race_m, "Condition x sex", s, df_)
    rx <- rx[!is.na(rx$term) & grepl("^cond.*:sexMen$", rx$term), , drop = FALSE]
    if (nrow(rx)) {
      rx <- invert(rx)
      rx$contrast <- paste0("Sex PR (men/women): osteoporosis relative to ",
                            sub(":sexMen$", "", sub("^cond", "", rx$term)))
      a1x[[length(a1x) + 1]] <- rx
    }
  }
}
a1m <- do.call(rbind, a1m); a1x <- if (length(a1x)) do.call(rbind, a1x) else data.frame()
write_out(a1m, "crosscond_A1_prevalence_ratios.csv")
write_out(a1x, "crosscond_A1_sex_specificity.csv")

# -----------------------------------------------------------------------------
# A2. Within-person paired comparison
# -----------------------------------------------------------------------------
say("")
say("--- A2: within-person paired comparison ---")
stack_A2 <- function(df_, conds) {
  ks <- which(COND$condition %in% conds)
  elig <- lapply(ks, function(k) base$analytic %in% TRUE & base$op %in% 1 & !is.na(base$told_op) &
                   base[[has_var(k, df_)]] %in% 1 & !is.na(base[[COND$told[k]]]))
  kcount <- Reduce(`+`, elig)
  out <- list()
  for (j in seq_along(ks)) for (tp in c("Comparison", "Osteoporosis")) {
    d <- base[, KEEP]
    d$pair_cond <- COND$condition[ks[j]]
    d$type      <- tp
    d$told      <- if (tp == "Osteoporosis") base$told_op else base[[COND$told[ks[j]]]]
    d$in_pair   <- elig[[j]]
    d$w_pair    <- base$w_pooled / pmax(kcount, 1)   # each person counts once in total
    out[[length(out) + 1]] <- d
  }
  st <- do.call(rbind, out)
  st$type <- factor(st$type, levels = c("Comparison", "Osteoporosis"))
  st
}
run_A2 <- function(df_, conds, set_label, samples = SAMPLES, sexes = SEXES, interaction = TRUE) {
  P  <- mk_design(stack_A2(df_, conds), "w_pair")
  vp <- P$variables
  pct <- list(); mod <- list()
  for (s in samples) {
    dom <- in_sample(vp, s) & vp$in_pair
    for (sx in sexes) {
      persons <- length(unique(vp$SEQN[dom & sex_cond(vp, sx)]))
      for (tp in levels(vp$type))
        pct[[length(pct) + 1]] <- cbind(data.frame(pairs = set_label, sample = s, definition = df_, sex = sx,
                                                   report = tp, persons = persons),
                                        est_cell(sub_des(P, dom & sex_cond(vp, sx) & vp$type == tp), "told"))
      r <- fit_model(sub_des(P, dom & sex_cond(vp, sx)), told ~ type, paste0("Paired, ", sx), s, df_)
      r <- r[!is.na(r$term) & r$term == "typeOsteoporosis", , drop = FALSE]
      if (nrow(r)) { r$pairs <- set_label; r$sex <- sx; r$persons <- persons
                     r$contrast <- "Osteoporosis report vs comparison report"; mod[[length(mod) + 1]] <- r }
    }
    if (interaction) {
      rx <- fit_model(sub_des(P, dom), told ~ type * sex, "Paired, report type x sex", s, df_)
      rx <- rx[!is.na(rx$term) & rx$term == "typeOsteoporosis:sexMen", , drop = FALSE]
      if (nrow(rx)) { rx$pairs <- set_label; rx$sex <- "Men vs women"; rx$persons <- NA
                      rx$contrast <- "Ratio of PRs (men/women): osteoporosis shortfall larger in men if < 1"
                      mod[[length(mod) + 1]] <- rx }
    }
  }
  list(pct = do.call(rbind, pct), mod = if (length(mod)) do.call(rbind, mod) else data.frame())
}
a2p <- list(); a2m <- list()
for (df_ in DEFS) {
  r <- run_A2(df_, COMPS, "All three conditions")
  a2p[[length(a2p) + 1]] <- r$pct; a2m[[length(a2m) + 1]] <- r$mod
}
# condition-specific pairs (primary definition and sample): hypertension both sexes;
# diabetes and cholesterol in women only (fewer than 30 men each)
r <- run_A2(PRIMARY_DEF, "Hypertension", "Hypertension only", samples = PRIMARY_SAMPLE)
a2p[[length(a2p) + 1]] <- r$pct; a2m[[length(a2m) + 1]] <- r$mod
for (cn in c("Diabetes", "High cholesterol")) {
  r <- run_A2(PRIMARY_DEF, cn, paste(cn, "only"), samples = PRIMARY_SAMPLE, sexes = "Women", interaction = FALSE)
  a2p[[length(a2p) + 1]] <- r$pct; a2m[[length(a2m) + 1]] <- r$mod
}
a2p <- do.call(rbind, a2p); a2p <- a2p[, setdiff(names(a2p), "pop_events_millions")]
a2m <- do.call(rbind, a2m)
write_out(a2p, "crosscond_A2_paired_percentages.csv")
write_out(a2m, "crosscond_A2_prevalence_ratios.csv")

# -----------------------------------------------------------------------------
# 4. Results summary and decision rules (as fixed in the plan)
# -----------------------------------------------------------------------------
say("")
say("=================== RESULTS (primary: health care contact, measured definitions) ===================")
say("A1 - % reporting a prior diagnosis, by condition:")
for (sx in c("Women", "Men")) for (cn in COND$condition) {
  r <- a1p[a1p$sample == PRIMARY_SAMPLE & a1p$definition == PRIMARY_DEF & a1p$sex == sx & a1p$condition == cn, ][1, ]
  say(sprintf("  %-5s %-17s %s", sx, cn, fmt_pct(r)))
}
say("A1 - PR, osteoporosis vs each condition (adjusted):")
for (sx in SEXES) for (cn in COMPS) {
  r <- a1m[a1m$sample == PRIMARY_SAMPLE & a1m$definition == PRIMARY_DEF & a1m$sex == sx &
             a1m$contrast == paste0("Osteoporosis vs ", cn), ]
  say(sprintf("  %-5s vs %-17s %s", sx, cn, fmt_pr(r)))
}
say("A2 - paired (adults with osteoporosis + >= 1 comparison condition), all three conditions:")
for (sx in SEXES) {
  ro <- a2p[a2p$pairs == "All three conditions" & a2p$sample == PRIMARY_SAMPLE & a2p$definition == PRIMARY_DEF &
              a2p$sex == sx & a2p$report == "Osteoporosis", ][1, ]
  rc <- a2p[a2p$pairs == "All three conditions" & a2p$sample == PRIMARY_SAMPLE & a2p$definition == PRIMARY_DEF &
              a2p$sex == sx & a2p$report == "Comparison", ][1, ]
  rm <- a2m[a2m$pairs == "All three conditions" & a2m$sample == PRIMARY_SAMPLE & a2m$definition == PRIMARY_DEF &
              a2m$sex == sx, ]
  say(sprintf("  %-5s persons=%s  osteoporosis %s | comparison %s | PR %s", sx, ro$persons,
              fmt_pct(ro), fmt_pct(rc), fmt_pr(rm)))
}
gx <- a2m[a2m$pairs == "All three conditions" & a2m$sample == PRIMARY_SAMPLE & a2m$definition == PRIMARY_DEF &
            a2m$sex == "Men vs women", ]
say("  Ratio of PRs (men/women): ", fmt_pr(gx))

say("")
say("Decision rules (fixed 2026-10-16):")
g2a <- sapply(COMPS, function(cn) {
  r <- a1m[a1m$sample == PRIMARY_SAMPLE & a1m$definition == PRIMARY_DEF & a1m$sex == "Both" &
             a1m$contrast == paste0("Osteoporosis vs ", cn), ]
  nrow(r) > 0 && !is.na(r$ucl[1]) && r$ucl[1] < 1
})
r2 <- a2m[a2m$pairs == "All three conditions" & a2m$sample == PRIMARY_SAMPLE & a2m$definition == PRIMARY_DEF &
            a2m$sex == "Both", ]
g2b <- nrow(r2) > 0 && !is.na(r2$ucl[1]) && r2$ucl[1] < 1
g3  <- nrow(gx) > 0 && !is.na(gx$ucl[1]) && gx$ucl[1] < 1
say(sprintf("  Gate 2 (osteoporosis reported less often): A1 upper limit < 1 for %s; A2 upper limit < 1: %s -> %s",
            paste(paste0(COMPS, "=", ifelse(g2a, "yes", "no")), collapse = ", "),
            ifelse(g2b, "yes", "no"), ifelse(all(g2a) && g2b, "MET", "NOT MET")))
say(sprintf("  Gate 3 (sex gap specific to osteoporosis): A2 ratio of PRs upper limit < 1 -> %s",
            ifelse(g3, "MET", "NOT MET (report as no evidence the shortfall differed by sex)")))
say("")
say("Sensitivity results (medication-inclusive; all adults; condition-specific pairs) are in the CSV files.")
say("Report only percentages whose nchs_flag and nchs_flag_complement are 'Reliable'.")
say("Done. Files in ", normalizePath(OUT_DIR))
