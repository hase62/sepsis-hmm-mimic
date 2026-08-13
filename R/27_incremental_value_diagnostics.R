
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
if (!exists("MIMIC_ACTIVE_ANNOT_MASTER", inherits = TRUE)) {
  source(code_path("17_analysis_config.R"))
}

if (!file.exists(code_path("17_analysis_config.R"))) stop("Missing shared analysis configuration: 17_analysis_config.R", call. = FALSE)
source(code_path("17_analysis_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

select    <- dplyr::select
filter    <- dplyr::filter
rename    <- dplyr::rename
mutate    <- dplyr::mutate
summarise <- dplyr::summarise
summarize <- dplyr::summarise
group_by  <- dplyr::group_by
ungroup   <- dplyr::ungroup
arrange   <- dplyr::arrange
left_join <- dplyr::left_join
bind_rows <- dplyr::bind_rows
all_of    <- tidyselect::all_of
any_of    <- tidyselect::any_of

BASE_PREFIX <- Sys.getenv("MIMIC_INCREMENTAL_BASE_PREFIX", "incremental_value_admission_sofa")
OUT_PREFIX  <- Sys.getenv("MIMIC_INCREMENTAL_DIAG_PREFIX", "incremental_value_diagnostics")
SELECTED_DAYS <- c(1L, 3L, 5L, 7L, 10L, 14L, 21L, 28L)
CALIBRATION_DAYS <- c(1L, 3L, 7L)
LONGSTAY_DAY <- 21L
MIN_EVENTS <- 10L
MIN_NONEVENTS <- 10L
BOOTSTRAP_N <- suppressWarnings(as.integer(Sys.getenv("MIMIC_INCREMENTAL_BOOTSTRAP_N", "500")))
if (!is.finite(BOOTSTRAP_N) || BOOTSTRAP_N < 0L) BOOTSTRAP_N <- 500L
BOOTSTRAP_SEED <- suppressWarnings(as.integer(Sys.getenv("MIMIC_INCREMENTAL_BOOTSTRAP_SEED", "20260610")))
if (!is.finite(BOOTSTRAP_SEED)) BOOTSTRAP_SEED <- 20260610L
BOOTSTRAP_INCLUDE_CURRENT_SOFA <- tolower(Sys.getenv("MIMIC_INCREMENTAL_BOOTSTRAP_CURRENT_SOFA", "false")) %in% c("1", "true", "yes", "y")

`%||%` <- function(a, b) if (!is.null(a)) a else b
num <- function(x) suppressWarnings(as.numeric(x))

.stopif <- function(ok, msg) if (!ok) stop(msg, call. = FALSE)
.write_both <- function(df, stem) {
  readr::write_csv(df, paste0(stem, ".csv"))
  readr::write_tsv(df, paste0(stem, ".tsv"))
  invisible(df)
}
.first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit) || !nzchar(hit)) return(NA_character_)
  hit
}
.pick_col_name <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

auc_rank <- function(y, score) {
  y <- as.integer(y)
  score <- as.numeric(score)
  ok <- is.finite(score) & !is.na(y)
  y <- y[ok]
  score <- score[ok]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(score, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fit_glm_safe <- function(formula, data) {
  tryCatch(
    suppressWarnings(stats::glm(formula, data = data, family = stats::binomial())),
    error = function(e) NULL
  )
}

pred_safe <- function(fit, newdata = NULL) {
  if (is.null(fit)) return(rep(NA_real_, if (is.null(newdata)) 0L else nrow(newdata)))
  tryCatch(stats::predict(fit, newdata = newdata, type = "response"), error = function(e) rep(NA_real_, if (is.null(newdata)) length(stats::fitted(fit)) else nrow(newdata)))
}

clip_prob <- function(p, eps = 1e-6) pmin(pmax(as.numeric(p), eps), 1 - eps)
logit <- function(p) log(clip_prob(p) / (1 - clip_prob(p)))

calibration_metrics <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  ok <- !is.na(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  n <- length(y)
  events <- sum(y == 1L)
  nonevents <- sum(y == 0L)
  if (n < 20L || events < 5L || nonevents < 5L) {
    return(tibble(n = n, events = events, brier = NA_real_, calibration_intercept = NA_real_, calibration_slope = NA_real_))
  }
  lp <- logit(p)
  brier <- mean((y - p)^2)
  intercept <- tryCatch(coef(suppressWarnings(glm(y ~ 1, offset = lp, family = binomial())))[1], error = function(e) NA_real_)
  slope <- tryCatch(coef(suppressWarnings(glm(y ~ lp, family = binomial())))[2], error = function(e) NA_real_)
  tibble(n = n, events = events, brier = brier, calibration_intercept = as.numeric(intercept), calibration_slope = as.numeric(slope))
}

extract_model_metrics <- function(fit, model_id, outcome, day, n, events, severity_adjustment_label, sofa_variable) {
  if (is.null(fit)) {
    return(tibble(
      outcome = outcome, day = day, model_id = model_id, n = n, events = events,
      event_rate = ifelse(n > 0, events / n, NA_real_), aic = NA_real_, bic = NA_real_,
      logLik = NA_real_, auc = NA_real_, severity_adjustment = severity_adjustment_label,
      sofa_variable = sofa_variable
    ))
  }
  p <- pred_safe(fit)
  y <- stats::model.response(stats::model.frame(fit))
  tibble(
    outcome = outcome,
    day = day,
    model_id = model_id,
    n = n,
    events = events,
    event_rate = ifelse(n > 0, events / n, NA_real_),
    aic = stats::AIC(fit),
    bic = stats::BIC(fit),
    logLik = as.numeric(stats::logLik(fit)),
    auc = auc_rank(y, p),
    severity_adjustment = severity_adjustment_label,
    sofa_variable = sofa_variable
  )
}

lrt_safe <- function(reduced, full, outcome, day, comparison, severity_adjustment_label, sofa_variable) {
  if (is.null(reduced) || is.null(full)) {
    return(tibble(outcome = outcome, day = day, comparison = comparison, df = NA_real_, chisq = NA_real_, p_value = NA_real_, severity_adjustment = severity_adjustment_label, sofa_variable = sofa_variable))
  }
  out <- tryCatch(stats::anova(reduced, full, test = "Chisq"), error = function(e) NULL)
  if (is.null(out) || nrow(out) < 2) {
    return(tibble(outcome = outcome, day = day, comparison = comparison, df = NA_real_, chisq = NA_real_, p_value = NA_real_, severity_adjustment = severity_adjustment_label, sofa_variable = sofa_variable))
  }
  tibble(
    outcome = outcome,
    day = day,
    comparison = comparison,
    df = suppressWarnings(as.numeric(out$Df[2])),
    chisq = suppressWarnings(as.numeric(out$Deviance[2])),
    p_value = suppressWarnings(as.numeric(out$`Pr(>Chi)`[2])),
    severity_adjustment = severity_adjustment_label,
    sofa_variable = sofa_variable
  )
}

make_current_day_sofa <- function(combined_df, row_ids) {
  .stopif(exists("combined", inherits = TRUE), "Object `combined` is required to construct current-day available-component SOFA surrogate.")
  df <- as.data.frame(combined_df)
  .stopif(max(row_ids, na.rm = TRUE) <= nrow(df), "row_id exceeds nrow(combined).")
  x <- df[row_ids, , drop = FALSE]

  getv <- function(cands) {
    nm <- .pick_col_name(x, cands)
    if (is.na(nm)) return(rep(NA_real_, nrow(x)))
    num(x[[nm]])
  }

  PF <- getv(c("PF", "PaO2_FiO2", "pO2_FO2I"))
  if (all(!is.finite(PF))) {
    pao2 <- getv(c("PaO2", "pO2", "PO2"))
    fio2 <- getv(c("FiO2", "FIO2", "fio2"))
    fio2_frac <- ifelse(is.finite(fio2) & fio2 > 1.5, fio2 / 100, fio2)
    PF <- ifelse(is.finite(pao2) & is.finite(fio2_frac) & fio2_frac > 0, pao2 / fio2_frac, NA_real_)
  }
  plt <- getv(c("PlateletCount", "Platelets", "platelets"))
  bili <- getv(c("TotalBilirubin", "bilirubin", "Bilirubin"))
  cr <- getv(c("Creatinine", "creatinine"))
  map <- getv(c("MAP", "MeanBP", "MEAN"))

  sofa_resp <- dplyr::case_when(!is.finite(PF) ~ NA_real_, PF < 100 ~ 3, PF < 200 ~ 2, PF < 300 ~ 1, TRUE ~ 0)
  sofa_coag <- dplyr::case_when(!is.finite(plt) ~ NA_real_, plt < 20 ~ 4, plt < 50 ~ 3, plt < 100 ~ 2, plt < 150 ~ 1, TRUE ~ 0)
  sofa_liver <- dplyr::case_when(!is.finite(bili) ~ NA_real_, bili >= 12 ~ 4, bili >= 6 ~ 3, bili >= 2 ~ 2, bili >= 1.2 ~ 1, TRUE ~ 0)
  sofa_renal <- dplyr::case_when(!is.finite(cr) ~ NA_real_, cr >= 5 ~ 4, cr >= 3.5 ~ 3, cr >= 2 ~ 2, cr >= 1.2 ~ 1, TRUE ~ 0)
  sofa_cv <- dplyr::case_when(!is.finite(map) ~ NA_real_, map < 70 ~ 1, TRUE ~ 0)

  n_components <- rowSums(!is.na(cbind(sofa_resp, sofa_coag, sofa_liver, sofa_renal, sofa_cv)))
  sofa_current <- rowSums(cbind(sofa_resp, sofa_coag, sofa_liver, sofa_renal, sofa_cv), na.rm = TRUE)
  sofa_current[n_components == 0] <- NA_real_

  tibble(
    current_sofa_surrogate = sofa_current,
    current_sofa_n_components = n_components,
    current_sofa_resp = sofa_resp,
    current_sofa_coag = sofa_coag,
    current_sofa_liver = sofa_liver,
    current_sofa_renal = sofa_renal,
    current_sofa_cv = sofa_cv
  )
}

prepare_day_data <- function(df, outcome_nm, d, sofa_variable) {
  df %>%
    filter(day == d) %>%
    select(age, sex, all_of(sofa_variable), current_subphenotype, all_of(outcome_nm)) %>%
    rename(y = all_of(outcome_nm), sofa_value = all_of(sofa_variable)) %>%
    filter(!is.na(y), !is.na(age), !is.na(sex), !is.na(sofa_value), !is.na(current_subphenotype)) %>%
    mutate(
      y = as.integer(y),
      sex = factor(sex),
      current_subphenotype = factor(current_subphenotype)
    ) %>%
    droplevels()
}

fit_model_set <- function(dat) {
  list(
    m0 = fit_glm_safe(y ~ age + sex, dat),
    m1 = fit_glm_safe(y ~ age + sex + sofa_value, dat),
    m2 = fit_glm_safe(y ~ age + sex + sofa_value + current_subphenotype, dat)
  )
}

run_model_grid <- function(df, sofa_variable, severity_adjustment_label) {
  metrics <- list()
  lrts <- list()
  idx <- 0L

  for (d in SELECTED_DAYS) {
    for (outcome_nm in c("death_3d", "death_7d", "longstay21")) {
      if (outcome_nm == "longstay21" && d >= LONGSTAY_DAY) next

      dat <- prepare_day_data(df, outcome_nm, d, sofa_variable)
      n <- nrow(dat)
      events <- sum(dat$y == 1L, na.rm = TRUE)
      nonevents <- sum(dat$y == 0L, na.rm = TRUE)

      if (n < (MIN_EVENTS + MIN_NONEVENTS) || events < MIN_EVENTS || nonevents < MIN_NONEVENTS || nlevels(dat$current_subphenotype) < 2) {
        idx <- idx + 1L
        metrics[[idx]] <- tibble(
          outcome = outcome_nm, day = d,
          model_id = c("M0_age_sex", "M1_age_sex_sofa", "M2_age_sex_sofa_hmm"),
          n = n, events = events, event_rate = ifelse(n > 0, events / n, NA_real_),
          aic = NA_real_, bic = NA_real_, logLik = NA_real_, auc = NA_real_,
          severity_adjustment = severity_adjustment_label, sofa_variable = sofa_variable,
          skipped_reason = sprintf("Insufficient data/events: n=%d, events=%d, nonevents=%d, levels=%d", n, events, nonevents, nlevels(dat$current_subphenotype))
        )
        next
      }

      fits <- fit_model_set(dat)
      idx <- idx + 1L
      metrics[[idx]] <- bind_rows(
        extract_model_metrics(fits$m0, "M0_age_sex", outcome_nm, d, n, events, severity_adjustment_label, sofa_variable),
        extract_model_metrics(fits$m1, "M1_age_sex_sofa", outcome_nm, d, n, events, severity_adjustment_label, sofa_variable),
        extract_model_metrics(fits$m2, "M2_age_sex_sofa_hmm", outcome_nm, d, n, events, severity_adjustment_label, sofa_variable)
      ) %>% mutate(skipped_reason = NA_character_)

      lrts[[length(lrts) + 1L]] <- bind_rows(
        lrt_safe(fits$m0, fits$m1, outcome_nm, d, "M1_vs_M0_add_SOFA", severity_adjustment_label, sofa_variable),
        lrt_safe(fits$m1, fits$m2, outcome_nm, d, "M2_vs_M1_add_HMM_subphenotype", severity_adjustment_label, sofa_variable)
      )
    }
  }

  list(metrics = bind_rows(metrics), lrt = bind_rows(lrts))
}

make_hierarchy_table <- function(metrics, lrt) {
  wide <- metrics %>%
    filter(model_id %in% c("M0_age_sex", "M1_age_sex_sofa", "M2_age_sex_sofa_hmm")) %>%
    select(outcome, day, severity_adjustment, sofa_variable, model_id, n, events, event_rate, aic, bic, auc) %>%
    pivot_wider(names_from = model_id, values_from = c(aic, bic, auc), names_sep = "__") %>%
    mutate(
      delta_auc_M2_minus_M1 = auc__M2_age_sex_sofa_hmm - auc__M1_age_sex_sofa,
      delta_auc_M1_minus_M0 = auc__M1_age_sex_sofa - auc__M0_age_sex,
      delta_aic_M2_minus_M1 = aic__M2_age_sex_sofa_hmm - aic__M1_age_sex_sofa,
      delta_bic_M2_minus_M1 = bic__M2_age_sex_sofa_hmm - bic__M1_age_sex_sofa
    )

  lrt_hmm <- lrt %>%
    filter(comparison == "M2_vs_M1_add_HMM_subphenotype") %>%
    select(outcome, day, severity_adjustment, sofa_variable, hmm_lrt_df = df, hmm_lrt_chisq = chisq, hmm_lrt_p = p_value)

  wide %>%
    left_join(lrt_hmm, by = c("outcome", "day", "severity_adjustment", "sofa_variable")) %>%
    arrange(severity_adjustment, outcome, day)
}

bootstrap_auc_delta <- function(df, sofa_variable, severity_adjustment_label) {
  if (BOOTSTRAP_N == 0L) {
    return(tibble())
  }
  set.seed(BOOTSTRAP_SEED)
  out <- list()
  k <- 0L
  for (d in SELECTED_DAYS) {
    for (outcome_nm in c("death_3d", "death_7d", "longstay21")) {
      if (outcome_nm == "longstay21" && d >= LONGSTAY_DAY) next
      dat <- prepare_day_data(df, outcome_nm, d, sofa_variable)
      n <- nrow(dat)
      events <- sum(dat$y == 1L, na.rm = TRUE)
      nonevents <- sum(dat$y == 0L, na.rm = TRUE)
      if (n < (MIN_EVENTS + MIN_NONEVENTS) || events < MIN_EVENTS || nonevents < MIN_NONEVENTS || nlevels(dat$current_subphenotype) < 2) next

      fits <- fit_model_set(dat)
      p1 <- pred_safe(fits$m1)
      p2 <- pred_safe(fits$m2)
      point_delta <- auc_rank(dat$y, p2) - auc_rank(dat$y, p1)

      vals <- rep(NA_real_, BOOTSTRAP_N)
      for (b in seq_len(BOOTSTRAP_N)) {
        ii <- sample.int(n, size = n, replace = TRUE)
        db <- dat[ii, , drop = FALSE] %>% droplevels()
        if (sum(db$y == 1L) < 2L || sum(db$y == 0L) < 2L || nlevels(db$current_subphenotype) < 2) next
        fb <- fit_model_set(db)
        if (is.null(fb$m1) || is.null(fb$m2)) next
        vals[b] <- auc_rank(db$y, pred_safe(fb$m2)) - auc_rank(db$y, pred_safe(fb$m1))
      }
      vals_ok <- vals[is.finite(vals)]
      k <- k + 1L
      out[[k]] <- tibble(
        outcome = outcome_nm,
        day = d,
        severity_adjustment = severity_adjustment_label,
        sofa_variable = sofa_variable,
        n = n,
        events = events,
        bootstrap_n_requested = BOOTSTRAP_N,
        bootstrap_n_valid = length(vals_ok),
        delta_auc_M2_minus_M1 = point_delta,
        delta_auc_bootstrap_mean = if (length(vals_ok)) mean(vals_ok) else NA_real_,
        delta_auc_ci_low = if (length(vals_ok) >= 20L) as.numeric(stats::quantile(vals_ok, 0.025, na.rm = TRUE)) else NA_real_,
        delta_auc_ci_high = if (length(vals_ok) >= 20L) as.numeric(stats::quantile(vals_ok, 0.975, na.rm = TRUE)) else NA_real_
      )
    }
  }
  bind_rows(out)
}

run_calibration <- function(df, sofa_variable, severity_adjustment_label) {
  out <- list()
  k <- 0L
  for (d in CALIBRATION_DAYS) {
    for (outcome_nm in c("death_3d", "death_7d")) {
      dat <- prepare_day_data(df, outcome_nm, d, sofa_variable)
      n <- nrow(dat)
      events <- sum(dat$y == 1L, na.rm = TRUE)
      nonevents <- sum(dat$y == 0L, na.rm = TRUE)
      if (n < (MIN_EVENTS + MIN_NONEVENTS) || events < MIN_EVENTS || nonevents < MIN_NONEVENTS || nlevels(dat$current_subphenotype) < 2) next
      fits <- fit_model_set(dat)
      for (mid in c("m1", "m2")) {
        model_id <- if (mid == "m1") "M1_age_sex_sofa" else "M2_age_sex_sofa_hmm"
        p <- pred_safe(fits[[mid]])
        cm <- calibration_metrics(dat$y, p)
        k <- k + 1L
        out[[k]] <- cm %>%
          mutate(outcome = outcome_nm, day = d, model_id = model_id, severity_adjustment = severity_adjustment_label, sofa_variable = sofa_variable) %>%
          select(outcome, day, model_id, severity_adjustment, sofa_variable, n, events, brier, calibration_intercept, calibration_slope)
      }
    }
  }
  bind_rows(out)
}

row_fp <- paste0(BASE_PREFIX, "_rowlevel_selected_days.csv")
.stopif(file.exists(row_fp), paste("Required input not found:", row_fp, "Run state_outcomes_8 first."))
row_df <- readr::read_csv(row_fp, show_col_types = FALSE) %>%
  mutate(
    day = as.integer(day),
    age = num(age),
    sofa_admission = num(sofa_admission),
    sex = factor(sex),
    current_subphenotype = factor(current_subphenotype),
    death_3d = as.integer(death_3d),
    death_7d = as.integer(death_7d),
    longstay21 = as.integer(longstay21)
  )

if ("A1 low burden" %in% levels(row_df$current_subphenotype)) {
  row_df$current_subphenotype <- stats::relevel(row_df$current_subphenotype, ref = "A1 low burden")
}

primary_sofa_label <- "baseline_sofa"
if ("sofa_source" %in% names(row_df)) {
  src <- unique(as.character(row_df$sofa_source[!is.na(row_df$sofa_source) & nzchar(row_df$sofa_source)]))
  if (length(src) >= 1L) primary_sofa_label <- src[[1]]
}
message("Primary baseline severity source: ", primary_sofa_label)

current_sofa_available <- FALSE
if (exists("combined", inherits = TRUE) && "row_id" %in% names(row_df)) {
  message("Constructing current-day available-component SOFA surrogate from combined[row_id, ].")
  cur <- make_current_day_sofa(get("combined", inherits = TRUE), row_df$row_id)
  row_df <- bind_cols(row_df, cur)
  current_sofa_available <- sum(is.finite(row_df$current_sofa_surrogate)) > 10 && length(unique(row_df$current_sofa_surrogate[is.finite(row_df$current_sofa_surrogate)])) > 1
} else {
  warning("combined or row_id not available; current-day SOFA surrogate sensitivity will be skipped.")
  row_df$current_sofa_surrogate <- NA_real_
  row_df$current_sofa_n_components <- NA_real_
}

.write_both(row_df, paste0(OUT_PREFIX, "_rowlevel_selected_days_with_current_sofa"))

primary <- run_model_grid(row_df, "sofa_admission", primary_sofa_label)
primary_hierarchy <- make_hierarchy_table(primary$metrics, primary$lrt)

boot_primary <- bootstrap_auc_delta(row_df, "sofa_admission", primary_sofa_label)
if (nrow(boot_primary) > 0) {
  primary_hierarchy <- primary_hierarchy %>%
    left_join(
      boot_primary %>% select(outcome, day, severity_adjustment, sofa_variable, bootstrap_n_valid, delta_auc_ci_low, delta_auc_ci_high),
      by = c("outcome", "day", "severity_adjustment", "sofa_variable")
    )
}

cal_primary <- run_calibration(row_df, "sofa_admission", primary_sofa_label)

.write_both(primary$metrics, paste0(OUT_PREFIX, "_primary_model_metrics"))
.write_both(primary$lrt, paste0(OUT_PREFIX, "_primary_lrt_summary"))
.write_both(primary_hierarchy, paste0(OUT_PREFIX, "_primary_riskset_model_hierarchy"))
.write_both(boot_primary, paste0(OUT_PREFIX, "_primary_auc_delta_bootstrap_ci"))
.write_both(cal_primary, paste0(OUT_PREFIX, "_primary_calibration_early_days"))

current_hierarchy <- tibble()
boot_current <- tibble()
cal_current <- tibble()
if (current_sofa_available) {
  message("Running current-day available-component SOFA surrogate sensitivity.")
  curres <- run_model_grid(row_df, "current_sofa_surrogate", "current_day_available_component_sofa_surrogate")
  current_hierarchy <- make_hierarchy_table(curres$metrics, curres$lrt)
  if (BOOTSTRAP_INCLUDE_CURRENT_SOFA) {
    boot_current <- bootstrap_auc_delta(row_df, "current_sofa_surrogate", "current_day_available_component_sofa_surrogate")
    if (nrow(boot_current) > 0) {
      current_hierarchy <- current_hierarchy %>%
        left_join(
          boot_current %>% select(outcome, day, severity_adjustment, sofa_variable, bootstrap_n_valid, delta_auc_ci_low, delta_auc_ci_high),
          by = c("outcome", "day", "severity_adjustment", "sofa_variable")
        )
    }
  }
  cal_current <- run_calibration(row_df, "current_sofa_surrogate", "current_day_available_component_sofa_surrogate")

  .write_both(curres$metrics, paste0(OUT_PREFIX, "_currentday_sofa_model_metrics"))
  .write_both(curres$lrt, paste0(OUT_PREFIX, "_currentday_sofa_lrt_summary"))
  .write_both(current_hierarchy, paste0(OUT_PREFIX, "_currentday_sofa_riskset_model_hierarchy"))
  if (nrow(boot_current) > 0L) {
    .write_both(boot_current, paste0(OUT_PREFIX, "_currentday_sofa_auc_delta_bootstrap_ci"))
  } else {
    unlink(c(paste0(OUT_PREFIX, "_currentday_sofa_auc_delta_bootstrap_ci.csv"),
             paste0(OUT_PREFIX, "_currentday_sofa_auc_delta_bootstrap_ci.tsv")), force = TRUE)
  }
  .write_both(cal_current, paste0(OUT_PREFIX, "_currentday_sofa_calibration_early_days"))
} else {
  warning("Current-day SOFA surrogate was not sufficiently available; sensitivity outputs will be empty.")
  .write_both(current_hierarchy, paste0(OUT_PREFIX, "_currentday_sofa_riskset_model_hierarchy"))
}

combined_hierarchy <- bind_rows(primary_hierarchy, current_hierarchy) %>%
  mutate(
    interpretation_window = case_when(
      day <= 14 ~ "primary_interpretation_window_day1_14",
      TRUE ~ "late_day_exploratory_due_to_risk_set_depletion"
    )
  )
.write_both(combined_hierarchy, paste0(OUT_PREFIX, "_combined_riskset_model_hierarchy"))

combined_calibration <- bind_rows(cal_primary, cal_current)
.write_both(combined_calibration, paste0(OUT_PREFIX, "_combined_calibration_early_days"))

manifest <- tibble(
  file = c(
    paste0(OUT_PREFIX, "_rowlevel_selected_days_with_current_sofa.csv"),
    paste0(OUT_PREFIX, "_primary_riskset_model_hierarchy.csv"),
    paste0(OUT_PREFIX, "_primary_auc_delta_bootstrap_ci.csv"),
    paste0(OUT_PREFIX, "_primary_calibration_early_days.csv"),
    paste0(OUT_PREFIX, "_currentday_sofa_riskset_model_hierarchy.csv"),
    paste0(OUT_PREFIX, "_combined_riskset_model_hierarchy.csv"),
    paste0(OUT_PREFIX, "_combined_calibration_early_days.csv")
  ),
  description = c(
    "Selected-day row-level dataset with current-day available-component SOFA surrogate components.",
    "Primary day-1 SOFA-surrogate risk-set/event-count and Model 0/1/2 AUROC hierarchy table, with HMM LRT and bootstrap CI when available.",
    "Bootstrap confidence intervals for Model 2 minus Model 1 AUROC increment under day-1 SOFA-surrogate adjustment.",
    "Early-day calibration intercept, calibration slope, and Brier score for Model 1 and Model 2 under day-1 SOFA-surrogate adjustment.",
    "Sensitivity hierarchy table using current-day available-component SOFA surrogate instead of day-1 surrogate.",
    "Combined primary and current-day SOFA-surrogate hierarchy table with late-day exploratory flag.",
    "Combined early-day calibration diagnostics for primary and current-day SOFA-surrogate analyses."
  ),
  bootstrap_n_requested = BOOTSTRAP_N,
  bootstrap_seed = BOOTSTRAP_SEED,
  current_day_sofa_sensitivity_available = current_sofa_available,
  note = "Diagnostics are descriptive. They do not constitute a calibrated prospective bedside prediction model."
)
.write_both(manifest, paste0(OUT_PREFIX, "_manifest"))

message("Done: incremental-value diagnostics written with prefix: ", OUT_PREFIX)
message("Bootstrap N requested: ", BOOTSTRAP_N)
message("Current-day SOFA sensitivity available: ", current_sofa_available)
