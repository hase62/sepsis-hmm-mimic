
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

select     <- dplyr::select
filter     <- dplyr::filter
rename     <- dplyr::rename
mutate     <- dplyr::mutate
summarise  <- dplyr::summarise
summarize  <- dplyr::summarise
group_by   <- dplyr::group_by
ungroup    <- dplyr::ungroup
arrange    <- dplyr::arrange
left_join  <- dplyr::left_join
inner_join <- dplyr::inner_join
bind_rows  <- dplyr::bind_rows
all_of     <- tidyselect::all_of
any_of     <- tidyselect::any_of

OUT_PREFIX <- Sys.getenv(
  "MIMIC_CURRENTDAY_SOFA_OUT_PREFIX",
  "incremental_value_current_sofa"
)

BASE_ROWLEVEL_CANDIDATES <- c(
  Sys.getenv("MIMIC_INCREMENTAL_ROWLEVEL", ""),
  "incremental_value_admission_sofa_rowlevel_selected_days.csv",
  "incremental_value_diagnostics_rowlevel_selected_days_with_current_sofa.csv"
)

SOFA_JOINED_CANDIDATES <- c(
  Sys.getenv("MIMIC_PATIENTDAY_SOFA_JOINED", ""),
  "patient_day_sofa_organsupport_joined.csv"
)

SOFA_COMPONENT_CANDIDATES <- c(
  Sys.getenv("MIMIC_PATIENTDAY_SOFA_COMPONENTS", ""),
  "patient_day_sofa_components.csv"
)

SUPPORT_CANDIDATES <- c(
  Sys.getenv("MIMIC_PATIENTDAY_ORGANSUPPORT", ""),
  "patient_day_organsupport.csv"
)

SELECTED_DAYS <- c(1L, 3L, 5L, 7L, 10L, 14L, 21L, 28L)
CALIBRATION_DAYS <- c(1L, 3L, 7L)
LONGSTAY_DAY <- 21L
MIN_EVENTS <- suppressWarnings(as.integer(Sys.getenv("MIMIC_CURRENTDAY_SOFA_MIN_EVENTS", "10")))
MIN_NONEVENTS <- suppressWarnings(as.integer(Sys.getenv("MIMIC_CURRENTDAY_SOFA_MIN_NONEVENTS", "10")))
MIN_AVAILABLE_COMPONENTS <- suppressWarnings(as.integer(Sys.getenv("MIMIC_CURRENTDAY_SOFA_MIN_COMPONENTS", "4")))
BOOTSTRAP_N <- suppressWarnings(as.integer(Sys.getenv("MIMIC_CURRENTDAY_SOFA_BOOTSTRAP_N", "500")))
BOOTSTRAP_SEED <- suppressWarnings(as.integer(Sys.getenv("MIMIC_CURRENTDAY_SOFA_BOOTSTRAP_SEED", "20260612")))

if (!is.finite(MIN_EVENTS) || MIN_EVENTS < 1L) MIN_EVENTS <- 10L
if (!is.finite(MIN_NONEVENTS) || MIN_NONEVENTS < 1L) MIN_NONEVENTS <- 10L
if (!is.finite(MIN_AVAILABLE_COMPONENTS) || MIN_AVAILABLE_COMPONENTS < 1L) MIN_AVAILABLE_COMPONENTS <- 4L
if (!is.finite(BOOTSTRAP_N) || BOOTSTRAP_N < 0L) BOOTSTRAP_N <- 500L
if (!is.finite(BOOTSTRAP_SEED)) BOOTSTRAP_SEED <- 20260612L

`%||%` <- function(a, b) if (!is.null(a)) a else b
num <- function(x) suppressWarnings(as.numeric(x))
as_bool <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  if (is.numeric(x)) return(ifelse(is.na(x), FALSE, x != 0))
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y")
}
.stopif <- function(ok, msg) if (!ok) stop(msg, call. = FALSE)

.write_both <- function(df, stem) {
  readr::write_csv(df, paste0(stem, ".csv"))
  readr::write_tsv(df, paste0(stem, ".tsv"))
  invisible(df)
}

.first_existing <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit) || !nzchar(hit)) return(NA_character_)
  hit
}

.clean_subphenotype <- function(x) {
  x <- as.character(x)
  x <- dplyr::case_when(
    x %in% c("A1 low burden", "A1 low burden", "A1 low burden") ~ "A1 low burden",
    x %in% c("A2 mild burden", "A2 mild burden") ~ "A2 mild burden",
    x %in% c("B1 renal mild", "B1 renal mild") ~ "B1 renal mild",
    x %in% c("B2 renal severe", "B2 renal severe") ~ "B2 renal severe",
    x %in% c("G respiratory", "G respiratory") ~ "G respiratory",
    x %in% c("D1 lactate/shock", "D1 lactate/shock") ~ "D1 lactate/shock",
    x %in% c("D2 hepatic", "D2 hepatic") ~ "D2 hepatic",
    x %in% c("D3 shock/coagulation-overlap", "D3 shock/coagulation-overlap", "D3 shock/coag overlap") ~ "D3 shock/coag overlap",
    is.na(x) | !nzchar(x) ~ NA_character_,
    TRUE ~ x
  )
  x
}

.subphenotype_levels <- c(
  "A1 low burden",
  "A2 mild burden",
  "B1 renal mild",
  "B2 renal severe",
  "G respiratory",
  "D1 lactate/shock",
  "D2 hepatic",
  "D3 shock/coag overlap"
)

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
  tryCatch(
    stats::predict(fit, newdata = newdata, type = "response"),
    error = function(e) rep(NA_real_, if (is.null(newdata)) length(stats::fitted(fit)) else nrow(newdata))
  )
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

extract_model_metrics <- function(fit, model_id, outcome, day, n, events, severity_adjustment, severity_variable) {
  if (is.null(fit)) {
    return(tibble(
      outcome = outcome, day = day, model_id = model_id, n = n, events = events,
      event_rate = ifelse(n > 0, events / n, NA_real_), aic = NA_real_, bic = NA_real_,
      logLik = NA_real_, auc = NA_real_, severity_adjustment = severity_adjustment,
      severity_variable = severity_variable
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
    severity_adjustment = severity_adjustment,
    severity_variable = severity_variable
  )
}

lrt_safe <- function(reduced, full, outcome, day, comparison, severity_adjustment, severity_variable) {
  if (is.null(reduced) || is.null(full)) {
    return(tibble(outcome = outcome, day = day, comparison = comparison, df = NA_real_, chisq = NA_real_, p_value = NA_real_, severity_adjustment = severity_adjustment, severity_variable = severity_variable))
  }
  out <- tryCatch(stats::anova(reduced, full, test = "Chisq"), error = function(e) NULL)
  if (is.null(out) || nrow(out) < 2) {
    return(tibble(outcome = outcome, day = day, comparison = comparison, df = NA_real_, chisq = NA_real_, p_value = NA_real_, severity_adjustment = severity_adjustment, severity_variable = severity_variable))
  }
  tibble(
    outcome = outcome,
    day = day,
    comparison = comparison,
    df = suppressWarnings(as.numeric(out$Df[2])),
    chisq = suppressWarnings(as.numeric(out$Deviance[2])),
    p_value = suppressWarnings(as.numeric(out$`Pr(>Chi)`[2])),
    severity_adjustment = severity_adjustment,
    severity_variable = severity_variable
  )
}

coef_table_m2 <- function(fit, outcome, day, n, events, severity_adjustment, severity_variable) {
  if (is.null(fit)) return(tibble())
  sm <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(sm)) return(tibble())
  cf <- as.data.frame(sm)
  cf$term <- rownames(cf)
  rownames(cf) <- NULL
  names(cf)[names(cf) == "Estimate"] <- "estimate"
  names(cf)[names(cf) == "Std. Error"] <- "std_error"
  names(cf)[names(cf) == "z value"] <- "z_value"
  names(cf)[names(cf) == "Pr(>|z|)"] <- "p_value"
  cf %>%
    as_tibble() %>%
    mutate(
      outcome = outcome,
      day = day,
      n = n,
      events = events,
      odds_ratio = exp(estimate),
      ci_low = exp(estimate - 1.96 * std_error),
      ci_high = exp(estimate + 1.96 * std_error),
      severity_adjustment = severity_adjustment,
      severity_variable = severity_variable
    ) %>%
    select(outcome, day, term, estimate, std_error, z_value, p_value,
           odds_ratio, ci_low, ci_high, n, events, severity_adjustment, severity_variable)
}

make_formula <- function(outcome_nm, covariates) {
  stats::as.formula(paste(outcome_nm, "~", paste(covariates, collapse = " + ")))
}

prepare_day_data <- function(df, outcome_nm, d, spec) {
  base_cols <- c("age", "sex", "current_subphenotype", outcome_nm)
  if (identical(spec$type, "total")) {
    cols <- c(base_cols, spec$var)
    dat <- df %>%
      filter(day == d) %>%
      select(all_of(cols), any_of("sofa_component_available_n")) %>%
      rename(y = all_of(outcome_nm), severity_total = all_of(spec$var)) %>%
      filter(!is.na(y), !is.na(age), !is.na(sex), !is.na(severity_total), !is.na(current_subphenotype))
    if (!is.null(spec$min_components) && "sofa_component_available_n" %in% names(dat)) {
      dat <- dat %>% filter(!is.na(sofa_component_available_n) & sofa_component_available_n >= spec$min_components)
    }
  } else if (identical(spec$type, "components")) {
    cols <- c(base_cols, spec$vars)
    dat <- df %>%
      filter(day == d) %>%
      select(all_of(cols), any_of(c("sofa_total_complete6", "sofa_component_available_n", "sofa_complete6_available_flag"))) %>%
      rename(y = all_of(outcome_nm)) %>%
      filter(!is.na(y), !is.na(age), !is.na(sex), !is.na(current_subphenotype)) %>%
      filter(dplyr::if_all(all_of(spec$vars), ~ !is.na(.x)))

    if ("sofa_complete6_available_flag" %in% names(dat)) {
      dat <- dat %>% filter(as_bool(sofa_complete6_available_flag))
    } else {
      if ("sofa_total_complete6" %in% names(dat)) {
        dat <- dat %>% filter(is.finite(as.numeric(sofa_total_complete6)))
      }
      if ("sofa_component_available_n" %in% names(dat)) {
        dat <- dat %>% filter(!is.na(sofa_component_available_n) & sofa_component_available_n >= 6)
      }
    }
  } else {
    stop("Unknown severity spec type: ", spec$type, call. = FALSE)
  }

  dat %>%
    mutate(
      y = as.integer(y),
      age = num(age),
      sex = factor(sex),
      current_subphenotype = factor(as.character(current_subphenotype), levels = .subphenotype_levels)
    ) %>%
    filter(!is.na(age), !is.na(sex), !is.na(current_subphenotype)) %>%
    droplevels()
}

fit_model_set <- function(dat, spec) {
  if (identical(spec$type, "total")) {
    m1_covars <- c("age", "sex", "severity_total")
  } else {
    m1_covars <- c("age", "sex", spec$vars)
  }
  list(
    m0 = fit_glm_safe(y ~ age + sex, dat),
    m1 = fit_glm_safe(make_formula("y", m1_covars), dat),
    m2 = fit_glm_safe(make_formula("y", c(m1_covars, "current_subphenotype")), dat)
  )
}

run_model_grid <- function(df, spec) {
  metrics <- list()
  lrts <- list()
  coefs <- list()
  idx <- 0L

  for (d in SELECTED_DAYS) {
    for (outcome_nm in c("death_3d", "death_7d", "longstay21")) {
      if (outcome_nm == "longstay21" && d >= LONGSTAY_DAY) next

      dat <- prepare_day_data(df, outcome_nm, d, spec)
      n <- nrow(dat)
      events <- sum(dat$y == 1L, na.rm = TRUE)
      nonevents <- sum(dat$y == 0L, na.rm = TRUE)
      severity_variable <- if (identical(spec$type, "total")) spec$var else paste(spec$vars, collapse = "+")

      if (n < (MIN_EVENTS + MIN_NONEVENTS) || events < MIN_EVENTS || nonevents < MIN_NONEVENTS || nlevels(dat$current_subphenotype) < 2) {
        idx <- idx + 1L
        metrics[[idx]] <- tibble(
          outcome = outcome_nm,
          day = d,
          model_id = c("M0_age_sex", "M1_age_sex_currentday_severity", "M2_age_sex_currentday_severity_HMM"),
          n = n,
          events = events,
          event_rate = ifelse(n > 0, events / n, NA_real_),
          aic = NA_real_,
          bic = NA_real_,
          logLik = NA_real_,
          auc = NA_real_,
          severity_adjustment = spec$label,
          severity_variable = severity_variable,
          skipped_reason = sprintf("Insufficient usable rows/events or HMM levels: n=%d, events=%d, nonevents=%d, levels=%d", n, events, nonevents, nlevels(dat$current_subphenotype))
        )
        next
      }

      fits <- fit_model_set(dat, spec)
      idx <- idx + 1L
      metrics[[idx]] <- bind_rows(
        extract_model_metrics(fits$m0, "M0_age_sex", outcome_nm, d, n, events, spec$label, severity_variable),
        extract_model_metrics(fits$m1, "M1_age_sex_currentday_severity", outcome_nm, d, n, events, spec$label, severity_variable),
        extract_model_metrics(fits$m2, "M2_age_sex_currentday_severity_HMM", outcome_nm, d, n, events, spec$label, severity_variable)
      ) %>% mutate(skipped_reason = NA_character_)

      lrts[[length(lrts) + 1L]] <- bind_rows(
        lrt_safe(fits$m0, fits$m1, outcome_nm, d, "M1_vs_M0_add_currentday_severity", spec$label, severity_variable),
        lrt_safe(fits$m1, fits$m2, outcome_nm, d, "M2_vs_M1_add_HMM_subphenotype", spec$label, severity_variable)
      )

      coefs[[length(coefs) + 1L]] <- coef_table_m2(fits$m2, outcome_nm, d, n, events, spec$label, severity_variable)
    }
  }

  list(metrics = bind_rows(metrics), lrt = bind_rows(lrts), coefs = bind_rows(coefs))
}

make_hierarchy_table <- function(metrics, lrt) {
  wide <- metrics %>%
    filter(model_id %in% c("M0_age_sex", "M1_age_sex_currentday_severity", "M2_age_sex_currentday_severity_HMM")) %>%
    select(outcome, day, severity_adjustment, severity_variable, model_id, n, events, event_rate, aic, bic, auc) %>%
    pivot_wider(names_from = model_id, values_from = c(aic, bic, auc), names_sep = "__") %>%
    mutate(
      delta_auc_M2_minus_M1 = auc__M2_age_sex_currentday_severity_HMM - auc__M1_age_sex_currentday_severity,
      delta_auc_M1_minus_M0 = auc__M1_age_sex_currentday_severity - auc__M0_age_sex,
      delta_aic_M2_minus_M1 = aic__M2_age_sex_currentday_severity_HMM - aic__M1_age_sex_currentday_severity,
      delta_bic_M2_minus_M1 = bic__M2_age_sex_currentday_severity_HMM - bic__M1_age_sex_currentday_severity
    )

  lrt_hmm <- lrt %>%
    filter(comparison == "M2_vs_M1_add_HMM_subphenotype") %>%
    select(outcome, day, severity_adjustment, severity_variable, hmm_lrt_df = df, hmm_lrt_chisq = chisq, hmm_lrt_p = p_value)

  wide %>%
    left_join(lrt_hmm, by = c("outcome", "day", "severity_adjustment", "severity_variable")) %>%
    mutate(
      interpretation_window = case_when(
        day <= 14 ~ "primary_interpretation_window_day1_14",
        TRUE ~ "late_day_exploratory_due_to_risk_set_depletion"
      )
    ) %>%
    arrange(severity_adjustment, outcome, day)
}

bootstrap_auc_delta <- function(df, spec) {
  if (BOOTSTRAP_N == 0L) {
    return(tibble(
      outcome = character(),
      day = integer(),
      severity_adjustment = character(),
      severity_variable = character(),
      n = integer(),
      events = integer(),
      bootstrap_n_requested = integer(),
      bootstrap_n_valid = integer(),
      delta_auc_M2_minus_M1 = double(),
      delta_auc_bootstrap_mean = double(),
      delta_auc_ci_low = double(),
      delta_auc_ci_high = double()
    ))
  }
  set.seed(BOOTSTRAP_SEED)
  out <- list()
  k <- 0L
  severity_variable <- if (identical(spec$type, "total")) spec$var else paste(spec$vars, collapse = "+")

  for (d in SELECTED_DAYS) {
    for (outcome_nm in c("death_3d", "death_7d", "longstay21")) {
      if (outcome_nm == "longstay21" && d >= LONGSTAY_DAY) next
      dat <- prepare_day_data(df, outcome_nm, d, spec)
      n <- nrow(dat)
      events <- sum(dat$y == 1L, na.rm = TRUE)
      nonevents <- sum(dat$y == 0L, na.rm = TRUE)
      if (n < (MIN_EVENTS + MIN_NONEVENTS) || events < MIN_EVENTS || nonevents < MIN_NONEVENTS || nlevels(dat$current_subphenotype) < 2) next

      fits <- fit_model_set(dat, spec)
      if (is.null(fits$m1) || is.null(fits$m2)) next
      point_delta <- auc_rank(dat$y, pred_safe(fits$m2)) - auc_rank(dat$y, pred_safe(fits$m1))

      vals <- rep(NA_real_, BOOTSTRAP_N)
      for (b in seq_len(BOOTSTRAP_N)) {
        ii <- sample.int(n, size = n, replace = TRUE)
        db <- dat[ii, , drop = FALSE] %>% droplevels()
        if (sum(db$y == 1L) < 2L || sum(db$y == 0L) < 2L || nlevels(db$current_subphenotype) < 2) next
        fb <- fit_model_set(db, spec)
        if (is.null(fb$m1) || is.null(fb$m2)) next
        vals[b] <- auc_rank(db$y, pred_safe(fb$m2)) - auc_rank(db$y, pred_safe(fb$m1))
      }
      vals_ok <- vals[is.finite(vals)]
      k <- k + 1L
      out[[k]] <- tibble(
        outcome = outcome_nm,
        day = d,
        severity_adjustment = spec$label,
        severity_variable = severity_variable,
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

run_calibration <- function(df, spec) {
  out <- list()
  k <- 0L
  severity_variable <- if (identical(spec$type, "total")) spec$var else paste(spec$vars, collapse = "+")
  for (d in CALIBRATION_DAYS) {
    for (outcome_nm in c("death_3d", "death_7d")) {
      dat <- prepare_day_data(df, outcome_nm, d, spec)
      n <- nrow(dat)
      events <- sum(dat$y == 1L, na.rm = TRUE)
      nonevents <- sum(dat$y == 0L, na.rm = TRUE)
      if (n < (MIN_EVENTS + MIN_NONEVENTS) || events < MIN_EVENTS || nonevents < MIN_NONEVENTS || nlevels(dat$current_subphenotype) < 2) next
      fits <- fit_model_set(dat, spec)
      for (mid in c("m1", "m2")) {
        model_id <- if (mid == "m1") "M1_age_sex_currentday_severity" else "M2_age_sex_currentday_severity_HMM"
        p <- pred_safe(fits[[mid]])
        cm <- calibration_metrics(dat$y, p)
        k <- k + 1L
        out[[k]] <- cm %>%
          mutate(outcome = outcome_nm, day = d, model_id = model_id, severity_adjustment = spec$label, severity_variable = severity_variable) %>%
          select(outcome, day, model_id, severity_adjustment, severity_variable, n, events, brier, calibration_intercept, calibration_slope)
      }
    }
  }
  bind_rows(out)
}

row_fp <- .first_existing(BASE_ROWLEVEL_CANDIDATES)
.stopif(!is.na(row_fp), paste("No row-level selected-day incremental-value input found among:", paste(BASE_ROWLEVEL_CANDIDATES, collapse = ", ")))
message("[current-day SOFA] Reading row-level HMM/outcome input: ", row_fp)

row_df <- readr::read_csv(row_fp, show_col_types = FALSE)
if (!"row_id" %in% names(row_df)) row_df$row_id <- NA_integer_
row_df <- row_df %>%
  mutate(
    row_id = as.integer(row_id),
    stay_id = as.character(stay_id),
    day = as.integer(day),
    age = num(age),
    sex = factor(sex),
    current_subphenotype = factor(.clean_subphenotype(current_subphenotype), levels = .subphenotype_levels),
    death_3d = as.integer(death_3d),
    death_7d = as.integer(death_7d),
    longstay21 = as.integer(longstay21)
  )

if ("A1 low burden" %in% levels(row_df$current_subphenotype)) {
  row_df$current_subphenotype <- stats::relevel(row_df$current_subphenotype, ref = "A1 low burden")
}

joined_fp <- .first_existing(SOFA_JOINED_CANDIDATES)
sofa_component_fp <- .first_existing(SOFA_COMPONENT_CANDIDATES)
support_fp <- .first_existing(SUPPORT_CANDIDATES)

if (!is.na(joined_fp)) {
  message("[current-day SOFA] Reading patient-day SOFA/support joined input: ", joined_fp)
  sofa_support <- readr::read_csv(joined_fp, show_col_types = FALSE)
} else {
  .stopif(!is.na(sofa_component_fp), paste("SOFA component file not found among:", paste(SOFA_COMPONENT_CANDIDATES, collapse = ", ")))
  message("[current-day SOFA] Reading patient-day SOFA component input: ", sofa_component_fp)
  sofa_support <- readr::read_csv(sofa_component_fp, show_col_types = FALSE)
  if (!is.na(support_fp)) {
    message("[current-day SOFA] Reading patient-day organ-support input: ", support_fp)
    support <- readr::read_csv(support_fp, show_col_types = FALSE)
    support_keep <- support %>%
      select(any_of(c(
        "hmm_row_index", "stay_id", "subject_id", "hadm_id", "charttime_year_day", "lab_time",
        "peep_max", "peep_median", "fio2_max", "fio2_median", "vasopressor_any",
        "vasopressor_sofa_any", "rrt_any", "crrt_any", "hemodialysis_any", "ecmo_any",
        "ecmo_flow_max", "vent_mode_chart_any", "vent_type_chart_any", "o2_device_chart_any",
        "dialysis_indicator_input", "iabp_indicator_input", "ecmo_indicator_input"
      )))
    by_cols <- intersect(c("hmm_row_index", "stay_id", "lab_time"), names(sofa_support))
    by_cols <- intersect(by_cols, names(support_keep))
    if ("hmm_row_index" %in% by_cols) by_cols <- "hmm_row_index"
    .stopif(length(by_cols) > 0, "Could not determine join key between SOFA components and organ-support table.")
    sofa_support <- sofa_support %>% left_join(support_keep, by = by_cols, suffix = c("", "_support"))
  }
}

.stopif("hmm_row_index" %in% names(sofa_support) || all(c("stay_id", "lab_time") %in% names(sofa_support)),
        "Patient-day SOFA/support input must contain hmm_row_index or stay_id + lab_time.")

sofa_keep_cols <- c(
  "hmm_row_index", "stay_id", "subject_id", "hadm_id", "charttime_year_day", "lab_time",
  "sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns", "sofa_renal",
  "sofa_total_available_component", "sofa_total_complete6", "sofa_component_available_n",
  "resp_available", "coag_available", "liver_available", "cv_available", "cns_available", "renal_available",
  "pfr_min", "plt_min", "bili_max", "crea_max", "map_min", "gcs_min", "uo_24h", "ventilation_active",
  "peep_max", "peep_median", "fio2_max", "fio2_median", "vasopressor_any", "vasopressor_sofa_any",
  "rrt_any", "crrt_any", "hemodialysis_any", "ecmo_any", "ecmo_flow_max",
  "dialysis_indicator_input", "iabp_indicator_input", "ecmo_indicator_input",
  "vent_mode_chart_any", "vent_type_chart_any", "o2_device_chart_any"
)
sofa_support <- sofa_support %>%
  select(any_of(sofa_keep_cols))
if (!"hmm_row_index" %in% names(sofa_support)) sofa_support$hmm_row_index <- NA_integer_
if (!"stay_id" %in% names(sofa_support)) sofa_support$stay_id <- NA_character_
if (!"lab_time" %in% names(sofa_support)) sofa_support$lab_time <- NA_integer_
sofa_support <- sofa_support %>%
  mutate(
    hmm_row_index = as.integer(hmm_row_index),
    stay_id = as.character(stay_id),
    lab_time = as.integer(lab_time)
  )

if ("row_id" %in% names(row_df) && "hmm_row_index" %in% names(sofa_support) && any(!is.na(row_df$row_id))) {
  join_df <- sofa_support %>% rename(row_id = hmm_row_index)
  analysis_df <- row_df %>% left_join(join_df, by = "row_id", suffix = c("", "_sofa"))
  join_key_used <- "row_id == hmm_row_index"
} else {
  join_df <- sofa_support %>% rename(day = lab_time)
  analysis_df <- row_df %>% left_join(join_df, by = c("stay_id", "day"), suffix = c("", "_sofa"))
  join_key_used <- "stay_id + day == stay_id + lab_time"
}

analysis_df <- analysis_df %>%
  mutate(
    sofa_resp = num(sofa_resp),
    sofa_coag = num(sofa_coag),
    sofa_liver = num(sofa_liver),
    sofa_cv = num(sofa_cv),
    sofa_cns = num(sofa_cns),
    sofa_renal = num(sofa_renal),
    sofa_total_available_component = num(sofa_total_available_component),
    sofa_total_complete6 = num(sofa_total_complete6),
    sofa_component_available_n = num(sofa_component_available_n)
  )

analysis_df$resp_available_flag <- if ("resp_available" %in% names(analysis_df)) as_bool(analysis_df$resp_available) else is.finite(analysis_df$pfr_min)
analysis_df$coag_available_flag <- if ("coag_available" %in% names(analysis_df)) as_bool(analysis_df$coag_available) else is.finite(analysis_df$plt_min)
analysis_df$liver_available_flag <- if ("liver_available" %in% names(analysis_df)) as_bool(analysis_df$liver_available) else is.finite(analysis_df$bili_max)
analysis_df$cv_available_flag <- if ("cv_available" %in% names(analysis_df)) as_bool(analysis_df$cv_available) else is.finite(analysis_df$map_min)
analysis_df$cns_available_flag <- if ("cns_available" %in% names(analysis_df)) as_bool(analysis_df$cns_available) else is.finite(analysis_df$gcs_min)
analysis_df$renal_available_flag <- if ("renal_available" %in% names(analysis_df)) as_bool(analysis_df$renal_available) else is.finite(analysis_df$crea_max) | is.finite(analysis_df$uo_24h)
analysis_df$sofa_complete6_available_flag <- is.finite(analysis_df$sofa_total_complete6) &
  is.finite(analysis_df$sofa_component_available_n) &
  analysis_df$sofa_component_available_n >= 6 &
  analysis_df$resp_available_flag &
  analysis_df$coag_available_flag &
  analysis_df$liver_available_flag &
  analysis_df$cv_available_flag &
  analysis_df$cns_available_flag &
  analysis_df$renal_available_flag

qc_by_day <- analysis_df %>%
  group_by(day) %>%
  summarise(
    rows = n(),
    stays = n_distinct(stay_id),
    current_subphenotype_available_pct = 100 * mean(!is.na(current_subphenotype)),
    currentday_available_component_sofa_pct = 100 * mean(is.finite(sofa_total_available_component)),
    currentday_available_component_sofa_ge_min_components_pct = 100 * mean(is.finite(sofa_total_available_component) & sofa_component_available_n >= MIN_AVAILABLE_COMPONENTS, na.rm = TRUE),
    currentday_complete6_sofa_pct = 100 * mean(sofa_complete6_available_flag, na.rm = TRUE),
    respiratory_component_available_pct = 100 * mean(resp_available_flag, na.rm = TRUE),
    liver_component_available_pct = 100 * mean(liver_available_flag, na.rm = TRUE),
    ventilation_active_pct = if ("ventilation_active" %in% names(analysis_df)) 100 * mean(as.logical(ventilation_active), na.rm = TRUE) else NA_real_,
    vasopressor_any_pct = if ("vasopressor_any" %in% names(analysis_df)) 100 * mean(as.logical(vasopressor_any), na.rm = TRUE) else NA_real_,
    rrt_any_pct = if ("rrt_any" %in% names(analysis_df)) 100 * mean(as.logical(rrt_any), na.rm = TRUE) else NA_real_,
    ecmo_any_pct = if ("ecmo_any" %in% names(analysis_df)) 100 * mean(as.logical(ecmo_any), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>% arrange(day)

qc_overall <- tibble(
  metric = c(
    "rowlevel_input_file", "patientday_join_key", "rows", "unique_stays",
    "selected_days", "min_available_components_for_available_sofa_analysis",
    "available_component_sofa_rows_pct", "available_component_sofa_ge_min_components_rows_pct",
    "complete6_sofa_rows_pct", "resp_component_available_rows_pct", "liver_component_available_rows_pct"
  ),
  value = c(
    row_fp, join_key_used, as.character(nrow(analysis_df)), as.character(n_distinct(analysis_df$stay_id)),
    paste(SELECTED_DAYS, collapse = ","), as.character(MIN_AVAILABLE_COMPONENTS),
    sprintf("%.3f", 100 * mean(is.finite(analysis_df$sofa_total_available_component))),
    sprintf("%.3f", 100 * mean(is.finite(analysis_df$sofa_total_available_component) & analysis_df$sofa_component_available_n >= MIN_AVAILABLE_COMPONENTS, na.rm = TRUE)),
    sprintf("%.3f", 100 * mean(analysis_df$sofa_complete6_available_flag, na.rm = TRUE)),
    sprintf("%.3f", 100 * mean(analysis_df$resp_available_flag, na.rm = TRUE)),
    sprintf("%.3f", 100 * mean(analysis_df$liver_available_flag, na.rm = TRUE))
  )
)

.write_both(qc_by_day, paste0(OUT_PREFIX, "_QC_by_day"))
.write_both(qc_overall, paste0(OUT_PREFIX, "_QC_overall"))
.write_both(analysis_df, paste0(OUT_PREFIX, "_rowlevel_selected_days_with_currentday_SOFA_organsupport"))

severity_specs <- list(
  list(
    id = "available_component_total_ge_min",
    label = paste0("current-day available-component SOFA total (>=", MIN_AVAILABLE_COMPONENTS, " components)"),
    type = "total",
    var = "sofa_total_available_component",
    min_components = MIN_AVAILABLE_COMPONENTS
  ),
  list(
    id = "complete6_total",
    label = "current-day complete 6-component SOFA total",
    type = "total",
    var = "sofa_total_complete6",
    min_components = 6,
    complete6_required = TRUE
  ),
  list(
    id = "complete6_components",
    label = "current-day complete 6-component SOFA components",
    type = "components",
    vars = c("sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns", "sofa_renal")
  )
)

severity_plot_levels <- vapply(severity_specs, `[[`, character(1), "label")

for (sp in severity_specs) {
  req <- if (identical(sp$type, "total")) sp$var else sp$vars
  missing_req <- setdiff(req, names(analysis_df))
  .stopif(length(missing_req) == 0, paste("Missing required severity columns for", sp$label, ":", paste(missing_req, collapse = ", ")))
}

message("Running current-day SOFA-adjusted model grids...")
all_results <- purrr::map(severity_specs, ~ run_model_grid(analysis_df, .x))

model_metrics <- bind_rows(purrr::map(all_results, "metrics")) %>%
  group_by(outcome, day, severity_adjustment, severity_variable) %>%
  mutate(
    aic_M1 = aic[model_id == "M1_age_sex_currentday_severity"][1],
    bic_M1 = bic[model_id == "M1_age_sex_currentday_severity"][1],
    auc_M1 = auc[model_id == "M1_age_sex_currentday_severity"][1],
    delta_aic_vs_M1 = aic - aic_M1,
    delta_bic_vs_M1 = bic - bic_M1,
    delta_auc_vs_M1 = auc - auc_M1
  ) %>%
  ungroup() %>%
  select(-aic_M1, -bic_M1, -auc_M1)

lrt_summary <- bind_rows(purrr::map(all_results, "lrt"))
model2_coefficients <- bind_rows(purrr::map(all_results, "coefs"))
riskset_model_hierarchy <- make_hierarchy_table(model_metrics, lrt_summary)

.write_both(model_metrics, paste0(OUT_PREFIX, "_model_metrics"))
.write_both(lrt_summary, paste0(OUT_PREFIX, "_lrt_summary"))
.write_both(model2_coefficients, paste0(OUT_PREFIX, "_model2_coefficients"))
.write_both(riskset_model_hierarchy, paste0(OUT_PREFIX, "_riskset_model_hierarchy"))

message("[current-day SOFA] Core non-bootstrap outputs written:")
message("  - ", paste0(OUT_PREFIX, "_QC_overall.csv"))
message("  - ", paste0(OUT_PREFIX, "_QC_by_day.csv"))
message("  - ", paste0(OUT_PREFIX, "_model_metrics.csv"))
message("  - ", paste0(OUT_PREFIX, "_lrt_summary.csv"))
message("  - ", paste0(OUT_PREFIX, "_model2_coefficients.csv"))
message("  - ", paste0(OUT_PREFIX, "_riskset_model_hierarchy.csv"))

message("[current-day SOFA] Running AUROC bootstrap CIs; BOOTSTRAP_N=", BOOTSTRAP_N)
boot_list <- purrr::map(severity_specs, ~ bootstrap_auc_delta(analysis_df, .x))
auc_bootstrap <- bind_rows(boot_list)
.write_both(auc_bootstrap, paste0(OUT_PREFIX, "_auc_delta_bootstrap_ci"))

calibration_list <- purrr::map(severity_specs, ~ run_calibration(analysis_df, .x))
calibration_early_days <- bind_rows(calibration_list)
.write_both(calibration_early_days, paste0(OUT_PREFIX, "_calibration_early_days"))

reporting_summary <- riskset_model_hierarchy %>%
  left_join(
    auc_bootstrap %>%
      select(outcome, day, severity_adjustment, severity_variable, bootstrap_n_valid, delta_auc_ci_low, delta_auc_ci_high),
    by = c("outcome", "day", "severity_adjustment", "severity_variable")
  ) %>%
  mutate(
    primary_sensitivity_flag = case_when(
      str_detect(severity_adjustment, "available-component") ~ "primary_currentday_available_component_SOFA_sensitivity",
      str_detect(severity_adjustment, "complete 6-component SOFA total") ~ "complete6_formal_SOFA_total_sensitivity",
      str_detect(severity_adjustment, "complete.*SOFA components") ~ "complete6_component_adjusted_sensitivity",
      TRUE ~ "other"
    ),
    model_comparison_interpretation = "secondary descriptive model-comparison; not prospective calibrated prediction"
  ) %>%
  arrange(primary_sensitivity_flag, outcome, day)
.write_both(reporting_summary, paste0(OUT_PREFIX, "_reporting_summary"))

support_cols <- intersect(
  c("ventilation_active", "peep_max", "fio2_max", "vasopressor_any", "vasopressor_sofa_any", "rrt_any", "crrt_any", "hemodialysis_any", "ecmo_any", "uo_24h"),
  names(analysis_df)
)
if (length(support_cols) > 0) {
  support_summary <- analysis_df %>%
    filter(day %in% SELECTED_DAYS, !is.na(current_subphenotype)) %>%
    group_by(day, current_subphenotype) %>%
    summarise(
      n = n(),
      ventilation_active_pct = if ("ventilation_active" %in% support_cols) 100 * mean(as.logical(ventilation_active), na.rm = TRUE) else NA_real_,
      peep_max_median = if ("peep_max" %in% support_cols) suppressWarnings(stats::median(num(peep_max), na.rm = TRUE)) else NA_real_,
      fio2_max_median = if ("fio2_max" %in% support_cols) suppressWarnings(stats::median(num(fio2_max), na.rm = TRUE)) else NA_real_,
      vasopressor_any_pct = if ("vasopressor_any" %in% support_cols) 100 * mean(as.logical(vasopressor_any), na.rm = TRUE) else NA_real_,
      rrt_any_pct = if ("rrt_any" %in% support_cols) 100 * mean(as.logical(rrt_any), na.rm = TRUE) else NA_real_,
      ecmo_any_pct = if ("ecmo_any" %in% support_cols) 100 * mean(as.logical(ecmo_any), na.rm = TRUE) else NA_real_,
      urine_output_24h_median = if ("uo_24h" %in% support_cols) suppressWarnings(stats::median(num(uo_24h), na.rm = TRUE)) else NA_real_,
      .groups = "drop"
    ) %>%
    arrange(day, current_subphenotype)
  .write_both(support_summary, paste0(OUT_PREFIX, "_organsupport_descriptive_by_day_subphenotype"))
}

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))

  plot_auc <- reporting_summary %>%
    filter(is.finite(delta_auc_M2_minus_M1)) %>%
    mutate(
      outcome = factor(outcome, levels = c("death_3d", "death_7d", "longstay21")),
      severity_adjustment = factor(severity_adjustment, levels = severity_plot_levels)
    )
  if (nrow(plot_auc) > 0) {
    p_auc <- ggplot(plot_auc, aes(x = day, y = delta_auc_M2_minus_M1, group = severity_adjustment)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_line() +
      geom_point(size = 1.7) +
      facet_grid(outcome ~ severity_adjustment, scales = "free_y") +
      labs(
        title = "AUROC increment after adding current HMM subphenotype",
        subtitle = "Model 2: age + sex + current-day severity + HMM; Model 1: age + sex + current-day severity",
        x = "ICU day",
        y = "AUROC difference: Model 2 - Model 1"
      ) +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank(), strip.text.x = element_text(size = 7))
    ggsave(paste0(OUT_PREFIX, "_delta_auc_plot.pdf"), p_auc, width = 11, height = 6.5, bg = "white")
    ggsave(paste0(OUT_PREFIX, "_delta_auc_plot.png"), p_auc, width = 11, height = 6.5, dpi = 300, bg = "white")
  }

  plot_lrt <- lrt_summary %>%
    filter(comparison == "M2_vs_M1_add_HMM_subphenotype", is.finite(p_value)) %>%
    mutate(
      minus_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
      outcome = factor(outcome, levels = c("death_3d", "death_7d", "longstay21")),
      severity_adjustment = factor(severity_adjustment, levels = severity_plot_levels)
    )
  if (nrow(plot_lrt) > 0) {
    p_lrt <- ggplot(plot_lrt, aes(x = day, y = minus_log10_p, group = severity_adjustment)) +
      geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed") +
      geom_line() +
      geom_point(size = 1.7) +
      facet_grid(outcome ~ severity_adjustment, scales = "free_y") +
      labs(
        title = "Likelihood-ratio evidence for adding current HMM subphenotype",
        subtitle = "Comparison: age + sex + current-day severity + HMM versus age + sex + current-day severity",
        x = "ICU day",
        y = "-log10(P value)"
      ) +
      theme_bw(base_size = 9) +
      theme(panel.grid.minor = element_blank(), strip.text.x = element_text(size = 7))
    ggsave(paste0(OUT_PREFIX, "_hmm_lrt_plot.pdf"), p_lrt, width = 11, height = 6.5, bg = "white")
    ggsave(paste0(OUT_PREFIX, "_hmm_lrt_plot.png"), p_lrt, width = 11, height = 6.5, dpi = 300, bg = "white")
  }
}

manifest <- tibble(
  file = c(
    paste0(OUT_PREFIX, "_rowlevel_selected_days_with_currentday_SOFA_organsupport.csv"),
    paste0(OUT_PREFIX, "_QC_overall.csv"),
    paste0(OUT_PREFIX, "_QC_by_day.csv"),
    paste0(OUT_PREFIX, "_model_metrics.csv"),
    paste0(OUT_PREFIX, "_lrt_summary.csv"),
    paste0(OUT_PREFIX, "_model2_coefficients.csv"),
    paste0(OUT_PREFIX, "_riskset_model_hierarchy.csv"),
    paste0(OUT_PREFIX, "_auc_delta_bootstrap_ci.csv"),
    paste0(OUT_PREFIX, "_calibration_early_days.csv"),
    paste0(OUT_PREFIX, "_reporting_summary.csv"),
    paste0(OUT_PREFIX, "_organsupport_descriptive_by_day_subphenotype.csv"),
    paste0(OUT_PREFIX, "_delta_auc_plot.pdf"),
    paste0(OUT_PREFIX, "_hmm_lrt_plot.pdf")
  ),
  description = c(
    "Selected-day row-level dataset with current-day patient-day SOFA components and organ-support descriptors joined to HMM subphenotype/outcome rows.",
    "Overall current-day SOFA join and component-availability QC.",
    "Selected-day current-day SOFA join and component-availability QC by ICU day.",
    "AIC, BIC, log-likelihood, and AUROC for Model 0, Model 1, and Model 2 under each current-day severity adjustment.",
    "Likelihood-ratio tests for adding current-day severity and for adding HMM subphenotype to age/sex/current-day severity.",
    "Model 2 logistic-regression coefficients and Wald odds ratios.",
    "Risk-set/event-count and Model 0/1/2 hierarchy summary with HMM LRT for each severity adjustment.",
    "Bootstrap confidence intervals for Model 2 minus Model 1 AUROC increment.",
    "Early-day apparent calibration and Brier diagnostics for Model 1 and Model 2.",
    "Summary of current-day SOFA-adjusted HMM information.",
    "Descriptive organ-support burden by selected ICU day and current HMM subphenotype.",
    "Plot of AUROC increment after adding current HMM subphenotype to current-day severity model.",
    "Plot of likelihood-ratio evidence for adding current HMM subphenotype to current-day severity model."
  ),
  rowlevel_input_file = row_fp,
  sofa_joined_input_file = ifelse(is.na(joined_fp), NA_character_, joined_fp),
  sofa_component_input_file = ifelse(is.na(sofa_component_fp), NA_character_, sofa_component_fp),
  support_input_file = ifelse(is.na(support_fp), NA_character_, support_fp),
  join_key_used = join_key_used,
  selected_days = paste(SELECTED_DAYS, collapse = ","),
  min_available_components = MIN_AVAILABLE_COMPONENTS,
  bootstrap_n_requested = BOOTSTRAP_N,
  bootstrap_seed = BOOTSTRAP_SEED,
  note = "Secondary descriptive model-comparison sensitivity. Current-day available-component SOFA is not necessarily complete formal SOFA; complete6 total and component outputs are enforced complete-case sensitivities."
)
.write_both(manifest, paste0(OUT_PREFIX, "_manifest"))

message("[current-day SOFA] Done.")
message("[current-day SOFA] Output prefix: ", OUT_PREFIX)
message("[current-day SOFA] Row-level input: ", row_fp)
message("[current-day SOFA] Join key used: ", join_key_used)
print(qc_overall)
