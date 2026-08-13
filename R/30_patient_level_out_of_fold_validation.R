
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
if (!exists("MIMIC_ACTIVE_ANNOT_MASTER", inherits = TRUE)) {
  source(code_path("17_analysis_config.R"))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (!is.null(a)) a else b
num <- function(x) suppressWarnings(as.numeric(x))
.stopif <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

first_existing <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)][1]
  if (!length(hit) || is.na(hit) || !nzchar(hit)) return(NA_character_)
  hit
}

valid_col_name <- function(x) {
  length(x) == 1L && !is.na(x) && nzchar(x)
}

first_present_col <- function(candidates, available_names) {
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[candidates %in% available_names]
  if (!length(hit)) return(NA_character_)
  hit[1]
}

write_both <- function(df, stem) {
  readr::write_csv(df, paste0(stem, ".csv"))
  readr::write_tsv(df, paste0(stem, ".tsv"))
  invisible(df)
}

OUT_PREFIX <- Sys.getenv(
  "SEPSIS_HMM_OOF_OUT_PREFIX",
  "patient_level_out_of_fold_validation"
)

BASE_ROWLEVEL_CANDIDATES <- c(
  Sys.getenv("SEPSIS_HMM_OOF_ROWLEVEL", ""),
  "incremental_value_diagnostics_rowlevel_selected_days_with_current_sofa.csv",
  "incremental_value_admission_sofa_rowlevel_selected_days.csv"
)

FILTERING_ROWLEVEL_CANDIDATES <- c(
  Sys.getenv("SEPSIS_HMM_OOF_FILTERING_ROWLEVEL", ""),
  "filtering_sensitivity_dominant_rowlevel_selected_days.csv"
)

SOFA_JOINED_CANDIDATES <- c(
  Sys.getenv("SEPSIS_HMM_OOF_SOFA_JOINED", ""),
  "patient_day_sofa_organsupport_joined.csv",
  "patient_day_sofa_components.csv"
)

SELECTED_DAYS <- c(1L, 3L, 5L, 7L, 10L, 14L, 21L, 28L)
PRIMARY_DAYS <- 1:14
OUTCOMES <- c("death_3d", "death_7d", "longstay21")
ASSIGNMENTS <- c("smoothing", "filtering")
SEX_LEVELS <- NULL

N_FOLDS <- suppressWarnings(as.integer(Sys.getenv("SEPSIS_HMM_OOF_FOLDS", "5")))
N_REPEATS <- suppressWarnings(as.integer(Sys.getenv("SEPSIS_HMM_OOF_REPEATS", "10")))
CV_SEED <- suppressWarnings(as.integer(Sys.getenv("SEPSIS_HMM_OOF_SEED", "20260722")))
BOOTSTRAP_N <- suppressWarnings(as.integer(Sys.getenv("SEPSIS_HMM_OOF_BOOTSTRAP_N", "500")))
BOOTSTRAP_SEED <- suppressWarnings(as.integer(Sys.getenv("SEPSIS_HMM_OOF_BOOTSTRAP_SEED", "20260723")))
MIN_EVENTS <- suppressWarnings(as.integer(Sys.getenv("SEPSIS_HMM_OOF_MIN_EVENTS", "10")))
MIN_NONEVENTS <- suppressWarnings(as.integer(Sys.getenv("SEPSIS_HMM_OOF_MIN_NONEVENTS", "10")))
MIN_COMPONENTS <- suppressWarnings(as.integer(Sys.getenv("SEPSIS_HMM_OOF_MIN_COMPONENTS", "4")))

if (!is.finite(N_FOLDS) || N_FOLDS < 2L) N_FOLDS <- 5L
if (!is.finite(N_REPEATS) || N_REPEATS < 1L) N_REPEATS <- 10L
if (!is.finite(CV_SEED)) CV_SEED <- 20260722L
if (!is.finite(BOOTSTRAP_N) || BOOTSTRAP_N < 0L) BOOTSTRAP_N <- 500L
if (!is.finite(BOOTSTRAP_SEED)) BOOTSTRAP_SEED <- 20260723L
if (!is.finite(MIN_EVENTS) || MIN_EVENTS < 1L) MIN_EVENTS <- 10L
if (!is.finite(MIN_NONEVENTS) || MIN_NONEVENTS < 1L) MIN_NONEVENTS <- 10L
if (!is.finite(MIN_COMPONENTS) || MIN_COMPONENTS < 1L) MIN_COMPONENTS <- 4L

auc_rank <- function(y, p) {
  ok <- is.finite(p) & !is.na(y)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 < 1L || n0 < 1L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

brier_score <- function(y, p) {
  ok <- is.finite(p) & !is.na(y)
  if (!any(ok)) return(NA_real_)
  mean((as.numeric(p[ok]) - as.numeric(y[ok]))^2)
}

calibration_metrics <- function(y, p) {
  ok <- is.finite(p) & !is.na(y)
  y <- as.integer(y[ok])
  p <- pmin(pmax(as.numeric(p[ok]), 1e-6), 1 - 1e-6)
  if (length(y) < 20L || length(unique(y)) < 2L) {
    return(tibble(calibration_intercept = NA_real_, calibration_slope = NA_real_))
  }
  lp <- qlogis(p)
  int_fit <- try(glm(y ~ 1, family = binomial(), offset = lp), silent = TRUE)
  slope_fit <- try(glm(y ~ lp, family = binomial()), silent = TRUE)
  tibble(
    calibration_intercept = if (inherits(int_fit, "try-error")) NA_real_ else unname(coef(int_fit)[1]),
    calibration_slope = if (inherits(slope_fit, "try-error") || length(coef(slope_fit)) < 2L) NA_real_ else unname(coef(slope_fit)[2])
  )
}

metric_row <- function(y, p1, p2) {
  cal1 <- calibration_metrics(y, p1)
  cal2 <- calibration_metrics(y, p2)
  tibble(
    n = length(y),
    events = sum(y == 1L, na.rm = TRUE),
    nonevents = sum(y == 0L, na.rm = TRUE),
    event_rate = mean(y == 1L, na.rm = TRUE),
    auc_M1 = auc_rank(y, p1),
    auc_M2 = auc_rank(y, p2),
    delta_auc_M2_minus_M1 = auc_rank(y, p2) - auc_rank(y, p1),
    brier_M1 = brier_score(y, p1),
    brier_M2 = brier_score(y, p2),
    delta_brier_M2_minus_M1 = brier_score(y, p2) - brier_score(y, p1),
    calibration_intercept_M1 = cal1$calibration_intercept,
    calibration_slope_M1 = cal1$calibration_slope,
    calibration_intercept_M2 = cal2$calibration_intercept,
    calibration_slope_M2 = cal2$calibration_slope
  )
}

clean_subphenotype <- function(x) {
  x <- as.character(x)
  x <- sub("^A1 low burden$", "A1 low burden", x)
  x <- sub("^A2 mild burden$", "A2 mild burden", x)
  x <- sub("^B1 renal mild$", "B1 renal mild", x)
  x <- sub("^B2 renal severe$", "B2 renal severe", x)
  x <- sub("^G respiratory$", "G respiratory", x)
  x <- sub("^D1 lactate/shock$", "D1 lactate/shock", x)
  x <- sub("^D2 hepatic$", "D2 hepatic", x)
  x <- sub("^D3 shock/coagulation-overlap$", "D3 shock/coag overlap", x)
  x
}

SUBPHENOTYPE_LEVELS <- c(
  "A1 low burden", "A2 mild burden", "B1 renal mild",
  "B2 renal severe", "G respiratory", "D1 lactate/shock",
  "D2 hepatic", "D3 shock/coag overlap"
)

make_design <- function(df, severity_col, phenotype_col, include_hmm) {
  dat <- df %>%
    transmute(
      age = num(age),
      sex = factor(as.character(sex), levels = SEX_LEVELS),
      severity = num(.data[[severity_col]]),
      phenotype = factor(clean_subphenotype(.data[[phenotype_col]]), levels = SUBPHENOTYPE_LEVELS)
    )
  f <- if (include_hmm) ~ age + sex + severity + phenotype else ~ age + sex + severity
  stats::model.matrix(f, data = dat)
}

fit_predict_glm <- function(x_train, y_train, x_test) {
  fit <- try(
    stats::glm.fit(x = x_train, y = y_train, family = stats::binomial()),
    silent = TRUE
  )
  if (inherits(fit, "try-error") || is.null(fit$coefficients)) {
    return(rep(mean(y_train, na.rm = TRUE), nrow(x_test)))
  }
  beta <- fit$coefficients
  beta[!is.finite(beta)] <- 0
  eta <- drop(x_test %*% beta)
  p <- plogis(eta)
  p[!is.finite(p)] <- mean(y_train, na.rm = TRUE)
  pmin(pmax(p, 1e-6), 1 - 1e-6)
}

make_stratified_folds <- function(y, k, seed) {
  set.seed(seed)
  fold <- integer(length(y))
  for (cls in sort(unique(y))) {
    idx <- which(y == cls)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

run_repeated_cv <- function(df, outcome_col, severity_col, phenotype_col,
                            nfolds, nrepeats, seed) {
  y <- as.integer(df[[outcome_col]])
  x1 <- make_design(df, severity_col, phenotype_col, include_hmm = FALSE)
  x2 <- make_design(df, severity_col, phenotype_col, include_hmm = TRUE)

  pred1_mat <- matrix(NA_real_, nrow = nrow(df), ncol = nrepeats)
  pred2_mat <- matrix(NA_real_, nrow = nrow(df), ncol = nrepeats)
  repeat_metrics <- vector("list", nrepeats)

  for (r in seq_len(nrepeats)) {
    folds <- make_stratified_folds(y, nfolds, seed + r * 1009L)
    p1 <- rep(NA_real_, length(y))
    p2 <- rep(NA_real_, length(y))

    for (f in seq_len(nfolds)) {
      test <- folds == f
      train <- !test
      if (!any(test) || sum(train) < 20L || length(unique(y[train])) < 2L) next
      p1[test] <- fit_predict_glm(x1[train, , drop = FALSE], y[train], x1[test, , drop = FALSE])
      p2[test] <- fit_predict_glm(x2[train, , drop = FALSE], y[train], x2[test, , drop = FALSE])
    }

    pred1_mat[, r] <- p1
    pred2_mat[, r] <- p2
    repeat_metrics[[r]] <- metric_row(y, p1, p2) %>% mutate(repeat_id = r)
  }

  p1_mean <- rowMeans(pred1_mat, na.rm = TRUE)
  p2_mean <- rowMeans(pred2_mat, na.rm = TRUE)
  p1_mean[!is.finite(p1_mean)] <- NA_real_
  p2_mean[!is.finite(p2_mean)] <- NA_real_

  list(
    repeat_metrics = bind_rows(repeat_metrics),
    prediction = tibble(
      stay_id = as.character(df$stay_id),
      day = as.integer(df$day),
      y = y,
      p_M1 = p1_mean,
      p_M2 = p2_mean
    ),
    pooled_metrics = metric_row(y, p1_mean, p2_mean)
  )
}

bootstrap_paired_metrics <- function(pred_df, nboot, seed) {
  if (nboot <= 0L) return(tibble())
  set.seed(seed)
  n <- nrow(pred_df)
  out <- vector("list", nboot)
  for (b in seq_len(nboot)) {
    ii <- sample.int(n, n, replace = TRUE)
    z <- pred_df[ii, , drop = FALSE]
    if (length(unique(z$y)) < 2L) next
    out[[b]] <- metric_row(z$y, z$p_M1, z$p_M2) %>% mutate(bootstrap = b)
  }
  bind_rows(out)
}

q025 <- function(x) if (sum(is.finite(x)) >= 20L) as.numeric(quantile(x, 0.025, na.rm = TRUE)) else NA_real_
q975 <- function(x) if (sum(is.finite(x)) >= 20L) as.numeric(quantile(x, 0.975, na.rm = TRUE)) else NA_real_

row_fp <- first_existing(BASE_ROWLEVEL_CANDIDATES)
.stopif(!is.na(row_fp), paste0(
  "No row-level incremental-value file found. Tried: ",
  paste(BASE_ROWLEVEL_CANDIDATES, collapse = ", ")
))
message("[OOF] Reading row-level data: ", row_fp)

row_df <- readr::read_csv(row_fp, show_col_types = FALSE)
required_base <- c(
  "stay_id", "day", "age", "sex", "sofa_admission",
  "current_subphenotype", "death_3d", "death_7d", "longstay21"
)
missing_base <- setdiff(required_base, names(row_df))
.stopif(length(missing_base) == 0L,
        paste("Row-level input is missing:", paste(missing_base, collapse = ", ")))

row_df <- row_df %>%
  mutate(
    row_id = if ("row_id" %in% names(.)) as.integer(row_id) else NA_integer_,
    stay_id = as.character(stay_id),
    day = as.integer(day),
    age = num(age),
    sex = as.character(sex),
    sofa_admission = num(sofa_admission),
    current_subphenotype = clean_subphenotype(current_subphenotype),
    death_3d = as.integer(death_3d),
    death_7d = as.integer(death_7d),
    longstay21 = as.integer(longstay21)
  ) %>%
  filter(day %in% SELECTED_DAYS)

SEX_LEVELS <- sort(unique(as.character(row_df$sex[!is.na(row_df$sex)])))
.stopif(length(SEX_LEVELS) >= 2L, "Fewer than two sex levels were available for model fitting.")

filter_fp <- first_existing(FILTERING_ROWLEVEL_CANDIDATES)
.stopif(!is.na(filter_fp), paste0(
  "Filtering assignment file not found. Tried: ",
  paste(FILTERING_ROWLEVEL_CANDIDATES, collapse = ", ")
))
message("[OOF] Reading filtering assignments: ", filter_fp)

filter_df <- readr::read_csv(filter_fp, show_col_types = FALSE) %>%
  transmute(
    row_id_filter = if ("row_id" %in% names(.)) as.integer(row_id) else NA_integer_,
    stay_id = as.character(stay_id),
    day = if ("lab_time" %in% names(.)) as.integer(lab_time) else as.integer(day),
    filtering_subphenotype = clean_subphenotype(dominant_filtering)
  )

if (any(is.finite(row_df$row_id)) && any(is.finite(filter_df$row_id_filter))) {
  row_df <- row_df %>%
    left_join(filter_df %>% select(row_id_filter, filtering_subphenotype),
              by = c("row_id" = "row_id_filter"))
} else {
  row_df <- row_df %>%
    left_join(filter_df %>% select(stay_id, day, filtering_subphenotype),
              by = c("stay_id", "day"))
}

.stopif(sum(!is.na(row_df$filtering_subphenotype)) > 0L,
        "No filtering assignments could be joined to the outcome-model rows.")

current_severity_col <- NA_character_
current_severity_label <- NA_character_
current_component_col <- NA_character_

sofa_fp <- first_existing(SOFA_JOINED_CANDIDATES)
if (!is.na(sofa_fp)) {
  message("[OOF] Reading patient-day SOFA input: ", sofa_fp)
  sofa_df <- readr::read_csv(sofa_fp, show_col_types = FALSE)

  id_candidates <- c("hmm_row_index", "row_id")
  id_hit <- first_present_col(id_candidates, names(sofa_df))
  severity_candidates <- c(
    Sys.getenv("SEPSIS_HMM_OOF_CURRENT_SOFA_COLUMN", ""),
    "sofa_total_available_component", "available_component_sofa",
    "current_day_available_component_sofa", "sofa_total", "total_sofa"
  )
  severity_hit <- first_present_col(severity_candidates, names(sofa_df))
  component_candidates <- c(
    Sys.getenv("SEPSIS_HMM_OOF_CURRENT_SOFA_NCOMP_COLUMN", ""),
    "sofa_component_available_n",
    "sofa_n_components", "n_sofa_components", "available_component_count",
    "n_components_available", "current_sofa_n_components"
  )
  component_hit <- first_present_col(component_candidates, names(sofa_df))

  message(
    "[OOF] SOFA columns: join=", ifelse(valid_col_name(id_hit), id_hit, "not found"),
    "; total=", ifelse(valid_col_name(severity_hit), severity_hit, "not found"),
    "; n_components=", ifelse(valid_col_name(component_hit), component_hit, "not found")
  )

  if (valid_col_name(id_hit) && valid_col_name(severity_hit) && any(is.finite(row_df$row_id))) {
    add <- tibble(
      row_id_join = as.integer(sofa_df[[id_hit]]),
      current_sofa_exact = num(sofa_df[[severity_hit]]),
      current_sofa_exact_ncomp = if (valid_col_name(component_hit)) {
        num(sofa_df[[component_hit]])
      } else {
        rep(NA_real_, nrow(sofa_df))
      }
    ) %>%
      distinct(row_id_join, .keep_all = TRUE)
    row_df <- row_df %>% left_join(add, by = c("row_id" = "row_id_join"))
    current_severity_col <- "current_sofa_exact"
    current_component_col <- "current_sofa_exact_ncomp"
    current_severity_label <- paste0("current_day_available_component_SOFA:", severity_hit)
  }
}

if (is.na(current_severity_col) && "current_sofa_surrogate" %in% names(row_df)) {
  row_df$current_sofa_surrogate <- num(row_df$current_sofa_surrogate)
  if (!"current_sofa_n_components" %in% names(row_df)) row_df$current_sofa_n_components <- NA_real_
  row_df$current_sofa_n_components <- num(row_df$current_sofa_n_components)
  current_severity_col <- "current_sofa_surrogate"
  current_component_col <- "current_sofa_n_components"
  current_severity_label <- "current_day_available_component_SOFA_surrogate"
}

severity_specs <- list(
  admission_total_SOFA = list(
    col = "sofa_admission",
    label = if ("sofa_source" %in% names(row_df)) {
      z <- unique(na.omit(as.character(row_df$sofa_source)))
      if (length(z)) paste0("admission_SOFA:", z[1]) else "admission_total_SOFA"
    } else "admission_total_SOFA",
    ncomp_col = NA_character_
  )
)
if (!is.na(current_severity_col)) {
  severity_specs$current_day_available_component_SOFA <- list(
    col = current_severity_col,
    label = current_severity_label,
    ncomp_col = current_component_col
  )
}

all_repeat <- list()
all_pooled <- list()
all_boot <- list()
all_pred <- list()
run_index <- 0L

for (assignment in ASSIGNMENTS) {
  phenotype_col <- if (assignment == "smoothing") "current_subphenotype" else "filtering_subphenotype"
  for (severity_name in names(severity_specs)) {
    spec <- severity_specs[[severity_name]]
    for (outcome in OUTCOMES) {
      for (d in SELECTED_DAYS) {
        dat <- row_df %>%
          filter(day == d) %>%
          mutate(
            phenotype_use = clean_subphenotype(.data[[phenotype_col]]),
            severity_use = num(.data[[spec$col]]),
            outcome_use = as.integer(.data[[outcome]])
          )

        if (!is.na(spec$ncomp_col) && spec$ncomp_col %in% names(dat)) {
          ncomp <- num(dat[[spec$ncomp_col]])
          dat <- dat[is.na(ncomp) | ncomp >= MIN_COMPONENTS, , drop = FALSE]
        }

        dat <- dat %>%
          filter(
            is.finite(age), !is.na(sex), is.finite(severity_use),
            phenotype_use %in% SUBPHENOTYPE_LEVELS,
            outcome_use %in% c(0L, 1L)
          ) %>%
          mutate(
            severity_for_model = severity_use,
            phenotype_for_model = phenotype_use,
            outcome_for_model = outcome_use
          )

        n_events <- sum(dat$outcome_for_model == 1L)
        n_nonevents <- sum(dat$outcome_for_model == 0L)
        n_levels <- n_distinct(dat$phenotype_for_model)
        if (n_events < MIN_EVENTS || n_nonevents < MIN_NONEVENTS || n_levels < 2L) {
          message(sprintf(
            "[OOF] Skip assignment=%s severity=%s outcome=%s day=%d (N=%d, events=%d, nonevents=%d, phenotypes=%d)",
            assignment, severity_name, outcome, d, nrow(dat), n_events, n_nonevents, n_levels
          ))
          next
        }

        run_index <- run_index + 1L
        message(sprintf(
          "[OOF] Run %d: assignment=%s severity=%s outcome=%s day=%d N=%d",
          run_index, assignment, severity_name, outcome, d, nrow(dat)
        ))

        cv_dat <- dat %>%
          transmute(
            stay_id, day, age, sex,
            severity_for_model,
            phenotype_for_model,
            outcome_for_model
          )
        names(cv_dat)[names(cv_dat) == "severity_for_model"] <- "severity_model"
        names(cv_dat)[names(cv_dat) == "phenotype_for_model"] <- "phenotype_model"
        names(cv_dat)[names(cv_dat) == "outcome_for_model"] <- outcome

        cv <- run_repeated_cv(
          cv_dat,
          outcome_col = outcome,
          severity_col = "severity_model",
          phenotype_col = "phenotype_model",
          nfolds = N_FOLDS,
          nrepeats = N_REPEATS,
          seed = CV_SEED + run_index * 100003L
        )

        apparent_x1 <- make_design(cv_dat, "severity_model", "phenotype_model", FALSE)
        apparent_x2 <- make_design(cv_dat, "severity_model", "phenotype_model", TRUE)
        y_now <- as.integer(cv_dat[[outcome]])
        apparent_p1 <- fit_predict_glm(apparent_x1, y_now, apparent_x1)
        apparent_p2 <- fit_predict_glm(apparent_x2, y_now, apparent_x2)
        apparent <- metric_row(y_now, apparent_p1, apparent_p2)

        boot <- bootstrap_paired_metrics(
          cv$prediction,
          nboot = BOOTSTRAP_N,
          seed = BOOTSTRAP_SEED + run_index * 1009L
        )

        meta <- tibble(
          assignment = assignment,
          severity_spec = severity_name,
          severity_label = spec$label,
          outcome = outcome,
          day = d,
          interpretation_window = ifelse(d %in% PRIMARY_DAYS, "primary_day_1_14", "exploratory_after_day_14")
        )

        all_repeat[[run_index]] <- bind_cols(meta[rep(1, nrow(cv$repeat_metrics)), ], cv$repeat_metrics)

        pooled <- bind_cols(meta, cv$pooled_metrics) %>%
          mutate(
            apparent_auc_M1 = apparent$auc_M1,
            apparent_auc_M2 = apparent$auc_M2,
            apparent_delta_auc_M2_minus_M1 = apparent$delta_auc_M2_minus_M1,
            apparent_brier_M1 = apparent$brier_M1,
            apparent_brier_M2 = apparent$brier_M2,
            apparent_delta_brier_M2_minus_M1 = apparent$delta_brier_M2_minus_M1,
            delta_auc_ci_low = if (nrow(boot)) q025(boot$delta_auc_M2_minus_M1) else NA_real_,
            delta_auc_ci_high = if (nrow(boot)) q975(boot$delta_auc_M2_minus_M1) else NA_real_,
            delta_brier_ci_low = if (nrow(boot)) q025(boot$delta_brier_M2_minus_M1) else NA_real_,
            delta_brier_ci_high = if (nrow(boot)) q975(boot$delta_brier_M2_minus_M1) else NA_real_,
            bootstrap_n_valid = if (nrow(boot)) sum(is.finite(boot$delta_auc_M2_minus_M1)) else 0L,
            folds = N_FOLDS,
            repeats = N_REPEATS
          )
        all_pooled[[run_index]] <- pooled

        if (nrow(boot)) all_boot[[run_index]] <- bind_cols(meta[rep(1, nrow(boot)), ], boot)
        all_pred[[run_index]] <- cv$prediction %>%
          mutate(
            assignment = assignment,
            severity_spec = severity_name,
            severity_label = spec$label,
            outcome = outcome,
            analysis_day = d
          )
      }
    }
  }
}

repeat_tbl <- bind_rows(all_repeat)
summary_tbl <- bind_rows(all_pooled)
bootstrap_tbl <- bind_rows(all_boot)
prediction_tbl <- bind_rows(all_pred)

.stopif(nrow(summary_tbl) > 0L, "No out-of-fold validation model could be evaluated.")

summary_tbl <- summary_tbl %>%
  mutate(
    delta_auc_95ci = ifelse(
      is.finite(delta_auc_ci_low) & is.finite(delta_auc_ci_high),
      sprintf("%.3f (%.3f to %.3f)", delta_auc_M2_minus_M1, delta_auc_ci_low, delta_auc_ci_high),
      sprintf("%.3f", delta_auc_M2_minus_M1)
    ),
    delta_brier_95ci = ifelse(
      is.finite(delta_brier_ci_low) & is.finite(delta_brier_ci_high),
      sprintf("%.4f (%.4f to %.4f)", delta_brier_M2_minus_M1, delta_brier_ci_low, delta_brier_ci_high),
      sprintf("%.4f", delta_brier_M2_minus_M1)
    )
  ) %>%
  arrange(severity_spec, assignment, factor(outcome, levels = OUTCOMES), day)

write_both(repeat_tbl, paste0(OUT_PREFIX, "_repeat_metrics"))
write_both(summary_tbl, paste0(OUT_PREFIX, "_summary"))
write_both(bootstrap_tbl, paste0(OUT_PREFIX, "_bootstrap_metrics"))
write_both(prediction_tbl, paste0(OUT_PREFIX, "_oof_predictions"))

key_day1 <- summary_tbl %>%
  filter(day == 1L) %>%
  select(
    assignment, severity_spec, severity_label, outcome, day,
    n, events, event_rate,
    auc_M1, auc_M2, delta_auc_M2_minus_M1,
    delta_auc_ci_low, delta_auc_ci_high, delta_auc_95ci,
    brier_M1, brier_M2, delta_brier_M2_minus_M1,
    delta_brier_ci_low, delta_brier_ci_high, delta_brier_95ci,
    calibration_intercept_M1, calibration_slope_M1,
    calibration_intercept_M2, calibration_slope_M2,
    apparent_delta_auc_M2_minus_M1,
    folds, repeats, bootstrap_n_valid
  )
write_both(key_day1, paste0(OUT_PREFIX, "_day1_key_results"))

validation_table <- summary_tbl %>%
  filter(day %in% PRIMARY_DAYS) %>%
  transmute(
    assignment,
    severity_adjustment = severity_label,
    outcome,
    ICU_day = day,
    active_N = n,
    event_N = events,
    AUROC_Model_1 = auc_M1,
    AUROC_Model_2_plus_HMM = auc_M2,
    delta_AUROC = delta_auc_M2_minus_M1,
    delta_AUROC_95CI_low = delta_auc_ci_low,
    delta_AUROC_95CI_high = delta_auc_ci_high,
    delta_AUROC_with_95CI = delta_auc_95ci,
    Brier_Model_1 = brier_M1,
    Brier_Model_2_plus_HMM = brier_M2,
    delta_Brier = delta_brier_M2_minus_M1,
    calibration_intercept_Model_2 = calibration_intercept_M2,
    calibration_slope_Model_2 = calibration_slope_M2,
    validation_scope = "Repeated patient-level CV conditional on the fixed HMM representation"
  )
write_both(validation_table, paste0(OUT_PREFIX, "_table"))

plot_df <- summary_tbl %>%
  filter(day %in% SELECTED_DAYS) %>%
  mutate(
    outcome_label = recode(
      outcome,
      death_3d = "Death within 3 days",
      death_7d = "Death within 7 days",
      longstay21 = "ICU LOS >=21 days"
    ),
    assignment_label = recode(
      assignment,
      smoothing = "Smoothed assignment",
      filtering = "Filtering assignment"
    ),
    severity_label_plot = recode(
      severity_spec,
      admission_total_SOFA = "Admission SOFA",
      current_day_available_component_SOFA = "Current-day available-component SOFA"
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x = day,
    y = delta_auc_M2_minus_M1,
    ymin = delta_auc_ci_low,
    ymax = delta_auc_ci_high,
    color = assignment_label,
    linetype = assignment_label,
    group = assignment_label
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_ribbon(aes(fill = assignment_label), alpha = 0.10, color = NA, show.legend = FALSE) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  facet_grid(severity_label_plot ~ outcome_label, scales = "free_y") +
  scale_x_continuous(breaks = SELECTED_DAYS) +
  labs(
    x = "ICU day",
    y = "Cross-validated AUROC increment (Model 2 - Model 1)",
    color = NULL,
    linetype = NULL,
    title = "Patient-level out-of-fold validation"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(size = 9)
  )

ggsave(
  paste0(OUT_PREFIX, "_auc_increment.pdf"),
  p, width = 10.5, height = 6.5, device = cairo_pdf, bg = "white"
)
ggsave(
  paste0(OUT_PREFIX, "_auc_increment.png"),
  p, width = 10.5, height = 6.5, dpi = 600, bg = "white"
)

manifest <- tibble::tribble(
  ~output, ~purpose,
  paste0(OUT_PREFIX, "_summary.csv"), "Main pooled out-of-fold metrics and paired-bootstrap confidence intervals.",
  paste0(OUT_PREFIX, "_repeat_metrics.csv"), "Metrics for each repeated patient-level cross-validation split.",
  paste0(OUT_PREFIX, "_bootstrap_metrics.csv"), "Paired patient-level bootstrap replicates of pooled out-of-fold predictions.",
  paste0(OUT_PREFIX, "_oof_predictions.csv"), "Averaged out-of-fold predictions for audit and reproducibility.",
  paste0(OUT_PREFIX, "_day1_key_results.csv"), "Day-1 key estimates.",
  paste0(OUT_PREFIX, "_table.csv"), "Conditional out-of-fold validation table for ICU days 1-14.",
  paste0(OUT_PREFIX, "_auc_increment.pdf"), "Cross-validated AUROC increments for smoothed and filtering assignments."
) %>%
  mutate(
    input_rowlevel = row_fp,
    input_filtering = filter_fp,
    input_currentday_sofa = ifelse(is.na(sofa_fp), "not used", sofa_fp),
    folds = N_FOLDS,
    repeats = N_REPEATS,
    cv_seed = CV_SEED,
    bootstrap_n = BOOTSTRAP_N,
    bootstrap_seed = BOOTSTRAP_SEED,
    validation_scope = "Outcome-model validation conditional on the fixed HMM representation"
  )
write_both(manifest, paste0(OUT_PREFIX, "_manifest"))

message("[OOF validation] Done. Main output: ", OUT_PREFIX, "_summary.csv")
message("[OOF validation] Table: ", OUT_PREFIX, "_table.csv")
message("[OOF validation] Figure: ", OUT_PREFIX, "_auc_increment.pdf")
