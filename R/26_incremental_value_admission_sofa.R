
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
inner_join <- dplyr::inner_join
bind_rows <- dplyr::bind_rows
all_of    <- tidyselect::all_of
any_of    <- tidyselect::any_of

OUT_PREFIX <- "incremental_value_admission_sofa"
ANNO_CANDIDATES <- c(
  MIMIC_ACTIVE_ANNOT_MASTER,
  MIMIC_ACTIVE_ANNOT_MASTER
)
SOFA_FILE_CANDIDATES <- c(
  "mimic_iv_sepsis3_cohort_final.csv.gz",
  "mimic_iv_sepsis3_cohort_final.csv"
)
SELECTED_DAYS <- c(1L, 3L, 5L, 7L, 10L, 14L, 21L, 28L)
LONGSTAY_DAY <- 21L
MIN_EVENTS <- 10L
MIN_NONEVENTS <- 10L

`%||%` <- function(a, b) if (!is.null(a)) a else b
num <- function(x) suppressWarnings(as.numeric(x))
.stopif <- function(ok, msg) if (!ok) stop(msg, call. = FALSE)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit) || !nzchar(hit)) return(NA_character_)
  hit
}

.pick_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  rep(NA, nrow(df))
}

.pick_col_name <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

write_both <- function(df, stem) {
  readr::write_csv(df, paste0(stem, ".csv"))
  readr::write_tsv(df, paste0(stem, ".tsv"))
  invisible(df)
}

state_id_from_name <- function(x) as.integer(stringr::str_remove(as.character(x), "^S"))

preferred_subphenotype_levels <- function(x = character(0)) {
  preferred <- c(
    "A1 low burden",
    "A2 mild burden",
    "B1 renal mild",
    "B2 renal severe",
    "G respiratory",
    "D1 lactate/shock",
    "D2 hepatic",
    "D3 shock/coagulation-overlap",
    "Other/NA"
  )
  x <- x[!is.na(x) & nzchar(x)]
  c(intersect(preferred, unique(x)), setdiff(unique(x), preferred))
}

make_safe_level <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "Other/NA"
  x
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

extract_model_metrics <- function(fit, model_id, outcome, day, n, events, sofa_source) {
  if (is.null(fit)) {
    return(tibble(
      outcome = outcome, day = day, model_id = model_id, n = n, events = events,
      event_rate = ifelse(n > 0, events / n, NA_real_), aic = NA_real_, bic = NA_real_,
      logLik = NA_real_, auc = NA_real_, sofa_source = sofa_source
    ))
  }
  pred <- tryCatch(stats::predict(fit, type = "response"), error = function(e) rep(NA_real_, nrow(model.frame(fit))))
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
    auc = auc_rank(y, pred),
    sofa_source = sofa_source
  )
}

lrt_safe <- function(reduced, full, outcome, day, comparison, sofa_source) {
  if (is.null(reduced) || is.null(full)) {
    return(tibble(outcome = outcome, day = day, comparison = comparison,
                  df = NA_real_, chisq = NA_real_, p_value = NA_real_, sofa_source = sofa_source))
  }
  out <- tryCatch(stats::anova(reduced, full, test = "Chisq"), error = function(e) NULL)
  if (is.null(out) || nrow(out) < 2) {
    return(tibble(outcome = outcome, day = day, comparison = comparison,
                  df = NA_real_, chisq = NA_real_, p_value = NA_real_, sofa_source = sofa_source))
  }
  tibble(
    outcome = outcome,
    day = day,
    comparison = comparison,
    df = suppressWarnings(as.numeric(out$Df[2])),
    chisq = suppressWarnings(as.numeric(out$Deviance[2])),
    p_value = suppressWarnings(as.numeric(out$`Pr(>Chi)`[2])),
    sofa_source = sofa_source
  )
}

coef_table_m2 <- function(fit, outcome, day, n, events, sofa_source) {
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
      sofa_source = sofa_source
    ) %>%
    select(outcome, day, term, estimate, std_error, z_value, p_value,
           odds_ratio, ci_low, ci_high, n, events, sofa_source)
}

needed <- c("post_smooth_full", "combined", "ids", "icustays_intubated_first_patient_admission")
miss <- needed[!vapply(needed, exists, logical(1), inherits = TRUE)]
.stopif(length(miss) == 0, paste("Missing required objects:", paste(miss, collapse = ", ")))

post_smooth_full <- as.matrix(get("post_smooth_full", envir = .GlobalEnv))
combined <- get("combined", envir = .GlobalEnv)
ids <- as.character(get("ids", envir = .GlobalEnv))
icu <- get("icustays_intubated_first_patient_admission", envir = .GlobalEnv) %>%
  as.data.frame() %>%
  mutate(stay_id = as.character(stay_id))

.stopif(nrow(post_smooth_full) == length(ids), "nrow(post_smooth_full) must equal length(ids).")
.stopif(nrow(post_smooth_full) == nrow(combined), "nrow(post_smooth_full) must equal nrow(combined).")
.stopif("lab_time" %in% names(combined), "combined must contain lab_time.")

anno_file <- first_existing(ANNO_CANDIDATES)
.stopif(!is.na(anno_file), paste("No annotation file found among:", paste(ANNO_CANDIDATES, collapse = ", ")))
anno_tbl <- readr::read_csv(anno_file, show_col_types = FALSE)
.stopif(all(c("state_id", "subphenotype") %in% names(anno_tbl)), "Annotation file must contain state_id and subphenotype.")

state_names <- colnames(post_smooth_full)
if (is.null(state_names) || any(is.na(state_names)) || any(!nzchar(state_names))) {
  state_names <- paste0("S", seq_len(ncol(post_smooth_full)))
  colnames(post_smooth_full) <- state_names
}
state_map <- tibble(state = state_names, state_id = state_id_from_name(state_names)) %>%
  left_join(anno_tbl %>% select(state_id, subphenotype), by = "state_id") %>%
  mutate(subphenotype = make_safe_level(subphenotype))
subph_levels <- preferred_subphenotype_levels(state_map$subphenotype)

subph_prob <- matrix(0, nrow = nrow(post_smooth_full), ncol = length(subph_levels),
                     dimnames = list(NULL, subph_levels))
for (sp in subph_levels) {
  st <- state_map$state[state_map$subphenotype == sp]
  st <- intersect(st, colnames(post_smooth_full))
  if (length(st) > 0) subph_prob[, sp] <- rowSums(post_smooth_full[, st, drop = FALSE], na.rm = TRUE)
}
dom_subph <- colnames(subph_prob)[max.col(subph_prob, ties.method = "first")]

age_vec_by_stay <- num(.pick_col(icu, c("age", "Age", "anchor_age")))
sex_raw_by_stay <- .pick_col(icu, c("sex", "gender", "Gender", "sex_2"))
sex_factor_by_stay <- as.factor(as.character(sex_raw_by_stay))

intime <- .pick_col(icu, c("intime", "icu_intime", "ICU_intime", "icu_in_time", "in_time", "admittime", "admit_time"))
outtime <- .pick_col(icu, c("outtime", "icu_outtime", "ICU_outtime", "out_time", "dischtime", "discharge_time"))
dod <- .pick_col(icu, c("dod", "death_time", "deathtime", "DOD"))

icu_time_tbl <- tibble(
  stay_id = as.character(icu$stay_id),
  age = age_vec_by_stay,
  sex = sex_factor_by_stay,
  intime = as.POSIXct(intime),
  outtime = as.POSIXct(outtime),
  dod = as.POSIXct(dod)
) %>%
  mutate(
    death_days = ifelse(!is.na(dod) & !is.na(intime), as.numeric(difftime(dod, intime, units = "days")), NA_real_),
    death_days = ifelse(!is.na(death_days) & death_days < 0, NA_real_, death_days),
    los_days = ifelse(!is.na(outtime) & !is.na(intime), as.numeric(difftime(outtime, intime, units = "days")), NA_real_),
    longstay21 = as.integer(!is.na(los_days) & los_days >= LONGSTAY_DAY)
  )

find_admission_sofa <- function(icu_df, combined_df, ids_vec) {
  sofa_candidates <- c(
    "admission_sofa", "sofa_admission", "SOFA_admission", "sofa_day1", "SOFA_day1",
    "first_sofa", "SOFA_first", "max_sofa_24h", "sofa_24hours", "SOFA_24hours",
    "sofa", "SOFA", "sofa_score", "SOFA_score", "total_sofa", "Total_SOFA"
  )

  nm <- .pick_col_name(icu_df, sofa_candidates)
  if (!is.na(nm)) {
    v <- num(icu_df[[nm]])
    if (sum(is.finite(v)) > 10 && length(unique(v[is.finite(v)])) > 1) {
      return(tibble(stay_id = as.character(icu_df$stay_id), sofa_admission = v,
                    sofa_source = paste0("icustays_intubated_first_patient_admission$", nm)))
    }
  }

  fp <- first_existing(SOFA_FILE_CANDIDATES)
  if (!is.na(fp)) {
    message("Trying to read SOFA from file: ", fp)
    sofa_file <- tryCatch(readr::read_csv(fp, show_col_types = FALSE), error = function(e) NULL)
    if (!is.null(sofa_file) && "stay_id" %in% names(sofa_file)) {
      nm2 <- .pick_col_name(sofa_file, sofa_candidates)
      if (!is.na(nm2)) {
        tmp <- sofa_file %>%
          transmute(stay_id = as.character(stay_id), sofa_admission = num(.data[[nm2]])) %>%
          group_by(stay_id) %>%
          summarise(sofa_admission = suppressWarnings(min(sofa_admission, na.rm = TRUE)), .groups = "drop") %>%
          mutate(sofa_admission = ifelse(is.infinite(sofa_admission), NA_real_, sofa_admission),
                 sofa_source = paste0(fp, "$", nm2))
        if (sum(is.finite(tmp$sofa_admission)) > 10 && length(unique(tmp$sofa_admission[is.finite(tmp$sofa_admission)])) > 1) {
          return(tmp)
        }
      }
    }
  }

  day_vec <- suppressWarnings(as.integer(combined_df$lab_time))
  df1 <- as.data.frame(combined_df[day_vec == 1L, , drop = FALSE])
  id1 <- as.character(ids_vec[day_vec == 1L])
  if (nrow(df1) > 0) {
    pick_num <- function(cands) {
      nm3 <- .pick_col_name(df1, cands)
      if (is.na(nm3)) return(rep(NA_real_, nrow(df1)))
      num(df1[[nm3]])
    }
    PF <- pick_num(c("PF", "PaO2_FiO2", "pO2_FO2I"))
    plt <- pick_num(c("PlateletCount", "Platelets", "platelets"))
    bili <- pick_num(c("TotalBilirubin", "bilirubin", "Bilirubin"))
    cr <- pick_num(c("Creatinine", "creatinine"))
    map <- pick_num(c("MAP", "MeanBP", "MEAN"))

    plausible <- sum(is.finite(PF) & PF > 20 & PF < 800) > 10 ||
      sum(is.finite(plt) & plt > 1 & plt < 1000) > 10 ||
      sum(is.finite(cr) & cr > 0.1 & cr < 20) > 10

    if (plausible) {
      sofa_resp <- case_when(!is.finite(PF) ~ NA_real_, PF < 100 ~ 3, PF < 200 ~ 2, PF < 300 ~ 1, TRUE ~ 0)
      sofa_coag <- case_when(!is.finite(plt) ~ NA_real_, plt < 20 ~ 4, plt < 50 ~ 3, plt < 100 ~ 2, plt < 150 ~ 1, TRUE ~ 0)
      sofa_liver <- case_when(!is.finite(bili) ~ NA_real_, bili >= 12 ~ 4, bili >= 6 ~ 3, bili >= 2 ~ 2, bili >= 1.2 ~ 1, TRUE ~ 0)
      sofa_renal <- case_when(!is.finite(cr) ~ NA_real_, cr >= 5 ~ 4, cr >= 3.5 ~ 3, cr >= 2 ~ 2, cr >= 1.2 ~ 1, TRUE ~ 0)
      sofa_cv <- case_when(!is.finite(map) ~ NA_real_, map < 70 ~ 1, TRUE ~ 0)

      tmp <- tibble(
        stay_id = id1,
        sofa_resp = sofa_resp,
        sofa_coag = sofa_coag,
        sofa_liver = sofa_liver,
        sofa_renal = sofa_renal,
        sofa_cv = sofa_cv
      ) %>%
        mutate(
          n_components = rowSums(!is.na(select(., sofa_resp, sofa_coag, sofa_liver, sofa_renal, sofa_cv))),
          sofa_admission = rowSums(select(., sofa_resp, sofa_coag, sofa_liver, sofa_renal, sofa_cv), na.rm = TRUE),
          sofa_admission = ifelse(n_components == 0, NA_real_, sofa_admission)
        ) %>%
        group_by(stay_id) %>%
        summarise(sofa_admission = suppressWarnings(max(sofa_admission, na.rm = TRUE)),
                  n_components = max(n_components, na.rm = TRUE), .groups = "drop") %>%
        mutate(sofa_admission = ifelse(is.infinite(sofa_admission), NA_real_, sofa_admission),
               sofa_source = "day1_available_component_sofa_surrogate")
      if (sum(is.finite(tmp$sofa_admission)) > 10 && length(unique(tmp$sofa_admission[is.finite(tmp$sofa_admission)])) > 1) {
        warning("Using an available-component day-1 SOFA surrogate because a recorded admission SOFA variable was not found.")
        return(tmp %>% select(stay_id, sofa_admission, sofa_source))
      }
    }
  }

  tibble(stay_id = as.character(icu_df$stay_id), sofa_admission = NA_real_, sofa_source = "not_found")
}

sofa_tbl <- find_admission_sofa(icu, combined, ids)
sofa_source_use <- sofa_tbl$sofa_source[which(!is.na(sofa_tbl$sofa_source) & sofa_tbl$sofa_source != "not_found")][1] %||% "not_found"
if (is.na(sofa_source_use)) sofa_source_use <- "not_found"

if (identical(sofa_source_use, "not_found")) {
  stop("Admission SOFA could not be found or constructed. Add a SOFA column to icustays_intubated_first_patient_admission or provide a Sepsis-3 cohort file with stay_id and SOFA.", call. = FALSE)
}
message("Using admission SOFA source: ", sofa_source_use)

analysis_df <- tibble(
  row_id = seq_along(ids),
  stay_id = ids,
  day = suppressWarnings(as.integer(combined$lab_time)),
  current_subphenotype = dom_subph
) %>%
  left_join(icu_time_tbl %>% select(stay_id, age, sex, death_days, los_days, longstay21), by = "stay_id") %>%
  left_join(sofa_tbl %>% select(stay_id, sofa_admission), by = "stay_id") %>%
  mutate(
    death_3d = as.integer(!is.na(death_days) & death_days <= (day + 2L)),
    death_7d = as.integer(!is.na(death_days) & death_days <= (day + 6L)),
    current_subphenotype = factor(current_subphenotype, levels = subph_levels),
    sex = factor(sex),
    sofa_admission = num(sofa_admission)
  )

if ("A1 low burden" %in% levels(analysis_df$current_subphenotype)) {
  analysis_df$current_subphenotype <- stats::relevel(analysis_df$current_subphenotype, ref = "A1 low burden")
} else {
  ref_sp <- analysis_df %>% count(current_subphenotype, sort = TRUE) %>% filter(!is.na(current_subphenotype)) %>% slice(1) %>% pull(current_subphenotype) %>% as.character()
  if (length(ref_sp) == 1 && nzchar(ref_sp)) analysis_df$current_subphenotype <- stats::relevel(analysis_df$current_subphenotype, ref = ref_sp)
}

write_both(
  analysis_df %>% filter(day %in% SELECTED_DAYS) %>% select(row_id, stay_id, day, age, sex, sofa_admission, current_subphenotype, death_3d, death_7d, longstay21) %>%
    mutate(sofa_source = sofa_source_use),
  paste0(OUT_PREFIX, "_rowlevel_selected_days")
)

outcomes <- tibble(
  outcome = c("death_3d", "death_7d", "longstay21"),
  label = c("3-day death", "7-day death", paste0("ICU long stay >= ", LONGSTAY_DAY, " days"))
)

metrics_list <- list()
lrt_list <- list()
coef_list <- list()
model_index <- 0L

for (d in SELECTED_DAYS) {
  for (outcome_nm in outcomes$outcome) {

    if (outcome_nm == "longstay21" && d >= LONGSTAY_DAY) next

    dat <- analysis_df %>%
      filter(day == d) %>%
      select(age, sex, sofa_admission, current_subphenotype, all_of(outcome_nm)) %>%
      rename(y = all_of(outcome_nm)) %>%
      filter(!is.na(y), !is.na(age), !is.na(sex), !is.na(sofa_admission), !is.na(current_subphenotype)) %>%
      mutate(y = as.integer(y)) %>%
      droplevels()

    n <- nrow(dat)
    events <- sum(dat$y == 1L, na.rm = TRUE)
    nonevents <- sum(dat$y == 0L, na.rm = TRUE)

    if (n < (MIN_EVENTS + MIN_NONEVENTS) || events < MIN_EVENTS || nonevents < MIN_NONEVENTS || nlevels(dat$current_subphenotype) < 2) {
      model_index <- model_index + 1L
      metrics_list[[model_index]] <- tibble(
        outcome = outcome_nm, day = d, model_id = c("M0_age_sex", "M1_age_sex_sofa", "M2_age_sex_sofa_hmm"),
        n = n, events = events, event_rate = ifelse(n > 0, events / n, NA_real_),
        aic = NA_real_, bic = NA_real_, logLik = NA_real_, auc = NA_real_, sofa_source = sofa_source_use,
        skipped_reason = sprintf("Insufficient usable rows/events or subphenotype levels: n=%d, events=%d, nonevents=%d, levels=%d", n, events, nonevents, nlevels(dat$current_subphenotype))
      )
      next
    }

    m0 <- fit_glm_safe(y ~ age + sex, dat)
    m1 <- fit_glm_safe(y ~ age + sex + sofa_admission, dat)
    m2 <- fit_glm_safe(y ~ age + sex + sofa_admission + current_subphenotype, dat)

    model_index <- model_index + 1L
    metrics_list[[model_index]] <- bind_rows(
      extract_model_metrics(m0, "M0_age_sex", outcome_nm, d, n, events, sofa_source_use),
      extract_model_metrics(m1, "M1_age_sex_sofa", outcome_nm, d, n, events, sofa_source_use),
      extract_model_metrics(m2, "M2_age_sex_sofa_hmm", outcome_nm, d, n, events, sofa_source_use)
    ) %>% mutate(skipped_reason = NA_character_)

    lrt_list[[length(lrt_list) + 1L]] <- bind_rows(
      lrt_safe(m0, m1, outcome_nm, d, "M1_vs_M0_add_SOFA", sofa_source_use),
      lrt_safe(m1, m2, outcome_nm, d, "M2_vs_M1_add_HMM_subphenotype", sofa_source_use)
    )

    coef_list[[length(coef_list) + 1L]] <- coef_table_m2(m2, outcome_nm, d, n, events, sofa_source_use)
  }
}

model_metrics <- bind_rows(metrics_list) %>%
  group_by(outcome, day) %>%
  mutate(
    aic_M1 = aic[model_id == "M1_age_sex_sofa"][1],
    bic_M1 = bic[model_id == "M1_age_sex_sofa"][1],
    auc_M1 = auc[model_id == "M1_age_sex_sofa"][1],
    delta_aic_vs_M1 = aic - aic_M1,
    delta_bic_vs_M1 = bic - bic_M1,
    delta_auc_vs_M1 = auc - auc_M1
  ) %>%
  ungroup() %>%
  select(-aic_M1, -bic_M1, -auc_M1)

lrt_summary <- bind_rows(lrt_list)
hmm_terms <- bind_rows(coef_list)

hmm_lrt <- lrt_summary %>%
  filter(comparison == "M2_vs_M1_add_HMM_subphenotype")

m1m2_wide <- model_metrics %>%
  filter(model_id %in% c("M1_age_sex_sofa", "M2_age_sex_sofa_hmm")) %>%
  select(outcome, day, model_id, n, events, event_rate, aic, bic, auc, sofa_source) %>%
  pivot_wider(
    names_from = model_id,
    values_from = c(aic, bic, auc),
    names_sep = "__"
  ) %>%
  mutate(
    delta_aic_M2_minus_M1 = aic__M2_age_sex_sofa_hmm - aic__M1_age_sex_sofa,
    delta_bic_M2_minus_M1 = bic__M2_age_sex_sofa_hmm - bic__M1_age_sex_sofa,
    delta_auc_M2_minus_M1 = auc__M2_age_sex_sofa_hmm - auc__M1_age_sex_sofa
  )

reporting_summary <- m1m2_wide %>%
  left_join(hmm_lrt %>% select(outcome, day, hmm_lrt_chisq = chisq, hmm_lrt_df = df, hmm_lrt_p = p_value), by = c("outcome", "day")) %>%
  arrange(outcome, day)

target_terms <- hmm_terms %>%
  filter(str_detect(term, "^current_subphenotype")) %>%
  mutate(
    subphenotype_term = str_replace(term, "^current_subphenotype", ""),
    target_signal = case_when(
      str_detect(subphenotype_term, "^D1\\b") ~ "D1 lactate/shock",
      str_detect(subphenotype_term, "^B2\\b") ~ "B2 renal severe",
      str_detect(subphenotype_term, "^G\\b") ~ "G respiratory",
      str_detect(subphenotype_term, "^D3\\b") ~ "D3 shock/coagulation-overlap",
      TRUE ~ "Other HMM subphenotype"
    )
  ) %>%
  arrange(outcome, day, p_value)

write_both(model_metrics, paste0(OUT_PREFIX, "_model_metrics"))
write_both(lrt_summary, paste0(OUT_PREFIX, "_lrt_summary"))
write_both(hmm_terms, paste0(OUT_PREFIX, "_model2_coefficients"))
write_both(target_terms, paste0(OUT_PREFIX, "_target_hmm_terms"))
write_both(reporting_summary, paste0(OUT_PREFIX, "_reporting_summary"))

manifest <- tibble(
  file = c(
    paste0(OUT_PREFIX, "_rowlevel_selected_days.csv"),
    paste0(OUT_PREFIX, "_model_metrics.csv"),
    paste0(OUT_PREFIX, "_lrt_summary.csv"),
    paste0(OUT_PREFIX, "_model2_coefficients.csv"),
    paste0(OUT_PREFIX, "_target_hmm_terms.csv"),
    paste0(OUT_PREFIX, "_reporting_summary.csv"),
    paste0(OUT_PREFIX, "_delta_auc_plot.pdf"),
    paste0(OUT_PREFIX, "_hmm_lrt_plot.pdf")
  ),
  description = c(
    "Row-level selected-day analysis dataset with age, sex, admission SOFA, current HMM subphenotype, and outcomes.",
    "AIC, BIC, log-likelihood, and AUROC for Model 0, Model 1, and Model 2.",
    "Likelihood-ratio tests for adding SOFA and for adding HMM subphenotype to age/sex/SOFA.",
    "Model 2 logistic-regression coefficients and Wald odds ratios.",
    "HMM subphenotype coefficient summary.",
    "Summary of Model 1 versus Model 2 incremental value.",
    "Plot of AUROC increment for adding HMM subphenotype to age/sex/SOFA.",
    "Plot of likelihood-ratio evidence for adding HMM subphenotype."
  ),
  sofa_source = sofa_source_use,
  annotation_file = anno_file,
  longstay_day = LONGSTAY_DAY
)
write_both(manifest, paste0(OUT_PREFIX, "_manifest"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))

  plot_auc <- reporting_summary %>%
    mutate(outcome = factor(outcome, levels = c("death_3d", "death_7d", "longstay21"))) %>%
    filter(is.finite(delta_auc_M2_minus_M1))

  if (nrow(plot_auc) > 0) {
    p_auc <- ggplot(plot_auc, aes(x = day, y = delta_auc_M2_minus_M1, group = outcome)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_line() +
      geom_point(size = 2) +
      facet_wrap(~ outcome, scales = "free_y") +
      labs(
        title = "Incremental AUROC after adding current HMM subphenotype",
        subtitle = "Model 2: age + sex + admission SOFA + current HMM subphenotype; Model 1: age + sex + admission SOFA",
        x = "ICU day",
        y = "AUROC difference: Model 2 - Model 1"
      ) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())
    ggsave(paste0(OUT_PREFIX, "_delta_auc_plot.pdf"), p_auc, width = 8, height = 5, bg = "white")
    ggsave(paste0(OUT_PREFIX, "_delta_auc_plot.png"), p_auc, width = 8, height = 5, dpi = 300, bg = "white")
  }

  plot_lrt <- hmm_lrt %>%
    mutate(
      minus_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
      outcome = factor(outcome, levels = c("death_3d", "death_7d", "longstay21"))
    ) %>%
    filter(is.finite(minus_log10_p))

  if (nrow(plot_lrt) > 0) {
    p_lrt <- ggplot(plot_lrt, aes(x = day, y = minus_log10_p, group = outcome)) +
      geom_hline(yintercept = -log10(0.05), linewidth = 0.3, linetype = "dashed") +
      geom_line() +
      geom_point(size = 2) +
      facet_wrap(~ outcome, scales = "free_y") +
      labs(
        title = "Likelihood-ratio evidence for adding current HMM subphenotype",
        subtitle = "Comparison: age + sex + admission SOFA + HMM subphenotype versus age + sex + admission SOFA",
        x = "ICU day",
        y = "-log10(P value)"
      ) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())
    ggsave(paste0(OUT_PREFIX, "_hmm_lrt_plot.pdf"), p_lrt, width = 8, height = 5, bg = "white")
    ggsave(paste0(OUT_PREFIX, "_hmm_lrt_plot.png"), p_lrt, width = 8, height = 5, dpi = 300, bg = "white")
  }
}

message("Done: incremental-value analysis written with prefix: ", OUT_PREFIX)
message("SOFA source: ", sofa_source_use)
