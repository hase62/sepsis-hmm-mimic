
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
  library(purrr)
  library(stringr)
  library(ggplot2)
  library(depmixS4)
})

OUT_PREFIX <- "filtering_sensitivity"
ANNO_CANDIDATES <- c(
  MIMIC_ACTIVE_ANNOT_MASTER,
  MIMIC_ACTIVE_ANNOT_MASTER
)
SELECTED_DAYS <- c(1L, 3L, 5L, 7L, 10L, 14L, 21L, 28L)
LONGSTAY_DAY <- if (exists("LONGSTAY_DAY", inherits = TRUE)) get("LONGSTAY_DAY", inherits = TRUE) else 21L

`%||%` <- function(a, b) if (!is.null(a)) a else b
num <- function(x) suppressWarnings(as.numeric(x))

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

.stop_missing <- function(nms) {
  miss <- nms[!vapply(nms, exists, logical(1), inherits = TRUE)]
  if (length(miss) > 0) {
    stop("Missing required object(s): ", paste(miss, collapse = ", "), call. = FALSE)
  }
}

.as_prob_mat <- function(x, K, label = "posterior") {
  m <- as.matrix(x)
  storage.mode(m) <- "double"

  if (ncol(m) == K) {
    out <- m
  } else {
    s_cols <- grep("^(S|St)[0-9]+$", colnames(m), value = TRUE)
    if (length(s_cols) >= K) {
      out <- m[, s_cols[seq_len(K)], drop = FALSE]
    } else if (ncol(m) > K) {
      out <- m[, (ncol(m) - K + 1):ncol(m), drop = FALSE]
    } else {
      stop(label, " has fewer columns than the number of states (K=", K, ").", call. = FALSE)
    }
  }

  colnames(out) <- paste0("S", seq_len(K))
  out
}

.state_id_from_name <- function(x) as.integer(stringr::str_remove(as.character(x), "^S"))

.preferred_subphenotype_levels <- function(x = character(0)) {
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
  c(intersect(preferred, x), setdiff(unique(x), preferred))
}

.pick_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  rep(NA, nrow(df))
}

.weighted_rate <- function(flag, w) {
  flag <- num(flag)
  w <- num(w)
  ok <- is.finite(flag) & is.finite(w) & w >= 0
  if (!any(ok) || sum(w[ok], na.rm = TRUE) <= 0) return(NA_real_)
  sum(flag[ok] * w[ok], na.rm = TRUE) / sum(w[ok], na.rm = TRUE)
}

.cohen_kappa <- function(a, b) {
  keep <- !is.na(a) & !is.na(b)
  a <- as.character(a[keep])
  b <- as.character(b[keep])
  if (!length(a)) return(NA_real_)
  lv <- sort(unique(c(a, b)))
  tab <- table(factor(a, levels = lv), factor(b, levels = lv))
  n <- sum(tab)
  if (n == 0) return(NA_real_)
  po <- sum(diag(tab)) / n
  pe <- sum(rowSums(tab) * colSums(tab)) / (n^2)
  if (!is.finite(pe) || pe >= 1) return(NA_real_)
  (po - pe) / (1 - pe)
}

.aggregate_prob_to_subphenotype <- function(post_mat, state_map, subphenotype_levels) {
  K <- ncol(post_mat)
  state_names <- colnames(post_mat) %||% paste0("S", seq_len(K))
  sm <- tibble(state = state_names, state_id = .state_id_from_name(state_names)) %>%
    left_join(state_map %>% select(state_id, subphenotype), by = "state_id") %>%
    mutate(subphenotype = if_else(is.na(subphenotype) | subphenotype == "", "Other/NA", subphenotype))

  out <- matrix(0, nrow = nrow(post_mat), ncol = length(subphenotype_levels))
  colnames(out) <- subphenotype_levels

  for (sp in subphenotype_levels) {
    cols <- which(sm$subphenotype == sp)
    if (length(cols) > 0) {
      out[, sp] <- rowSums(post_mat[, cols, drop = FALSE], na.rm = TRUE)
    }
  }
  out
}

.dominant_subphenotype <- function(prob_mat) {
  out <- colnames(prob_mat)[max.col(prob_mat, ties.method = "first")]
  out[rowSums(prob_mat, na.rm = TRUE) <= 0] <- NA_character_
  out
}

mimic_promote_final_model(expected_k = MIMIC_FINAL_STATE_DIMENSION, allow_rds = TRUE)

.stop_missing(c(
  "fitted_model", "response_formulas", "imputed_data_scale_df",
  "ids", "combined", "icustays_intubated_first_patient_admission"
))

fitted_model <- get("fitted_model", inherits = TRUE)
response_formulas <- get("response_formulas", inherits = TRUE)
imputed_data_scale_df <- get("imputed_data_scale_df", inherits = TRUE)
ids <- as.character(get("ids", inherits = TRUE))
combined <- get("combined", inherits = TRUE)
icu <- get("icustays_intubated_first_patient_admission", inherits = TRUE) %>%
  mutate(stay_id = as.character(stay_id))

K <- length(fitted_model@response)
if (K != as.integer(MIMIC_FINAL_STATE_DIMENSION)) stop("Expected final HMM K=", MIMIC_FINAL_STATE_DIMENSION, ", found K=", K, call. = FALSE)
state_names <- paste0("S", seq_len(K))

anno_file <- .first_existing(ANNO_CANDIDATES)
if (is.na(anno_file)) {
  stop("No annotation file found. Tried: ", paste(ANNO_CANDIDATES, collapse = ", "), call. = FALSE)
}

anno_tbl <- readr::read_csv(anno_file, show_col_types = FALSE) %>%
  mutate(state_id = as.integer(state_id))

if (!"subphenotype" %in% names(anno_tbl)) {
  stop("Annotation table must contain 'subphenotype'.", call. = FALSE)
}

state_map <- tibble(state = state_names, state_id = seq_len(K)) %>%
  left_join(anno_tbl %>% select(state_id, subphenotype), by = "state_id") %>%
  mutate(subphenotype = if_else(is.na(subphenotype) | subphenotype == "", "Other/NA", subphenotype))

subphenotype_levels <- .preferred_subphenotype_levels(state_map$subphenotype)

if (!exists("transition_rhs", inherits = TRUE)) {
  message("transition_rhs not found; using '~ age + sex' as fallback.")
  transition_rhs <- "~ age + sex"
} else {
  transition_rhs <- get("transition_rhs", inherits = TRUE)
}

rle_ids <- rle(ids)
if (length(unique(ids)) != length(rle_ids$lengths)) {
  warning("ids are not contiguous by stay. depmix ntimes may be invalid unless rows are already ordered by stay.")
}
ntimes_vec <- as.numeric(rle_ids$lengths)

mod_fixed <- depmix(
  response   = response_formulas,
  data       = imputed_data_scale_df,
  nstates    = K,
  family     = replicate(length(response_formulas), gaussian(), simplify = FALSE),
  transition = as.formula(transition_rhs),
  ntimes     = ntimes_vec
)
mod_fixed <- setpars(mod_fixed, getpars(fitted_model))

message("Computing filtering posterior from the fixed fitted HMM ...")
post_filtering_raw <- posterior(mod_fixed, type = "filtering")
post_filtering <- .as_prob_mat(post_filtering_raw, K, label = "filtering posterior")

if (exists("post_smooth_full", inherits = TRUE)) {
  post_smoothing <- .as_prob_mat(get("post_smooth_full", inherits = TRUE), K, label = "smoothing posterior")
} else {
  message("post_smooth_full not found; computing smoothing posterior from the fixed fitted HMM ...")
  post_smoothing_raw <- posterior(mod_fixed, type = "smoothing")
  post_smoothing <- .as_prob_mat(post_smoothing_raw, K, label = "smoothing posterior")
}

stopifnot(nrow(post_filtering) == length(ids))
stopifnot(nrow(post_smoothing) == length(ids))

saveRDS(post_filtering, paste0(OUT_PREFIX, "_post_filtering_stateprob.rds"))
saveRDS(post_smoothing, paste0(OUT_PREFIX, "_post_smoothing_stateprob_used.rds"))

prob_filter_sp <- .aggregate_prob_to_subphenotype(post_filtering, state_map, subphenotype_levels)
prob_smooth_sp <- .aggregate_prob_to_subphenotype(post_smoothing, state_map, subphenotype_levels)

saveRDS(prob_filter_sp, paste0(OUT_PREFIX, "_post_filtering_subphenotypeprob.rds"))
saveRDS(prob_smooth_sp, paste0(OUT_PREFIX, "_post_smoothing_subphenotypeprob_used.rds"))

lab_day <- num(combined$lab_time)
if (length(lab_day) != length(ids)) {
  stop("combined$lab_time length does not match ids length.", call. = FALSE)
}

dom_filter <- .dominant_subphenotype(prob_filter_sp)
dom_smooth <- .dominant_subphenotype(prob_smooth_sp)

row_df <- tibble(
  row_id = seq_along(ids),
  stay_id = ids,
  lab_time = as.integer(lab_day),
  dominant_filtering = dom_filter,
  dominant_smoothing = dom_smooth,
  same_dominant = dominant_filtering == dominant_smoothing
)

death_ids <- icu %>% filter(!is.na(dod)) %>% pull(stay_id)
death_flag <- as.integer(ids %in% death_ids)

icu_intime <- .pick_col(icu, c("intime", "icu_intime", "ICU_intime", "icu_in_time", "in_time", "admittime", "admit_time"))
outtime <- .pick_col(icu, c("outtime", "icu_outtime", "ICU_outtime", "out_time", "dischtime", "discharge_time"))

death_day_tbl <- icu %>%
  transmute(stay_id = as.character(stay_id), dod = dod, icu_intime = icu_intime) %>%
  mutate(
    death_days = ifelse(!is.na(dod) & !is.na(icu_intime), as.numeric(difftime(dod, icu_intime, units = "days")), NA_real_),
    death_days = ifelse(!is.na(death_days) & death_days < 0, NA_real_, death_days)
  ) %>%
  select(stay_id, death_days)

death_days_vec <- death_day_tbl$death_days[match(ids, death_day_tbl$stay_id)]
death3_flag <- as.integer(!is.na(death_days_vec) & death_days_vec <= (lab_day + 2L))
death7_flag <- as.integer(!is.na(death_days_vec) & death_days_vec <= (lab_day + 6L))

los_tbl <- icu %>%
  transmute(stay_id = as.character(stay_id), icu_intime = icu_intime, outtime = outtime) %>%
  mutate(los_days = ifelse(!is.na(outtime) & !is.na(icu_intime), as.numeric(difftime(outtime, icu_intime, units = "days")), NA_real_)) %>%
  select(stay_id, los_days)
los_days_vec <- los_tbl$los_days[match(ids, los_tbl$stay_id)]
longstay_flag <- as.integer(!is.na(los_days_vec) & los_days_vec >= LONGSTAY_DAY)

row_df <- row_df %>%
  mutate(
    death_eventual = death_flag,
    death_3d = death3_flag,
    death_7d = death7_flag,
    longstay_ge21d = longstay_flag
  )

.write_both(
  row_df %>% filter(lab_time %in% SELECTED_DAYS),
  paste0(OUT_PREFIX, "_dominant_rowlevel_selected_days")
)

agreement_overall <- tibble(
  n_rows = nrow(row_df),
  agreement = mean(row_df$same_dominant, na.rm = TRUE),
  kappa = .cohen_kappa(row_df$dominant_filtering, row_df$dominant_smoothing)
)
.write_both(agreement_overall, paste0(OUT_PREFIX, "_agreement_overall"))

agreement_by_day <- row_df %>%
  filter(!is.na(lab_time)) %>%
  group_by(day = lab_time) %>%
  summarise(
    n_rows = n(),
    n_stays = n_distinct(stay_id),
    agreement = mean(same_dominant, na.rm = TRUE),
    kappa = .cohen_kappa(dominant_filtering, dominant_smoothing),
    .groups = "drop"
  ) %>%
  arrange(day)
.write_both(agreement_by_day, paste0(OUT_PREFIX, "_agreement_by_day"))

confusion_overall <- row_df %>%
  count(dominant_smoothing, dominant_filtering, name = "n") %>%
  group_by(dominant_smoothing) %>%
  mutate(row_prop = n / sum(n)) %>%
  ungroup() %>%
  arrange(dominant_smoothing, desc(n))
.write_both(confusion_overall, paste0(OUT_PREFIX, "_confusion_overall"))

make_occ <- function(prob_sp, method_label) {
  purrr::map_dfr(SELECTED_DAYS, function(d) {
    mask <- !is.na(lab_day) & lab_day == d
    if (!any(mask)) {
      return(tibble(day = d, subphenotype = subphenotype_levels, occupancy = NA_real_, n_rows = 0L, n_stays = 0L, method = method_label))
    }
    tibble(
      day = d,
      subphenotype = subphenotype_levels,
      occupancy = colMeans(prob_sp[mask, subphenotype_levels, drop = FALSE], na.rm = TRUE),
      n_rows = sum(mask),
      n_stays = n_distinct(ids[mask]),
      method = method_label
    )
  })
}

occupancy_compare <- bind_rows(
  make_occ(prob_smooth_sp, "smoothing"),
  make_occ(prob_filter_sp, "filtering")
) %>%
  mutate(subphenotype = factor(subphenotype, levels = subphenotype_levels))

.write_both(occupancy_compare, paste0(OUT_PREFIX, "_occupancy_compare_selected_days"))

make_outcome <- function(prob_sp, method_label) {
  purrr::map_dfr(SELECTED_DAYS, function(d) {
    mask <- !is.na(lab_day) & lab_day == d
    if (!any(mask)) {
      return(tibble(
        day = d, subphenotype = subphenotype_levels,
        N_weighted = NA_real_, p_death_eventual = NA_real_, p_death_3d = NA_real_,
        p_death_7d = NA_real_, p_longstay_ge21d = NA_real_,
        n_rows = 0L, n_stays = 0L, method = method_label
      ))
    }

    purrr::map_dfr(subphenotype_levels, function(sp) {
      w <- prob_sp[mask, sp]
      tibble(
        day = d,
        subphenotype = sp,
        N_weighted = sum(w, na.rm = TRUE),
        p_death_eventual = .weighted_rate(death_flag[mask], w),
        p_death_3d = .weighted_rate(death3_flag[mask], w),
        p_death_7d = .weighted_rate(death7_flag[mask], w),
        p_longstay_ge21d = .weighted_rate(longstay_flag[mask], w),
        n_rows = sum(mask),
        n_stays = n_distinct(ids[mask]),
        method = method_label
      )
    })
  })
}

outcome_compare <- bind_rows(
  make_outcome(prob_smooth_sp, "smoothing"),
  make_outcome(prob_filter_sp, "filtering")
) %>%
  mutate(subphenotype = factor(subphenotype, levels = subphenotype_levels))

.write_both(outcome_compare, paste0(OUT_PREFIX, "_outcome_compare_selected_days"))
.write_both(
  outcome_compare %>% filter(method == "filtering"),
  paste0(OUT_PREFIX, "_outcome_filtering_selected_days")
)

outcome_delta <- outcome_compare %>%
  select(method, day, subphenotype, p_death_3d, p_death_7d, p_longstay_ge21d) %>%
  pivot_wider(names_from = method, values_from = c(p_death_3d, p_death_7d, p_longstay_ge21d)) %>%
  mutate(
    delta_p_death_3d_filter_minus_smooth = p_death_3d_filtering - p_death_3d_smoothing,
    delta_p_death_7d_filter_minus_smooth = p_death_7d_filtering - p_death_7d_smoothing,
    delta_p_longstay_filter_minus_smooth = p_longstay_ge21d_filtering - p_longstay_ge21d_smoothing
  )
.write_both(outcome_delta, paste0(OUT_PREFIX, "_outcome_delta_selected_days"))

p_agree <- agreement_by_day %>%
  filter(day <= max(SELECTED_DAYS, na.rm = TRUE)) %>%
  ggplot(aes(x = day, y = agreement)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_x_continuous(breaks = SELECTED_DAYS) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "ICU day",
    y = "Dominant subphenotype agreement",
    title = "Filtering vs smoothing posterior: dominant-subphenotype agreement",
    subtitle = "Filtering uses information up to the current day; smoothing uses the full retrospective trajectory."
  ) +
  theme_bw(base_size = 12)

ggsave(paste0(OUT_PREFIX, "_agreement_by_day.pdf"), p_agree, width = 7, height = 4.5, device = cairo_pdf)
ggsave(paste0(OUT_PREFIX, "_agreement_by_day.png"), p_agree, width = 7, height = 4.5, dpi = 600, bg = "white")

p_occ <- occupancy_compare %>%
  filter(day %in% SELECTED_DAYS) %>%
  ggplot(aes(x = factor(day), y = occupancy, group = method, linetype = method, shape = method)) +
  geom_line(aes(group = method), linewidth = 0.7) +
  geom_point(size = 1.5) +
  facet_wrap(~ subphenotype, scales = "free_y") +
  labs(
    x = "ICU day",
    y = "Mean posterior occupancy",
    title = "Filtering vs smoothing posterior: subphenotype occupancy"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUT_PREFIX, "_occupancy_compare_selected_days.pdf"), p_occ, width = 10, height = 7, device = cairo_pdf)
ggsave(paste0(OUT_PREFIX, "_occupancy_compare_selected_days.png"), p_occ, width = 10, height = 7, dpi = 600, bg = "white")

outcome_plot <- outcome_compare %>%
  filter(day %in% SELECTED_DAYS) %>%
  select(method, day, subphenotype, p_death_3d, p_death_7d, p_longstay_ge21d) %>%
  pivot_longer(
    cols = c(p_death_3d, p_death_7d, p_longstay_ge21d),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      p_death_3d = "3-day death",
      p_death_7d = "7-day death",
      p_longstay_ge21d = paste0("Long stay >=", LONGSTAY_DAY, "d")
    )
  )

p_out <- outcome_plot %>%
  ggplot(aes(x = factor(day), y = value, group = method, linetype = method, shape = method)) +
  geom_line(linewidth = 0.7, na.rm = TRUE) +
  geom_point(size = 1.4, na.rm = TRUE) +
  facet_grid(metric ~ subphenotype, scales = "free_y") +
  labs(
    x = "ICU day",
    y = "Posterior-weighted descriptive rate",
    title = "Filtering vs smoothing posterior: descriptive outcome-associated profiles"
  ) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUT_PREFIX, "_outcome_compare_selected_days.pdf"), p_out, width = 13, height = 8, device = cairo_pdf)
ggsave(paste0(OUT_PREFIX, "_outcome_compare_selected_days.png"), p_out, width = 13, height = 8, dpi = 600, bg = "white")

summary_selected <- agreement_by_day %>%
  filter(day %in% SELECTED_DAYS) %>%
  transmute(
    day,
    n_stays,
    agreement = round(agreement, 3),
    kappa = round(kappa, 3)
  )
.write_both(summary_selected, paste0(OUT_PREFIX, "_agreement_selected_days"))

manifest <- tibble(
  output = c(
    "post_filtering_stateprob.rds",
    "post_filtering_subphenotypeprob.rds",
    "agreement_overall",
    "agreement_by_day",
    "confusion_overall",
    "dominant_rowlevel_selected_days",
    "occupancy_compare_selected_days",
    "outcome_compare_selected_days",
    "outcome_delta_selected_days",
    "agreement_by_day figure",
    "occupancy_compare_selected_days figure",
    "outcome_compare_selected_days figure"
  ),
  file = c(
    paste0(OUT_PREFIX, "_post_filtering_stateprob.rds"),
    paste0(OUT_PREFIX, "_post_filtering_subphenotypeprob.rds"),
    paste0(OUT_PREFIX, "_agreement_overall.csv/.tsv"),
    paste0(OUT_PREFIX, "_agreement_by_day.csv/.tsv"),
    paste0(OUT_PREFIX, "_confusion_overall.csv/.tsv"),
    paste0(OUT_PREFIX, "_dominant_rowlevel_selected_days.csv/.tsv"),
    paste0(OUT_PREFIX, "_occupancy_compare_selected_days.csv/.tsv"),
    paste0(OUT_PREFIX, "_outcome_compare_selected_days.csv/.tsv"),
    paste0(OUT_PREFIX, "_outcome_delta_selected_days.csv/.tsv"),
    paste0(OUT_PREFIX, "_agreement_by_day.pdf/.png"),
    paste0(OUT_PREFIX, "_occupancy_compare_selected_days.pdf/.png"),
    paste0(OUT_PREFIX, "_outcome_compare_selected_days.pdf/.png")
  ),
  description = c(
    "State-level filtering posterior probabilities from the fixed fitted HMM.",
    "Subphenotype-level filtering posterior probabilities after aggregation by annotation.",
    "Overall dominant-subphenotype agreement between smoothing and filtering.",
    "Day-wise dominant-subphenotype agreement and Cohen kappa.",
    "Confusion table: smoothing dominant subphenotype by filtering dominant subphenotype.",
    "Row-level dominant assignments and outcome flags for selected days.",
    "Selected-day posterior occupancy comparison between smoothing and filtering.",
    "Selected-day posterior-weighted outcome-associated profile comparison.",
    "Filtering-minus-smoothing differences for selected outcome-associated rates.",
    "Line plot of day-wise dominant-subphenotype agreement.",
    "Faceted line plot of selected-day occupancy under smoothing vs filtering.",
    "Faceted line plot of selected-day outcome-associated profiles under smoothing vs filtering."
  )
)
.write_both(manifest, paste0(OUT_PREFIX, "_manifest"))

message("Done: filtering-posterior sensitivity outputs written with prefix: ", OUT_PREFIX)
