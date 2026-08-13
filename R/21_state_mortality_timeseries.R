
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

ANNO_FILE <- MIMIC_ACTIVE_ANNOT_MASTER
MAX_DAY <- 28L

.check_exists <- function(nm) {
  if (!exists(nm, envir = .GlobalEnv)) {
    stop(sprintf("Object '%s' not found. Run 14_load_final_model.R before this script.", nm))
  }
  get(nm, envir = .GlobalEnv)
}

.pick_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  rep(NA, nrow(df))
}

post_smooth_full <- .check_exists("post_smooth_full")
combined <- .check_exists("combined")
ids <- .check_exists("ids")
icu <- .check_exists("icustays_intubated_first_patient_admission")

stopifnot(nrow(post_smooth_full) == length(ids))
stopifnot(nrow(post_smooth_full) == nrow(combined))
if (!("lab_time" %in% names(combined))) stop("combined must contain 'lab_time'.")

K <- ncol(post_smooth_full)
state_levels <- colnames(post_smooth_full)
if (is.null(state_levels)) state_levels <- paste0("S", seq_len(K))

if (!file.exists(ANNO_FILE)) {
  stop("Annotation file not found: ", ANNO_FILE)
}

anno_tbl <- readr::read_csv(ANNO_FILE, show_col_types = FALSE)
required_anno <- c("state_id", "subphenotype", "assigned_subtype", "tags_main", "orthogonal_types", "ABGD")
missing_anno <- setdiff(required_anno, names(anno_tbl))
if (length(missing_anno) > 0) {
  stop("Annotation file is missing required columns: ", paste(missing_anno, collapse = ", "))
}
anno_tbl <- anno_tbl %>%
  transmute(
    state_id = as.integer(state_id),
    subphenotype,
    assigned_subtype,
    tags = tags_main,
    union_Orthogonal = orthogonal_types,
    union_ABGD = ABGD
  )

icu <- icu %>% mutate(stay_id = as.character(stay_id))
ids <- as.character(ids)
death_ids <- icu %>% filter(!is.na(dod)) %>% pull(stay_id)
death_flag <- as.integer(ids %in% death_ids)

icu_intime <- .pick_col(
  icu,
  c("intime", "icu_intime", "ICU_intime", "icu_in_time", "in_time", "admittime", "admit_time")
)

death_day_tbl <- icu %>%
  transmute(
    stay_id = as.character(stay_id),
    dod = dod,
    icu_intime = icu_intime
  ) %>%
  mutate(
    death_days = ifelse(
      !is.na(dod) & !is.na(icu_intime),
      as.numeric(difftime(dod, icu_intime, units = "days")),
      NA_real_
    ),
    death_days = ifelse(!is.na(death_days) & death_days < 0, NA_real_, death_days)
  ) %>%
  select(stay_id, death_days)

death_days_vec <- death_day_tbl$death_days[match(ids, death_day_tbl$stay_id)]

lab_day <- suppressWarnings(as.integer(combined$lab_time))
death3_flag <- as.integer(!is.na(death_days_vec) & death_days_vec <= (lab_day + 2L))
death7_flag <- as.integer(!is.na(death_days_vec) & death_days_vec <= (lab_day + 6L))

days <- seq_len(MAX_DAY)

state_day_tbls <- purrr::map(days, function(d) {
  mask_d <- lab_day == d
  if (!any(mask_d, na.rm = TRUE)) {
    return(tibble(
      state = state_levels,
      !!paste0("N_state_any_d", d) := 0,
      !!paste0("p_death_eventual_d", d) := NA_real_,
      !!paste0("p_death_3d_d", d) := NA_real_,
      !!paste0("p_death_7d_d", d) := NA_real_
    ))
  }

  post_d <- post_smooth_full[mask_d, , drop = FALSE]

  N_any <- colSums(post_d)
  N_evt <- colSums(post_d * death_flag[mask_d])
  N_3d <- colSums(post_d * death3_flag[mask_d])
  N_7d <- colSums(post_d * death7_flag[mask_d])

  tibble(
    state = state_levels,
    !!paste0("N_state_any_d", d) := as.numeric(N_any),
    !!paste0("p_death_eventual_d", d) := ifelse(N_any > 0, as.numeric(N_evt / N_any), NA_real_),
    !!paste0("p_death_3d_d", d) := ifelse(N_any > 0, as.numeric(N_3d / N_any), NA_real_),
    !!paste0("p_death_7d_d", d) := ifelse(N_any > 0, as.numeric(N_7d / N_any), NA_real_)
  )
})

state_day_wide <- purrr::reduce(
  state_day_tbls,
  .init = tibble(state = state_levels),
  .f = ~ left_join(.x, .y, by = "state")
) %>%
  mutate(state_id = as.integer(str_remove(state, "^S"))) %>%
  left_join(anno_tbl, by = "state_id") %>%
  relocate(state, state_id, assigned_subtype, subphenotype, tags, union_Orthogonal, union_ABGD)

p_cols <- names(state_day_wide)[grepl("^p_death_", names(state_day_wide))]
state_day_wide_round <- state_day_wide %>%
  mutate(across(all_of(p_cols), ~ round(.x, 3)))

write_csv(state_day_wide_round, "state_mortality_table_state_pdeath_timeseries_all_days.csv")

state_day_long <- state_day_wide %>%
  pivot_longer(
    cols = matches("^(N_state_any_d|p_death_eventual_d|p_death_3d_d|p_death_7d_d)\\d+$"),
    names_to = c(".value", "day"),
    names_pattern = "^(N_state_any_d|p_death_eventual_d|p_death_3d_d|p_death_7d_d)(\\d+)$"
  ) %>%
  mutate(day = as.integer(day)) %>%
  arrange(subphenotype, state_id, day)

write_csv(state_day_long, "state_mortality_table_state_pdeath_timeseries_all_days_long.csv")

N_wide <- state_day_wide %>%
  select(state, state_id, assigned_subtype, subphenotype, tags, union_Orthogonal, union_ABGD, starts_with("N_state_any_d"))
write_csv(N_wide, "state_mortality_table_state_weight_timeseries_all_days.csv")

message("Done: state-by-day mortality tables were written.")
