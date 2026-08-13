
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
all_of    <- tidyselect::all_of
any_of    <- tidyselect::any_of

OUT_PREFIX <- "day21_destinations"
DAY_DEST <- 21L
ANNO_CANDIDATES <- c(
  MIMIC_ACTIVE_ANNOT_MASTER,
  MIMIC_ACTIVE_ANNOT_MASTER,
  MIMIC_ACTIVE_ANNOT_MASTER
)

.stopif <- function(ok, msg) if (!ok) stop(msg, call. = FALSE)
.pick_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(rep(NA, nrow(df)))
  df[[hit[1]]]
}
first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) NA_character_ else hit[1]
}
state_id_from_name <- function(x) as.integer(sub("^S", "", as.character(x)))

format_subphenotype <- function(x) {
  dplyr::case_when(
    x %in% c("A1 low burden", "A1 low burden") ~ "A1 low burden",
    x %in% c("A2 mild burden", "A2 mild burden") ~ "A2 mild burden",
    x %in% c("B1 renal mild", "B1 renal mild") ~ "B1 renal mild",
    x %in% c("B2 renal severe", "B2 renal severe") ~ "B2 renal severe",
    x %in% c("G respiratory", "G respiratory") ~ "G respiratory",
    x %in% c("D1 lactate/shock", "D1 lactate/shock") ~ "D1 lactate/shock",
    x %in% c("D2 hepatic", "D2 hepatic") ~ "D2 hepatic",
    x %in% c("D3 shock/coagulation-overlap", "D3 shock/coagulation-overlap") ~ "D3 shock/coagulation-overlap",
    TRUE ~ as.character(x)
  )
}
subph_levels <- c(
  "A1 low burden", "A2 mild burden", "B1 renal mild", "B2 renal severe",
  "G respiratory", "D1 lactate/shock", "D2 hepatic", "D3 shock/coagulation-overlap"
)
destination_levels <- c(
  "ICU discharge alive before day 21",
  "Death recorded within 21 days of ICU admission",
  "Active in ICU at day 21",
  "Not observed / unavailable"
)

.stopif(exists("post_smooth_full", envir = .GlobalEnv), "post_smooth_full not found")
.stopif(exists("combined", envir = .GlobalEnv), "combined not found")
.stopif(exists("ids", envir = .GlobalEnv), "ids not found")
.stopif(exists("icustays_intubated_first_patient_admission", envir = .GlobalEnv), "icustays_intubated_first_patient_admission not found")

post <- get("post_smooth_full", envir = .GlobalEnv)
combined <- get("combined", envir = .GlobalEnv)
ids <- as.character(get("ids", envir = .GlobalEnv))
icu <- get("icustays_intubated_first_patient_admission", envir = .GlobalEnv) %>%
  mutate(stay_id = as.character(stay_id))

.stopif(nrow(post) == length(ids), "nrow(post_smooth_full) must equal length(ids)")
.stopif(nrow(combined) == length(ids), "nrow(combined) must equal length(ids)")
.stopif(!is.null(colnames(post)), "post_smooth_full must have state columns named S1, S2, ...")

anno_file <- first_existing(ANNO_CANDIDATES)
.stopif(!is.na(anno_file), "State annotation file not found")
anno <- readr::read_csv(anno_file, show_col_types = FALSE)
.stopif(all(c("state_id", "subphenotype") %in% names(anno)), "Annotation file must contain state_id and subphenotype")

state_map <- tibble(
  state = colnames(post),
  state_id = state_id_from_name(colnames(post))
) %>%
  left_join(anno %>% select(state_id, subphenotype), by = "state_id") %>%
  mutate(subphenotype_label = format_subphenotype(subphenotype))

max_state_col <- max.col(post, ties.method = "first")
dom_state <- colnames(post)[max_state_col]
dom_tbl <- tibble(
  row_id = seq_along(ids),
  stay_id = ids,
  lab_day = suppressWarnings(as.integer(combined$lab_time)),
  dominant_state = dom_state,
  state_id = state_id_from_name(dom_state)
) %>%
  left_join(state_map %>% select(state_id, subphenotype_label), by = "state_id") %>%
  mutate(subphenotype_label = factor(subphenotype_label, levels = subph_levels))

admission_tbl <- dom_tbl %>%
  arrange(stay_id, lab_day, row_id) %>%
  group_by(stay_id) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    stay_id,
    admission_day = lab_day,
    admission_subphenotype = as.character(subphenotype_label),
    admission_state_id = state_id,
    admission_state = dominant_state
  )

icu_intime <- .pick_col(icu, c("intime", "icu_intime", "ICU_intime", "icu_in_time", "in_time", "admittime", "admit_time"))
icu_outtime <- .pick_col(icu, c("outtime", "icu_outtime", "ICU_outtime", "out_time", "dischtime", "discharge_time"))

status_tbl <- icu %>%
  mutate(
    icu_intime = icu_intime,
    icu_outtime = icu_outtime,
    death_days = ifelse(!is.na(dod) & !is.na(icu_intime), as.numeric(difftime(dod, icu_intime, units = "days")), NA_real_),
    los_days = ifelse(!is.na(icu_outtime) & !is.na(icu_intime), as.numeric(difftime(icu_outtime, icu_intime, units = "days")), NA_real_),
    death_days = ifelse(!is.na(death_days) & death_days < 0, NA_real_, death_days),
    los_days = ifelse(!is.na(los_days) & los_days < 0, NA_real_, los_days),
    day21_destination = dplyr::case_when(
      !is.na(death_days) & death_days < DAY_DEST ~ "Death recorded within 21 days of ICU admission",
      is.na(death_days) & !is.na(los_days) & los_days < DAY_DEST ~ "ICU discharge alive before day 21",
      !is.na(death_days) & !is.na(los_days) & los_days < DAY_DEST & death_days >= DAY_DEST ~ "ICU discharge alive before day 21",
      !is.na(los_days) & los_days >= DAY_DEST ~ "Active in ICU at day 21",
      TRUE ~ "Not observed / unavailable"
    ),
    day21_destination = factor(day21_destination, levels = destination_levels)
  ) %>%
  select(stay_id, death_days, los_days, day21_destination)

rowlevel <- admission_tbl %>%
  left_join(status_tbl, by = "stay_id") %>%
  mutate(
    admission_subphenotype = factor(admission_subphenotype, levels = subph_levels),
    day21_destination = factor(as.character(day21_destination), levels = destination_levels)
  )

summary_tbl <- rowlevel %>%
  count(admission_subphenotype, day21_destination, name = "n") %>%
  group_by(admission_subphenotype) %>%
  mutate(
    total_n = sum(n),
    proportion = ifelse(total_n > 0, n / total_n, NA_real_),
    percent = 100 * proportion
  ) %>%
  ungroup() %>%
  arrange(admission_subphenotype, day21_destination)

wide_tbl <- summary_tbl %>%
  select(admission_subphenotype, day21_destination, n, percent, total_n) %>%
  mutate(destination_safe = make.names(as.character(day21_destination))) %>%
  select(-day21_destination) %>%
  pivot_wider(
    names_from = destination_safe,
    values_from = c(n, percent),
    values_fill = 0
  ) %>%
  arrange(admission_subphenotype)

readr::write_csv(rowlevel, paste0(OUT_PREFIX, "_rowlevel.csv"))
readr::write_csv(summary_tbl, paste0(OUT_PREFIX, "_summary_long.csv"))
readr::write_csv(wide_tbl, paste0(OUT_PREFIX, "_summary_wide.csv"))

manifest <- tibble(
  key = c("out_prefix", "day_destination", "annotation_file", "n_patient_stays", "note"),
  value = c(OUT_PREFIX, as.character(DAY_DEST), anno_file, as.character(n_distinct(rowlevel$stay_id)),
            "Day-21 admission-anchored course by admission-day dominant HMM subphenotype.")
)
readr::write_csv(manifest, paste0(OUT_PREFIX, "_manifest.csv"))

message("Done: day-21 destination outputs written with prefix: ", OUT_PREFIX)
