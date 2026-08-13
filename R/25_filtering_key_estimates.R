
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
group_by  <- dplyr::group_by
ungroup   <- dplyr::ungroup
arrange   <- dplyr::arrange
bind_rows <- dplyr::bind_rows

OUT_PREFIX <- Sys.getenv("MIMIC_FILTERING_KEY_PREFIX", "filtering_key_estimates")
BASE_PREFIX_CANDIDATES <- c(
  Sys.getenv("MIMIC_FILTERING_BASE_PREFIX", ""),
  "filtering_sensitivity",
  "filtering_sensitivity"
)
BASE_PREFIX_CANDIDATES <- BASE_PREFIX_CANDIDATES[nzchar(BASE_PREFIX_CANDIDATES)]
SELECTED_DAYS <- c(1L, 3L, 5L, 7L, 10L, 14L, 21L, 28L)

.write_both <- function(df, stem) {
  readr::write_csv(df, paste0(stem, ".csv"))
  readr::write_tsv(df, paste0(stem, ".tsv"))
  invisible(df)
}

find_prefix <- function(prefixes) {
  for (p in prefixes) {
    if (file.exists(paste0(p, "_agreement_overall.csv")) && file.exists(paste0(p, "_outcome_compare_selected_days.csv"))) return(p)
  }
  NA_character_
}

base_prefix <- find_prefix(BASE_PREFIX_CANDIDATES)
if (is.na(base_prefix)) {
  stop("Could not find filtering sensitivity outputs. Tried prefixes: ", paste(BASE_PREFIX_CANDIDATES, collapse = ", "), call. = FALSE)
}
message("Using filtering sensitivity prefix: ", base_prefix)

agreement_overall <- readr::read_csv(paste0(base_prefix, "_agreement_overall.csv"), show_col_types = FALSE)
agreement_by_day <- readr::read_csv(paste0(base_prefix, "_agreement_by_day.csv"), show_col_types = FALSE)
outcome_compare <- readr::read_csv(paste0(base_prefix, "_outcome_compare_selected_days.csv"), show_col_types = FALSE)

pick_signal <- function(sp_regex, metric_col, signal_label, days = SELECTED_DAYS) {
  outcome_compare %>%
    filter(day %in% days, str_detect(as.character(subphenotype), sp_regex)) %>%
    select(method, day, subphenotype, value = all_of(metric_col), N_weighted, n_stays) %>%
    pivot_wider(names_from = method, values_from = c(value, N_weighted), names_sep = "__") %>%
    mutate(
      signal = signal_label,
      metric = metric_col,
      smoothing_estimate = value__smoothing,
      filtering_estimate = value__filtering,
      absolute_difference_filtering_minus_smoothing = filtering_estimate - smoothing_estimate
    ) %>%
    select(signal, day, subphenotype, metric, smoothing_estimate, filtering_estimate,
           absolute_difference_filtering_minus_smoothing,
           N_weighted__smoothing, N_weighted__filtering, n_stays)
}

key_long <- bind_rows(
  pick_signal("^D1\\b", "p_death_3d", "D1 3-day death"),
  pick_signal("^D1\\b", "p_death_7d", "D1 7-day death"),
  pick_signal("^B2\\b", "p_longstay_ge21d", "B2 long stay >=21d"),
  pick_signal("^G\\b", "p_longstay_ge21d", "G respiratory long stay >=21d")
) %>%
  arrange(signal, day)

compact_days <- c(1L, 3L, 7L, 14L)
key_compact <- key_long %>%
  filter(day %in% compact_days) %>%
  mutate(
    smoothing_percent = 100 * smoothing_estimate,
    filtering_percent = 100 * filtering_estimate,
    absolute_difference_percentage_points = 100 * absolute_difference_filtering_minus_smoothing
  ) %>%
  select(signal, day, subphenotype, metric, smoothing_percent, filtering_percent,
         absolute_difference_percentage_points, N_weighted__smoothing, N_weighted__filtering, n_stays)

selected_agree <- agreement_by_day %>%
  filter(day %in% SELECTED_DAYS)

overview <- tibble(
  quantity = c(
    "overall_dominant_subphenotype_agreement",
    "overall_cohen_kappa",
    "day1_dominant_subphenotype_agreement",
    "day1_cohen_kappa",
    "minimum_selected_day_agreement",
    "minimum_selected_day_kappa",
    "day_of_minimum_selected_day_agreement"
  ),
  estimate = c(
    agreement_overall$agreement[1],
    agreement_overall$kappa[1],
    selected_agree$agreement[selected_agree$day == 1L][1],
    selected_agree$kappa[selected_agree$day == 1L][1],
    min(selected_agree$agreement, na.rm = TRUE),
    min(selected_agree$kappa, na.rm = TRUE),
    selected_agree$day[which.min(selected_agree$agreement)][1]
  ),
  interpretation = c(
    "Overall dominant-subphenotype agreement between smoothing and filtering assignments.",
    "Overall Cohen kappa between smoothing and filtering dominant assignments.",
    "Day 1 agreement, expected to be the most sensitive to retrospective smoothing because subsequent observations are not yet available to filtering.",
    "Day 1 Cohen kappa.",
    "Lowest agreement among selected reporting days.",
    "Lowest kappa among selected reporting days.",
    "Selected ICU day at which agreement was lowest."
  )
)

.write_both(overview, paste0(OUT_PREFIX, "_overview"))
.write_both(key_long, paste0(OUT_PREFIX, "_key_estimates_long"))
.write_both(key_compact, paste0(OUT_PREFIX, "_key_estimates_compact"))

manifest <- tibble(
  file = c(
    paste0(OUT_PREFIX, "_overview.csv"),
    paste0(OUT_PREFIX, "_key_estimates_long.csv"),
    paste0(OUT_PREFIX, "_key_estimates_compact.csv")
  ),
  description = c(
    "Overall and selected-day smoothing/filtering agreement summary, including day-1 agreement and selected-day minima.",
    "Long-form smoothing versus filtering estimates for D1 short-term mortality and B2/G respiratory long-stay signals across selected ICU days.",
    "Compact smoothing versus filtering key estimates at ICU days 1, 3, 7, and 14."
  ),
  source_prefix = base_prefix,
  note = "Filtering summaries support robustness to prospective-compatible assignment but do not validate real-time clinical prediction."
)
.write_both(manifest, paste0(OUT_PREFIX, "_manifest"))

message("Done: filtering key estimate tables written with prefix: ", OUT_PREFIX)
