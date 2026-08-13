
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
if (!exists("MIMIC_ACTIVE_ANNOT_MASTER", inherits = TRUE)) {
  source(code_path("17_analysis_config.R"))
}

if (!file.exists(code_path("17_analysis_config.R"))) stop("Missing shared analysis configuration: 17_analysis_config.R", call. = FALSE)
source(code_path("17_analysis_config.R"))


suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

state_master <- tibble::as_tibble(mimic_read_active_annotation())

select <- dplyr::select
filter <- dplyr::filter
mutate <- dplyr::mutate
rename <- dplyr::rename
arrange <- dplyr::arrange
summarise <- dplyr::summarise
group_by <- dplyr::group_by
ungroup <- dplyr::ungroup
left_join <- dplyr::left_join
bind_rows <- dplyr::bind_rows
all_of <- tidyselect::all_of

OUT_PREFIX <- Sys.getenv("MIMIC_STATE_DOMINANCE_OUT_PREFIX", unset = "state_dominance_qc")
OUT_DIR <- Sys.getenv("MIMIC_STATE_DOMINANCE_OUT_DIR", unset = ".")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_path <- function(suffix) file.path(OUT_DIR, paste0(OUT_PREFIX, suffix))

subphenotype_order <- c(
  "A1 low burden",
  "A2 mild burden",
  "B1 renal mild",
  "B2 renal severe",
  "G respiratory",
  "D1 lactate/shock",
  "D2 hepatic",
  "D3 shock/coagulation-overlap"
)

normalize_subphenotype <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    str_detect(x, "^A1\\b") ~ "A1 low burden",
    str_detect(x, "^A2\\b") ~ "A2 mild burden",
    str_detect(x, "^B1\\b") ~ "B1 renal mild",
    str_detect(x, "^B2\\b") ~ "B2 renal severe",
    str_detect(x, "^G\\b") ~ "G respiratory",
    str_detect(x, "^D1\\b") ~ "D1 lactate/shock",
    str_detect(x, "^D2\\b") ~ "D2 hepatic",
    str_detect(x, "^D3\\b") ~ "D3 shock/coagulation-overlap",
    TRUE ~ x
  )
}

find_annotation_table <- function() {
  if (exists("state_master", inherits = TRUE)) {
    sm <- get("state_master", inherits = TRUE)
    if (is.data.frame(sm) && all(c("state_id", "subphenotype") %in% names(sm))) {
      message("Using annotation object: state_master")
      return(as_tibble(sm))
    }
  }
  candidates <- c(MIMIC_ACTIVE_ANNOT_MASTER)
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0) {
    candidates <- list.files(pattern = "state_master.*\\.csv$", ignore.case = TRUE, full.names = FALSE)
  }
  if (length(candidates) == 0) stop("State annotation table was not found.", call. = FALSE)
  message("Using annotation file: ", candidates[1])
  readr::read_csv(candidates[1], show_col_types = FALSE)
}

ann <- find_annotation_table()
if (!all(c("state_id", "subphenotype") %in% names(ann))) {
  stop("Annotation table must contain columns: state_id, subphenotype", call. = FALSE)
}
ann_map <- ann %>%
  transmute(
    state_id = as.integer(state_id),
    subphenotype_raw = as.character(subphenotype),
    subphenotype = normalize_subphenotype(subphenotype_raw)
  ) %>%
  distinct(state_id, .keep_all = TRUE)

if (!exists("post_smooth_full", inherits = TRUE)) {
  stop("Required object post_smooth_full was not found. Run the HMM state analysis first.", call. = FALSE)
}
post <- get("post_smooth_full", inherits = TRUE)
if (!is.matrix(post) && !is.data.frame(post)) stop("post_smooth_full must be a matrix or data frame.", call. = FALSE)
post <- as.matrix(post)
state_cols <- colnames(post)
if (is.null(state_cols)) state_cols <- paste0("S", seq_len(ncol(post)))
state_id <- suppressWarnings(as.integer(gsub("^S", "", state_cols)))
if (anyNA(state_id)) {
  state_id <- seq_len(ncol(post))
  warning("Could not parse state ids from posterior column names; using column order as state_id.")
}

state_weight_tbl <- tibble(
  state_id = state_id,
  state_column = state_cols,
  posterior_weight = as.numeric(colSums(post, na.rm = TRUE)),
  mean_posterior_occupancy = as.numeric(colMeans(post, na.rm = TRUE))
) %>%
  left_join(ann_map, by = "state_id")
if (any(is.na(state_weight_tbl$subphenotype))) {
  missing_states <- paste(state_weight_tbl$state_id[is.na(state_weight_tbl$subphenotype)], collapse = ", ")
  stop("Some posterior states were not found in annotation table: ", missing_states, call. = FALSE)
}

total_weight <- sum(state_weight_tbl$posterior_weight, na.rm = TRUE)
state_contrib <- state_weight_tbl %>%
  group_by(subphenotype) %>%
  mutate(
    subphenotype_posterior_weight = sum(posterior_weight, na.rm = TRUE),
    state_share_within_subphenotype = ifelse(subphenotype_posterior_weight > 0, posterior_weight / subphenotype_posterior_weight, NA_real_),
    state_occupancy_share_all_states = ifelse(total_weight > 0, posterior_weight / total_weight, NA_real_),
    state_rank_within_subphenotype = rank(-state_share_within_subphenotype, ties.method = "first")
  ) %>%
  ungroup() %>%
  mutate(
    subphenotype = factor(subphenotype, levels = subphenotype_order)
  ) %>%
  arrange(subphenotype, state_rank_within_subphenotype) %>%
  mutate(subphenotype = as.character(subphenotype))

summary_tbl <- state_contrib %>%
  group_by(subphenotype) %>%
  summarise(
    n_states = dplyr::n(),
    posterior_weight = sum(posterior_weight, na.rm = TRUE),
    subphenotype_occupancy_share_all_states = first(subphenotype_posterior_weight) / total_weight,
    largest_state_id = state_id[which.max(state_share_within_subphenotype)],
    largest_state_share_within_subphenotype = max(state_share_within_subphenotype, na.rm = TRUE),
    top2_state_ids = paste(state_id[order(-state_share_within_subphenotype)][seq_len(min(2, dplyr::n()))], collapse = ";"),
    top2_state_share_within_subphenotype = sum(sort(state_share_within_subphenotype, decreasing = TRUE)[seq_len(min(2, dplyr::n()))], na.rm = TRUE),
    top3_state_ids = paste(state_id[order(-state_share_within_subphenotype)][seq_len(min(3, dplyr::n()))], collapse = ";"),
    top3_state_share_within_subphenotype = sum(sort(state_share_within_subphenotype, decreasing = TRUE)[seq_len(min(3, dplyr::n()))], na.rm = TRUE),
    single_state_exceeds_50pct = largest_state_share_within_subphenotype > 0.50,
    top2_exceeds_75pct = top2_state_share_within_subphenotype > 0.75,
    .groups = "drop"
  ) %>%
  mutate(
    subphenotype = factor(subphenotype, levels = subphenotype_order),
    dominance_interpretation = dplyr::case_when(
      single_state_exceeds_50pct ~ "Largest latent state accounts for >50% of this subphenotype's posterior occupancy; interpret with state-level checks.",
      top2_exceeds_75pct ~ "No single state exceeds 50%, but the top two states account for >75%; interpret with state-level checks.",
      TRUE ~ "Posterior occupancy is distributed across multiple latent states."
    )
  ) %>%
  arrange(subphenotype) %>%
  mutate(subphenotype = as.character(subphenotype))

reporting_summary <- summary_tbl %>%
  transmute(
    subphenotype,
    number_of_latent_states = n_states,
    subphenotype_posterior_occupancy_percent = 100 * subphenotype_occupancy_share_all_states,
    largest_state_id,
    largest_state_share_percent_within_subphenotype = 100 * largest_state_share_within_subphenotype,
    top2_state_ids,
    top2_state_share_percent_within_subphenotype = 100 * top2_state_share_within_subphenotype,
    top3_state_ids,
    top3_state_share_percent_within_subphenotype = 100 * top3_state_share_within_subphenotype,
    single_state_exceeds_50pct,
    top2_exceeds_75pct,
    dominance_interpretation
  )

manifest <- tibble(
  key = c(
    "script", "out_prefix", "n_rows_posterior", "n_states_posterior", "total_posterior_weight",
    "annotation_columns", "interpretation"
  ),
  value = c(
    "23_state_dominance_qc.R",
    OUT_PREFIX,
    as.character(nrow(post)),
    as.character(ncol(post)),
    format(total_weight, digits = 12),
    "state_id, subphenotype",
    "Lightweight dominance QC based on smoothed posterior state mass. It evaluates whether aggregated subphenotype occupancy is concentrated in one or two latent states; it does not refit the HMM or perform causal/treatment-response analysis."
  )
)

readr::write_csv(state_contrib, out_path("_state_contributions.csv"))
readr::write_csv(summary_tbl, out_path("_summary.csv"))
readr::write_csv(reporting_summary, out_path("_reporting_summary.csv"))
readr::write_csv(manifest, out_path("_manifest.csv"))

message("Done: subphenotype state dominance QC written with prefix: ", OUT_PREFIX)
message("Outputs in: ", normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE))
