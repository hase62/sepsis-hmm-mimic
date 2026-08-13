
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tibble)
})

get_env_path <- function(name, default) {
  x <- Sys.getenv(name, unset = default)
  if (!nzchar(x)) default else x
}

paths <- list(
  cohort_before = get_env_path(
    "MIMIC_PRE_HMM_COHORT_FILE",
    "icustays_intubated_first_patient_admission.csv"
  ),
  primary_matrix_before = get_env_path(
    "MIMIC_PRE_HMM_PRIMARY_MATRIX_FILE",
    "mimiciv.combined.wide.csv"
  ),
  analytical_day_by_stay = get_env_path(
    "MIMIC_ANALYTICAL_DAY_BY_STAY_FILE",
    "analytical_day_by_stay.csv"
  ),
  analytical_day_detailed = get_env_path(
    "MIMIC_ANALYTICAL_DAY_DETAILED_FILE",
    "analytical_day_detailed_summary.csv"
  ),
  final_cohort = get_env_path(
    "MIMIC_HMM_COHORT_FILE",
    "icustays_intubated_first_patient_admission.HMM_eligible.csv"
  ),
  final_primary_matrix = get_env_path(
    "MIMIC_HMM_PRIMARY_MATRIX_FILE",
    "mimiciv.combined.wide.HMM_raw.csv"
  ),
  qc_dir = get_env_path("MIMIC_HMM_INPUT_QC_DIR", "hmm_input_qc")
)

dir.create(paths$qc_dir, recursive = TRUE, showWarnings = FALSE)
final_paths <- list(
  eligibility = file.path(paths$qc_dir, "eligibility_by_stay.csv"),
  excluded = file.path(paths$qc_dir, "excluded_stays.csv"),
  exclusion_summary = file.path(paths$qc_dir, "exclusion_summary.csv"),
  flow = file.path(paths$qc_dir, "cohort_flow.csv"),
  cohort_summary = file.path(paths$qc_dir, "final_cohort_summary.csv"),
  analytical_day_qc = file.path(paths$qc_dir, "analytical_day_post_exclusion_qc.csv"),
  matrix_qc = file.path(paths$qc_dir, "input_matrix_qc.csv"),
  manifest = file.path(paths$qc_dir, "hmm_input_manifest.rds")
)

required_files <- c(
  paths$cohort_before,
  paths$primary_matrix_before,
  paths$analytical_day_by_stay,
  paths$analytical_day_detailed
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing required input file(s):\n",
    paste0("  - ", missing_files, collapse = "\n"),
    "\nRun the cohort, analytical-matrix, and analytical-day steps first."
  )
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  readr::write_csv(x, tmp, na = "")
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(tmp, path)) {
    unlink(tmp, force = TRUE)
    stop("Could not write output: ", path)
  }
  invisible(path)
}

atomic_fwrite <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  data.table::fwrite(x, tmp, na = "")
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(tmp, path)) {
    unlink(tmp, force = TRUE)
    stop("Could not write output: ", path)
  }
  invisible(path)
}

file_signature <- function(path, md5 = FALSE) {
  if (!file.exists(path)) {
    return(list(path = normalizePath(path, mustWork = FALSE), exists = FALSE))
  }
  z <- file.info(path)
  list(
    path = normalizePath(path, mustWork = TRUE),
    exists = TRUE,
    size = unname(z$size),
    mtime = format(z$mtime, "%Y-%m-%d %H:%M:%OS6 %Z"),
    md5 = if (md5) unname(tools::md5sum(path)) else NULL
  )
}

id_md5 <- function(ids) {
  tf <- tempfile(fileext = ".txt")
  on.exit(unlink(tf), add = TRUE)
  writeLines(sort(unique(as.character(ids))), tf)
  unname(tools::md5sum(tf))
}

cohort <- readr::read_csv(paths$cohort_before, show_col_types = FALSE, progress = FALSE)
if (!all(c("stay_id", "age") %in% names(cohort))) {
  stop("Cohort metadata must contain stay_id and age.")
}
cohort <- cohort %>%
  mutate(
    stay_id = as.character(stay_id),
    age = suppressWarnings(as.numeric(age))
  )
if (anyDuplicated(cohort$stay_id)) stop("Cohort metadata contains duplicated stay_id values.")
if (nrow(cohort) == 0L) stop("Cohort metadata is empty.")

day_qc <- readr::read_csv(paths$analytical_day_by_stay, show_col_types = FALSE, progress = FALSE)
required_qc_cols <- c(
  "stay_id", "n_union_days", "n_included_lab_days",
  "n_excluded_vital_or_proc_days", "n_before_first_lab",
  "n_between_lab_days", "n_after_last_lab"
)
missing_qc_cols <- setdiff(required_qc_cols, names(day_qc))
if (length(missing_qc_cols) > 0L) {
  stop("Analytical-day QC lacks: ", paste(missing_qc_cols, collapse = ", "))
}
day_qc <- day_qc %>% mutate(stay_id = as.character(stay_id))
if (anyDuplicated(day_qc$stay_id)) stop("Analytical-day QC contains duplicated stay_id values.")

primary <- data.table::fread(paths$primary_matrix_before, showProgress = FALSE)
if (!"stay_id" %in% names(primary)) stop("Primary raw matrix lacks stay_id.")
primary[, stay_id := as.character(stay_id)]
primary_counts <- primary[, .(n_matrix_analytical_days = .N), by = stay_id]

eligibility <- cohort %>%
  left_join(day_qc, by = "stay_id") %>%
  left_join(as.data.frame(primary_counts), by = "stay_id") %>%
  mutate(
    across(
      any_of(c(
        "n_union_days", "n_included_lab_days", "n_excluded_vital_or_proc_days",
        "n_before_first_lab", "n_between_lab_days", "n_after_last_lab",
        "n_matrix_analytical_days"
      )),
      ~ dplyr::coalesce(as.integer(.x), 0L)
    ),
    lab_anchor_eligible = n_included_lab_days >= 1L,
    age_eligible = is.finite(age) & age >= 18 & age < 80,
    sequence_length_eligible = n_matrix_analytical_days >= 2L,
    exclusion_reason_code = case_when(
      !lab_anchor_eligible ~ "NO_LAB_ANCHORED_ANALYTICAL_DAY",
      !is.finite(age) ~ "MISSING_AGE",
      age < 18 ~ "AGE_LT_18",
      age >= 80 ~ "AGE_GE_80",
      n_matrix_analytical_days == 0L ~ "NO_PRIMARY_MATRIX_ROW",
      n_matrix_analytical_days == 1L ~ "ONLY_ONE_ANALYTICAL_DAY",
      TRUE ~ NA_character_
    ),
    exclusion_stage = case_when(
      !lab_anchor_eligible ~ "1_laboratory_anchor",
      !age_eligible ~ "2_age",
      !sequence_length_eligible ~ "3_sequence_length",
      TRUE ~ "4_included"
    ),
    exclusion_reason = recode(
      exclusion_reason_code,
      NO_LAB_ANCHORED_ANALYTICAL_DAY = "No laboratory-anchored analytical day",
      MISSING_AGE = "Age missing",
      AGE_LT_18 = "Age younger than 18 years",
      AGE_GE_80 = "Age 80 years or older",
      NO_PRIMARY_MATRIX_ROW = "No analytical-day row in the primary raw matrix",
      ONLY_ONE_ANALYTICAL_DAY = "Only one analytical day",
      .default = NA_character_
    ),
    included_in_hmm_raw_cohort = is.na(exclusion_reason_code)
  )

extra_qc_ids <- setdiff(day_qc$stay_id, cohort$stay_id)
extra_matrix_ids <- setdiff(primary_counts$stay_id, cohort$stay_id)
if (length(extra_qc_ids) > 0L || length(extra_matrix_ids) > 0L) {
  stop("Cohort mismatch between metadata, analytical-day QC, and primary matrix.")
}

final_ids <- eligibility$stay_id[eligibility$included_in_hmm_raw_cohort]
if (length(final_ids) == 0L) stop("No stay met final pre-HMM eligibility.")

final_cohort <- cohort %>% filter(stay_id %in% final_ids)
final_primary <- primary[stay_id %in% final_ids]
final_primary_counts <- final_primary[, .N, by = stay_id]
if (!setequal(final_cohort$stay_id, final_primary_counts$stay_id)) {
  stop("Final cohort and final primary matrix contain different stay sets.")
}
if (any(final_primary_counts$N < 2L)) stop("A final HMM stay has fewer than two analytical-day rows.")

sort_columns <- intersect(c("stay_id", "charttime_year_day", "lab_time"), names(final_primary))
if (length(sort_columns) > 0L) data.table::setorderv(final_primary, sort_columns)
atomic_write_csv(final_cohort, paths$final_cohort)
atomic_fwrite(final_primary, paths$final_primary_matrix)

n_subjects <- function(x) {
  if (!"subject_id" %in% names(x)) return(NA_integer_)
  dplyr::n_distinct(x$subject_id[!is.na(x$subject_id)])
}

stage1_enter <- eligibility
stage1_excluded <- stage1_enter %>% filter(!lab_anchor_eligible)
stage1_remain <- stage1_enter %>% filter(lab_anchor_eligible)
stage2_excluded <- stage1_remain %>% filter(!age_eligible)
stage2_remain <- stage1_remain %>% filter(age_eligible)
stage3_excluded <- stage2_remain %>% filter(!sequence_length_eligible)
stage3_remain <- stage2_remain %>% filter(sequence_length_eligible)

flow_row <- function(order, code, label, entering, excluded_now, remaining, reason) {
  tibble(
    stage_order = order,
    stage_code = code,
    stage = label,
    n_stays_entering = n_distinct(entering$stay_id),
    n_subjects_entering = n_subjects(entering),
    n_stays_excluded = n_distinct(excluded_now$stay_id),
    n_subjects_excluded = n_subjects(excluded_now),
    exclusion_reason = reason,
    n_stays_remaining = n_distinct(remaining$stay_id),
    n_subjects_remaining = n_subjects(remaining)
  )
}

flow <- bind_rows(
  flow_row(
    1L, "LAB_ANCHOR", "Laboratory-anchored analytical-day eligibility",
    stage1_enter, stage1_excluded, stage1_remain,
    "No laboratory-anchored analytical day"
  ),
  flow_row(
    2L, "AGE", "Prespecified age eligibility",
    stage1_remain, stage2_excluded, stage2_remain,
    "Age missing, younger than 18 years, or 80 years or older"
  ),
  flow_row(
    3L, "SEQUENCE_LENGTH", "Minimum HMM sequence length",
    stage2_remain, stage3_excluded, stage3_remain,
    "Fewer than two analytical days"
  )
)

excluded <- eligibility %>%
  filter(!included_in_hmm_raw_cohort) %>%
  relocate(
    stay_id,
    any_of(c("subject_id", "hadm_id", "age")),
    exclusion_stage,
    exclusion_reason_code,
    exclusion_reason
  )

exclusion_summary <- excluded %>%
  group_by(exclusion_stage, exclusion_reason_code, exclusion_reason) %>%
  summarise(
    n_stays = n_distinct(stay_id),
    n_subjects = if ("subject_id" %in% names(excluded)) {
      n_distinct(subject_id[!is.na(subject_id)])
    } else {
      NA_integer_
    },
    n_union_days = sum(n_union_days, na.rm = TRUE),
    n_primary_matrix_days = sum(n_matrix_analytical_days, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(exclusion_stage, exclusion_reason_code)

cohort_summary <- tibble(
  definition = c("Pre-HMM candidate cohort", "Final HMM raw cohort before missing-data completion and normalization"),
  n_stays = c(n_distinct(cohort$stay_id), n_distinct(final_cohort$stay_id)),
  n_subjects = c(n_subjects(cohort), n_subjects(final_cohort)),
  n_primary_matrix_rows = c(nrow(primary), nrow(final_primary)),
  minimum_sequence_length = c(min(primary_counts$n_matrix_analytical_days), min(final_primary_counts$N)),
  maximum_sequence_length = c(max(primary_counts$n_matrix_analytical_days), max(final_primary_counts$N))
)

analytical_day_qc <- bind_rows(
  eligibility %>%
    summarise(
      population = "Pre-HMM candidate cohort",
      n_stays = n_distinct(stay_id),
      n_union_days = sum(n_union_days, na.rm = TRUE),
      n_lab_anchored_days = sum(n_included_lab_days, na.rm = TRUE),
      n_non_lab_days = sum(n_excluded_vital_or_proc_days, na.rm = TRUE),
      n_before_first_lab = sum(n_before_first_lab, na.rm = TRUE),
      n_between_lab_days = sum(n_between_lab_days, na.rm = TRUE),
      n_after_last_lab = sum(n_after_last_lab, na.rm = TRUE)
    ),
  eligibility %>%
    filter(included_in_hmm_raw_cohort) %>%
    summarise(
      population = "Final HMM raw cohort",
      n_stays = n_distinct(stay_id),
      n_union_days = sum(n_union_days, na.rm = TRUE),
      n_lab_anchored_days = sum(n_included_lab_days, na.rm = TRUE),
      n_non_lab_days = sum(n_excluded_vital_or_proc_days, na.rm = TRUE),
      n_before_first_lab = sum(n_before_first_lab, na.rm = TRUE),
      n_between_lab_days = sum(n_between_lab_days, na.rm = TRUE),
      n_after_last_lab = sum(n_after_last_lab, na.rm = TRUE)
    )
) %>%
  mutate(
    proportion_union_days_retained = if_else(
      n_union_days > 0L,
      n_lab_anchored_days / n_union_days,
      NA_real_
    )
  )

matrix_qc <- tibble(
  primary_source_file = paths$primary_matrix_before,
  primary_output_file = paths$final_primary_matrix,
  n_primary_rows_before = nrow(primary),
  n_primary_stays_before = uniqueN(primary$stay_id),
  n_primary_rows_final = nrow(final_primary),
  n_primary_stays_final = uniqueN(final_primary$stay_id),
  final_stay_set_matches_metadata = setequal(final_cohort$stay_id, unique(final_primary$stay_id)),
  minimum_final_sequence_length = min(final_primary_counts$N),
  maximum_final_sequence_length = max(final_primary_counts$N)
)

atomic_write_csv(eligibility, final_paths$eligibility)
atomic_write_csv(excluded, final_paths$excluded)
atomic_write_csv(exclusion_summary, final_paths$exclusion_summary)
atomic_write_csv(flow, final_paths$flow)
atomic_write_csv(cohort_summary, final_paths$cohort_summary)
atomic_write_csv(analytical_day_qc, final_paths$analytical_day_qc)
atomic_write_csv(matrix_qc, final_paths$matrix_qc)

manifest <- list(
  version = "hmm_input_v1",
  criteria = list(
    laboratory_anchor = "n_included_lab_days >= 1",
    age = "age >= 18 & age < 80",
    sequence_length = "n_matrix_analytical_days >= 2"
  ),
  input_files = list(
    cohort = file_signature(paths$cohort_before, md5 = TRUE),
    primary_matrix = file_signature(paths$primary_matrix_before),
    analytical_day_by_stay = file_signature(paths$analytical_day_by_stay, md5 = TRUE),
    analytical_day_detailed = file_signature(paths$analytical_day_detailed, md5 = TRUE)
  ),
  output_files = list(
    cohort = file_signature(paths$final_cohort, md5 = TRUE),
    primary_matrix = file_signature(paths$final_primary_matrix),
    flow = file_signature(final_paths$flow, md5 = TRUE),
    eligibility = file_signature(final_paths$eligibility, md5 = TRUE)
  ),
  n_stays_before = n_distinct(cohort$stay_id),
  n_stays_excluded_no_lab_anchor = n_distinct(stage1_excluded$stay_id),
  n_stays_excluded_age = n_distinct(stage2_excluded$stay_id),
  n_stays_excluded_sequence_length = n_distinct(stage3_excluded$stay_id),
  n_stays_final = n_distinct(final_cohort$stay_id),
  final_stay_ids_md5 = id_md5(final_ids),
  final_cohort_file = paths$final_cohort,
  final_primary_matrix_file = paths$final_primary_matrix
)
saveRDS(manifest, final_paths$manifest, version = 3)

cat("Final HMM input completed.\n")
cat("Candidate stays            :", n_distinct(cohort$stay_id), "\n")
cat("Excluded: no lab anchor    :", n_distinct(stage1_excluded$stay_id), "\n")
cat("Excluded: age              :", n_distinct(stage2_excluded$stay_id), "\n")
cat("Excluded: sequence <2 days :", n_distinct(stage3_excluded$stay_id), "\n")
cat("Final HMM cohort           :", n_distinct(final_cohort$stay_id), "\n")
cat("Final analytical rows      :", nrow(final_primary), "\n")
cat("QC directory               :", paths$qc_dir, "\n")
