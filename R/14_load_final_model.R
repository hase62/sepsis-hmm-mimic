
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
if (!file.exists("model_config.R")) stop("model_config.R not found. Run 13_select_final_model.R first.", call. = FALSE)
source("model_config.R")
suppressPackageStartupMessages({
  library(depmixS4)
  library(readr)
})

fitted_model <- readRDS(MIMIC_FINAL_MODEL_RDS)
prep <- readRDS(MIMIC_FINAL_PREPROCESSING_RDS)
if (is.null(prep$fullcohort_cache) || !file.exists(prep$fullcohort_cache)) {
  stop("Full-cohort preprocessing cache is unavailable: ", prep$fullcohort_cache)
}
bundle <- readRDS(prep$fullcohort_cache)

response_formulas <- prep$response_formulas
bests <- prep$normalizers
imputed_data_scale_glm_list <- prep$age_sex_reference_glms
imputed_data_scale_df <- as.data.frame(bundle$hmm_data)
imputed_data <- as.data.frame(bundle$imputed_data)
combined <- as.data.frame(bundle$combined)
ids <- as.character(bundle$pat_ids)
if (!"stay_id" %in% names(combined)) combined$stay_id <- ids
if (!"lab_time" %in% names(combined)) combined$lab_time <- bundle$metadata$lab_time

icustays_intubated_first_patient_admission <- readr::read_csv(
  "icustays_intubated_first_patient_admission.csv",
  show_col_types = FALSE
)
icustays_intubated_first_patient_admission$stay_id <-
  as.character(icustays_intubated_first_patient_admission$stay_id)

post_raw <- depmixS4::posterior(fitted_model, type = "smoothing")
expected_state_cols <- paste0("S", seq_len(MIMIC_FINAL_STATE_DIMENSION))
if (all(expected_state_cols %in% colnames(post_raw))) {
  post_smooth_full <- as.matrix(post_raw[, expected_state_cols, drop = FALSE])
} else if (ncol(post_raw) == MIMIC_FINAL_STATE_DIMENSION) {
  post_smooth_full <- as.matrix(post_raw)
  colnames(post_smooth_full) <- expected_state_cols
} else if (
  ncol(post_raw) == MIMIC_FINAL_STATE_DIMENSION + 1L &&
  any(colnames(post_raw) %in% c("state", "State"))
) {
  label_col <- which(colnames(post_raw) %in% c("state", "State"))[1]
  post_smooth_full <- as.matrix(post_raw[, -label_col, drop = FALSE])
  colnames(post_smooth_full) <- expected_state_cols
} else {
  stop("Unexpected posterior output dimensions.")
}

stopifnot(
  nrow(post_smooth_full) == length(ids),
  nrow(combined) == length(ids),
  nrow(imputed_data_scale_df) == length(ids),
  all(is.finite(post_smooth_full)),
  all(abs(rowSums(post_smooth_full) - 1) < 1e-6)
)

message("Loaded final HMM environment: ", MIMIC_FINAL_MODEL_RDS)
message("Rows: ", length(ids), "; stays: ", length(unique(ids)),
        "; states: ", MIMIC_FINAL_STATE_DIMENSION)
