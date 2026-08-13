
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript 13_select_final_model.R <dimension> [n_starts=20]")
}
dimension <- as.integer(args[[1L]])
n_starts <- if (length(args) >= 2L) as.integer(args[[2L]]) else 20L

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})
source(code_path("06_preprocessing.R"))

out_dir <- Sys.getenv("MIMIC_FULLCOHORT_START_DIR", unset = "full_cohort_starts")
metrics_files <- file.path(out_dir, sprintf("dim%03d_start%02d.metrics.csv", dimension, seq_len(n_starts)))
missing_files <- metrics_files[!file.exists(metrics_files)]
if (length(missing_files) && !get_env_bool("MIMIC_ALLOW_PARTIAL_STARTS", FALSE)) {
  stop("Some requested start metrics are missing:\n", paste(missing_files, collapse = "\n"))
}
metrics_files <- metrics_files[file.exists(metrics_files)]
if (!length(metrics_files)) stop("No full-cohort start metrics were found.")

metrics <- dplyr::bind_rows(lapply(metrics_files, function(f) {
  x <- readr::read_csv(f, show_col_types = FALSE)
  x$metrics_file <- f
  x$model_file <- sub("\\.metrics\\.csv$", ".rds", f)
  x
}))
metrics$model_exists <- file.exists(metrics$model_file)
metrics$eligible <- metrics$valid %in% TRUE & metrics$model_exists & is.finite(metrics$logLik)
readr::write_csv(metrics, sprintf("full_cohort_K%d_start_comparison.csv", dimension))

eligible <- metrics %>% dplyr::filter(eligible) %>% dplyr::arrange(dplyr::desc(logLik))
if (!nrow(eligible)) stop("No valid full-cohort start is available.")
if (nrow(eligible) < 3L) warning("Fewer than three valid starts were available: ", nrow(eligible))

best <- eligible[1L, , drop = FALSE]
final_model <- sprintf("final_hmm_K%d.rds", dimension)
if (!file.copy(best$model_file, final_model, overwrite = TRUE)) {
  stop("Could not copy selected model to: ", final_model)
}

cache_path <- fullcohort_cache_path()
if (!file.exists(cache_path)) {
  stop(
    "Required full-cohort cache is missing: ", cache_path, "\n",
    "Run first: Rscript 11_prepare_full_cohort.R"
  )
}
bundle <- tryCatch(
  readRDS(cache_path),
  error = function(e) stop("Could not read full-cohort cache: ", conditionMessage(e))
)
if (!identical(bundle$version, MIMIC_PREPROCESSING_VERSION) ||
    !identical(bundle$stage, "final_fullcohort")) {
  stop("Full-cohort cache is incompatible with the current analysis request.")
}
if (anyNA(bundle$ntimes) || sum(bundle$ntimes) != nrow(bundle$hmm_data)) {
  stop("Full-cohort cache has invalid sequence lengths.")
}

transition_rhs <- "~ age + sex"
preprocessing_file <- sprintf(
  "final_preprocessing_K%d.rds",
  dimension
)
preprocessing_object <- list(
  version = bundle$version,
  normalizers = bundle$normalizers,
  age_sex_reference_glms = bundle$stable_reference_glms,
  stable_reference_glms = bundle$stable_reference_glms,
  stable_reference = bundle$stable_reference,
  response_formulas = bundle$response_formulas,
  transition_rhs = transition_rhs,
  state_dimension = dimension,
  selected_start_id = best$start_id,
  selected_start_seed = best$seed,
  full_ntimes = bundle$ntimes,
  observation_variable_names = bundle$observation_variable_names,
  fullcohort_cache = cache_path,
  created_at = Sys.time()
)
saveRDS(preprocessing_object, preprocessing_file)

fit_metrics <- data.frame(
  model_stage = "final_full_cohort_estimation",
  state_dimension = dimension,
  selected_start_id = best$start_id,
  selected_seed = best$seed,
  n_requested_starts = n_starts,
  n_valid_starts = nrow(eligible),
  n_patient_stays = length(bundle$ntimes),
  n_rows = sum(bundle$ntimes),
  logLik = best$logLik,
  AIC = best$AIC,
  BIC = best$BIC,
  minimum_state_occupancy = best$minimum_state_occupancy,
  model_rds = final_model,
  preprocessing_rds = preprocessing_file,
  stringsAsFactors = FALSE
)
readr::write_tsv(
  fit_metrics,
  sprintf("final_hmm_K%d_fit_metrics.tsv", dimension)
)

config_lines <- c(
  sprintf('MIMIC_FINAL_MODEL_RDS <- "%s"', final_model),
  sprintf('MIMIC_FINAL_PREPROCESSING_RDS <- "%s"', preprocessing_file),
  sprintf('MIMIC_FINAL_STATE_DIMENSION <- %dL', dimension),
  sprintf('MIMIC_FINAL_START_ID <- %dL', best$start_id),
  sprintf('MIMIC_FINAL_START_SEED <- %dL', best$seed)
)
writeLines(config_lines, "model_config.R")

cat("Selected full-cohort start:\n")
print(
  best[, c("start_id", "seed", "logLik", "AIC", "BIC", "minimum_state_occupancy")],
  row.names = FALSE
)
cat("Final model:", final_model, "\n")
cat("Model configuration: model_config.R\n")
