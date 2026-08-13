
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    paste0(
      "Usage: Rscript 09_fit_model_selection_candidate.R ",
      "<dimension> <iteration> [save_model_rds=0]"
    )
  )
}

dimension <- suppressWarnings(as.integer(args[[1L]]))
iteration <- suppressWarnings(as.integer(args[[2L]]))
save_model_rds <- if (length(args) >= 3L) {
  suppressWarnings(as.integer(args[[3L]])) == 1L
} else {
  identical(Sys.getenv("MIMIC_SAVE_SPLIT_MODEL_RDS", unset = "0"), "1")
}

if (!is.finite(dimension) || dimension < 2L) stop("Invalid dimension.")
if (!is.finite(iteration) || iteration < 1L) stop("Invalid iteration.")
if (is.na(save_model_rds)) stop("save_model_rds must be 0 or 1.")

source(code_path("06_preprocessing.R"))
source(code_path("07_hmm_fit_helpers.R"))

manifest <- read_hmm_finalization_manifest()
cat(sprintf(
  "Using finalized hmm_input cohort: %d stays; matrix=%s\n",
  manifest$n_stays_final,
  manifest$resolved_paths$matrix
))
cat(
  "Split fitted-model RDS: ",
  if (save_model_rds) "ENABLED" else "disabled (metrics only)",
  "\n",
  sep = ""
)

cache_path <- iteration_cache_path(iteration)
if (!file.exists(cache_path)) {
  stop(
    "Required iteration cache is missing: ", cache_path, "\n",
    "Run first: Rscript 08_prepare_model_selection_split.R ",
    iteration
  )
}

bundle <- tryCatch(
  readRDS(cache_path),
  error = function(e) {
    stop("Could not read required iteration cache: ", cache_path, "\n", conditionMessage(e))
  }
)

required_fields <- c(
  "version", "stage", "iteration",
  "train_ids", "val_ids", "train_df", "val_df",
  "train_ntimes", "val_ntimes", "response_formulas",
  "stable_reference", "observation_variable_names",
  "normalization_spec_version", "normalization_families"
)
missing_fields <- setdiff(required_fields, names(bundle))
if (length(missing_fields)) {
  stop("Iteration cache is incomplete; missing: ", paste(missing_fields, collapse = ", "))
}
if (!identical(bundle$version, MIMIC_PREPROCESSING_VERSION)) {
  stop("Iteration cache preprocessing version does not match the current code.")
}
if (!identical(bundle$normalization_spec_version, MIMIC_NORMALIZATION_SPEC_VERSION) ||
    !identical(
      bundle$normalization_families,
      normalization_family_map(bundle$observation_variable_names)
    )) {
  stop(
    "Iteration cache uses a different normalization-family specification. ",
    "Rerun 08_prepare_model_selection_split.R for this iteration."
  )
}
if (!identical(bundle$stage, "model_selection_iteration") ||
    !identical(as.integer(bundle$iteration), as.integer(iteration))) {
  stop("Iteration cache identity does not match the requested iteration.")
}
if (anyNA(bundle$train_ntimes) || anyNA(bundle$val_ntimes) ||
    any(!is.finite(bundle$train_ntimes)) || any(!is.finite(bundle$val_ntimes))) {
  stop("Iteration cache contains invalid sequence lengths.")
}
if (sum(bundle$train_ntimes) != nrow(bundle$train_df) ||
    sum(bundle$val_ntimes) != nrow(bundle$val_df)) {
  stop("Iteration cache sequence lengths do not match matrix row counts.")
}
if (length(intersect(as.character(bundle$train_ids), as.character(bundle$val_ids)))) {
  stop("Iteration cache has overlapping train and validation stay IDs.")
}

cat(sprintf(
  paste0(
    "Read iteration cache: iteration=%d; train=%d stays/%d rows; ",
    "validation=%d stays/%d rows\n"
  ),
  iteration,
  length(bundle$train_ids), nrow(bundle$train_df),
  length(bundle$val_ids), nrow(bundle$val_df)
))

transition_rhs <- "~ age + sex"
prefix <- sprintf("model_selection_K%03d_split%02d", dimension, iteration)

fit_checkpoint_path <- if (save_model_rds) paste0(prefix, ".rds") else NULL
fit_checkpoint_meta_path <- if (save_model_rds) {
  paste0(prefix, ".fit_checkpoint_meta.rds")
} else {
  NULL
}

cache_md5 <- unname(tools::md5sum(cache_path))
fit_seed <- iteration
result <- fit_and_eval_hmm_presplit(
  train_df = bundle$train_df,
  val_df = bundle$val_df,
  train_ntimes = bundle$train_ntimes,
  val_ntimes = bundle$val_ntimes,
  response_formulas = bundle$response_formulas,
  dimension = dimension,
  transition_rhs = transition_rhs,
  seed = fit_seed,
  maxit = 2000,
  tol = 1e-4,
  checkpoint_path = fit_checkpoint_path,
  checkpoint_meta_path = fit_checkpoint_meta_path,
  preprocessing_cache_md5 = cache_md5
)

write_scalar <- function(x, suffix) {
  write.table(
    as.numeric(x),
    paste0(prefix, suffix),
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
}

write_scalar(result$val_mse1, ".mse1.tsv")
write_scalar(result$val_mse2, ".mse2.tsv")
write_scalar(result$val_mse3, ".mse3.tsv")
write_scalar(result$train_logLik, ".logLikelihood.tsv")
write_scalar(result$train_BIC, ".BIC.tsv")
write_scalar(result$train_AIC, ".AIC.tsv")
write_scalar(result$train_ICL_approx, ".ICL_approx.tsv")

readr::write_csv(
  data.frame(
    dimension = dimension,
    iteration = iteration,
    fit_seed = fit_seed,
    save_model_rds = save_model_rds,
    fitted_model_rds = if (save_model_rds) fit_checkpoint_path else NA_character_,
    finalized_hmm_stays = manifest$n_stays_final,
    finalized_stay_ids_md5 = manifest$final_stay_ids_md5,
    n_train_stays = length(bundle$train_ids),
    n_validation_stays = length(bundle$val_ids),
    n_stable_reference = nrow(bundle$stable_reference),
    logLik = as.numeric(result$train_logLik),
    AIC = result$train_AIC,
    BIC = result$train_BIC,
    ICL_approx = result$train_ICL_approx,
    posterior_entropy = result$posterior_entropy,
    mse1 = result$val_mse1,
    mse2 = result$val_mse2,
    mse3 = result$val_mse3,
    preprocessing_cache = cache_path,
    stringsAsFactors = FALSE
  ),
  paste0(prefix, ".metrics.csv")
)

cat(sprintf(
  paste0(
    "Completed K=%d iteration=%d: MSE1=%.6f MSE2=%.6f ",
    "MSE3=%.6f logLik=%.2f; model_rds=%s\n"
  ),
  dimension,
  iteration,
  result$val_mse1,
  result$val_mse2,
  result$val_mse3,
  as.numeric(result$train_logLik),
  if (save_model_rds) fit_checkpoint_path else "not_saved"
))

rm(result, bundle)
invisible(gc())
