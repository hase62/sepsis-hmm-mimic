
if (!file.exists("model_config.R")) stop("model_config.R not found.", call. = FALSE)
source("model_config.R")
MIMIC_ACTIVE_ANNOT_VER <- "final"
MIMIC_ACTIVE_ANNOT_MASTER <- "state_annotation.csv"
MIMIC_ACTIVE_REFERENCE_PROFILES <- "state_reference_profiles.csv"
MIMIC_ACTIVE_THRESHOLDS <- "annotation_thresholds_used.csv"
MIMIC_ACTIVE_K <- as.integer(MIMIC_FINAL_STATE_DIMENSION)
MIMIC_FINAL_HMM_RDS <- MIMIC_FINAL_MODEL_RDS

mimic_read_active_annotation <- function(path = MIMIC_ACTIVE_ANNOT_MASTER, expected_k = MIMIC_ACTIVE_K) {
  if (!file.exists(path)) stop("State annotation file not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("state_id", "subphenotype") %in% names(x))) stop("State annotation lacks state_id/subphenotype.", call. = FALSE)
  x$state_id <- suppressWarnings(as.integer(x$state_id))
  if (!identical(sort(unique(x$state_id)), seq_len(as.integer(expected_k)))) stop("State annotation does not match the fitted state dimension.", call. = FALSE)
  x
}

mimic_assert_annotation_matches <- function(annotation_tbl, path = MIMIC_ACTIVE_ANNOT_MASTER, expected_k = MIMIC_ACTIVE_K, object_name = "annotation object") {
  active <- mimic_read_active_annotation(path, expected_k)
  current <- annotation_tbl[, c("state_id", "subphenotype"), drop = FALSE]
  current$state_id <- as.integer(current$state_id)
  current <- current[!duplicated(current$state_id), , drop = FALSE]
  active <- active[!duplicated(active$state_id), c("state_id", "subphenotype"), drop = FALSE]
  current <- current[order(current$state_id), , drop = FALSE]
  active <- active[order(active$state_id), , drop = FALSE]
  if (!identical(current$state_id, active$state_id) || !identical(as.character(current$subphenotype), as.character(active$subphenotype))) stop(object_name, " is inconsistent with state_annotation.csv.", call. = FALSE)
  invisible(TRUE)
}

mimic_promote_final_model <- function(expected_k = MIMIC_ACTIVE_K, allow_rds = TRUE) {
  model <- if (exists("final_fitted_model", envir = .GlobalEnv, inherits = FALSE)) get("final_fitted_model", envir = .GlobalEnv) else if (exists("fitted_model", envir = .GlobalEnv, inherits = FALSE)) get("fitted_model", envir = .GlobalEnv) else if (allow_rds && file.exists(MIMIC_FINAL_MODEL_RDS)) readRDS(MIMIC_FINAL_MODEL_RDS) else NULL
  if (is.null(model)) stop("Final HMM is unavailable.", call. = FALSE)
  k <- tryCatch(length(model@response), error = function(e) NA_integer_)
  if (!is.finite(k) || as.integer(k) != as.integer(expected_k)) stop("Fitted HMM state dimension is inconsistent with model_config.R.", call. = FALSE)
  assign("fitted_model", model, envir = .GlobalEnv)
  invisible(model)
}

for (.f in c(MIMIC_ACTIVE_ANNOT_MASTER, MIMIC_ACTIVE_REFERENCE_PROFILES, MIMIC_ACTIVE_THRESHOLDS)) if (!file.exists(.f)) stop("Required analysis input not found: ", .f, call. = FALSE)
rm(.f)
