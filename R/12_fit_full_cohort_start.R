
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript 12_fit_full_cohort_start.R ",
    "<dimension> <start_id>"
  )
}

dimension <- suppressWarnings(as.integer(args[[1L]]))
start_id <- suppressWarnings(as.integer(args[[2L]]))

if (!is.finite(dimension) || dimension < 2L) stop("Invalid dimension.")
if (!is.finite(start_id) || start_id < 1L) stop("Invalid start_id.")

suppressPackageStartupMessages(library(depmixS4))
source(code_path("06_preprocessing.R"))
source(code_path("07_hmm_fit_helpers.R"))

cache_path <- fullcohort_cache_path()
if (!file.exists(cache_path)) {
  stop(
    "Required full-cohort cache is missing: ", cache_path, "\n",
    "Run first: Rscript 11_prepare_full_cohort.R"
  )
}

bundle <- tryCatch(
  readRDS(cache_path),
  error = function(e) {
    stop("Could not read required full-cohort cache: ", cache_path, "\n", conditionMessage(e))
  }
)

required_fields <- c(
  "version", "stage", "hmm_data", "ntimes",
  "response_formulas", "observation_variable_names",
  "normalization_spec_version", "normalization_families"
)
missing_fields <- setdiff(required_fields, names(bundle))
if (length(missing_fields)) {
  stop("Full-cohort cache is incomplete; missing: ", paste(missing_fields, collapse = ", "))
}
if (!identical(bundle$version, MIMIC_PREPROCESSING_VERSION)) {
  stop("Full-cohort cache preprocessing version does not match the current code.")
}
if (!identical(bundle$normalization_spec_version, MIMIC_NORMALIZATION_SPEC_VERSION) ||
    !identical(
      bundle$normalization_families,
      normalization_family_map(bundle$observation_variable_names)
    )) {
  stop(
    "Full-cohort cache uses a different normalization-family specification. ",
    "Rerun 11_prepare_full_cohort.R."
  )
}
if (!identical(bundle$stage, "final_fullcohort")) {
  stop("Full-cohort cache has an unexpected stage identifier.")
}
if (anyNA(bundle$ntimes) || any(!is.finite(bundle$ntimes))) {
  stop("Full-cohort cache contains invalid sequence lengths.")
}
if (sum(bundle$ntimes) != nrow(bundle$hmm_data)) {
  stop("Full-cohort cache sequence lengths do not match HMM matrix rows.")
}

transition_rhs <- "~ age + sex"
seed <- 310000L + start_id
out_dir <- Sys.getenv(
  "MIMIC_FULLCOHORT_START_DIR",
  unset = "full_cohort_starts"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, sprintf("dim%03d_start%02d", dimension, start_id))
model_file <- paste0(prefix, ".rds")
metrics_file <- paste0(prefix, ".metrics.csv")
failed_file <- paste0(prefix, ".failed")

cat(sprintf(
  "Read full-cohort cache: rows=%d stays=%d; fitting K=%d start=%d seed=%d\n",
  nrow(bundle$hmm_data), length(bundle$ntimes), dimension, start_id, seed
))

set.seed(seed)
model <- depmixS4::depmix(
  response = bundle$response_formulas,
  family = replicate(
    length(bundle$response_formulas),
    stats::gaussian(),
    simplify = FALSE
  ),
  nstates = dimension,
  transition = stats::as.formula(transition_rhs),
  data = bundle$hmm_data,
  ntimes = bundle$ntimes
)

fit_error <- NULL
fitted <- tryCatch(
  depmixS4::fit(
    model,
    emcontrol = depmixS4::em.control(
      maxit = 2000,
      tol = 1e-4,
      random.start = TRUE
    ),
    verbose = TRUE
  ),
  error = function(e) {
    fit_error <<- conditionMessage(e)
    NULL
  }
)

if (is.null(fitted)) {
  readr::write_csv(
    data.frame(
      dimension = dimension,
      start_id = start_id,
      seed = seed,
      valid = FALSE,
      error = fit_error,
      stringsAsFactors = FALSE
    ),
    metrics_file
  )
  writeLines(fit_error %||% "Unknown fitting error", failed_file)
  quit(save = "no", status = 1L)
}

ll <- depmix_logLik_safe(fitted)
ic <- depmix_information_criteria(ll)

post_error <- NULL
post_prob <- tryCatch(
  extract_posterior_probabilities(
    depmixS4::posterior(fitted, type = "smoothing"),
    dimension
  ),
  error = function(e) {
    post_error <<- conditionMessage(e)
    NULL
  }
)
occupancy <- if (is.null(post_prob)) rep(NA_real_, dimension) else colMeans(post_prob)

emission_sds <- unlist(lapply(seq_len(dimension), function(k) {
  unlist(lapply(fitted@response[[k]], function(r) {
    s <- tryCatch(as.numeric(r@parameters$sd), error = function(e) NA_real_)
    if (!length(s)) NA_real_ else s[[1L]]
  }))
}))

finite_parameters <- all(is.finite(depmixS4::getpars(fitted)))
minimum_emission_sd <- if (any(is.finite(emission_sds))) {
  min(emission_sds[is.finite(emission_sds)])
} else {
  NA_real_
}
minimum_state_occupancy <- if (any(is.finite(occupancy))) {
  min(occupancy[is.finite(occupancy)])
} else {
  NA_real_
}
valid <- is.finite(as.numeric(ll)) &&
  finite_parameters &&
  !is.null(post_prob) &&
  all(is.finite(occupancy)) &&
  !(is.finite(minimum_emission_sd) && minimum_emission_sd < 1e-8)

saveRDS(fitted, model_file)
readr::write_csv(
  data.frame(
    dimension = dimension,
    start_id = start_id,
    seed = seed,
    valid = valid,
    logLik = as.numeric(ll),
    logLik_per_row = as.numeric(ll) / nrow(bundle$hmm_data),
    AIC = ic$AIC,
    BIC = ic$BIC,
    minimum_state_occupancy = minimum_state_occupancy,
    number_states_below_1e4 = sum(occupancy < 1e-4, na.rm = TRUE),
    minimum_emission_sd = minimum_emission_sd,
    finite_parameters = finite_parameters,
    posterior_error = post_error %||% "",
    preprocessing_cache = cache_path,
    stringsAsFactors = FALSE
  ),
  metrics_file
)
if (file.exists(failed_file)) unlink(failed_file)

cat(sprintf(
  "Completed full-cohort start %d: valid=%s logLik=%.2f\n",
  start_id, valid, as.numeric(ll)
))
