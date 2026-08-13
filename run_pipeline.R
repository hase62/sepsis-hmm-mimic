args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop(paste(
    "Usage: Rscript run_pipeline.R <command> <MIMIC_ROOT> <OUTPUT_DIR> [command arguments]",
    "Commands:",
    "  check-inputs",
    "  prepare-data",
    "  prepare-split <iteration>",
    "  fit-selection <K> <iteration>",
    "  select-k",
    "  prepare-final",
    "  fit-final <K> <start_id>",
    "  select-final <K> [n_starts]",
    "  annotate",
    "  analyze",
    sep = "\n"
  ), call. = FALSE)
}
command <- args[[1L]]
repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
code_dir <- file.path(repo, "R")
resource_dir <- file.path(repo, "resources")
out <- normalizePath(args[[3L]], winslash = "/", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(
  MIMIC_ROOT = normalizePath(args[[2L]], winslash = "/", mustWork = TRUE),
  SEPSIS_HMM_CODE_DIR = code_dir,
  SEPSIS_HMM_RESOURCE_DIR = resource_dir,
  SEPSIS_HMM_OUTPUT_DIR = out,
  MIMIC_K_MIN = "5",
  MIMIC_K_MAX = "65"
)
setwd(out)
rscript <- file.path(R.home("bin"), "Rscript")
run <- function(script, script_args = character()) {
  status <- system2(rscript, c(file.path(code_dir, script), script_args))
  if (!identical(as.integer(status), 0L)) stop("Failed: ", script, call. = FALSE)
}

if (command == "check-inputs") {
  status <- system2(rscript, c(file.path(repo, "check_inputs.R"), Sys.getenv("MIMIC_ROOT")))
  if (!identical(as.integer(status), 0L)) stop("Input check failed.", call. = FALSE)
} else if (command == "prepare-data") {
  status <- system2(rscript, c(file.path(repo, "check_inputs.R"), Sys.getenv("MIMIC_ROOT")))
  if (!identical(as.integer(status), 0L)) stop("Input check failed.", call. = FALSE)
  run("01_extract_sepsis_cohort.R")
  run("03_build_analytical_matrix.R")
  run("04_analytical_day_qc.R")
  run("05_finalize_hmm_input.R")
} else if (command == "prepare-split") {
  if (length(args) < 4L) stop("prepare-split requires <iteration>.")
  run("08_prepare_model_selection_split.R", args[[4L]])
} else if (command == "fit-selection") {
  if (length(args) < 5L) stop("fit-selection requires <K> <iteration>.")
  run("09_fit_model_selection_candidate.R", c(args[[4L]], args[[5L]]))
} else if (command == "select-k") {
  run("10_select_state_dimension.R", "20")
} else if (command == "prepare-final") {
  run("11_prepare_full_cohort.R")
} else if (command == "fit-final") {
  if (length(args) < 5L) stop("fit-final requires <K> <start_id>.")
  run("12_fit_full_cohort_start.R", c(args[[4L]], args[[5L]]))
} else if (command == "select-final") {
  if (length(args) < 4L) stop("select-final requires <K> [n_starts].")
  n_starts <- if (length(args) >= 5L) args[[5L]] else "20"
  run("13_select_final_model.R", c(args[[4L]], n_starts))
} else if (command == "annotate") {
  run("19_run_annotation.R")
} else if (command == "analyze") {
  run("34_run_downstream_analyses.R")
} else {
  stop("Unknown command: ", command, call. = FALSE)
}
