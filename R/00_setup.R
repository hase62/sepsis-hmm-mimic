code_dir <- Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R"))
resource_dir <- Sys.getenv("SEPSIS_HMM_RESOURCE_DIR", unset = file.path(getwd(), "resources"))
output_dir <- Sys.getenv("SEPSIS_HMM_OUTPUT_DIR", unset = getwd())
mimic_root <- Sys.getenv("MIMIC_ROOT", unset = "")

if (!nzchar(mimic_root)) stop("Set MIMIC_ROOT to the MIMIC-IV v3.1 directory containing hosp/ and icu/.", call. = FALSE)
if (!dir.exists(file.path(mimic_root, "hosp")) || !dir.exists(file.path(mimic_root, "icu"))) {
  stop("MIMIC_ROOT must contain hosp/ and icu/.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mimic_path <- function(...) file.path(mimic_root, ...)
resource_path <- function(...) file.path(resource_dir, ...)
code_path <- function(...) file.path(code_dir, ...)
source(file.path(code_dir, "00_packages.R"))

if (!nzchar(Sys.getenv("MIMIC_HMM_VARIABLE_CONFIG", unset = ""))) {
  Sys.setenv(MIMIC_HMM_VARIABLE_CONFIG = resource_path("observation_variables.csv"))
}
if (!nzchar(Sys.getenv("MIMIC_HMM_NORMALIZATION_CONFIG", unset = ""))) {
  Sys.setenv(MIMIC_HMM_NORMALIZATION_CONFIG = resource_path("normalization_families.csv"))
}
if (!nzchar(Sys.getenv("MIMIC_ANNOT_THRESHOLD_FORMAL", unset = ""))) {
  Sys.setenv(MIMIC_ANNOT_THRESHOLD_FORMAL = resource_path("annotation_thresholds.csv"))
}
