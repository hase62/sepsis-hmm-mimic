args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
mimic_root <- if (length(args) >= 1L) args[[1L]] else Sys.getenv("MIMIC_ROOT", unset = "")
if (!nzchar(mimic_root)) stop("Usage: Rscript check_inputs.R <MIMIC_ROOT>", call. = FALSE)
mimic_root <- normalizePath(mimic_root, winslash = "/", mustWork = TRUE)
resource_dir <- file.path(repo, "resources")

required_raw <- c(
  "hosp/admissions.csv.gz",
  "hosp/labevents.csv.gz",
  "hosp/microbiologyevents.csv.gz",
  "hosp/omr.csv.gz",
  "hosp/patients.csv.gz",
  "hosp/prescriptions.csv.gz",
  "icu/chartevents.csv.gz",
  "icu/d_items.csv.gz",
  "icu/icustays.csv.gz",
  "icu/inputevents.csv.gz",
  "icu/outputevents.csv.gz",
  "icu/procedureevents.csv.gz"
)
missing_raw <- required_raw[!file.exists(file.path(mimic_root, required_raw))]
if (length(missing_raw)) {
  stop("Missing MIMIC-IV files:\n", paste0("  - ", missing_raw, collapse = "\n"), call. = FALSE)
}

required_resources <- c(
  "antibiotics_list.txt",
  "mimic_item_map.csv",
  "observation_variables.csv",
  "normalization_families.csv",
  "annotation_thresholds.csv"
)
missing_resources <- required_resources[!file.exists(file.path(resource_dir, required_resources))]
if (length(missing_resources)) {
  stop("Missing repository resources:\n", paste0("  - resources/", missing_resources, collapse = "\n"), call. = FALSE)
}

antibiotics <- trimws(readLines(file.path(resource_dir, "antibiotics_list.txt"), warn = FALSE))
antibiotics <- antibiotics[nzchar(antibiotics)]
if (!length(antibiotics)) stop("resources/antibiotics_list.txt is empty.", call. = FALSE)

item_map <- read.csv(file.path(resource_dir, "mimic_item_map.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("itemid", "matched_test") %in% names(item_map))) {
  stop("resources/mimic_item_map.csv must contain itemid and matched_test.", call. = FALSE)
}
if (!nrow(item_map) || any(!nzchar(trimws(as.character(item_map$matched_test))))) {
  stop("resources/mimic_item_map.csv contains empty mappings.", call. = FALSE)
}

obs <- read.csv(file.path(resource_dir, "observation_variables.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("variable_order", "variable") %in% names(obs))) {
  stop("resources/observation_variables.csv must contain variable_order and variable.", call. = FALSE)
}
obs <- obs[order(as.integer(obs$variable_order)), , drop = FALSE]
if (nrow(obs) != 38L || anyDuplicated(obs$variable) || anyDuplicated(obs$variable_order)) {
  stop("resources/observation_variables.csv must define 38 uniquely ordered variables.", call. = FALSE)
}

norm <- read.csv(file.path(resource_dir, "normalization_families.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("variable_order", "variable", "normalization_family") %in% names(norm))) {
  stop("resources/normalization_families.csv must contain variable_order, variable, and normalization_family.", call. = FALSE)
}
norm <- norm[order(as.integer(norm$variable_order)), , drop = FALSE]
if (!identical(as.character(norm$variable), as.character(obs$variable))) {
  stop("Normalization variables must match observation_variables.csv in the same order.", call. = FALSE)
}

thr <- read.csv(file.path(resource_dir, "annotation_thresholds.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("threshold_name", "threshold") %in% names(thr))) {
  stop("resources/annotation_thresholds.csv must contain threshold_name and threshold.", call. = FALSE)
}
if (!nrow(thr) || anyDuplicated(thr$threshold_name) || any(!is.finite(suppressWarnings(as.numeric(thr$threshold))))) {
  stop("resources/annotation_thresholds.csv contains invalid threshold rows.", call. = FALSE)
}

cat("Input check passed.\n")
cat("MIMIC-IV root:", mimic_root, "\n")
cat("Observation variables:", nrow(obs), "\n")
cat("Antibiotic terms:", length(antibiotics), "\n")
cat("Item mappings:", nrow(item_map), "\n")
cat("Annotation thresholds:", nrow(thr), "\n")
