
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
suppressPackageStartupMessages({
  library(bestNormalize)
  library(data.table)
  library(dplyr)
  library(mice)
  library(readr)
  library(tidyr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

MIMIC_BASE_PREPROCESSING_VERSION <- "sepsis_hmm_base_preprocessing_v1"
MIMIC_PREPROCESSING_VERSION <- "sepsis_hmm_preprocessing_v1"

get_env_bool <- function(name, default = FALSE) {
  x <- Sys.getenv(name, unset = if (default) "1" else "0")
  tolower(x) %in% c("1", "true", "yes", "y")
}

get_env_num <- function(name, default) {
  x <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (!is.finite(x)) default else x
}

hmm_final_input_paths <- function() {
  list(
    cohort = Sys.getenv(
      "MIMIC_HMM_COHORT_FILE",
      unset = "icustays_intubated_first_patient_admission.HMM_eligible.csv"
    ),
    matrix = Sys.getenv(
      "MIMIC_HMM_PRIMARY_MATRIX_FILE",
      unset = "mimiciv.combined.wide.HMM_raw.csv"
    ),
    manifest = Sys.getenv(
      "MIMIC_HMM_FINALIZATION_MANIFEST",
      unset = file.path("hmm_input_qc", "hmm_input_manifest.rds")
    )
  )
}

read_hmm_finalization_manifest <- function() {
  paths <- hmm_final_input_paths()
  missing <- c(paths$cohort, paths$matrix, paths$manifest)
  missing <- missing[!file.exists(missing)]
  if (length(missing) > 0L) {
    stop(
      "Final pre-HMM input has not been created. Missing:\n",
      paste0("  - ", missing, collapse = "\n"),
      "\nRun 05_finalize_hmm_input.R before HMM preprocessing."
    )
  }
  manifest <- tryCatch(readRDS(paths$manifest), error = function(e) NULL)
  if (!is.list(manifest) || is.null(manifest$final_stay_ids_md5)) {
    stop("Invalid HMM finalization manifest: ", paths$manifest)
  }
  manifest$resolved_paths <- paths
  manifest
}

file_signature <- function(path) {
  if (!file.exists(path)) return(list(path = normalizePath(path, mustWork = FALSE), exists = FALSE))
  z <- file.info(path)
  list(
    path = normalizePath(path, mustWork = TRUE),
    exists = TRUE,
    size = unname(z$size),
    mtime = format(z$mtime, "%Y-%m-%d %H:%M:%OS6 %Z")
  )
}

object_md5 <- function(x) {
  tf <- tempfile(fileext = ".rds")
  on.exit(unlink(tf), add = TRUE)
  saveRDS(x, tf, version = 3)
  unname(tools::md5sum(tf))
}

atomic_save_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(object, tmp, version = 3)
  if (file.exists(path)) unlink(path, force = TRUE)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Could not atomically move RDS into place: ", path)
  }
  invisible(path)
}

with_cache_lock <- function(lock_dir, target_file, code,
                            wait_seconds = 30,
                            timeout_hours = 24,
                            stale_hours = 48) {
  start <- Sys.time()
  repeat {
    if (dir.create(lock_dir, recursive = FALSE, showWarnings = FALSE)) break

    info <- suppressWarnings(file.info(lock_dir))
    if (nrow(info) == 1L && is.finite(as.numeric(info$mtime))) {
      age_hours <- as.numeric(difftime(Sys.time(), info$mtime, units = "hours"))
      if (is.finite(age_hours) && age_hours > stale_hours) {
        warning("Removing stale cache lock: ", lock_dir)
        unlink(lock_dir, recursive = TRUE, force = TRUE)
        next
      }
    }

    if (as.numeric(difftime(Sys.time(), start, units = "hours")) > timeout_hours) {
      stop("Timed out waiting for cache lock: ", lock_dir)
    }
    message("Another job is preparing the shared cache; waiting: ", target_file)
    Sys.sleep(wait_seconds)
  }

  on.exit(unlink(lock_dir, recursive = TRUE, force = TRUE), add = TRUE)
  force(code)
  invisible(TRUE)
}

read_frozen_hmm_variables <- function(path = Sys.getenv(
  "MIMIC_HMM_VARIABLE_CONFIG",
  unset = "observation_variables.csv"
)) {
  if (!file.exists(path)) {
    stop(
      "Frozen HMM-variable file not found: ", path, "\n",
      "Bundle resources/observation_variables.csv with the repository."
    )
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("variable_order", "variable") %in% names(x))) {
    stop("Frozen variable file must contain variable_order and variable columns: ", path)
  }
  x <- x[order(x$variable_order), , drop = FALSE]
  vars <- as.character(x$variable)
  if (!length(vars) || anyNA(vars) || any(!nzchar(vars)) || anyDuplicated(vars)) {
    stop("Frozen HMM-variable file contains invalid variable names.")
  }
  list(vars = vars, path = path, md5 = unname(tools::md5sum(path)))
}

read_frozen_normalization_families <- function(path = Sys.getenv(
  "MIMIC_HMM_NORMALIZATION_CONFIG",
  unset = "normalization_families.csv"
)) {
  if (!file.exists(path)) {
    stop("Frozen normalization-family file not found: ", path)
  }
  x <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("variable_order", "variable", "normalization_family")
  missing <- setdiff(required, names(x))
  if (length(missing)) {
    stop(
      "Frozen normalization-family file lacks: ",
      paste(missing, collapse = ", "), ": ", path
    )
  }
  x <- x[order(x$variable_order), required, drop = FALSE]
  x$variable <- as.character(x$variable)
  x$normalization_family <- as.character(x$normalization_family)
  allowed <- c("orderNorm", "lambert_s", "log10", "sqrt", "arcsinh", "zscore")
  if (!nrow(x) || anyNA(x$variable) || any(!nzchar(x$variable)) ||
      anyDuplicated(x$variable)) {
    stop("Frozen normalization-family file contains invalid variable names: ", path)
  }
  bad <- setdiff(unique(x$normalization_family), allowed)
  if (length(bad)) {
    stop("Unknown frozen normalization families: ", paste(bad, collapse = ", "))
  }
  families <- x$normalization_family
  names(families) <- x$variable
  list(
    table = x,
    families = families,
    path = path,
    md5 = unname(tools::md5sum(path))
  )
}

MIMIC_NORMALIZATION_CONFIG <- read_frozen_normalization_families()
MIMIC_NORMALIZATION_FAMILIES <- MIMIC_NORMALIZATION_CONFIG$families
MIMIC_NORMALIZATION_SPEC_VERSION <- paste0(
  "normalization-family-v1-md5-",
  MIMIC_NORMALIZATION_CONFIG$md5
)

normalization_family_map <- function(variables = names(MIMIC_NORMALIZATION_FAMILIES)) {
  variables <- as.character(variables)
  missing <- setdiff(variables, names(MIMIC_NORMALIZATION_FAMILIES))
  extra <- setdiff(names(MIMIC_NORMALIZATION_FAMILIES), variables)
  if (length(missing)) {
    stop("No frozen normalization family for: ", paste(missing, collapse = ", "))
  }
  if (length(variables) == length(MIMIC_NORMALIZATION_FAMILIES) && length(extra)) {
    stop("Normalization-family file contains variables outside the frozen HMM set: ",
         paste(extra, collapse = ", "))
  }
  families <- unname(MIMIC_NORMALIZATION_FAMILIES[variables])
  names(families) <- variables
  families
}

qm_align <- function(primary, secondary) {
  primary <- primary[is.finite(primary)]
  secondary <- secondary[is.finite(secondary)]
  if (length(primary) < 10L || length(secondary) < 10L) return(function(z) z)
  F_sec <- stats::ecdf(secondary)
  Q_pri <- function(p) as.numeric(stats::quantile(primary, probs = p, type = 7, na.rm = TRUE))
  function(z) Q_pri(F_sec(z))
}

lm_align <- function(primary, secondary) {
  df <- data.frame(y = primary, x = secondary)
  df <- df[is.finite(df$y) & is.finite(df$x), , drop = FALSE]
  if (nrow(df) < 10L) return(function(z) z)
  fit <- stats::lm(y ~ x, data = df)
  function(z) as.numeric(stats::coef(fit)[1] + stats::coef(fit)[2] * z)
}

harmonize_and_fuse_one <- function(df, spec,
                                   prefer = c("quantile", "lm"),
                                   min_overlap = 30,
                                   min_cor = 0.5,
                                   clip_q = c(0.001, 0.999)) {
  prefer <- match.arg(prefer)
  p <- spec$primary
  s <- spec$secondaries %||% character(0)
  out <- spec$out %||% p
  present <- c(p, s)[c(p, s) %in% names(df)]
  if (!length(present)) {
    df[[out]] <- NA_real_
    return(list(df = df, info = data.frame()))
  }
  used_primary <- if (p %in% names(df)) p else present[1]
  secondaries <- setdiff(present, used_primary)
  P <- as.numeric(df[[used_primary]])
  fused <- P
  info_rows <- list()

  for (sec in secondaries) {
    S <- as.numeric(df[[sec]])
    idx <- is.finite(P) & is.finite(S)
    n_both <- sum(idx)
    r_both <- if (n_both >= 2L) suppressWarnings(stats::cor(P[idx], S[idx])) else NA_real_
    if (n_both >= min_overlap && is.finite(r_both) && abs(r_both) >= min_cor) {
      aligner <- if (prefer == "quantile") qm_align(P[idx], S[idx]) else lm_align(P[idx], S[idx])
      method <- prefer
    } else if (n_both >= 10L) {
      aligner <- lm_align(P[idx], S[idx])
      method <- "lm_fallback"
    } else {
      aligner <- function(z) z
      method <- "identity"
    }
    S_aligned <- aligner(S)
    if (sum(is.finite(P)) >= 50L) {
      lo <- as.numeric(stats::quantile(P, clip_q[1], na.rm = TRUE))
      hi <- as.numeric(stats::quantile(P, clip_q[2], na.rm = TRUE))
      S_aligned <- pmin(pmax(S_aligned, lo), hi)
    }
    fused <- dplyr::coalesce(fused, S_aligned)
    info_rows[[length(info_rows) + 1L]] <- data.frame(
      out = out,
      used_primary = used_primary,
      used_secondary = sec,
      overlap = n_both,
      cor = r_both,
      method = method,
      stringsAsFactors = FALSE
    )
  }
  df[[out]] <- fused
  list(df = df, info = dplyr::bind_rows(info_rows))
}

make_spec_if_present <- function(cols, primary, secondaries, out) {
  present <- intersect(c(primary, secondaries), cols)
  if (!length(present)) return(NULL)
  list(
    primary = if (primary %in% present) primary else present[1],
    secondaries = setdiff(present, primary),
    out = out
  )
}

build_specs_from_columns <- function(cols) {
  specs <- list(
    make_spec_if_present(cols, "Lactate", c("cLac."), "Lactate_final"),
    make_spec_if_present(cols, "PH", c("cPH"), "pH_final"),
    make_spec_if_present(cols, "CO2", c("cCO2"), "CO2_final"),
    make_spec_if_present(cols, "pO2", c("pO2(T)", "pO2_FO2I"), "PaO2_final"),
    make_spec_if_present(cols, "Sodium", c("cNa."), "Sodium_final"),
    make_spec_if_present(cols, "Chlor", c("Chlor_whole", "cCl."), "Chloride_final"),
    make_spec_if_present(cols, "Glucose", c("Glucose_whole", "cGlu."), "Glucose_final"),
    make_spec_if_present(cols, "Creatinine", c("Creatinine_whole", "cCrea"), "Creatinine_final"),
    make_spec_if_present(cols, "Hemoglobin", c("cHem."), "Hemoglobin_final"),
    make_spec_if_present(cols, "MEAN", c("ART.M."), "MAP_final"),
    make_spec_if_present(cols, "SYS.", c("ART.S."), "SYS_final"),
    make_spec_if_present(cols, "DIAS.", c("ART.D."), "DIAS_final")
  )
  Filter(Negate(is.null), specs)
}

harmonize_and_fuse_many <- function(df, specs, ...) {
  infos <- list()
  for (i in seq_along(specs)) {
    z <- harmonize_and_fuse_one(df, specs[[i]], ...)
    df <- z$df
    infos[[i]] <- z$info
  }
  list(df = df, info = dplyr::bind_rows(infos))
}

clean_after_fusion_minloss <- function(df, specs) {
  for (sp in specs) {
    out <- sp$out %||% sp$primary
    if (!out %in% names(df)) next
    base <- sub("_final$", "", out)
    if (base != out) {
      if (base %in% names(df)) {
        df[[base]] <- dplyr::coalesce(df[[base]], df[[out]])
        df[[out]] <- NULL
      } else {
        names(df)[match(out, names(df))] <- base
      }
      out <- base
    }
    drop_cols <- intersect(setdiff(c(sp$primary, sp$secondaries), out), names(df))
    if (length(drop_cols)) df <- df[, setdiff(names(df), drop_cols), drop = FALSE]
    meta_cols <- names(df)[startsWith(names(df), paste0(out, "_")) & names(df) != out]
    final_meta <- names(df)[grepl("_final_", names(df), fixed = TRUE)]
    drop_meta <- unique(c(meta_cols, final_meta))
    if (length(drop_meta)) df <- df[, setdiff(names(df), drop_meta), drop = FALSE]
  }
  names(df) <- sub("_final$", "", names(df))
  df
}

fill_short_gaps <- function(x, max_run = 1, time = seq_along(x)) {
  n <- length(x)
  if (!n || !anyNA(x)) return(x)
  r <- rle(is.na(x))
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L
  for (i in seq_along(r$values)) {
    if (!r$values[i]) next
    s <- starts[i]
    e <- ends[i]
    if ((e - s + 1L) > max_run) next
    left_ok <- s > 1L && is.finite(x[s - 1L])
    right_ok <- e < n && is.finite(x[e + 1L])
    if (!(left_ok && right_ok)) next
    x[s:e] <- stats::approx(
      x = time[c(s - 1L, e + 1L)],
      y = x[c(s - 1L, e + 1L)],
      xout = time[s:e]
    )$y
  }
  x
}

fit_one_frozen_normalizer <- function(x_train, family, variable) {
  x_train <- as.numeric(x_train)
  if (!length(x_train) || any(!is.finite(x_train))) {
    stop("Non-finite training values before normalization: ", variable)
  }
  if (length(unique(x_train)) < 2L || !is.finite(stats::sd(x_train)) ||
      stats::sd(x_train) <= 0) {
    stop("Training variable is constant and cannot be normalized: ", variable)
  }

  obj <- switch(
    family,
    zscore = list(
      type = "zscore",
      family = "zscore",
      mean = mean(x_train),
      sd = stats::sd(x_train)
    ),
    orderNorm = bestNormalize::orderNorm(
      x_train,
      standardize = TRUE,
      warn = FALSE
    ),
    lambert_s = bestNormalize::lambert(
      x_train,
      type = "s",
      standardize = TRUE,
      warn = FALSE
    ),
    log10 = bestNormalize::log_x(
      x_train,
      b = 10,
      standardize = TRUE,
      warn = FALSE
    ),
    sqrt = bestNormalize::sqrt_x(
      x_train,
      standardize = TRUE
    ),
    arcsinh = bestNormalize::arcsinh_x(
      x_train,
      standardize = TRUE
    ),
    stop("Unknown frozen normalization family for ", variable, ": ", family)
  )

  if (!identical(family, "zscore")) {
    attr(obj, "mimic_normalization_family") <- family
    z_train <- if (inherits(obj, "orderNorm")) {
      suppressWarnings(as.numeric(stats::predict(obj, newdata = x_train, warn = FALSE)))
    } else {
      suppressWarnings(as.numeric(stats::predict(obj, newdata = x_train)))
    }
    if (length(z_train) != length(x_train) || any(!is.finite(z_train))) {
      stop(
        "Frozen normalizer produced non-finite training values: ", variable,
        " [", family, "]"
      )
    }
  }
  obj
}

fit_normalizers <- function(df, train_rows) {
  stopifnot(all(vapply(df, is.numeric, logical(1))))
  if (length(train_rows) != nrow(df) || anyNA(train_rows)) {
    stop("Invalid train_rows for normalization.")
  }
  families <- normalization_family_map(names(df))
  out <- setNames(vector("list", ncol(df)), names(df))
  for (v in names(df)) {
    out[[v]] <- fit_one_frozen_normalizer(
      x_train = df[[v]][train_rows],
      family = families[[v]],
      variable = v
    )
  }
  out
}

normalizer_family <- function(x) {
  if (!is.null(x$type) && identical(x$type, "zscore")) return("zscore")
  fam <- attr(x, "mimic_normalization_family")
  if (!is.null(fam)) return(as.character(fam))
  class(x)[1L]
}

transform_with_normalizers <- function(df, norms) {
  missing <- setdiff(names(norms), names(df))
  if (length(missing)) {
    stop("Normalization input is missing variables: ", paste(missing, collapse = ", "))
  }
  z <- lapply(names(norms), function(v) {
    b <- norms[[v]]
    x <- df[[v]]
    out <- if (!is.null(b$type) && identical(b$type, "zscore")) {
      (x - b$mean) / b$sd
    } else if (inherits(b, "orderNorm")) {
      suppressWarnings(as.numeric(stats::predict(b, newdata = x, warn = FALSE)))
    } else {
      suppressWarnings(as.numeric(stats::predict(b, newdata = x)))
    }
    bad <- which(!is.finite(out))
    if (length(bad)) {
      stop(
        "Frozen normalization produced non-finite values: ", v,
        " [", normalizer_family(b), "]; n=", length(bad),
        "; first rows=", paste(utils::head(bad, 10L), collapse = ",")
      )
    }
    out
  })
  names(z) <- names(norms)
  as.data.frame(z, check.names = FALSE)
}

clip_ranges <- function(df) {
  ranges <- list(
    HR = c(20, 300), MAP = c(20, 200), SYS = c(40, 300), DIAS = c(20, 200),
    SpO2 = c(30, 100), RR = c(5, 120), Temperature = c(25, 45),
    Glucose = c(20, 2000), Sodium = c(110, 190), Chloride = c(60, 150),
    UreaNitrogen = c(2, 300), Creatinine = c(0.1, 20), HCO3 = c(5, 60),
    Potassium = c(1, 12), Magnesium = c(0.2, 10), Phosphorous = c(0.5, 50),
    AnionGap = c(0, 50), FiO2 = c(0, 100), WBC = c(1, 500),
    Hemoglobin = c(2, 25), Hematocrit = c(10, 70), PlateletCount = c(5, 2000),
    RBC = c(0.5, 10), MCV = c(50, 200), MCH = c(15, 50), MCHC = c(15, 70),
    pH = c(6.8, 7.8), pCO2 = c(10, 120), PaO2 = c(20, 500), CO2 = c(5, 60),
    cBase = c(-30, 30), Lactate = c(0, 20), ALT = c(0, 10000), AST = c(0, 10000),
    TotalBilirubin = c(0, 50), INR = c(0.7, 10), PTT = c(10, 500),
    CalciumIonized = c(0.2, 2.5)
  )
  for (v in intersect(names(ranges), names(df))) {
    lo <- ranges[[v]][1]
    hi <- ranges[[v]][2]
    x <- suppressWarnings(as.numeric(df[[v]]))
    x[x < lo | x > hi] <- NA_real_
    df[[v]] <- x
  }
  df
}

standardize_fio2 <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- ifelse(is.finite(x) & x > 1.5, x / 100, x)
  x[!is.finite(x) | x < 0.21 | x > 1] <- NA_real_
  x
}

derive_pf <- function(df) {
  if (!all(c("PaO2", "FiO2") %in% names(df))) stop("PaO2 and FiO2 are required to derive PF.")
  fio2 <- standardize_fio2(df$FiO2)
  df$PF <- as.numeric(df$PaO2) / fio2
  df
}

build_base_bundle <- function(variable_config = Sys.getenv(
                                "MIMIC_HMM_VARIABLE_CONFIG",
                                unset = "observation_variables.csv"
                              )) {
  frozen <- read_frozen_hmm_variables(variable_config)
  cohort_finalization <- read_hmm_finalization_manifest()
  final_paths <- cohort_finalization$resolved_paths
  input_file <- final_paths$matrix
  cohort_file <- final_paths$cohort

  if (!file.exists(input_file)) stop("Input matrix not found: ", input_file)
  if (!file.exists(cohort_file)) stop("Finalized HMM cohort metadata not found: ", cohort_file)

  message("Reading base matrix: ", input_file)
  combined <- readr::read_csv(input_file, show_col_types = FALSE, progress = FALSE)
  combined$stay_id <- as.character(combined$stay_id)

  cohort <- readr::read_csv(cohort_file, show_col_types = FALSE, progress = FALSE)
  cohort$stay_id <- as.character(cohort$stay_id)
  cohort$age <- suppressWarnings(as.numeric(cohort$age))

  if (any(!is.finite(cohort$age) | cohort$age < 18 | cohort$age >= 80)) {
    stop("Finalized HMM cohort contains an age-ineligible stay.")
  }
  if (anyDuplicated(cohort$stay_id)) {
    stop("Finalized HMM cohort contains duplicated stay_id values.")
  }

  count_per_stay <- table(combined$stay_id)
  cohort_ids <- sort(unique(cohort$stay_id))
  matrix_ids <- sort(names(count_per_stay))
  extra_matrix_ids <- setdiff(matrix_ids, cohort_ids)
  missing_matrix_ids <- setdiff(cohort_ids, matrix_ids)
  short_sequence_ids <- names(count_per_stay[count_per_stay < 2L])
  if (length(extra_matrix_ids) || length(missing_matrix_ids) || length(short_sequence_ids)) {
    stop(
      "hmm_input outputs are inconsistent: ",
      length(extra_matrix_ids), " extra matrix stay(s), ",
      length(missing_matrix_ids), " missing matrix stay(s), and ",
      length(short_sequence_ids), " stay(s) with fewer than two rows."
    )
  }
  if (all(c("Hematocrit", "cHct.") %in% names(combined))) {
    idx <- is.finite(combined$Hematocrit) & is.finite(combined$`cHct.`)
    if (sum(idx) >= 30L) {
      aligner <- qm_align(combined$Hematocrit[idx], combined$`cHct.`[idx])
      h2 <- aligner(combined$`cHct.`)
      combined$Hematocrit <- dplyr::coalesce(combined$Hematocrit, h2)
    }
  }
  if ("Potassium" %in% names(combined)) {
    for (sec in intersect(c("Potassium_whole", "cK."), names(combined))) {
      idx <- is.finite(combined$Potassium) & is.finite(combined[[sec]])
      if (sum(idx) >= 30L) {
        aligner <- qm_align(combined$Potassium[idx], combined[[sec]][idx])
        combined$Potassium <- dplyr::coalesce(combined$Potassium, aligner(combined[[sec]]))
      }
    }
  }

  specs <- build_specs_from_columns(names(combined))
  harm <- harmonize_and_fuse_many(
    combined, specs,
    prefer = "quantile", min_overlap = 30, min_cor = 0.5,
    clip_q = c(0.001, 0.999)
  )
  combined <- clean_after_fusion_minloss(harm$df, specs)

  if ("INR.PT." %in% names(combined) && !"INR" %in% names(combined)) {
    names(combined)[names(combined) == "INR.PT."] <- "INR"
  }
  if ("INR(PT)" %in% names(combined) && !"INR" %in% names(combined)) {
    names(combined)[names(combined) == "INR(PT)"] <- "INR"
  }

  frozen_vars <- frozen$vars
  pre_mice_vars <- frozen_vars
  if ("PF" %in% pre_mice_vars) {
    pre_mice_vars <- unique(c(setdiff(pre_mice_vars, "PF"), "PaO2", "FiO2"))
  }
  stable_aux <- c(
    "MAP", "Lactate", "pH", "SpO2", "FiO2", "PaO2", "Creatinine",
    "TotalBilirubin", "PlateletCount", "INR", "HCO3",
    "Dialysis", "IABP", "ECMO"
  )
  required_pre <- unique(pre_mice_vars)
  missing_pre <- setdiff(required_pre, names(combined))
  if (length(missing_pre)) {
    stop("Required frozen/pre-MICE variables absent after rebuilding: ", paste(missing_pre, collapse = ", "))
  }

  metadata_cols <- intersect(c("stay_id", "hadm_id", "charttime_year_day", "lab_time"), names(combined))
  keep <- unique(c(metadata_cols, required_pre, intersect(stable_aux, names(combined))))
  combined <- combined %>% dplyr::select(dplyr::all_of(keep))
  combined <- clip_ranges(combined)

  if ("charttime_year_day" %in% names(combined)) {
    combined <- combined %>%
      dplyr::arrange(stay_id, charttime_year_day) %>%
      dplyr::group_by(stay_id) %>%
      dplyr::mutate(lab_time = as.numeric(charttime_year_day - min(charttime_year_day)) + 1) %>%
      dplyr::ungroup()
  } else if ("lab_time" %in% names(combined)) {
    combined <- combined %>% dplyr::arrange(stay_id, lab_time)
  } else {
    stop("Neither charttime_year_day nor lab_time is available.")
  }

  combined <- combined %>% dplyr::arrange(stay_id, lab_time)
  reference_raw <- combined

  imputation_source <- combined
  imputation_source <- imputation_source %>%
    dplyr::group_by(stay_id) %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(required_pre),
      ~ fill_short_gaps(.x, max_run = 1, time = lab_time)
    )) %>%
    dplyr::ungroup()

  cohort_used <- cohort %>% dplyr::filter(stay_id %in% unique(combined$stay_id))
  if (anyDuplicated(cohort_used$stay_id)) {
    cohort_used <- cohort_used %>%
      dplyr::arrange(stay_id) %>%
      dplyr::distinct(stay_id, .keep_all = TRUE)
  }

  list(
    version = MIMIC_BASE_PREPROCESSING_VERSION,
    input_file = input_file,
    input_signature = file_signature(input_file),
    cohort_signature = file_signature(cohort_file),
    cohort_finalization_manifest_signature = file_signature(
      final_paths$manifest
    ),
    variable_config = frozen,
    frozen_vars = frozen_vars,
    pre_mice_vars = required_pre,
    reference_raw = reference_raw,
    imputation_source = imputation_source,
    cohort = cohort_used,
    cohort_finalization = cohort_finalization,
    harmonization_info = harm$info,
    created_at = Sys.time()
  )
}

base_cache_path <- function() {
  cache_dir <- Sys.getenv("MIMIC_HMM_CACHE_DIR", unset = "preprocessing_cache")
  file.path(cache_dir, "base.rds")
}

base_signature_expected <- function(variable_config) {
  frozen <- read_frozen_hmm_variables(variable_config)
  cohort_finalization <- read_hmm_finalization_manifest()
  final_paths <- cohort_finalization$resolved_paths
  list(
    version = MIMIC_BASE_PREPROCESSING_VERSION,
    input = file_signature(final_paths$matrix),
    cohort = file_signature(final_paths$cohort),
    cohort_finalization_manifest = file_signature(final_paths$manifest),
    final_stay_ids_md5 = cohort_finalization$final_stay_ids_md5,
    variable_md5 = frozen$md5
  )
}

base_cache_compatible <- function(x, expected) {
  identical(x$version, expected$version) &&
    identical(x$input_signature, expected$input) &&
    identical(x$cohort_signature, expected$cohort) &&
    identical(
      x$cohort_finalization_manifest_signature,
      expected$cohort_finalization_manifest
    ) &&
    identical(x$cohort_finalization$final_stay_ids_md5, expected$final_stay_ids_md5) &&
    identical(x$variable_config$md5, expected$variable_md5)
}

ensure_base_bundle <- function(variable_config = Sys.getenv(
                               "MIMIC_HMM_VARIABLE_CONFIG",
                               unset = "observation_variables.csv"
                             ),
                             force = get_env_bool("MIMIC_FORCE_PREPROCESSING", FALSE)) {
  path <- base_cache_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  expected <- base_signature_expected(variable_config)
  if (!force && file.exists(path)) {
    old <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(old) && base_cache_compatible(old, expected)) return(old)
  }

  lock <- paste0(path, ".lock")
  with_cache_lock(lock, path, {
    old <- NULL
    if (!force && file.exists(path)) {
      old <- tryCatch(readRDS(path), error = function(e) NULL)
    }
    if (is.null(old) || !base_cache_compatible(old, expected)) {
      z <- build_base_bundle(variable_config = variable_config)
      atomic_save_rds(z, path)
    }
  })

  out <- tryCatch(readRDS(path), error = function(e) {
    stop("Could not read prepared base cache: ", path, "\n", conditionMessage(e))
  })
  if (!base_cache_compatible(out, expected)) {
    stop("Prepared base cache is incompatible: ", path)
  }
  out
}

mean_mice_completions <- function(mids_object, m = mids_object$m) {
  completed <- lapply(seq_len(m), function(i) as.data.frame(mice::complete(mids_object, action = i)))
  out <- completed[[1]]
  for (j in seq_along(out)) {
    vals <- vapply(completed, function(z) z[[j]], numeric(nrow(out)))
    out[[j]] <- rowMeans(vals)
  }
  out
}

make_observation_matrix <- function(imputed_pre, frozen_vars) {
  x <- derive_pf(imputed_pre)
  missing_final <- setdiff(frozen_vars, names(x))
  if (length(missing_final)) stop("Frozen HMM variables missing after PF derivation: ", paste(missing_final, collapse = ", "))
  x <- x[, frozen_vars, drop = FALSE]
  if (any(!vapply(x, is.numeric, logical(1)))) stop("Non-numeric observation variable after imputation.")
  x
}

row_metadata <- function(base) {
  m <- base$imputation_source %>% dplyr::select(stay_id, lab_time)
  cmeta <- base$cohort %>% dplyr::select(dplyr::any_of(c("stay_id", "age", "sex", "dod")))
  m <- m %>% dplyr::left_join(cmeta, by = "stay_id")
  m$sex <- ifelse(as.character(m$sex) == "M", 0, 1)
  if (any(!is.finite(m$age)) || any(m$age <= 0) || any(!is.finite(m$sex))) {
    stop("Invalid age/sex metadata after cohort join.")
  }
  m
}

normal_indicator <- function(x, lo = -Inf, hi = Inf) {
  ifelse(is.finite(x), as.numeric(x >= lo & x <= hi), NA_real_)
}

prepare_stable_scoring_table <- function(reference_raw) {
  x <- reference_raw
  n <- nrow(x)
  getv <- function(v) if (v %in% names(x)) suppressWarnings(as.numeric(x[[v]])) else rep(NA_real_, n)
  map <- getv("MAP")
  lac <- getv("Lactate")
  ph <- getv("pH")
  spo2 <- getv("SpO2")
  fio2 <- standardize_fio2(getv("FiO2"))
  pao2 <- getv("PaO2")
  pf <- pao2 / fio2
  cr <- getv("Creatinine")
  bili <- getv("TotalBilirubin")
  plt <- getv("PlateletCount")
  inr <- getv("INR")
  hco3 <- getv("HCO3")

  support_cols <- intersect(c("Dialysis", "IABP", "ECMO"), names(x))
  support_any <- if (length(support_cols)) {
    rowSums(as.data.frame(lapply(x[support_cols], function(z) as.numeric(z == 1))), na.rm = TRUE) > 0
  } else {
    rep(FALSE, n)
  }

  indicators <- data.frame(
    map_ok = normal_indicator(map, 65, 105),
    lactate_ok = normal_indicator(lac, -Inf, 2.0),
    ph_ok = normal_indicator(ph, 7.30, 7.50),
    spo2_ok = normal_indicator(spo2, 92, 100),
    fio2_ok = normal_indicator(fio2, 0.21, 0.40),
    pf_ok = normal_indicator(pf, 200, Inf),
    creatinine_ok = normal_indicator(cr, -Inf, 1.5),
    bilirubin_ok = normal_indicator(bili, -Inf, 2.0),
    platelet_ok = normal_indicator(plt, 100, Inf),
    inr_ok = normal_indicator(inr, -Inf, 1.5),
    hco3_ok = normal_indicator(hco3, 18, 30)
  )
  n_measured <- rowSums(!is.na(indicators))
  score <- rowSums(indicators, na.rm = TRUE) + 0.15 * n_measured - 6 * as.numeric(support_any)

  severe <- support_any |
    (is.finite(map) & map < 55) |
    (is.finite(lac) & lac > 4) |
    (is.finite(ph) & (ph < 7.20 | ph > 7.60)) |
    (is.finite(pf) & pf < 100) |
    (is.finite(cr) & cr > 3.5) |
    (is.finite(bili) & bili > 6) |
    (is.finite(plt) & plt < 50)

  tier1 <- !severe & !support_any & is.finite(map) & map >= 65 &
    is.finite(lac) & lac <= 2.0 & (!is.finite(fio2) | fio2 <= 0.40) & n_measured >= 5
  tier2 <- !severe & !support_any & is.finite(map) & map >= 65 &
    (!is.finite(lac) | lac <= 3.0) & (!is.finite(fio2) | fio2 <= 0.50) & n_measured >= 4
  tier3 <- !severe & !support_any & (!is.finite(map) | map >= 60) &
    (!is.finite(lac) | lac <= 4.0) & n_measured >= 3
  tier4 <- !severe & !support_any & n_measured >= 2

  data.frame(
    row_index = seq_len(n),
    stay_id = as.character(x$stay_id),
    lab_time = as.numeric(x$lab_time),
    stability_score = score,
    n_stability_measures = n_measured,
    support_any = support_any,
    severe_instability = severe,
    tier1 = tier1,
    tier2 = tier2,
    tier3 = tier3,
    tier4 = tier4,
    MAP = map, Lactate = lac, pH = ph, SpO2 = spo2, FiO2 = fio2, PF = pf,
    Creatinine = cr, TotalBilirubin = bili, PlateletCount = plt,
    INR = inr, HCO3 = hco3,
    stringsAsFactors = FALSE
  )
}

select_stable_reference <- function(base, allowed_ids, target_n) {
  score_tab <- prepare_stable_scoring_table(base$reference_raw)
  score_tab <- score_tab[score_tab$stay_id %in% allowed_ids, , drop = FALSE]
  target_n <- min(as.integer(target_n), length(unique(score_tab$stay_id)))
  if (target_n < 20L) stop("Stable-reference target is too small: ", target_n)

  tier_names <- c("tier1", "tier2", "tier3", "tier4")
  chosen_tier <- tail(tier_names, 1)
  candidate <- NULL
  for (tt in tier_names) {
    z <- score_tab[score_tab[[tt]], , drop = FALSE]
    if (length(unique(z$stay_id)) >= target_n) {
      candidate <- z
      chosen_tier <- tt
      break
    }
  }
  if (is.null(candidate)) candidate <- score_tab[score_tab$tier4, , drop = FALSE]

  best_per_patient <- candidate %>%
    dplyr::arrange(stay_id, dplyr::desc(stability_score), dplyr::desc(n_stability_measures), lab_time) %>%
    dplyr::group_by(stay_id) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()

  if (nrow(best_per_patient) < target_n) {
    already <- best_per_patient$stay_id
    fallback <- score_tab %>%
      dplyr::filter(!stay_id %in% already) %>%
      dplyr::arrange(stay_id, dplyr::desc(stability_score), dplyr::desc(n_stability_measures), lab_time) %>%
      dplyr::group_by(stay_id) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
    fallback$reference_selection_tier <- "fallback_best_available"
    best_per_patient$reference_selection_tier <- chosen_tier
    best_per_patient <- dplyr::bind_rows(best_per_patient, fallback)
  } else {
    best_per_patient$reference_selection_tier <- chosen_tier
  }

  selected <- best_per_patient %>%
    dplyr::arrange(dplyr::desc(stability_score), dplyr::desc(n_stability_measures), lab_time, stay_id) %>%
    dplyr::slice_head(n = target_n)

  selected$selected_for_reference <- TRUE
  selected
}

reference_selection_summary <- function(ref) {
  data.frame(
    n_reference_rows = nrow(ref),
    n_unique_stays = length(unique(ref$stay_id)),
    median_stability_score = median(ref$stability_score, na.rm = TRUE),
    median_n_stability_measures = median(ref$n_stability_measures, na.rm = TRUE),
    fallback_rows = sum(ref$reference_selection_tier == "fallback_best_available", na.rm = TRUE),
    selection_tiers = paste(names(table(ref$reference_selection_tier)), as.integer(table(ref$reference_selection_tier)), sep = ":", collapse = ";"),
    stringsAsFactors = FALSE
  )
}

fit_demographic_glms <- function(Z, reference_rows, meta) {
  reference_rows <- unique(as.integer(reference_rows))
  reference_rows <- reference_rows[reference_rows >= 1L & reference_rows <= nrow(Z)]
  if (length(reference_rows) < 50L) stop("Too few demographic reference rows: ", length(reference_rows))
  age_ref <- meta$age[reference_rows]
  sex_ref <- meta$sex[reference_rows]
  out <- vector("list", ncol(Z))
  names(out) <- names(Z)
  for (j in seq_along(Z)) {
    y <- Z[[j]][reference_rows]
    out[[j]] <- stats::glm(y ~ age_ref + log(age_ref) + sex_ref)
  }
  out
}

apply_demographic_glms <- function(Z, glms, meta) {
  out <- Z
  for (j in seq_along(Z)) {
    cf <- stats::coef(glms[[j]])
    if (length(cf) != 4L || any(!is.finite(cf))) stop("Invalid demographic GLM coefficients for: ", names(Z)[j])
    offset <- cf[1] + cf[2] * meta$age + cf[3] * log(meta$age) + cf[4] * meta$sex
    out[[j]] <- Z[[j]] - offset
  }
  out
}

add_hmm_covariates <- function(Z_corrected, meta) {
  data.frame(Z_corrected, age = meta$age, sex = meta$sex, check.names = FALSE)
}

sequence_lengths <- function(ids) {
  ids <- as.character(ids)
  rr <- rle(ids)
  if (anyDuplicated(rr$values)) stop("At least one stay_id occurs in non-contiguous blocks.")
  as.numeric(rr$lengths)
}

build_iteration_bundle <- function(iteration, train_frac = 0.8) {
  base <- ensure_base_bundle()
  meta <- row_metadata(base)
  all_ids <- unique(as.character(meta$stay_id))
  set.seed(as.integer(iteration))
  n_train <- floor(length(all_ids) * train_frac)
  train_ids <- sample(all_ids, n_train, replace = FALSE)
  val_ids <- setdiff(all_ids, train_ids)
  train_rows <- meta$stay_id %in% train_ids
  val_rows <- !train_rows

  mice_data <- as.data.frame(base$imputation_source[, base$pre_mice_vars, drop = FALSE])
  if (any(!vapply(mice_data, is.numeric, logical(1)))) stop("MICE input contains non-numeric columns.")
  set.seed(500000L + as.integer(iteration))
  mice_fit <- mice::mice(
    mice_data,
    m = 5,
    maxit = 5,
    method = "pmm",
    ignore = val_rows,
    seed = 500000L + as.integer(iteration),
    printFlag = FALSE
  )
  imputed_pre <- mean_mice_completions(mice_fit)
  observations <- make_observation_matrix(imputed_pre, base$frozen_vars)

  set.seed(600000L + as.integer(iteration))
  normalizers <- fit_normalizers(observations, train_rows = train_rows)
  Z_uncentered <- transform_with_normalizers(observations, normalizers)

  survivor_fraction <- mean(is.na(base$cohort$dod))
  target_n <- max(50L, round(survivor_fraction * length(train_ids)))
  stable_ref <- select_stable_reference(base, train_ids, target_n)
  stable_glms <- fit_demographic_glms(Z_uncentered, stable_ref$row_index, meta)
  Z_corrected <- apply_demographic_glms(Z_uncentered, stable_glms, meta)

  hmm_data <- add_hmm_covariates(Z_corrected, meta)
  response_formulas <- lapply(base$frozen_vars, function(v) stats::as.formula(paste(v, "~ 1")))
  names(response_formulas) <- base$frozen_vars

  train_df <- hmm_data[train_rows, , drop = FALSE]
  val_df <- hmm_data[val_rows, , drop = FALSE]
  train_pat_ids <- meta$stay_id[train_rows]
  val_pat_ids <- meta$stay_id[val_rows]

  list(
    version = MIMIC_PREPROCESSING_VERSION,
    stage = "model_selection_iteration",
    iteration = as.integer(iteration),
    train_frac = train_frac,
    base_signature = list(
      input = base$input_signature,
      cohort = base$cohort_signature,
      cohort_finalization_manifest = base$cohort_finalization_manifest_signature,
      variable_md5 = base$variable_config$md5
    ),
    train_ids = train_ids,
    val_ids = val_ids,
    train_df = train_df,
    val_df = val_df,
    train_pat_ids = train_pat_ids,
    val_pat_ids = val_pat_ids,
    train_ntimes = sequence_lengths(train_pat_ids),
    val_ntimes = sequence_lengths(val_pat_ids),
    response_formulas = response_formulas,
    observation_variable_names = base$frozen_vars,
    normalization_spec_version = MIMIC_NORMALIZATION_SPEC_VERSION,
    normalization_config_md5 = MIMIC_NORMALIZATION_CONFIG$md5,
    normalization_families = normalization_family_map(base$frozen_vars),
    normalizers = normalizers,
    stable_reference_glms = stable_glms,
    stable_reference = stable_ref,
    mice_method = mice_fit$method,
    mice_predictor_matrix = mice_fit$predictorMatrix,
    created_at = Sys.time()
  )
}

iteration_cache_path <- function(iteration) {
  cache_dir <- Sys.getenv("MIMIC_HMM_CACHE_DIR", unset = "preprocessing_cache")
  file.path(cache_dir, sprintf("iteration_%02d.rds", as.integer(iteration)))
}

iteration_cache_compatible <- function(x, iteration, base) {
  is.list(x) && identical(x$version, MIMIC_PREPROCESSING_VERSION) &&
    identical(x$stage, "model_selection_iteration") &&
    identical(x$normalization_spec_version, MIMIC_NORMALIZATION_SPEC_VERSION) &&
    identical(x$normalization_config_md5, MIMIC_NORMALIZATION_CONFIG$md5) &&
    identical(x$normalization_families, normalization_family_map(base$frozen_vars)) &&
    identical(x$iteration, as.integer(iteration)) &&
    identical(x$base_signature$input, base$input_signature) &&
    identical(x$base_signature$cohort, base$cohort_signature) &&
    identical(
      x$base_signature$cohort_finalization_manifest,
      base$cohort_finalization_manifest_signature
    ) &&
    identical(x$base_signature$variable_md5, base$variable_config$md5)
}

ensure_iteration_bundle <- function(iteration,
                                    force = get_env_bool("MIMIC_FORCE_PREPROCESSING", FALSE)) {
  base_current <- ensure_base_bundle(force = force)
  path <- iteration_cache_path(iteration)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!force && file.exists(path)) {
    old <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(old) && iteration_cache_compatible(old, iteration, base_current)) return(old)
  }

  x <- build_iteration_bundle(iteration = iteration)
  atomic_save_rds(x, path)
  prefix <- sub("\\.rds$", "", path)
  readr::write_csv(x$stable_reference, paste0(prefix, "_stable_reference.csv"))
  readr::write_csv(
    reference_selection_summary(x$stable_reference),
    paste0(prefix, "_stable_reference_summary.csv")
  )
  readr::write_csv(
    data.frame(stay_id = x$train_ids, split = "train"),
    paste0(prefix, "_train_ids.csv")
  )
  readr::write_csv(
    data.frame(stay_id = x$val_ids, split = "validation"),
    paste0(prefix, "_validation_ids.csv")
  )
  x
}

build_fullcohort_bundle <- function() {
  base <- ensure_base_bundle()
  meta <- row_metadata(base)
  all_rows <- rep(TRUE, nrow(meta))
  mice_data <- as.data.frame(base$imputation_source[, base$pre_mice_vars, drop = FALSE])
  set.seed(700000L)
  mice_fit <- mice::mice(
    mice_data,
    m = 5,
    maxit = 5,
    method = "pmm",
    seed = 700000L,
    printFlag = FALSE
  )
  imputed_pre <- mean_mice_completions(mice_fit)
  observations <- make_observation_matrix(imputed_pre, base$frozen_vars)

  set.seed(800000L)
  normalizers <- fit_normalizers(observations, train_rows = all_rows)
  Z_uncentered <- transform_with_normalizers(observations, normalizers)

  all_ids <- unique(meta$stay_id)
  target_n <- sum(is.na(base$cohort$dod))
  target_n <- min(target_n, length(all_ids))
  stable_ref <- select_stable_reference(base, all_ids, target_n)
  stable_glms <- fit_demographic_glms(Z_uncentered, stable_ref$row_index, meta)
  Z_corrected <- apply_demographic_glms(Z_uncentered, stable_glms, meta)

  hmm_data <- add_hmm_covariates(Z_corrected, meta)
  response_formulas <- lapply(base$frozen_vars, function(v) stats::as.formula(paste(v, "~ 1")))
  names(response_formulas) <- base$frozen_vars

  list(
    version = MIMIC_PREPROCESSING_VERSION,
    stage = "final_fullcohort",
    base_signature = list(
      input = base$input_signature,
      cohort = base$cohort_signature,
      cohort_finalization_manifest = base$cohort_finalization_manifest_signature,
      variable_md5 = base$variable_config$md5
    ),
    hmm_data = hmm_data,
    combined = as.data.frame(base$imputation_source),
    imputed_data = observations,
    pat_ids = meta$stay_id,
    ntimes = sequence_lengths(meta$stay_id),
    response_formulas = response_formulas,
    observation_variable_names = base$frozen_vars,
    normalization_spec_version = MIMIC_NORMALIZATION_SPEC_VERSION,
    normalization_config_md5 = MIMIC_NORMALIZATION_CONFIG$md5,
    normalization_families = normalization_family_map(base$frozen_vars),
    normalizers = normalizers,
    stable_reference_glms = stable_glms,
    stable_reference = stable_ref,
    transformed_uncentered = Z_uncentered,
    metadata = meta,
    mice_method = mice_fit$method,
    mice_predictor_matrix = mice_fit$predictorMatrix,
    created_at = Sys.time()
  )
}

fullcohort_cache_path <- function() {
  cache_dir <- Sys.getenv("MIMIC_HMM_CACHE_DIR", unset = "preprocessing_cache")
  file.path(cache_dir, "fullcohort.rds")
}

fullcohort_cache_compatible <- function(x, base) {
  is.list(x) && identical(x$version, MIMIC_PREPROCESSING_VERSION) &&
    identical(x$stage, "final_fullcohort") &&
    identical(x$normalization_spec_version, MIMIC_NORMALIZATION_SPEC_VERSION) &&
    identical(x$normalization_config_md5, MIMIC_NORMALIZATION_CONFIG$md5) &&
    identical(x$normalization_families, normalization_family_map(base$frozen_vars)) &&
    identical(x$base_signature$input, base$input_signature) &&
    identical(x$base_signature$cohort, base$cohort_signature) &&
    identical(
      x$base_signature$cohort_finalization_manifest,
      base$cohort_finalization_manifest_signature
    ) &&
    identical(x$base_signature$variable_md5, base$variable_config$md5)
}

ensure_fullcohort_bundle <- function(force = get_env_bool("MIMIC_FORCE_PREPROCESSING", FALSE)) {
  base_current <- ensure_base_bundle(force = force)
  path <- fullcohort_cache_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  if (!force && file.exists(path)) {
    old <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(old) && fullcohort_cache_compatible(old, base_current)) {
      return(old)
    }
  }

  x <- build_fullcohort_bundle()
  atomic_save_rds(x, path)
  prefix <- sub("\\.rds$", "", path)
  readr::write_csv(x$stable_reference, paste0(prefix, "_stable_reference.csv"))
  readr::write_csv(
    reference_selection_summary(x$stable_reference),
    paste0(prefix, "_stable_reference_summary.csv")
  )
  x
}
