
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))

.required_objects <- c(
  "fitted_model", "response_formulas", "bests",
  "imputed_data_scale_glm_list", "imputed_data", "imputed_data_scale_df"
)
.missing_objects <- .required_objects[
  !vapply(.required_objects, exists, logical(1), inherits = TRUE)
]
if (length(.missing_objects)) {
  stop("Missing required objects: ", paste(.missing_objects, collapse = ", "), call. = FALSE)
}

obs_names <- vapply(response_formulas, function(f) {
  vars <- all.vars(f)
  if (length(vars) != 1L) {
    stop("Unexpected response formula: ", paste(deparse(f), collapse = ""), call. = FALSE)
  }
  vars[[1L]]
}, character(1))

if (!length(obs_names) || anyDuplicated(obs_names)) {
  stop("Invalid or duplicated observation names.", call. = FALSE)
}

K <- length(fitted_model@response)
R <- length(fitted_model@response[[1L]])
if (K < 2L || R != length(obs_names)) {
  stop(
    "Model/formula dimension mismatch: K=", K,
    ", R=", R, ", formulas=", length(obs_names),
    call. = FALSE
  )
}

for (nm in c("bests", "imputed_data_scale_glm_list")) {
  obj <- get(nm, inherits = TRUE)
  if (is.null(names(obj)) || any(!nzchar(names(obj))) || anyDuplicated(names(obj))) {
    stop(nm, " must be uniquely and completely named.", call. = FALSE)
  }
  miss <- setdiff(obs_names, names(obj))
  extra <- setdiff(names(obj), obs_names)
  if (length(miss) || length(extra)) {
    stop(
      nm, " variable mismatch; missing=", paste(miss, collapse = ";"),
      "; extra=", paste(extra, collapse = ";"),
      call. = FALSE
    )
  }
}

imputed_data <- as.data.frame(imputed_data, check.names = FALSE)
imputed_data_scale_df <- as.data.frame(imputed_data_scale_df, check.names = FALSE)
if (!setequal(names(imputed_data), obs_names)) {
  stop("imputed_data variables do not match HMM response variables.", call. = FALSE)
}
imputed_data <- imputed_data[, obs_names, drop = FALSE]
if (any(!vapply(imputed_data, is.numeric, logical(1)))) {
  stop("imputed_data contains non-numeric observation columns.", call. = FALSE)
}
if (!all(c(obs_names, "age", "sex") %in% names(imputed_data_scale_df))) {
  stop("imputed_data_scale_df lacks HMM variables and/or age/sex.", call. = FALSE)
}
if (nrow(imputed_data) != nrow(imputed_data_scale_df)) {
  stop("Raw and HMM-input row counts differ.", call. = FALSE)
}

mu_norm <- matrix(
  NA_real_, nrow = K, ncol = R,
  dimnames = list(paste0("State", seq_len(K)), obs_names)
)
sd_norm <- matrix(
  NA_real_, nrow = K, ncol = R,
  dimnames = list(paste0("State", seq_len(K)), obs_names)
)

for (k in seq_len(K)) {
  if (length(fitted_model@response[[k]]) != R) {
    stop("Response count differs across states at state ", k, call. = FALSE)
  }
  for (j in seq_len(R)) {
    par_j <- fitted_model@response[[k]][[j]]@parameters
    mu_j <- suppressWarnings(as.numeric(par_j$coefficients))
    sd_j <- suppressWarnings(as.numeric(par_j$sd))
    if (
      length(mu_j) != 1L || length(sd_j) != 1L ||
      !is.finite(mu_j) || !is.finite(sd_j) || sd_j <= 0
    ) {
      stop(
        "Invalid Gaussian parameters at state=", k,
        ", variable=", obs_names[[j]],
        ": mean=", paste(mu_j, collapse = ","),
        ", sd=", paste(sd_j, collapse = ","),
        call. = FALSE
      )
    }
    mu_norm[k, j] <- mu_j
    sd_norm[k, j] <- sd_j
  }
}

mimic_normalizer_family <- function(obj) {
  if (!is.null(obj$type) && identical(obj$type, "zscore")) return("zscore")
  fam <- attr(obj, "mimic_normalization_family")
  if (!is.null(fam) && length(fam)) return(as.character(fam[[1L]]))
  class(obj)[[1L]]
}

mimic_forward_normalizer <- function(obj, x) {
  x <- suppressWarnings(as.numeric(x))
  if (!length(x) || any(!is.finite(x))) {
    stop("Non-finite raw input to normalizer.", call. = FALSE)
  }
  if (!is.null(obj$type) && identical(obj$type, "zscore")) {
    mu <- suppressWarnings(as.numeric(obj$mean))[[1L]]
    s <- suppressWarnings(as.numeric(obj$sd))[[1L]]
    if (!is.finite(mu) || !is.finite(s) || s <= 0) {
      stop("Invalid z-score normalizer.", call. = FALSE)
    }
    return((x - mu) / s)
  }
  out <- if (inherits(obj, "orderNorm")) {
    stats::predict(obj, newdata = x, warn = FALSE)
  } else {
    stats::predict(obj, newdata = x)
  }
  out <- suppressWarnings(as.numeric(out))
  if (length(out) != length(x) || any(!is.finite(out))) {
    stop(
      "Forward normalization failed for class ",
      paste(class(obj), collapse = "/"),
      call. = FALSE
    )
  }
  out
}

mimic_inverse_normalizer <- function(obj, z) {
  z <- suppressWarnings(as.numeric(z))
  if (!length(z) || any(!is.finite(z))) {
    stop("Non-finite normalized input to inverse normalizer.", call. = FALSE)
  }
  if (!is.null(obj$type) && identical(obj$type, "zscore")) {
    mu <- suppressWarnings(as.numeric(obj$mean))[[1L]]
    s <- suppressWarnings(as.numeric(obj$sd))[[1L]]
    if (!is.finite(mu) || !is.finite(s) || s <= 0) {
      stop("Invalid z-score normalizer.", call. = FALSE)
    }
    return(z * s + mu)
  }
  out <- if (inherits(obj, "orderNorm")) {
    stats::predict(obj, newdata = z, inverse = TRUE, warn = FALSE)
  } else {
    stats::predict(obj, newdata = z, inverse = TRUE)
  }
  out <- suppressWarnings(as.numeric(out))
  if (length(out) != length(z) || any(!is.finite(out))) {
    stop(
      "Inverse normalization failed for class ",
      paste(class(obj), collapse = "/"),
      call. = FALSE
    )
  }
  out
}

.accepted_glm_terms <- list(
  intercept = c("(Intercept)"),
  age = c("age_ref", "age", "age_sex$age"),
  logage = c("log(age_ref)", "log(age)", "`log(age)`", "log(age_sex$age)"),
  sex = c("sex_ref", "sex", "age_sex$sex")
)

mimic_validate_demographic_glm <- function(glm_obj, variable = "") {
  if (is.null(glm_obj)) stop("Missing demographic GLM for ", variable, call. = FALSE)
  cf <- stats::coef(glm_obj)
  nm <- names(cf)
  ok <- length(cf) == 4L && all(is.finite(cf)) && length(nm) == 4L &&
    nm[[1L]] %in% .accepted_glm_terms$intercept &&
    nm[[2L]] %in% .accepted_glm_terms$age &&
    nm[[3L]] %in% .accepted_glm_terms$logage &&
    nm[[4L]] %in% .accepted_glm_terms$sex
  if (!ok) {
    stop(
      "Unexpected demographic GLM contract for ", variable,
      ": terms=", paste(nm, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

mimic_demographic_offset <- function(glm_obj, age, sex, variable = "") {
  mimic_validate_demographic_glm(glm_obj, variable)
  cf <- stats::coef(glm_obj)
  age <- suppressWarnings(as.numeric(age))
  sex <- suppressWarnings(as.numeric(sex))
  if (
    length(age) != length(sex) && length(age) != 1L && length(sex) != 1L
  ) {
    stop("age/sex lengths are incompatible.", call. = FALSE)
  }
  if (any(!is.finite(age)) || any(age <= 0) || any(!is.finite(sex))) {
    stop("Invalid age/sex supplied for demographic offset.", call. = FALSE)
  }
  cf[[1L]] + cf[[2L]] * age + cf[[3L]] * log(age) + cf[[4L]] * sex
}

for (v in obs_names) {
  mimic_validate_demographic_glm(imputed_data_scale_glm_list[[v]], v)
}

normalizer_family_contract_qc <- data.frame(
  variable = obs_names,
  actual_family = vapply(bests[obs_names], mimic_normalizer_family, character(1)),
  expected_family = NA_character_,
  family_match = NA,
  stringsAsFactors = FALSE
)
if (exists("MIMIC_NORMALIZATION_FAMILIES", inherits = TRUE)) {
  expected_map <- get("MIMIC_NORMALIZATION_FAMILIES", inherits = TRUE)
  normalizer_family_contract_qc$expected_family <- as.character(expected_map[obs_names])
  normalizer_family_contract_qc$family_match <-
    normalizer_family_contract_qc$actual_family ==
    normalizer_family_contract_qc$expected_family
  if (anyNA(normalizer_family_contract_qc$expected_family) ||
      any(!normalizer_family_contract_qc$family_match)) {
    bad <- normalizer_family_contract_qc$variable[
      is.na(normalizer_family_contract_qc$family_match) |
        !normalizer_family_contract_qc$family_match
    ]
    stop(
      "Frozen normalizer-family contract failed for: ",
      paste(bad, collapse = ", "),
      call. = FALSE
    )
  }
}
utils::write.csv(
  normalizer_family_contract_qc,
  "state_profile_normalizer_qc.csv",
  row.names = FALSE
)

mimic_display_unit_spec <- do.call(rbind, lapply(obs_names, function(v) {
  x <- suppressWarnings(as.numeric(imputed_data[[v]]))
  x <- x[is.finite(x)]
  if (!length(x)) stop("No finite training values for ", v, call. = FALSE)
  q <- stats::quantile(x, probs = c(0.05, 0.50, 0.95), na.rm = TRUE, names = FALSE)
  model_unit <- "display"
  multiplier <- 1
  if (v %in% c("FiO2", "SpO2", "Hematocrit")) {
    if (is.finite(q[[3L]]) && q[[3L]] <= 1.5) {
      model_unit <- "fraction"
      multiplier <- 100
    } else {
      model_unit <- "percent"
    }
  } else if (v == "PlateletCount") {
    if (is.finite(q[[2L]]) && q[[2L]] > 1000) {
      model_unit <- "per_uL"
      multiplier <- 0.001
    } else {
      model_unit <- "10^3_per_uL"
    }
  }
  data.frame(
    variable = v,
    model_unit = model_unit,
    display_multiplier = multiplier,
    training_q05_model = q[[1L]],
    training_q50_model = q[[2L]],
    training_q95_model = q[[3L]],
    training_q05_display = q[[1L]] * multiplier,
    training_q50_display = q[[2L]] * multiplier,
    training_q95_display = q[[3L]] * multiplier,
    stringsAsFactors = FALSE
  )
}))
rownames(mimic_display_unit_spec) <- mimic_display_unit_spec$variable

mimic_to_display_unit <- function(variable, x) {
  if (!variable %in% rownames(mimic_display_unit_spec)) {
    stop("Unknown display-unit variable: ", variable, call. = FALSE)
  }
  suppressWarnings(as.numeric(x)) * mimic_display_unit_spec[variable, "display_multiplier"]
}

mimic_to_model_unit <- function(variable, x) {
  if (!variable %in% rownames(mimic_display_unit_spec)) {
    stop("Unknown display-unit variable: ", variable, call. = FALSE)
  }
  suppressWarnings(as.numeric(x)) / mimic_display_unit_spec[variable, "display_multiplier"]
}

utils::write.csv(
  mimic_display_unit_spec,
  "state_profile_display_units.csv",
  row.names = FALSE
)

age_vec <- suppressWarnings(as.numeric(imputed_data_scale_df$age))
sex_vec <- suppressWarnings(as.numeric(imputed_data_scale_df$sex))
if (any(!is.finite(age_vec)) || any(age_vec <= 0) || any(!is.finite(sex_vec))) {
  stop("Invalid age/sex in imputed_data_scale_df.", call. = FALSE)
}

preprocessing_contract_qc <- do.call(rbind, lapply(obs_names, function(v) {
  z_uncentered <- mimic_forward_normalizer(bests[[v]], imputed_data[[v]])
  offset <- mimic_demographic_offset(
    imputed_data_scale_glm_list[[v]], age_vec, sex_vec, variable = v
  )
  reconstructed <- z_uncentered - offset
  actual <- suppressWarnings(as.numeric(imputed_data_scale_df[[v]]))
  if (length(actual) != length(reconstructed) || any(!is.finite(actual))) {
    stop("Invalid HMM-input column for ", v, call. = FALSE)
  }
  err <- abs(reconstructed - actual)
  data.frame(
    variable = v,
    normalizer_family = mimic_normalizer_family(bests[[v]]),
    glm_terms = paste(names(stats::coef(imputed_data_scale_glm_list[[v]])), collapse = ";"),
    n_rows = length(err),
    max_absolute_hmm_input_error = max(err),
    median_absolute_hmm_input_error = stats::median(err),
    hmm_input_reproduction_ok = max(err) <= 1e-8,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  preprocessing_contract_qc,
  "state_profile_preprocessing_qc.csv",
  row.names = FALSE
)
if (any(!preprocessing_contract_qc$hmm_input_reproduction_ok)) {
  stop(
    "Raw-to-HMM preprocessing reproduction failed for: ",
    paste(
      preprocessing_contract_qc$variable[
        !preprocessing_contract_qc$hmm_input_reproduction_ok
      ],
      collapse = ", "
    ),
    call. = FALSE
  )
}

.roundtrip_tol <- suppressWarnings(as.numeric(Sys.getenv(
  "MIMIC_ANNOT_RAW_ROUNDTRIP_TOL", unset = "0.02"
)))
if (!is.finite(.roundtrip_tol) || .roundtrip_tol <= 0) .roundtrip_tol <- 0.02

roundtrip_qc <- do.call(rbind, lapply(obs_names, function(v) {
  x <- suppressWarnings(as.numeric(imputed_data[[v]]))
  x <- x[is.finite(x)]
  test_x <- unique(as.numeric(stats::quantile(
    x,
    probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99),
    na.rm = TRUE,
    names = FALSE
  )))
  z <- mimic_forward_normalizer(bests[[v]], test_x)
  back <- mimic_inverse_normalizer(bests[[v]], z)
  abs_err <- abs(back - test_x)
  ref_scale <- max(
    stats::IQR(x, na.rm = TRUE),
    abs(stats::median(x, na.rm = TRUE)) * 0.01,
    1e-8
  )
  scaled_err <- max(abs_err) / ref_scale
  data.frame(
    variable = v,
    family = mimic_normalizer_family(bests[[v]]),
    n_test = length(test_x),
    max_absolute_error_model = max(abs_err),
    max_scaled_error = scaled_err,
    tolerance = .roundtrip_tol,
    roundtrip_ok = scaled_err <= .roundtrip_tol,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  roundtrip_qc,
  "state_profile_inverse_transform_qc.csv",
  row.names = FALSE
)
if (any(!roundtrip_qc$roundtrip_ok)) {
  stop(
    "Normalizer round-trip QC failed for: ",
    paste(roundtrip_qc$variable[!roundtrip_qc$roundtrip_ok], collapse = ", "),
    call. = FALSE
  )
}

AGES <- c(55, 65, 75)
SEXS <- c(0, 1)
mean_profiles <- list()
low_profiles <- list()
high_profiles <- list()
halfwidth_profiles <- list()
sd_delta_profiles <- list()
coordinate_qc_rows <- list()
qci <- 1L

train_z_uncentered <- as.data.frame(lapply(obs_names, function(v) {
  mimic_forward_normalizer(bests[[v]], imputed_data[[v]])
}), check.names = FALSE)
names(train_z_uncentered) <- obs_names

.analytic_state_z_tol <- suppressWarnings(as.numeric(Sys.getenv(
  "MIMIC_ANNOT_ANALYTIC_STATE_Z_TOL", unset = "1e-5"
)))
if (!is.finite(.analytic_state_z_tol) || .analytic_state_z_tol <= 0) {
  .analytic_state_z_tol <- 1e-5
}
.orderNorm_warn_frac_sd <- suppressWarnings(as.numeric(Sys.getenv(
  "MIMIC_ANNOT_ORDERNORM_STATE_Z_WARN_FRAC_SD", unset = "0.10"
)))
if (!is.finite(.orderNorm_warn_frac_sd) || .orderNorm_warn_frac_sd <= 0) {
  .orderNorm_warn_frac_sd <- 0.10
}
.orderNorm_fail_frac_sd <- suppressWarnings(as.numeric(Sys.getenv(
  "MIMIC_ANNOT_ORDERNORM_STATE_Z_FAIL_FRAC_SD", unset = "0.50"
)))
if (!is.finite(.orderNorm_fail_frac_sd) ||
    .orderNorm_fail_frac_sd <= .orderNorm_warn_frac_sd) {
  .orderNorm_fail_frac_sd <- 0.50
}

for (age_ in AGES) {
  for (sex_ in SEXS) {
    key <- paste(age_, sex_, sep = "_")
    mean_mat <- low_mat <- high_mat <- half_mat <- sd_delta_mat <- matrix(
      NA_real_, nrow = K, ncol = R,
      dimnames = list(rownames(mu_norm), obs_names)
    )
    for (j in seq_along(obs_names)) {
      v <- obs_names[[j]]
      off <- mimic_demographic_offset(
        imputed_data_scale_glm_list[[v]], age_, sex_, variable = v
      )
      z_mean <- mu_norm[, j] + off
      z_low <- z_mean - sd_norm[, j]
      z_high <- z_mean + sd_norm[, j]

      raw_mean_model <- mimic_inverse_normalizer(bests[[v]], z_mean)
      raw_low_model <- mimic_inverse_normalizer(bests[[v]], z_low)
      raw_high_model <- mimic_inverse_normalizer(bests[[v]], z_high)

      z_back <- mimic_forward_normalizer(bests[[v]], raw_mean_model)
      z_error <- abs(z_back - z_mean)
      max_z_err <- max(z_error)
      error_frac_sd <- z_error / sd_norm[, j]
      max_error_frac_sd <- max(error_frac_sd)
      family_v <- mimic_normalizer_family(bests[[v]])

      display_mean <- mimic_to_display_unit(v, raw_mean_model)
      display_low <- mimic_to_display_unit(v, raw_low_model)
      display_high <- mimic_to_display_unit(v, raw_high_model)

      eps <- 1e-4
      raw_eps_plus <- mimic_inverse_normalizer(bests[[v]], z_mean + eps)
      raw_eps_minus <- mimic_inverse_normalizer(bests[[v]], z_mean - eps)
      dydz_display <- mimic_to_display_unit(
        v, (raw_eps_plus - raw_eps_minus) / (2 * eps)
      )
      display_sd_delta <- abs(dydz_display) * sd_norm[, j]

      if (any(display_low > display_mean) || any(display_mean > display_high)) {
        stop(
          "Non-monotone inverse interval for ", v,
          " at age=", age_, ", sex=", sex_,
          call. = FALSE
        )
      }

      mean_mat[, j] <- display_mean
      low_mat[, j] <- display_low
      high_mat[, j] <- display_high
      half_mat[, j] <- (display_high - display_low) / 2
      sd_delta_mat[, j] <- display_sd_delta

      tz <- suppressWarnings(as.numeric(train_z_uncentered[[v]]))
      tz <- tz[is.finite(tz)]
      z_min <- min(tz)
      z_max <- max(tz)
      in_support <- z_mean >= z_min & z_mean <= z_max
      max_in_support <- if (any(in_support)) max(z_error[in_support]) else NA_real_
      max_out_support <- if (any(!in_support)) max(z_error[!in_support]) else NA_real_
      ties_status <- if (!is.null(bests[[v]]$ties_status)) {
        suppressWarnings(as.integer(bests[[v]]$ties_status[[1L]]))
      } else {
        NA_integer_
      }

      hard_fail <- if (identical(family_v, "orderNorm")) {
        !is.finite(max_error_frac_sd)
      } else {
        !is.finite(max_z_err) || max_z_err > .analytic_state_z_tol
      }
      warn_flag <- identical(family_v, "orderNorm") &&
        is.finite(max_error_frac_sd) &&
        max_error_frac_sd > .orderNorm_warn_frac_sd
      coordinate_qc_rows[[qci]] <- data.frame(
        variable = v,
        normalizer_family = family_v,
        orderNorm_ties_status = ties_status,
        age = age_,
        sex = sex_,
        max_inverse_forward_z_error = max_z_err,
        max_inverse_forward_z_error_in_training_support = max_in_support,
        max_inverse_forward_z_error_outside_training_support = max_out_support,
        max_inverse_forward_error_fraction_emission_sd = max_error_frac_sd,
        median_inverse_forward_error_fraction_emission_sd = stats::median(error_frac_sd),
        analytic_z_tolerance = .analytic_state_z_tol,
        orderNorm_warn_fraction_emission_sd = .orderNorm_warn_frac_sd,
        orderNorm_fail_fraction_emission_sd = .orderNorm_fail_frac_sd,
        orderNorm_severe_approximation = identical(family_v, "orderNorm") &&
          is.finite(max_error_frac_sd) && max_error_frac_sd > .orderNorm_fail_frac_sd,
        inverse_forward_warning = warn_flag,
        inverse_forward_hard_fail = hard_fail,
        n_state_means_below_training_z = sum(z_mean < z_min),
        n_state_means_above_training_z = sum(z_mean > z_max),
        n_lower_bounds_below_training_z = sum(z_low < z_min),
        n_upper_bounds_above_training_z = sum(z_high > z_max),
        training_z_min = z_min,
        training_z_max = z_max,
        stringsAsFactors = FALSE
      )
      qci <- qci + 1L
    }

    mean_profiles[[key]] <- as.data.frame(mean_mat, check.names = FALSE)
    low_profiles[[key]] <- as.data.frame(low_mat, check.names = FALSE)
    high_profiles[[key]] <- as.data.frame(high_mat, check.names = FALSE)
    halfwidth_profiles[[key]] <- as.data.frame(half_mat, check.names = FALSE)
    sd_delta_profiles[[key]] <- as.data.frame(sd_delta_mat, check.names = FALSE)

    utils::write.table(
      mean_mat,
      paste0("state_profiles_mu_orig.", age_, ".", sex_, ".csv"),
      row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ","
    )
    utils::write.table(
      sd_delta_mat,
      paste0("state_profiles_sd_orig.", age_, ".", sex_, ".csv"),
      row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ","
    )
    utils::write.table(
      low_mat,
      paste0("state_profiles_interval_low_orig.", age_, ".", sex_, ".csv"),
      row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ","
    )
    utils::write.table(
      high_mat,
      paste0("state_profiles_interval_high_orig.", age_, ".", sex_, ".csv"),
      row.names = TRUE, col.names = TRUE, quote = FALSE, sep = ","
    )
  }
}

coordinate_qc <- do.call(rbind, coordinate_qc_rows)
utils::write.csv(
  coordinate_qc,
  "state_profile_coordinate_qc.csv",
  row.names = FALSE
)
if (any(coordinate_qc$inverse_forward_hard_fail)) {
  bad <- coordinate_qc[coordinate_qc$inverse_forward_hard_fail, , drop = FALSE]
  stop(
    "State-profile inverse/forward QC found material failures: ",
    paste0(
      bad$variable, "@", bad$age, "/", bad$sex,
      " (family=", bad$normalizer_family,
      ", max_z_error=", signif(bad$max_inverse_forward_z_error, 5),
      ", max_error/SD=", signif(bad$max_inverse_forward_error_fraction_emission_sd, 5),
      ")",
      collapse = "; "
    ),
    ". See state_profile_coordinate_qc.csv",
    call. = FALSE
  )
}
if (any(coordinate_qc$inverse_forward_warning)) {
  warning(
    "orderNorm state-profile inverse/forward approximation exceeded the warning threshold for ",
    paste(unique(coordinate_qc$variable[coordinate_qc$inverse_forward_warning]), collapse = ", "),
    ". This is diagnostic, not a preprocessing-contract failure; see state_coordinate_QC.",
    call. = FALSE
  )
}

.bind_profile_list <- function(lst) {
  out <- do.call(rbind, lapply(names(lst), function(key) {
    p <- strsplit(key, "_", fixed = TRUE)[[1L]]
    data.frame(
      state_id = seq_len(K),
      age = as.numeric(p[[1L]]),
      sex = as.numeric(p[[2L]]),
      lst[[key]],
      check.names = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

.median_profile <- function(df) {
  vars <- setdiff(names(df), c("state_id", "age", "sex"))
  out <- stats::aggregate(
    df[, vars, drop = FALSE],
    by = list(state_id = df$state_id),
    FUN = function(x) stats::median(as.numeric(x), na.rm = TRUE)
  )
  out[order(out$state_id), , drop = FALSE]
}

mean_all <- .bind_profile_list(mean_profiles)
low_all <- .bind_profile_list(low_profiles)
high_all <- .bind_profile_list(high_profiles)
half_all <- .bind_profile_list(halfwidth_profiles)
sd_delta_all <- .bind_profile_list(sd_delta_profiles)

representative_mean <- .median_profile(mean_all)
representative_low <- .median_profile(low_all)
representative_high <- .median_profile(high_all)
representative_half <- .median_profile(half_all)
representative_sd_delta <- .median_profile(sd_delta_all)

.display_guardrails <- list(
  MAP = c(20, 200),
  UreaNitrogen = c(2, 300),
  HCO3 = c(5, 60),
  PlateletCount = c(5, 2000),
  PF = c(20, 2000),
  FiO2 = c(15, 100),
  RR = c(5, 120),
  AST = c(0, 10000),
  ALT = c(0, 10000),
  INR = c(0.5, 10),
  Lactate = c(0, 20),
  Creatinine = c(0.05, 20)
)
profile_scale_qc <- do.call(rbind, lapply(names(.display_guardrails), function(v) {
  x <- suppressWarnings(as.numeric(representative_mean[[v]]))
  lim <- .display_guardrails[[v]]
  med <- stats::median(x, na.rm = TRUE)
  data.frame(
    variable = v,
    min_state_mean = min(x, na.rm = TRUE),
    median_state_mean = med,
    max_state_mean = max(x, na.rm = TRUE),
    guardrail_low = lim[[1L]],
    guardrail_high = lim[[2L]],
    scale_ok = is.finite(med) && med >= lim[[1L]] && med <= lim[[2L]],
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  profile_scale_qc,
  "state_profile_scale_qc.csv",
  row.names = FALSE
)
if (any(!profile_scale_qc$scale_ok)) {
  stop(
    "State profiles are not on display scale for: ",
    paste(profile_scale_qc$variable[!profile_scale_qc$scale_ok], collapse = ", "),
    call. = FALSE
  )
}

mimic_state_profile_bundle <- list(
  version = "state_profile_v1",
  state_ids = seq_len(K),
  observation_names = obs_names,
  ages = AGES,
  sexes = SEXS,
  mu_residualized_normalized = mu_norm,
  sd_residualized_normalized = sd_norm,
  mean_by_stratum = mean_profiles,
  interval_low_by_stratum = low_profiles,
  interval_high_by_stratum = high_profiles,
  interval_halfwidth_by_stratum = halfwidth_profiles,
  sd_delta_by_stratum = sd_delta_profiles,
  mean_all = mean_all,
  interval_low_all = low_all,
  interval_high_all = high_all,
  interval_halfwidth_all = half_all,
  sd_delta_all = sd_delta_all,
  representative_mean = representative_mean,
  representative_interval_low = representative_low,
  representative_interval_high = representative_high,
  representative_interval_halfwidth = representative_half,
  representative_sd_delta = representative_sd_delta,
  display_unit_spec = mimic_display_unit_spec,
  normalizer_family_contract_qc = normalizer_family_contract_qc,
  preprocessing_contract_qc = preprocessing_contract_qc,
  normalizer_roundtrip_qc = roundtrip_qc,
  coordinate_qc = coordinate_qc,
  profile_scale_qc = profile_scale_qc
)

saveRDS(
  mimic_state_profile_bundle,
  "state_profile_bundle.rds"
)
utils::write.csv(
  representative_mean,
  "state_representative_mean.csv",
  row.names = FALSE
)
utils::write.csv(
  representative_low,
  "state_representative_interval_low.csv",
  row.names = FALSE
)
utils::write.csv(
  representative_high,
  "state_representative_interval_high.csv",
  row.names = FALSE
)

if (
  requireNamespace("uwot", quietly = TRUE) &&
  requireNamespace("ggplot2", quietly = TRUE) &&
  requireNamespace("ggrepel", quietly = TRUE)
) {
  tryCatch({
    set.seed(123)
    emb_mu <- uwot::umap(
      mu_norm,
      n_neighbors = min(15, K - 1),
      metric = "euclidean",
      scale = FALSE
    )
    k_clust <- min(10, K)
    cl <- stats::kmeans(emb_mu, centers = k_clust, nstart = 20)
    df_mu <- data.frame(
      UMAP1 = emb_mu[, 1],
      UMAP2 = emb_mu[, 2],
      state = factor(rownames(mu_norm), levels = rownames(mu_norm)),
      cluster = factor(cl$cluster)
    )
    p_umap <- ggplot2::ggplot(
      df_mu,
      ggplot2::aes(UMAP1, UMAP2, color = cluster, label = state)
    ) +
      ggplot2::geom_point(size = 3) +
      ggrepel::geom_text_repel(show.legend = FALSE, size = 3.4) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = "UMAP of HMM state means (residualized normalized space)",
        subtitle = paste0("K=", K, ", clusters=", k_clust),
        color = "cluster"
      )
    ggplot2::ggsave(
      "state_profiles_umap_states_normalized.png",
      p_umap, width = 8, height = 6, dpi = 600
    )
    ggplot2::ggsave(
      "state_profiles_umap_states_normalized.pdf",
      p_umap, width = 8, height = 6
    )
  }, error = function(e) {
    warning("Optional state UMAP was skipped: ", conditionMessage(e), call. = FALSE)
  })
} else {
  message("[state profiles] optional UMAP skipped")
}

message(
  "[state profiles] completed: K=", K,
  ", observations=", R,
  "; raw-to-HMM contract verified"
)
