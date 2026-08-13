
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
extract_posterior_probabilities <- function(post, nstates) {
  expected <- paste0("S", seq_len(nstates))
  if (all(expected %in% colnames(post))) {
    ans <- as.matrix(post[, expected, drop = FALSE])
  } else if (ncol(post) == nstates) {
    ans <- as.matrix(post)
    colnames(ans) <- expected
  } else if (ncol(post) == nstates + 1L && any(colnames(post) %in% c("state", "State"))) {
    label_col <- which(colnames(post) %in% c("state", "State"))[1]
    ans <- as.matrix(post[, -label_col, drop = FALSE])
    colnames(ans) <- expected
  } else {
    stop("Unexpected posterior dimensions: ", nrow(post), " x ", ncol(post), "; expected ", nstates, " probability columns.")
  }
  if (any(!is.finite(ans))) stop("Posterior probabilities contain non-finite values.")
  ans
}

formula_lhs_names <- function(response_formulas) {
  vapply(response_formulas, function(f) as.character(stats::as.formula(f))[2], character(1))
}

depmix_logLik_safe <- function(object) {

  method <- methods::selectMethod(
    "logLik",
    signature = c(object = "depmix"),
    optional = TRUE
  )
  if (is.null(method)) {
    stop("The depmixS4 S4 logLik method for class 'depmix' is unavailable.")
  }

  ll <- method(object)
  ll_value <- as.numeric(ll)
  if (length(ll_value) != 1L || !is.finite(ll_value)) {
    stop("Non-finite depmixS4 log-likelihood.")
  }

  df <- attr(ll, "df")
  if (is.null(df) || !is.finite(df)) df <- object@npars
  nobs <- attr(ll, "nobs")
  if (is.null(nobs) || !is.finite(nobs)) nobs <- sum(object@ntimes)

  structure(
    ll_value,
    class = "logLik",
    df = as.numeric(df),
    nobs = as.numeric(nobs)
  )
}

depmix_information_criteria <- function(ll) {
  ll_value <- as.numeric(ll)
  df <- as.numeric(attr(ll, "df"))
  nobs <- as.numeric(attr(ll, "nobs"))
  if (!is.finite(df) || df <= 0) stop("Invalid logLik degrees of freedom.")
  if (!is.finite(nobs) || nobs <= 0) stop("Invalid logLik observation count.")
  list(
    AIC = -2 * ll_value + 2 * df,
    BIC = -2 * ll_value + log(nobs) * df
  )
}

checkpoint_metadata <- function(dimension, seed, transition_rhs,
                                obs_names, train_df, train_ntimes,
                                preprocessing_cache_md5) {
  list(
    checkpoint_version = 1L,
    dimension = as.integer(dimension),
    seed = as.integer(seed),
    transition_rhs = as.character(transition_rhs),
    observation_variables = as.character(obs_names),
    train_rows = as.integer(nrow(train_df)),
    train_ntimes = as.integer(train_ntimes),
    preprocessing_cache_md5 = as.character(preprocessing_cache_md5)
  )
}

fit_and_eval_hmm_presplit <- function(train_df,
                                      val_df,
                                      train_ntimes,
                                      val_ntimes,
                                      response_formulas,
                                      dimension,
                                      transition_rhs,
                                      seed,
                                      maxit = 2000,
                                      tol = 1e-4,
                                      checkpoint_path = NULL,
                                      checkpoint_meta_path = NULL,
                                      preprocessing_cache_md5 = NA_character_) {
  suppressPackageStartupMessages(library(depmixS4))
  set.seed(as.integer(seed))
  nresp <- length(response_formulas)
  obs_names <- formula_lhs_names(response_formulas)
  if (!all(obs_names %in% names(train_df)) || !all(obs_names %in% names(val_df))) {
    stop("Some response variables are absent from the train/validation matrices.")
  }
  if (sum(train_ntimes) != nrow(train_df) || sum(val_ntimes) != nrow(val_df)) {
    stop("Sequence lengths do not match train/validation rows.")
  }

  hmm_model <- depmixS4::depmix(
    response = response_formulas,
    family = replicate(nresp, stats::gaussian(), simplify = FALSE),
    nstates = dimension,
    transition = stats::as.formula(transition_rhs),
    data = train_df,
    ntimes = train_ntimes
  )

  expected_checkpoint_meta <- checkpoint_metadata(
    dimension = dimension,
    seed = seed,
    transition_rhs = transition_rhs,
    obs_names = obs_names,
    train_df = train_df,
    train_ntimes = train_ntimes,
    preprocessing_cache_md5 = preprocessing_cache_md5
  )

  use_checkpoint <- FALSE
  if (!is.null(checkpoint_path) && !is.null(checkpoint_meta_path) &&
      file.exists(checkpoint_path) && file.exists(checkpoint_meta_path)) {
    saved_meta <- tryCatch(readRDS(checkpoint_meta_path), error = function(e) NULL)
    if (identical(saved_meta, expected_checkpoint_meta)) {
      fitted_model <- readRDS(checkpoint_path)
      if (!methods::is(fitted_model, "depmix.fitted") ||
          fitted_model@nstates != dimension ||
          !identical(as.integer(fitted_model@ntimes), as.integer(train_ntimes))) {
        stop("Saved HMM checkpoint is structurally incompatible.")
      }
      use_checkpoint <- TRUE
      cat("Resuming from fitted-model checkpoint: ", checkpoint_path, "\n", sep = "")
    } else {
      warning("Ignoring an incompatible fitted-model checkpoint: ", checkpoint_path)
    }
  }

  if (!use_checkpoint) {
    fitted_model <- depmixS4::fit(
      hmm_model,
      emcontrol = depmixS4::em.control(maxit = maxit, tol = tol, random.start = TRUE),
      verbose = TRUE
    )

    if (!is.null(checkpoint_path) && !is.null(checkpoint_meta_path)) {
      saveRDS(fitted_model, checkpoint_path)
      saveRDS(expected_checkpoint_meta, checkpoint_meta_path)
      cat("Saved fitted-model checkpoint: ", checkpoint_path, "\n", sep = "")
    }
  }

  tr_coefs <- lapply(fitted_model@transition, function(ti) ti@parameters$coefficients)
  mu <- matrix(NA_real_, nrow = dimension, ncol = nresp,
               dimnames = list(paste0("State", seq_len(dimension)), obs_names))
  for (k in seq_len(dimension)) {
    for (j in seq_len(nresp)) {
      mu[k, j] <- as.numeric(fitted_model@response[[k]][[j]]@parameters$coefficients)[1]
    }
  }

  mod_val <- depmixS4::depmix(
    response = response_formulas,
    family = replicate(nresp, stats::gaussian(), simplify = FALSE),
    nstates = dimension,
    transition = stats::as.formula(transition_rhs),
    data = val_df,
    ntimes = val_ntimes
  )
  mod_val <- depmixS4::setpars(mod_val, depmixS4::getpars(fitted_model))
  phi_prob <- extract_posterior_probabilities(
    depmixS4::posterior(mod_val, type = "filtering"), dimension
  )

  cov_names <- all.vars(stats::as.formula(transition_rhs))
  cov_df <- cbind(Intercept = 1, val_df[, cov_names, drop = FALSE])
  mse <- rep(NA_real_, 3L)

  for (tau in seq_len(3L)) {
    errors <- list()
    idx <- 1L
    for (i in seq_along(val_ntimes)) {
      len <- val_ntimes[i]
      if (len <= tau) {
        idx <- idx + len
        next
      }
      phi_i <- phi_prob[idx:(idx + len - 1L), , drop = FALSE]
      df_i <- val_df[idx:(idx + len - 1L), , drop = FALSE]
      cov_i <- cov_df[idx:(idx + len - 1L), , drop = FALSE]
      errs_i <- numeric(len - tau)

      for (tt in seq_len(len - tau)) {
        phi_t <- phi_i[tt, , drop = FALSE]
        for (kk in seq_len(tau)) {
          transition_rows <- lapply(tr_coefs, function(eta) {
            linear <- as.numeric(t(eta) %*% as.numeric(cov_i[tt + kk - 1L, ]))
            ex <- exp(linear - max(linear))
            ex / sum(ex)
          })
          A <- do.call(rbind, transition_rows)
          phi_t <- phi_t %*% A
        }
        pred <- as.numeric(phi_t %*% mu)
        obs <- as.numeric(df_i[tt + tau, obs_names, drop = TRUE])
        errs_i[tt] <- sum((pred - obs)^2)
      }
      errors[[length(errors) + 1L]] <- errs_i
      idx <- idx + len
    }
    mse[tau] <- mean(unlist(errors), na.rm = TRUE)
  }

  ll <- depmix_logLik_safe(fitted_model)
  ic <- depmix_information_criteria(ll)
  post_prob <- extract_posterior_probabilities(
    depmixS4::posterior(fitted_model, type = "smoothing"), dimension
  )
  entropy <- -sum(post_prob * log(post_prob + 1e-12))

  icl_approx <- ic$BIC + 2 * entropy

  list(
    fitted_model = fitted_model,
    train_logLik = ll,
    train_AIC = ic$AIC,
    train_BIC = ic$BIC,
    train_ICL_approx = icl_approx,
    val_mse1 = mse[1],
    val_mse2 = mse[2],
    val_mse3 = mse[3],
    posterior_entropy = entropy
  )
}
