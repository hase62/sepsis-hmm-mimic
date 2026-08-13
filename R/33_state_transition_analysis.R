
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
if (!exists("MIMIC_ACTIVE_ANNOT_MASTER", inherits = TRUE)) {
  source(code_path("17_analysis_config.R"))
}

if (!file.exists(code_path("17_analysis_config.R"))) stop("Missing shared analysis configuration: 17_analysis_config.R", call. = FALSE)
source(code_path("17_analysis_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(tibble)
})

ANNO_FILE <- MIMIC_ACTIVE_ANNOT_MASTER
MODEL_RDS_CANDIDATES <- c(Sys.getenv("MIMIC_FINAL_HMM_RDS", unset = ""))
MODEL_RDS_CANDIDATES <- MODEL_RDS_CANDIDATES[nzchar(MODEL_RDS_CANDIDATES)]
AGE_GRID <- c(55, 65, 75)
SEX_GRID <- c(0, 1)
SELECTED_DAYS <- c(1L, 2L, 4L, 7L)
MAX_DAY <- 28L

`%||%` <- function(a, b) if (!is.null(a)) a else b
.stopif <- function(ok, msg) if (!ok) stop(msg, call. = FALSE)

.pick_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  rep(NA, nrow(df))
}

state_id_from_name <- function(x) as.integer(sub("^S", "", x))

preferred_subphenotype_levels <- function(x = character(0)) {
  preferred <- c(
    "A1 low burden",
    "A2 mild burden",
    "B1 renal mild",
    "B2 renal severe",
    "G respiratory",
    "D1 lactate/shock",
    "D2 hepatic",
    "D3 shock/coagulation-overlap",
    "Other/NA"
  )
  x <- x[!is.na(x) & nzchar(x)]
  c(intersect(preferred, x), setdiff(unique(x), preferred))
}

write_both <- function(df, stem) {
  readr::write_csv(df, paste0(stem, ".csv"))
  readr::write_tsv(df, paste0(stem, ".tsv"))
}

get_subphenotype_colors <- function(levels_needed) {
  cols <- NULL
  if (exists("subphenotype_colors", envir = .GlobalEnv)) {
    cols <- get("subphenotype_colors", envir = .GlobalEnv)
  } else if (exists("SUBPHENO_COLORS", envir = .GlobalEnv)) {
    cols <- get("SUBPHENO_COLORS", envir = .GlobalEnv)
  }
  if (!is.null(cols) && is.character(cols) && !is.null(names(cols))) {
    out <- cols[intersect(names(cols), levels_needed)]
  } else {
    base_cols <- grDevices::hcl.colors(max(3, length(levels_needed)), palette = "Dark 3")
    out <- setNames(base_cols[seq_along(levels_needed)], levels_needed)
  }
  miss <- setdiff(levels_needed, names(out))
  if (length(miss) > 0) {
    extra <- grDevices::hcl.colors(length(miss), palette = "Set 2")
    out <- c(out, setNames(extra, miss))
  }
  out[levels_needed]
}

if (file.exists("18_plot_palette.R")) {
  source(code_path("18_plot_palette.R"))
}

.stopif(file.exists(ANNO_FILE), paste("Annotation file not found:", ANNO_FILE))
anno_tbl <- readr::read_csv(ANNO_FILE, show_col_types = FALSE) %>%
  mutate(state_id = as.integer(state_id))

.stopif(all(c("state_id", "subphenotype") %in% names(anno_tbl)),
        "Annotation file must contain state_id and subphenotype.")

subphenotype_levels_order <- preferred_subphenotype_levels(unique(anno_tbl$subphenotype))
subphenotype_colors_use <- get_subphenotype_colors(subphenotype_levels_order)
status_levels <- c("Discharged", "Deceased", "Missing_on_ICU")
status_colors <- c(Discharged = "grey85", Deceased = "grey55", Missing_on_ICU = "grey95")
plot_levels_all <- c(subphenotype_levels_order, status_levels)
plot_fill_colors <- c(subphenotype_colors_use, status_colors)

build_state_subphenotype_map <- function(K, anno_tbl, subphenotype_col = "subphenotype") {
  base <- tibble(
    state_id = seq_len(K),
    state    = paste0("S", state_id)
  )

  df <- base %>%
    left_join(anno_tbl, by = "state_id")

  subphenotype_vec <- df[[subphenotype_col]]
  subphenotype_vec[is.na(subphenotype_vec) | subphenotype_vec == ""] <- "Other/NA"

  tibble(
    state_id = df$state_id,
    state    = df$state,
    subphenotype = subphenotype_vec
  )
}

mimic_promote_final_model(expected_k = MIMIC_FINAL_STATE_DIMENSION, allow_rds = TRUE)

.make_cov_vec <- function(cov_names, newdata) {
  cov_vec <- numeric(length(cov_names))
  names(cov_vec) <- cov_names
  for (nm in cov_names) {
    if (nm %in% c("(Intercept)", "Intercept", "1")) {
      cov_vec[nm] <- 1
    } else {
      cov_vec[nm] <- as.numeric(newdata[[nm]] %||% 0)
    }
  }
  cov_vec
}

.softmax <- function(x) {
  z <- x - max(x, na.rm = TRUE)
  ex <- exp(z)
  ex / sum(ex)
}

.row_or_col_states_logits <- function(eta, cov_vec, K) {
  if (!is.matrix(eta)) stop("transition coefficients are not a matrix.")
  rn <- rownames(eta)
  cn <- colnames(eta)
  looks_state <- function(v) !is.null(v) && length(v) >= K && all(grepl("^(S|St)\\d+$", v))

  if (looks_state(cn)) return(drop(t(eta) %*% cov_vec))
  if (looks_state(rn)) return(drop(eta %*% cov_vec))

  if (!is.null(cn) && length(cn) == length(cov_vec)) return(drop(eta %*% cov_vec))
  if (!is.null(rn) && length(rn) == K) return(drop(t(eta) %*% cov_vec))
  stop("Could not infer transition coefficient orientation from row/column names.")
}

get_transition_matrix <- function(
    fitted_model,
    newdata = list(age = 65, sex = 0, Dialysis = 0, IABP = 0, ECMO = 0)
) {
  K <- length(fitted_model@response)
  A <- matrix(NA_real_, nrow = K, ncol = K)

  eta1 <- fitted_model@transition[[1]]@parameters$coefficients
  if (!is.matrix(eta1)) stop("transition coefficients are not a matrix.")

  rn1 <- rownames(eta1)
  cn1 <- colnames(eta1)
  is_states <- function(v) !is.null(v) && length(v) >= K && all(grepl("^(S|St)\\d+$", v))

  if (!is.null(rn1) && !is_states(rn1)) {
    cov_names <- rn1
  } else if (!is.null(cn1) && !is_states(cn1)) {
    cov_names <- cn1
  } else {
    stop("Could not infer covariate names from transition coefficients.")
  }

  cov_vec <- .make_cov_vec(cov_names, newdata)
  for (from in seq_len(K)) {
    eta <- fitted_model@transition[[from]]@parameters$coefficients
    logits <- .row_or_col_states_logits(eta, cov_vec, K)
    A[from, ] <- .softmax(logits)
  }

  rownames(A) <- paste0("S", seq_len(K))
  colnames(A) <- paste0("S", seq_len(K))
  A
}

aggregate_transition_to_subphenotype <- function(A, state_subphenotype_map) {
  stopifnot(is.matrix(A))

  if (is.null(rownames(A))) rownames(A) <- paste0("S", seq_len(nrow(A)))
  if (is.null(colnames(A))) colnames(A) <- paste0("S", seq_len(ncol(A)))

  sm <- state_subphenotype_map %>%
    mutate(state = factor(state, levels = rownames(A))) %>%
    arrange(state)
  sm$state <- as.character(sm$state)

  missing_states <- setdiff(rownames(A), sm$state)
  if (length(missing_states) > 0) {
    sm <- bind_rows(
      sm,
      tibble(
        state_id = NA_integer_,
        state = missing_states,
        subphenotype = "Other/NA"
      )
    )
  }

  levels_use <- preferred_subphenotype_levels(unique(sm$subphenotype))
  G <- length(levels_use)
  M <- matrix(0, nrow = G, ncol = G, dimnames = list(levels_use, levels_use))

  for (g_from in levels_use) {
    from_states <- sm$state[sm$subphenotype == g_from]
    if (!length(from_states)) next
    A_from <- A[from_states, , drop = FALSE]

    for (g_to in levels_use) {
      to_states <- sm$state[sm$subphenotype == g_to]
      if (!length(to_states)) next
      prob_sum <- sum(A_from[, to_states, drop = FALSE])
      M[g_from, g_to] <- prob_sum / length(from_states)
    }
  }
  M
}

plot_state_transition_matrix <- function(A, title = "State-level transition matrix") {
  dfA <- A %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from") %>%
    tidyr::pivot_longer(cols = -from, names_to = "to", values_to = "prob") %>%
    mutate(
      from = factor(from, levels = rownames(A)),
      to   = factor(to,   levels = colnames(A))
    )

  ggplot(dfA, aes(x = to, y = from, fill = prob)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    coord_fixed() +
    labs(title = title, x = "To state", y = "From state", fill = "P") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

plot_subphenotype_transition_matrix <- function(M, title, out_png) {
  df <- M %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from_subphenotype") %>%
    tidyr::pivot_longer(cols = -from_subphenotype, names_to = "to_subphenotype", values_to = "prob") %>%
    mutate(
      from_subphenotype = factor(from_subphenotype, levels = rownames(M)),
      to_subphenotype = factor(to_subphenotype, levels = colnames(M))
    )

  p <- ggplot(df, aes(x = to_subphenotype, y = from_subphenotype, fill = prob)) +
    geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1)) +
    coord_fixed() +
    labs(title = title, x = "To subphenotype", y = "From subphenotype", fill = "P") +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), panel.grid = element_blank())

  ggsave(filename = out_png, plot = p, width = 7, height = 6, dpi = 600, bg = "white")
  ggsave(filename = sub("\\.png$", ".pdf", out_png), plot = p, width = 7, height = 6, bg = "white", device = cairo_pdf)
  invisible(p)
}

plot_subphenotype_transition_diff <- function(M1, M0, title, out_png) {
  stopifnot(identical(rownames(M1), rownames(M0)))
  stopifnot(identical(colnames(M1), colnames(M0)))
  D <- M1 - M0
  lim <- max(abs(D), na.rm = TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 0.05

  df <- D %>%
    as.data.frame() %>%
    tibble::rownames_to_column("from_subphenotype") %>%
    tidyr::pivot_longer(cols = -from_subphenotype, names_to = "to_subphenotype", values_to = "diff") %>%
    mutate(
      from_subphenotype = factor(from_subphenotype, levels = rownames(M1)),
      to_subphenotype = factor(to_subphenotype, levels = colnames(M1))
    )

  p <- ggplot(df, aes(x = to_subphenotype, y = from_subphenotype, fill = diff)) +
    geom_tile() +
    scale_fill_gradient2(limits = c(-lim, lim), midpoint = 0, name = "ΔP") +
    coord_fixed() +
    labs(title = title, x = "To subphenotype", y = "From subphenotype") +
    theme_bw(base_size = 13) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), panel.grid = element_blank())

  ggsave(filename = out_png, plot = p, width = 7, height = 6, dpi = 600, bg = "white")
  ggsave(filename = sub("\\.png$", ".pdf", out_png), plot = p, width = 7, height = 6, bg = "white", device = cairo_pdf)
  invisible(p)
}

compute_stationary_row <- function(M, tol = 1e-10, maxit = 10000) {
  stopifnot(is.matrix(M))
  if (is.null(rownames(M))) rownames(M) <- paste0("meta", seq_len(nrow(M)))
  G <- nrow(M)
  pi <- rep(1 / G, G)
  for (it in seq_len(maxit)) {
    pi_new <- as.numeric(pi %*% M)
    if (max(abs(pi_new - pi)) < tol) {
      pi <- pi_new
      break
    }
    pi <- pi_new
  }
  pi <- pmax(pi, 0)
  pi <- pi / sum(pi)
  names(pi) <- rownames(M)
  pi
}

fitted_model_loaded <- FALSE
if (!exists("fitted_model", envir = .GlobalEnv)) {
  rds_use <- MODEL_RDS_CANDIDATES[file.exists(MODEL_RDS_CANDIDATES)][1]
  if (!is.na(rds_use) && nzchar(rds_use)) {
    fitted_model <- readRDS(rds_use)
    fitted_model_loaded <- TRUE
  }
} else {
  fitted_model <- get("fitted_model", envir = .GlobalEnv)
  fitted_model_loaded <- TRUE
}

if (fitted_model_loaded) {
  K <- length(fitted_model@response)
  state_subphenotype_map <- build_state_subphenotype_map(K, anno_tbl, subphenotype_col = "subphenotype")

  cov_grid <- tidyr::expand_grid(age = AGE_GRID, sex = SEX_GRID) %>%
    mutate(Dialysis = 0, IABP = 0, ECMO = 0, scenario = sprintf("age%d_sex%d", age, sex))

  subphenotype_mat_list <- list()
  state_mat_list <- list()
  subphenotype_long_list <- list()

  for (i in seq_len(nrow(cov_grid))) {
    row_i <- cov_grid[i, ]
    cov_list <- as.list(row_i[c("age", "sex", "Dialysis", "IABP", "ECMO")])
    A_state <- get_transition_matrix(fitted_model, newdata = cov_list)
    M_subphenotype <- aggregate_transition_to_subphenotype(A_state, state_subphenotype_map)
    label <- row_i$scenario

    state_mat_list[[label]] <- A_state
    subphenotype_mat_list[[label]] <- M_subphenotype

    df_long <- M_subphenotype %>%
      as.data.frame() %>%
      tibble::rownames_to_column("from_subphenotype") %>%
      tidyr::pivot_longer(cols = -from_subphenotype, names_to = "to_subphenotype", values_to = "p_transition") %>%
      mutate(age = row_i$age, sex = row_i$sex, scenario = label)
    subphenotype_long_list[[label]] <- df_long

    readr::write_csv(
      A_state %>% as.data.frame() %>% tibble::rownames_to_column("from_state"),
      sprintf("transition_state_transition_matrix_%s.csv", label)
    )
    readr::write_csv(
      M_subphenotype %>% as.data.frame() %>% tibble::rownames_to_column("from_subphenotype"),
      sprintf("transition_subphenotype_transition_matrix_%s.csv", label)
    )
  }

  subphenotype_transition_long <- dplyr::bind_rows(subphenotype_long_list)
  readr::write_csv(subphenotype_transition_long, "transition_subphenotype_transition_prob_by_age_sex.csv")

  if (!is.null(subphenotype_mat_list[["age65_sex0"]])) {
    plot_subphenotype_transition_matrix(
      subphenotype_mat_list[["age65_sex0"]],
      title = "Subphenotype transitions @ age=65, sex=0 (male)",
      out_png = "transition_subphenotype_transition_heatmap_age65_sex0.png"
    )
  }

  if (!is.null(state_mat_list[["age65_sex0"]])) {
    p_state_65M <- plot_state_transition_matrix(
      state_mat_list[["age65_sex0"]],
      title = "State-level transitions @ age=65, sex=0 (male)"
    )
    ggsave("transition_state_transition_heatmap_age65_sex0.png", p_state_65M, width = 8, height = 7, dpi = 600, bg = "white")
    ggsave("transition_state_transition_heatmap_age65_sex0.pdf", p_state_65M, width = 8, height = 7, bg = "white", device = cairo_pdf)
  }

  if (all(c("age65_sex0", "age65_sex1") %in% names(subphenotype_mat_list))) {
    plot_subphenotype_transition_diff(
      M1 = subphenotype_mat_list[["age65_sex1"]],
      M0 = subphenotype_mat_list[["age65_sex0"]],
      title = "ΔP (female - male) @ age=65",
      out_png = "transition_subphenotype_transition_diff_sex_age65.png"
    )
  }

  if (all(c("age55_sex0", "age75_sex0") %in% names(subphenotype_mat_list))) {
    plot_subphenotype_transition_diff(
      M1 = subphenotype_mat_list[["age75_sex0"]],
      M0 = subphenotype_mat_list[["age55_sex0"]],
      title = "ΔP (age 75 - age 55) @ sex=0 (male)",
      out_png = "transition_subphenotype_transition_diff_age75_55_sex0.png"
    )
  }

  self_transitions <- purrr::imap_dfr(
    subphenotype_mat_list,
    ~ tibble(scenario = .y, from_subphenotype = rownames(.x), p_self = diag(.x))
  ) %>%
    left_join(cov_grid %>% select(scenario, age, sex), by = "scenario") %>%
    arrange(from_subphenotype, age, sex)
  readr::write_csv(self_transitions, "transition_subphenotype_self_transition_by_age_sex.csv")

  stationary_df <- purrr::imap_dfr(
    subphenotype_mat_list,
    ~ {
      pi <- compute_stationary_row(.x)
      tibble(scenario = .y, subphenotype = names(pi), pi_stationary = as.numeric(pi))
    }
  ) %>%
    left_join(cov_grid %>% select(scenario, age, sex), by = "scenario") %>%
    arrange(age, sex, subphenotype)
  readr::write_csv(stationary_df, "transition_subphenotype_stationary_distribution_by_age_sex.csv")

  if ("age65_sex0" %in% stationary_df$scenario) {
    df_plot <- stationary_df %>%
      filter(scenario == "age65_sex0") %>%
      mutate(subphenotype = factor(subphenotype, levels = subphenotype_levels_order))
    p_stat_65M <- ggplot(df_plot, aes(x = subphenotype, y = pi_stationary, fill = subphenotype)) +
      geom_col() +
      scale_fill_manual(values = subphenotype_colors_use, guide = "none") +
      labs(title = "Stationary distribution of subphenotypes @ age=65, sex=0 (male)", x = "Subphenotype", y = "Stationary probability") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave("transition_subphenotype_stationary_age65_sex0_barplot.png", p_stat_65M, width = 8, height = 5, dpi = 600, bg = "white")
    ggsave("transition_subphenotype_stationary_age65_sex0_barplot.pdf", p_stat_65M, width = 8, height = 5, device = cairo_pdf, bg = "white")
  }

  message("transition_analysis model-based transition section: done.")
} else {
  message("transition_analysis model-based transition section: skipped (fitted_model / RDS not found).")
}

CORE_TRAJ_DAYS <- c(1L, 2L, 4L, 7L)
EXTENDED_TRAJ_DAYS <- c(1L, 2L, 4L, 7L, 10L, 14L, 21L, 28L)
DAILY_TRAJ_DAYS <- seq_len(MAX_DAY)

.day_label <- function(d) {
  ifelse(as.integer(d) == 1L, "ICU admission", paste0("ICU day ", as.integer(d)))
}

.have_empirical <- all(vapply(
  c("post_smooth_full", "combined", "ids", "icustays_intubated_first_patient_admission"),
  exists,
  logical(1),
  envir = .GlobalEnv
))

if (.have_empirical) {
  post_smooth_full <- get("post_smooth_full", envir = .GlobalEnv)
  post_smooth_full <- as.matrix(post_smooth_full)
  storage.mode(post_smooth_full) <- "double"
  combined <- get("combined", envir = .GlobalEnv)
  ids <- as.character(get("ids", envir = .GlobalEnv))
  icu <- get("icustays_intubated_first_patient_admission", envir = .GlobalEnv) %>%
    mutate(stay_id = as.character(stay_id))

  .stopif(nrow(post_smooth_full) == length(ids), "length(ids) must equal nrow(post_smooth_full).")
  .stopif(nrow(post_smooth_full) == nrow(combined), "nrow(post_smooth_full) must equal nrow(combined).")
  .stopif("lab_time" %in% names(combined), "combined must contain lab_time.")

  state_levels <- colnames(post_smooth_full)
  if (is.null(state_levels)) state_levels <- paste0("S", seq_len(ncol(post_smooth_full)))
  state_id_tbl <- tibble(state = state_levels, state_id = state_id_from_name(state_levels))
  state_meta <- state_id_tbl %>% left_join(anno_tbl %>% select(state_id, subphenotype), by = "state_id")
  state_meta <- state_meta %>% mutate(subphenotype = ifelse(is.na(subphenotype) | subphenotype == "", "Other/NA", subphenotype))

  subph_levels_emp <- preferred_subphenotype_levels(unique(state_meta$subphenotype))
  G <- length(subph_levels_emp)
  K <- nrow(state_meta)
  M_state_to_subph <- matrix(0, nrow = K, ncol = G, dimnames = list(state_meta$state, subph_levels_emp))
  for (j in seq_len(K)) {
    g <- match(state_meta$subphenotype[j], subph_levels_emp)
    if (!is.na(g)) M_state_to_subph[j, g] <- 1
  }

  icu_intime <- .pick_col(icu, c("intime", "icu_intime", "ICU_intime", "icu_in_time", "in_time", "admittime", "admit_time"))
  outtime <- .pick_col(icu, c("outtime", "icu_outtime", "ICU_outtime", "out_time", "dischtime", "discharge_time"))
  icu_status_tbl <- icu %>%
    transmute(
      stay_id = as.character(stay_id),
      dod = dod,
      icu_intime = icu_intime,
      outtime = outtime
    ) %>%
    mutate(
      death_days = ifelse(!is.na(dod) & !is.na(icu_intime), as.numeric(difftime(dod, icu_intime, units = "days")), NA_real_),
      los_days = ifelse(!is.na(outtime) & !is.na(icu_intime), as.numeric(difftime(outtime, icu_intime, units = "days")), NA_real_),
      death_days = ifelse(!is.na(death_days) & death_days < 0, NA_real_, death_days),
      los_days = ifelse(!is.na(los_days) & los_days < 0, NA_real_, los_days)
    )

  assign_status_at_day <- function(stay_ids, day_num, dominant_tbl, status_tbl) {
    dom <- dominant_tbl %>% filter(day == day_num) %>% select(stay_id, dominant_subphenotype)
    out <- tibble(stay_id = stay_ids) %>%
      left_join(dom, by = "stay_id") %>%
      left_join(status_tbl %>% select(stay_id, death_days, los_days), by = "stay_id") %>%
      mutate(
        status = case_when(
          !is.na(dominant_subphenotype) ~ dominant_subphenotype,
          !is.na(death_days) & death_days < day_num ~ "Deceased",
          !is.na(los_days) & los_days < day_num ~ "Discharged",
          TRUE ~ "Missing_on_ICU"
        )
      ) %>%
      select(stay_id, status)
    out
  }

  compute_dominant_by_days <- function(days_use) {
    day_vec <- suppressWarnings(as.integer(combined$lab_time))
    stay_day_dominant_list <- vector("list", length(days_use))
    names(stay_day_dominant_list) <- paste0("day", days_use)

    for (d in days_use) {
      mask_d <- !is.na(day_vec) & day_vec == d
      if (!any(mask_d)) {
        stay_day_dominant_list[[paste0("day", d)]] <- tibble(
          stay_id = character(0),
          day = integer(0),
          dominant_subphenotype = character(0),
          p_dominant = numeric(0)
        )
        next
      }

      Wd <- post_smooth_full[mask_d, , drop = FALSE]
      rownames(Wd) <- ids[mask_d]
      Psub_row <- Wd %*% M_state_to_subph
      stay_sum <- rowsum(Psub_row, group = ids[mask_d], reorder = FALSE)
      stay_n <- as.numeric(table(ids[mask_d])[rownames(stay_sum)])
      stay_avg <- stay_sum / stay_n
      dom_idx <- max.col(stay_avg, ties.method = "first")
      dom_name <- colnames(stay_avg)[dom_idx]
      dom_prob <- stay_avg[cbind(seq_len(nrow(stay_avg)), dom_idx)]

      stay_day_dominant_list[[paste0("day", d)]] <- tibble(
        stay_id = rownames(stay_avg),
        day = d,
        dominant_subphenotype = dom_name,
        p_dominant = as.numeric(dom_prob)
      )
    }

    bind_rows(stay_day_dominant_list)
  }

  build_traj_wide <- function(stay_ids, days_use, dominant_tbl, status_tbl) {
    status_by_day_list <- lapply(days_use, function(d) {
      assign_status_at_day(stay_ids, d, dominant_tbl, status_tbl) %>%
        rename(!!paste0("day", d) := status)
    })
    wide <- Reduce(function(x, y) dplyr::left_join(x, y, by = "stay_id"), status_by_day_list)
    wide %>% mutate(admission_subphenotype = .data[[paste0("day", days_use[1])]])
  }

  build_traj_long <- function(traj_wide, days_use) {
    traj_wide %>%
      pivot_longer(cols = paste0("day", days_use), names_to = "timepoint", values_to = "status") %>%
      mutate(
        day = as.integer(readr::parse_number(timepoint)),
        time_label = .day_label(day),
        status = factor(status, levels = plot_levels_all)
      )
  }

  build_counts_by_day <- function(traj_long, active_levels) {
    counts_all <- traj_long %>%
      count(day, time_label, status, name = "n") %>%
      group_by(day) %>%
      mutate(prop_all_stays = ifelse(sum(n) > 0, n / sum(n), NA_real_)) %>%
      ungroup()

    active_totals <- counts_all %>%
      filter(as.character(status) %in% active_levels) %>%
      group_by(day) %>%
      summarise(n_active_stays = sum(n), .groups = "drop")

    counts_all %>%
      left_join(active_totals, by = "day") %>%
      mutate(
        prop_active_stays = ifelse(as.character(status) %in% active_levels & !is.na(n_active_stays) & n_active_stays > 0,
                                   n / n_active_stays, NA_real_)
      )
  }

  build_adjacent_transition_counts <- function(traj_wide, days_use) {
    adjacent_pairs <- tibble(day_from = days_use[-length(days_use)], day_to = days_use[-1])
    purrr::pmap_dfr(adjacent_pairs, function(day_from, day_to) {
      from_col <- paste0("day", day_from)
      to_col <- paste0("day", day_to)
      traj_wide %>%
        count(from = .data[[from_col]], to = .data[[to_col]], name = "n") %>%
        mutate(day_from = day_from, day_to = day_to, time_from = .day_label(day_from), time_to = .day_label(day_to))
    })
  }

  build_change_overall <- function(traj_wide, base_day, target_days, active_levels) {
    base_col <- paste0("day", base_day)
    purrr::map_dfr(target_days[target_days != base_day], function(d) {
      dcol <- paste0("day", d)
      tmp <- traj_wide %>%
        transmute(admission = .data[[base_col]], current = .data[[dcol]]) %>%
        filter(admission %in% active_levels, current %in% active_levels)
      tibble(
        day = d,
        time_label = .day_label(d),
        n_active = nrow(tmp),
        n_same = sum(tmp$admission == tmp$current, na.rm = TRUE),
        n_changed = sum(tmp$admission != tmp$current, na.rm = TRUE),
        p_changed = ifelse(n_active > 0, n_changed / n_active, NA_real_),
        p_same = ifelse(n_active > 0, n_same / n_active, NA_real_)
      )
    })
  }

  build_change_by_admission <- function(traj_wide, base_day, target_days, active_levels) {
    base_col <- paste0("day", base_day)
    purrr::map_dfr(target_days[target_days != base_day], function(d) {
      dcol <- paste0("day", d)
      traj_wide %>%
        transmute(admission = .data[[base_col]], current = .data[[dcol]]) %>%
        filter(admission %in% active_levels, current %in% active_levels) %>%
        group_by(admission) %>%
        summarise(
          day = d,
          time_label = .day_label(d),
          n_active = n(),
          n_same = sum(admission == current, na.rm = TRUE),
          n_changed = sum(admission != current, na.rm = TRUE),
          p_changed = ifelse(n_active > 0, n_changed / n_active, NA_real_),
          p_same = ifelse(n_active > 0, n_same / n_active, NA_real_),
          .groups = "drop"
        )
    }) %>%
      mutate(admission = factor(admission, levels = active_levels))
  }

  build_admission_destination_tables <- function(traj_wide, base_day, target_days, active_levels) {
    base_col <- paste0("day", base_day)
    base_df <- traj_wide %>%
      transmute(stay_id, admission = .data[[base_col]]) %>%
      filter(admission %in% active_levels)

    allstatus_tbl <- purrr::map_dfr(target_days, function(d) {
      dcol <- paste0("day", d)
      base_df %>%
        left_join(traj_wide %>% transmute(stay_id, destination = .data[[dcol]]), by = "stay_id") %>%
        count(admission, destination, name = "n") %>%
        group_by(admission) %>%
        mutate(
          day = d,
          time_label = .day_label(d),
          n_admission_total = sum(n),
          prop_within_admission_all = ifelse(n_admission_total > 0, n / n_admission_total, NA_real_)
        ) %>%
        ungroup()
    })

    active_tbl <- allstatus_tbl %>%
      filter(destination %in% active_levels) %>%
      group_by(day, time_label, admission) %>%
      mutate(
        n_active_destination = sum(n),
        prop_within_admission_active = ifelse(n_active_destination > 0, n / n_active_destination, NA_real_)
      ) %>%
      ungroup()

    list(allstatus = allstatus_tbl, active = active_tbl)
  }

  stay_day_dominant_daily <- compute_dominant_by_days(DAILY_TRAJ_DAYS)
  write_both(stay_day_dominant_daily, "transition_table_stay_day_dominant_subphenotype_daily")

  all_stays <- sort(unique(as.character(icu_status_tbl$stay_id)))
  if (!length(all_stays)) all_stays <- sort(unique(ids))

  traj_wide_selected <- build_traj_wide(all_stays, EXTENDED_TRAJ_DAYS, stay_day_dominant_daily, icu_status_tbl)
  traj_long_selected <- build_traj_long(traj_wide_selected, EXTENDED_TRAJ_DAYS)
  write_both(traj_wide_selected, "transition_table_subphenotype_trajectory_status_wide_extended_days")
  write_both(traj_long_selected, "transition_table_subphenotype_trajectory_status_long_extended_days")

  stay_day_dominant_selected <- stay_day_dominant_daily %>% filter(day %in% EXTENDED_TRAJ_DAYS)
  write_both(stay_day_dominant_selected, "transition_table_stay_day_dominant_subphenotype_extended_days")

  counts_by_day_selected <- build_counts_by_day(traj_long_selected, subph_levels_emp)
  write_both(counts_by_day_selected, "transition_table_subphenotype_trajectory_counts_by_day_extended_days")

  traj_wide_daily <- build_traj_wide(all_stays, DAILY_TRAJ_DAYS, stay_day_dominant_daily, icu_status_tbl)
  write_both(traj_wide_daily, "transition_table_subphenotype_trajectory_status_wide_daily")

  traj_long_daily <- purrr::map_dfr(DAILY_TRAJ_DAYS, function(d) {
    assign_status_at_day(all_stays, d, stay_day_dominant_daily, icu_status_tbl) %>%
      mutate(day = d, time_label = .day_label(d), status = factor(status, levels = plot_levels_all))
  })
  write_both(traj_long_daily, "transition_table_subphenotype_trajectory_status_long_daily")

  counts_by_day_daily <- build_counts_by_day(traj_long_daily, subph_levels_emp)
  write_both(counts_by_day_daily, "transition_table_subphenotype_trajectory_counts_by_day_daily")

  adjacent_daily_allstatus <- build_adjacent_transition_counts(traj_wide_daily, DAILY_TRAJ_DAYS) %>%
    group_by(day_from, time_from, from) %>%
    mutate(
      n_from_total = sum(n),
      prop_to_allstatus = ifelse(n_from_total > 0, n / n_from_total, NA_real_)
    ) %>%
    ungroup()
  write_both(adjacent_daily_allstatus, "transition_table_subphenotype_nextday_transition_allstatus_daily")

  adjacent_daily_active <- adjacent_daily_allstatus %>%
    filter(from %in% subph_levels_emp, to %in% subph_levels_emp) %>%
    group_by(day_from, time_from, from) %>%
    mutate(
      n_from_active = sum(n),
      prop_to_active = ifelse(n_from_active > 0, n / n_from_active, NA_real_)
    ) %>%
    ungroup() %>%
    mutate(
      from = factor(from, levels = subph_levels_emp),
      to = factor(to, levels = subph_levels_emp)
    )
  write_both(adjacent_daily_active, "transition_table_subphenotype_nextday_transition_active_daily")

  adjacent_daily_summary <- adjacent_daily_active %>%
    group_by(day_from, time_from, from) %>%
    summarise(
      n_from_active = unique(n_from_active)[1],
      p_stay_same = sum(prop_to_active[as.character(to) == as.character(from)], na.rm = TRUE),
      top_next_subphenotype = if (all(is.na(prop_to_active) | as.character(to) == as.character(from))) NA_character_ else as.character(to[which.max(ifelse(as.character(to) == as.character(from), -Inf, prop_to_active))]),
      p_top_next = if (all(is.na(prop_to_active) | as.character(to) == as.character(from))) NA_real_ else max(ifelse(as.character(to) == as.character(from), NA_real_, prop_to_active), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(from = factor(from, levels = subph_levels_emp))
  write_both(adjacent_daily_summary, "transition_table_subphenotype_nextday_transition_summary_daily")

  transition_pairs_df <- build_adjacent_transition_counts(traj_wide_selected, EXTENDED_TRAJ_DAYS)
  write_both(transition_pairs_df, "transition_table_subphenotype_transition_counts_adjacent_extended_days")

  overall_change_tbl <- build_change_overall(traj_wide_selected, EXTENDED_TRAJ_DAYS[1], EXTENDED_TRAJ_DAYS, subph_levels_emp)
  write_both(overall_change_tbl, "transition_table_subphenotype_transition_change_overall_extended_days")

  subtype_change_tbl <- build_change_by_admission(traj_wide_selected, EXTENDED_TRAJ_DAYS[1], EXTENDED_TRAJ_DAYS, subph_levels_emp)
  write_both(subtype_change_tbl, "transition_table_subphenotype_transition_change_by_admission_extended_days")

  admission_dest <- build_admission_destination_tables(traj_wide_selected, EXTENDED_TRAJ_DAYS[1], EXTENDED_TRAJ_DAYS, subph_levels_emp)
  write_both(admission_dest$allstatus, "transition_table_admission_to_day_destination_allstatus_extended_days")
  write_both(admission_dest$active, "transition_table_admission_to_day_destination_active_extended_days")

  write_both(
    admission_dest$allstatus %>% filter(day == 2L),
    "transition_table_admission_to_day2_destination_allstatus"
  )
  write_both(
    admission_dest$active %>% filter(day == 2L),
    "transition_table_admission_to_day2_destination_active"
  )

  traj_core <- traj_wide_selected %>%
    transmute(
      day1 = factor(.data[["day1"]], levels = plot_levels_all),
      day2 = factor(.data[["day2"]], levels = plot_levels_all),
      day4 = factor(.data[["day4"]], levels = plot_levels_all),
      day7 = factor(.data[["day7"]], levels = plot_levels_all)
    ) %>%
    count(day1, day2, day4, day7, name = "Freq")
  write_both(traj_core, "transition_table_subphenotype_trajectory_alluvial_counts_core_days")

  if (requireNamespace("ggalluvial", quietly = TRUE)) {
    p_alluvial <- ggplot(
      traj_core,
      aes(axis1 = day1, axis2 = day2, axis3 = day4, axis4 = day7, y = Freq)
    ) +
      ggalluvial::geom_alluvium(aes(fill = day1), width = 1/12, alpha = 0.75, knot.pos = 0.45) +
      ggalluvial::geom_stratum(width = 1/8, fill = "white", color = "grey35") +
      ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
      scale_x_discrete(
        limits = c("day1", "day2", "day4", "day7"),
        labels = c("ICU admission", "ICU day 2", "ICU day 4", "ICU day 7"),
        expand = c(0.08, 0.02)
      ) +
      scale_y_continuous(expand = c(0.01, 0.01)) +
      scale_fill_manual(values = plot_fill_colors, drop = FALSE, guide = "none") +
      labs(x = NULL, y = "Number of patient-stays", title = "Subphenotype trajectories across ICU admission / day 2 / day 4 / day 7") +
      theme_bw(base_size = 12) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

    ggsave("transition_figure_subphenotype_trajectory_alluvial_core_days.png", p_alluvial, width = 10, height = 6, dpi = 300, bg = "white")
    ggsave("transition_figure_subphenotype_trajectory_alluvial_core_days.pdf", p_alluvial, width = 10, height = 6, bg = "white", device = cairo_pdf)

    facet_dat <- traj_core %>% filter(as.character(day1) %in% subph_levels_emp)
    p_alluvial_facet <- ggplot(
      facet_dat,
      aes(axis1 = day1, axis2 = day2, axis3 = day4, axis4 = day7, y = Freq)
    ) +
      ggalluvial::geom_alluvium(aes(fill = day1), width = 1/12, alpha = 0.75, knot.pos = 0.45) +
      ggalluvial::geom_stratum(width = 1/8, fill = "white", color = "grey35") +
      ggplot2::geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6) +
      scale_x_discrete(
        limits = c("day1", "day2", "day4", "day7"),
        labels = c("ICU admission", "ICU day 2", "ICU day 4", "ICU day 7"),
        expand = c(0.08, 0.02)
      ) +
      scale_y_continuous(expand = c(0.01, 0.01)) +
      scale_fill_manual(values = plot_fill_colors, drop = FALSE, guide = "none") +
      facet_wrap(~ day1, ncol = 2, scales = "free_y") +
      labs(x = NULL, y = "Number of patient-stays", title = "Subphenotype trajectories stratified by admission subphenotype (core days)") +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

    ggsave("transition_figure_subphenotype_trajectory_alluvial_by_admission_core_days.png", p_alluvial_facet, width = 12, height = 8, dpi = 300, bg = "white")
    ggsave("transition_figure_subphenotype_trajectory_alluvial_by_admission_core_days.pdf", p_alluvial_facet, width = 12, height = 8, bg = "white", device = cairo_pdf)
  } else {
    message("ggalluvial is not installed; core alluvial figures were skipped. Tables were still written.")
  }

  p_stack_ext <- counts_by_day_selected %>%
    mutate(
      status = factor(status, levels = rev(plot_levels_all)),
      time_label = factor(time_label, levels = .day_label(EXTENDED_TRAJ_DAYS))
    ) %>%
    ggplot(aes(x = time_label, y = n, fill = status)) +
    geom_col(color = "grey20", linewidth = 0.2) +
    scale_fill_manual(values = plot_fill_colors, drop = FALSE) +
    labs(x = NULL, y = "Number of patient-stays", fill = NULL, title = "Subtype/status distribution across extended ICU days") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("transition_figure_subphenotype_trajectory_stackedbar_extended_days.png", p_stack_ext, width = 11, height = 6, dpi = 300, bg = "white")
  ggsave("transition_figure_subphenotype_trajectory_stackedbar_extended_days.pdf", p_stack_ext, width = 11, height = 6, bg = "white", device = cairo_pdf)

  p_line_daily <- counts_by_day_daily %>%
    filter(as.character(status) %in% subph_levels_emp) %>%
    mutate(status = factor(as.character(status), levels = subph_levels_emp)) %>%
    ggplot(aes(x = day, y = prop_active_stays, color = status)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.2) +
    scale_color_manual(values = subphenotype_colors_use, drop = FALSE) +
    scale_x_continuous(breaks = c(1, 2, 4, 7, 10, 14, 21, 28)) +
    labs(x = "ICU day", y = "Proportion among active ICU stays", color = NULL,
         title = "Daily subphenotype composition among active ICU stays") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave("transition_figure_subphenotype_active_proportion_by_day_daily.png", p_line_daily, width = 10, height = 5.5, dpi = 300, bg = "white")
  ggsave("transition_figure_subphenotype_active_proportion_by_day_daily.pdf", p_line_daily, width = 10, height = 5.5, bg = "white", device = cairo_pdf)

  p_next_heat <- adjacent_daily_active %>%
    mutate(
      from = factor(as.character(from), levels = subph_levels_emp),
      to = factor(as.character(to), levels = rev(subph_levels_emp))
    ) %>%
    ggplot(aes(x = day_from, y = to, fill = prop_to_active)) +
    geom_tile(color = "white") +
    facet_wrap(~ from, ncol = 3) +
    scale_fill_viridis_c(limits = c(0, 1), na.value = "grey95") +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 7, 10, 14, 21, 27)) +
    labs(x = "ICU day d", y = "Subtype at day d+1", fill = "Prop",
         title = "Among patients in each subtype on day d, where do they go on day d+1?") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("transition_figure_subphenotype_nextday_transition_heatmap_active_daily.png", p_next_heat, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave("transition_figure_subphenotype_nextday_transition_heatmap_active_daily.pdf", p_next_heat, width = 12, height = 8, bg = "white", device = cairo_pdf)

  p_next_stack <- adjacent_daily_active %>%
    mutate(from = factor(as.character(from), levels = subph_levels_emp)) %>%
    ggplot(aes(x = day_from, y = prop_to_active, fill = to)) +
    geom_col(color = "grey20", linewidth = 0.1) +
    facet_wrap(~ from, ncol = 3) +
    scale_fill_manual(values = subphenotype_colors_use, drop = FALSE) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 7, 10, 14, 21, 27)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "ICU day d", y = "Destination mix at day d+1", fill = "Subtype at day d+1",
         title = "Next-day destination mix for each current-day subtype") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave("transition_figure_subphenotype_nextday_transition_stacked_active_daily.png", p_next_stack, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave("transition_figure_subphenotype_nextday_transition_stacked_active_daily.pdf", p_next_stack, width = 12, height = 8, bg = "white", device = cairo_pdf)

  p_next_same <- adjacent_daily_summary %>%
    ggplot(aes(x = day_from, y = p_stay_same, color = from)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.1) +
    scale_color_manual(values = subphenotype_colors_use, drop = FALSE) +
    scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 7, 10, 14, 21, 27)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "ICU day d", y = "Remain in same subtype at day d+1", color = "Subtype at day d",
         title = "Daily next-day persistence by current subtype") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave("transition_figure_subphenotype_nextday_persistence_daily.png", p_next_same, width = 10, height = 5.5, dpi = 300, bg = "white")
  ggsave("transition_figure_subphenotype_nextday_persistence_daily.pdf", p_next_same, width = 10, height = 5.5, bg = "white", device = cairo_pdf)

  heat_active <- admission_dest$active %>%
    filter(day != EXTENDED_TRAJ_DAYS[1]) %>%
    mutate(
      admission = factor(admission, levels = subph_levels_emp),
      destination = factor(destination, levels = subph_levels_emp),
      time_label = factor(time_label, levels = .day_label(EXTENDED_TRAJ_DAYS[EXTENDED_TRAJ_DAYS != EXTENDED_TRAJ_DAYS[1]]))
    )

  p_heat <- ggplot(heat_active, aes(x = destination, y = admission, fill = prop_within_admission_active)) +
    geom_tile(color = "white") +
    geom_text(aes(label = ifelse(is.na(prop_within_admission_active), "", sprintf("%.2f", prop_within_admission_active))), size = 2.8) +
    facet_wrap(~ time_label, ncol = 3) +
    scale_fill_viridis_c(limits = c(0, 1), na.value = "grey95") +
    labs(x = "Destination subtype", y = "Admission subtype", fill = "Prop",
         title = "Where each admission subtype goes among active ICU stays") +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("transition_figure_admission_destination_heatmap_active_extended_days.png", p_heat, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave("transition_figure_admission_destination_heatmap_active_extended_days.pdf", p_heat, width = 12, height = 8, bg = "white", device = cairo_pdf)

  stack_by_adm <- admission_dest$allstatus %>%
    mutate(
      admission = factor(admission, levels = subph_levels_emp),
      destination = factor(destination, levels = plot_levels_all),
      time_label = factor(time_label, levels = .day_label(EXTENDED_TRAJ_DAYS))
    )
  p_stack_by_adm <- ggplot(stack_by_adm, aes(x = time_label, y = prop_within_admission_all, fill = destination)) +
    geom_col(color = "grey20", linewidth = 0.15) +
    facet_wrap(~ admission, ncol = 3) +
    scale_fill_manual(values = plot_fill_colors, drop = FALSE) +
    labs(x = NULL, y = "Proportion within admission subtype", fill = NULL,
         title = "Longer-horizon destination mix for each admission subtype") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave("transition_figure_admission_destination_stacked_allstatus_extended_days.png", p_stack_by_adm, width = 12, height = 8, dpi = 300, bg = "white")
  ggsave("transition_figure_admission_destination_stacked_allstatus_extended_days.pdf", p_stack_by_adm, width = 12, height = 8, bg = "white", device = cairo_pdf)

  p_same <- subtype_change_tbl %>%
    mutate(time_label = factor(time_label, levels = .day_label(EXTENDED_TRAJ_DAYS[EXTENDED_TRAJ_DAYS != EXTENDED_TRAJ_DAYS[1]]))) %>%
    ggplot(aes(x = day, y = p_same, color = admission)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.3) +
    scale_color_manual(values = subphenotype_colors_use, drop = FALSE) +
    scale_x_continuous(breaks = EXTENDED_TRAJ_DAYS[EXTENDED_TRAJ_DAYS != EXTENDED_TRAJ_DAYS[1]]) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "ICU day", y = "Same subtype as admission (among active stays)", color = "Admission subtype",
         title = "Subtype persistence from ICU admission across extended days") +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  ggsave("transition_figure_subphenotype_persistence_by_admission_extended_days.png", p_same, width = 10, height = 5.5, dpi = 300, bg = "white")
  ggsave("transition_figure_subphenotype_persistence_by_admission_extended_days.pdf", p_same, width = 10, height = 5.5, bg = "white", device = cairo_pdf)

  message("transition_analysis empirical trajectory section: done.")
} else {
  message("Empirical trajectory section skipped because the final-model objects are not loaded.")
}

message("state transition analysis: done.")
