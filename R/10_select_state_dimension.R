
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

n_iterations <- if (length(args) >= 1L) as.integer(args[[1L]]) else 20L
target_h <- as.integer(Sys.getenv("MIMIC_TARGET_H", unset = "3"))

if (!target_h %in% 1:3) {
  stop("MIMIC_TARGET_H must be one of 1, 2, or 3.")
}
if (!is.finite(n_iterations) || n_iterations < 2L) {
  stop("At least two split iterations are required to estimate an SE.")
}

data_dir <- Sys.getenv("MIMIC_MODEL_SELECTION_DIR", unset = ".")
out_dir <- Sys.getenv(
  "MIMIC_MODEL_SELECTION_OUT",
  unset = file.path(data_dir, "model_selection_outputs")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

parse_k_values <- function() {
  explicit <- trimws(Sys.getenv("MIMIC_K_VALUES", unset = ""))
  if (nzchar(explicit)) {
    vals <- suppressWarnings(as.integer(strsplit(explicit, ",", fixed = TRUE)[[1L]]))
    vals <- sort(unique(vals[is.finite(vals) & vals >= 2L]))
    if (!length(vals)) {
      stop("MIMIC_K_VALUES was supplied but no valid state numbers were parsed.")
    }
    return(vals)
  }

  k_min <- as.integer(Sys.getenv("MIMIC_K_MIN", unset = "5"))
  k_max <- as.integer(Sys.getenv("MIMIC_K_MAX", unset = "65"))
  if (!is.finite(k_min) || !is.finite(k_max) || k_min > k_max) {
    stop("Invalid MIMIC_K_MIN/MIMIC_K_MAX.")
  }
  seq.int(k_min, k_max)
}

expected_k <- parse_k_values()

metric_pattern <- "^model_selection_K[0-9]+_split[0-9]+\\.metrics\\.csv$"

files <- list.files(
  path = data_dir,
  pattern = metric_pattern,
  full.names = TRUE
)

if (!length(files)) {
  stop(
    "No metrics files found in ", normalizePath(data_dir, mustWork = FALSE),
    " using pattern: ", metric_pattern
  )
}

read_one <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  x$source_file <- basename(path)
  x
}

raw_all <- bind_rows(lapply(files, read_one))

required_columns <- c(
  "dimension", "iteration", "fit_seed",
  "finalized_hmm_stays", "finalized_stay_ids_md5",
  "n_train_stays", "n_validation_stays", "n_stable_reference",
  "logLik", "AIC", "BIC", "ICL_approx", "posterior_entropy",
  "mse1", "mse2", "mse3", "preprocessing_cache"
)

missing_columns <- setdiff(required_columns, names(raw_all))
if (length(missing_columns)) {
  stop(
    "Metrics files are missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

raw <- raw_all %>%
  filter(
    dimension %in% expected_k,
    iteration >= 1L,
    iteration <= n_iterations
  ) %>%
  arrange(dimension, iteration)

if (!nrow(raw)) {
  stop("No metrics rows remained after K/iteration filtering.")
}

duplicates <- raw %>%
  count(dimension, iteration, name = "n_records") %>%
  filter(n_records != 1L)

if (nrow(duplicates)) {
  write_csv(
    duplicates,
    file.path(out_dir, "duplicate_dimension_iteration_records.csv")
  )
  stop(
    "Duplicate K-iteration metrics records were found. ",
    "See duplicate_dimension_iteration_records.csv."
  )
}

expected_grid <- tidyr::crossing(
  dimension = expected_k,
  iteration = seq_len(n_iterations)
)

missing_jobs <- anti_join(
  expected_grid,
  raw %>% distinct(dimension, iteration),
  by = c("dimension", "iteration")
)

write_csv(
  missing_jobs,
  file.path(out_dir, "missing_dimension_iteration_jobs.csv")
)

identity_qc <- raw %>%
  summarise(
    n_unique_finalized_stays = n_distinct(finalized_hmm_stays),
    n_unique_stay_md5 = n_distinct(finalized_stay_ids_md5),
    n_unique_train_stays = n_distinct(n_train_stays),
    n_unique_validation_stays = n_distinct(n_validation_stays)
  )

write_csv(
  identity_qc,
  file.path(out_dir, "model_selection_identity_QC.csv")
)

if (
  identity_qc$n_unique_finalized_stays != 1L ||
  identity_qc$n_unique_stay_md5 != 1L
) {
  stop(
    "Cohort/model identity differs across metrics files. ",
    "See model_selection_identity_QC.csv."
  )
}

metric_columns <- c(
  "mse1", "mse2", "mse3",
  "logLik", "AIC", "BIC", "ICL_approx", "posterior_entropy"
)

nonfinite_qc <- raw %>%
  transmute(
    dimension,
    iteration,
    source_file,
    across(
      all_of(metric_columns),
      ~ !is.finite(.x),
      .names = "nonfinite_{.col}"
    )
  ) %>%
  filter(if_any(starts_with("nonfinite_"), identity))

write_csv(
  nonfinite_qc,
  file.path(out_dir, "nonfinite_model_selection_metrics.csv")
)

if (nrow(nonfinite_qc)) {
  warning(
    nrow(nonfinite_qc),
    " result rows contain at least one non-finite metric. ",
    "They remain in the raw export but are omitted metric-by-metric from summaries."
  )
}

summarise_metric <- function(x) {
  n_ok <- sum(is.finite(x))
  mean_x <- if (n_ok) mean(x[is.finite(x)]) else NA_real_
  median_x <- if (n_ok) median(x[is.finite(x)]) else NA_real_
  sd_x <- if (n_ok >= 2L) sd(x[is.finite(x)]) else NA_real_
  se_x <- if (n_ok >= 2L) sd_x / sqrt(n_ok) else NA_real_

  tibble(
    n = n_ok,
    mean = mean_x,
    median = median_x,
    sd = sd_x,
    se = se_x,
    min = if (n_ok) min(x[is.finite(x)]) else NA_real_,
    max = if (n_ok) max(x[is.finite(x)]) else NA_real_
  )
}

summary_long <- raw %>%
  select(dimension, iteration, all_of(metric_columns)) %>%
  pivot_longer(
    cols = all_of(metric_columns),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(dimension, metric) %>%
  group_modify(~ summarise_metric(.x$value)) %>%
  ungroup()

summary_wide <- summary_long %>%
  pivot_wider(
    names_from = metric,
    values_from = c(n, mean, median, sd, se, min, max),
    names_glue = "{metric}_{.value}"
  ) %>%
  arrange(dimension)

mse_long <- summary_long %>%
  filter(metric %in% c("mse1", "mse2", "mse3")) %>%
  mutate(
    horizon = as.integer(str_remove(metric, "^mse")),
    metric = factor(metric, levels = c("mse1", "mse2", "mse3"))
  ) %>%
  arrange(horizon, dimension)

one_se_for_horizon <- function(h) {
  z <- mse_long %>%
    filter(horizon == h, is.finite(mean), is.finite(se)) %>%
    arrange(dimension)

  if (!nrow(z)) {
    return(tibble())
  }

  best_index <- which.min(z$mean)
  best_k <- z$dimension[[best_index]]
  best_mean <- z$mean[[best_index]]
  best_se <- z$se[[best_index]]
  threshold <- best_mean + best_se

  eligible <- z %>%
    filter(mean <= threshold) %>%
    arrange(dimension)

  selected_k <- eligible$dimension[[1L]]

  tibble(
    horizon = h,
    minimum_MSE_dimension = best_k,
    minimum_MSE_mean = best_mean,
    SE_at_minimum = best_se,
    one_SE_threshold = threshold,
    one_SE_dimension = selected_k,
    n_eligible_dimensions = nrow(eligible),
    eligible_dimensions = paste(eligible$dimension, collapse = ";")
  )
}

one_se_candidates <- bind_rows(lapply(1:3, one_se_for_horizon))

if (!nrow(one_se_candidates)) {
  stop("No valid MSE horizon could be summarized.")
}

primary_recommendation <- one_se_candidates %>%
  filter(horizon == target_h)

if (!nrow(primary_recommendation)) {
  stop("The target horizon has no valid 1-SE result.")
}

eligible_by_horizon <- lapply(1:3, function(h) {
  x <- one_se_candidates %>% filter(horizon == h)
  if (!nrow(x)) return(integer())
  as.integer(strsplit(x$eligible_dimensions[[1L]], ";", fixed = TRUE)[[1L]])
})

consensus_eligible <- Reduce(intersect, eligible_by_horizon)

cross_horizon_diagnostic <- tibble(
  target_horizon = target_h,
  primary_one_SE_dimension = primary_recommendation$one_SE_dimension,
  primary_minimum_MSE_dimension = primary_recommendation$minimum_MSE_dimension,
  all_horizon_one_SE_dimensions = paste(
    paste0(
      "h", one_se_candidates$horizon,
      "=", one_se_candidates$one_SE_dimension
    ),
    collapse = ";"
  ),
  smallest_dimension_within_1SE_for_all_horizons =
    if (length(consensus_eligible)) min(consensus_eligible) else NA_integer_,
  all_horizon_eligible_dimensions =
    if (length(consensus_eligible)) paste(consensus_eligible, collapse = ";") else ""
)

support_at_primary <- summary_wide %>%
  filter(dimension == primary_recommendation$one_SE_dimension) %>%
  transmute(
    dimension,
    BIC_mean,
    ICL_approx_mean,
    posterior_entropy_mean,
    mse1_mean,
    mse2_mean,
    mse3_mean
  )

primary_recommendation <- primary_recommendation %>%
  left_join(support_at_primary, by = c("one_SE_dimension" = "dimension")) %>%
  mutate(
    selection_role =
      "Primary recommendation: smallest K within 1-SE of mean held-out MSE at target horizon"
  )

write_csv(
  raw,
  file.path(out_dir, "model_selection_all_metrics.csv")
)
write_csv(
  summary_long,
  file.path(out_dir, "model_selection_summary_long.csv")
)
write_csv(
  summary_wide,
  file.path(out_dir, "model_selection_summary_wide.csv")
)
write_csv(
  one_se_candidates,
  file.path(out_dir, "model_selection_oneSE_by_horizon.csv")
)
write_csv(
  primary_recommendation,
  file.path(out_dir, paste0(
    "model_selection_primary_h",
    target_h,
    "_recommendation.csv"
  ))
)
write_csv(
  cross_horizon_diagnostic,
  file.path(out_dir, "model_selection_cross_horizon_diagnostic.csv")
)

p_all_mse <- ggplot(
  mse_long,
  aes(x = dimension, y = mean)
) +
  geom_ribbon(
    aes(ymin = mean - se, ymax = mean + se),
    alpha = 0.15
  ) +
  geom_line() +
  geom_point(size = 1.4) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 1,
    labeller = as_labeller(c(
      mse1 = "1-step held-out MSE",
      mse2 = "2-step held-out MSE",
      mse3 = "3-step held-out MSE"
    ))
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Number of HMM states (K)",
    y = "Mean held-out MSE across patient-level splits",
    title = "Held-out prediction error by HMM state number",
    subtitle = "Ribbon: mean +/- one SE across patient-level splits"
  )

ggsave(
  file.path(out_dir, "model_selection_MSE_all_horizons.pdf"),
  p_all_mse,
  width = 9,
  height = 11
)

ggsave(
  file.path(out_dir, "model_selection_MSE_all_horizons.png"),
  p_all_mse,
  width = 9,
  height = 11,
  dpi = 300,
  bg = "white"
)

primary_curve <- mse_long %>%
  filter(horizon == target_h)

p_primary <- ggplot(
  primary_curve,
  aes(x = dimension, y = mean)
) +
  geom_ribbon(
    aes(ymin = mean - se, ymax = mean + se),
    alpha = 0.15
  ) +
  geom_line() +
  geom_point(size = 1.7) +
  geom_hline(
    yintercept = primary_recommendation$one_SE_threshold,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = primary_recommendation$one_SE_dimension,
    linetype = "dashed"
  ) +
  geom_point(
    data = primary_curve %>%
      filter(dimension == primary_recommendation$one_SE_dimension),
    size = 3.2,
    shape = 21,
    fill = "white",
    stroke = 1
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Number of HMM states (K)",
    y = paste0("Mean ", target_h, "-step held-out MSE"),
    title = paste0(
      "Primary model selection at horizon h=", target_h
    ),
    subtitle = paste0(
      "Horizontal: 1-SE threshold; vertical: smallest eligible K = ",
      primary_recommendation$one_SE_dimension
    )
  )

ggsave(
  file.path(
    out_dir,
    paste0("model_selection_primary_h", target_h, "_1SE.pdf")
  ),
  p_primary,
  width = 10,
  height = 6.5
)

ggsave(
  file.path(
    out_dir,
    paste0("model_selection_primary_h", target_h, "_1SE.png")
  ),
  p_primary,
  width = 10,
  height = 6.5,
  dpi = 300,
  bg = "white"
)

criteria_long <- summary_long %>%
  filter(metric %in% c("BIC", "ICL_approx"), is.finite(mean)) %>%
  group_by(metric) %>%
  mutate(delta_from_minimum = mean - min(mean)) %>%
  ungroup()

p_criteria <- ggplot(
  criteria_long,
  aes(x = dimension, y = delta_from_minimum)
) +
  geom_line() +
  geom_point(size = 1.5) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 1,
    labeller = as_labeller(c(
      BIC = "BIC",
      ICL_approx = "Approximate ICL"
    ))
  ) +
  geom_vline(
    xintercept = primary_recommendation$one_SE_dimension,
    linetype = "dashed"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Number of HMM states (K)",
    y = "Mean criterion minus its minimum",
    title = "Supporting information criteria",
    subtitle = "Dashed line: target-horizon 1-SE recommendation"
  )

ggsave(
  file.path(out_dir, "model_selection_BIC_ICL_support.pdf"),
  p_criteria,
  width = 9,
  height = 8
)

pareto_data <- summary_wide %>%
  transmute(
    dimension,
    target_MSE_mean = .data[[paste0("mse", target_h, "_mean")]],
    ICL_mean = ICL_approx_mean
  ) %>%
  filter(is.finite(target_MSE_mean), is.finite(ICL_mean))

p_pareto <- ggplot(
  pareto_data,
  aes(x = ICL_mean, y = target_MSE_mean)
) +
  geom_path() +
  geom_point(size = 2, alpha = 0.75) +
  geom_point(
    data = pareto_data %>%
      filter(dimension == primary_recommendation$one_SE_dimension),
    size = 4,
    shape = 21,
    fill = "white",
    stroke = 1.2
  ) +
  theme_minimal(base_size = 12) +
  labs(
    x = "Mean approximate ICL",
    y = paste0("Mean ", target_h, "-step held-out MSE"),
    title = "Held-out MSE versus approximate ICL",
    subtitle = "Outlined point: target-horizon 1-SE recommendation"
  )

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p_pareto <- p_pareto +
    ggrepel::geom_text_repel(
      aes(label = dimension),
      size = 3,
      max.overlaps = 30
    )
}

ggsave(
  file.path(out_dir, "model_selection_MSE_ICL_pareto.pdf"),
  p_pareto,
  width = 8,
  height = 6.5
)

manifest <- tibble(
  item = c(
    "created_at",
    "data_directory",
    "iterations_expected",
    "target_horizon",
    "expected_K_values",
    "metrics_files_read",
    "metrics_rows_used",
    "missing_jobs",
    "primary_rule",
    "primary_one_SE_dimension",
    "primary_minimum_MSE_dimension"
  ),
  value = c(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    normalizePath(data_dir, mustWork = FALSE),
    as.character(n_iterations),
    as.character(target_h),
    paste(expected_k, collapse = ";"),
    as.character(length(files)),
    as.character(nrow(raw)),
    as.character(nrow(missing_jobs)),
    "Mean target-horizon held-out MSE across patient-level splits; smallest K within 1-SE",
    as.character(primary_recommendation$one_SE_dimension),
    as.character(primary_recommendation$minimum_MSE_dimension)
  )
)

write_csv(
  manifest,
  file.path(out_dir, "model_selection_manifest.csv")
)

cat("\nModel-selection aggregation completed.\n")
cat("Metrics rows used: ", nrow(raw), "\n", sep = "")
cat("Missing K-iteration jobs: ", nrow(missing_jobs), "\n", sep = "")
cat("Primary horizon: h=", target_h, "\n", sep = "")
cat(
  "Minimum-mean-MSE K: ",
  primary_recommendation$minimum_MSE_dimension,
  "\n",
  sep = ""
)
cat(
  "Primary 1-SE recommendation: K=",
  primary_recommendation$one_SE_dimension,
  "\n",
  sep = ""
)
cat("Cross-horizon 1-SE candidates:\n")
print(one_se_candidates, row.names = FALSE)
cat(
  "\nSupporting h=1/h=2 curves and BIC/ICL summaries were also written.\n"
)
