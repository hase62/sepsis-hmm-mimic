
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
if (!exists("MIMIC_ACTIVE_ANNOT_MASTER", inherits = TRUE)) {
  source(code_path("17_analysis_config.R"))
}

if (!file.exists(code_path("17_analysis_config.R"))) stop("Missing shared analysis configuration: 17_analysis_config.R", call. = FALSE)
source(code_path("17_analysis_config.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

OUT_PREFIX <- Sys.getenv(
  "MIMIC_ORGANSUPPORT_QC_OUT_PREFIX",
  "organsupport_descriptive"
)

SELECTED_DAYS <- as.integer(strsplit(Sys.getenv("MIMIC_ORGANSUPPORT_SELECTED_DAYS", "1,3,5,7,10,14,21,28"), ",")[[1]])
PRIMARY_DAYS_FOR_PLOTS <- as.integer(strsplit(Sys.getenv("MIMIC_ORGANSUPPORT_PRIMARY_DAYS", "1,3,7,14"), ",")[[1]])

ROWLEVEL_CANDIDATES <- c(
  Sys.getenv("MIMIC_ORGANSUPPORT_ROWLEVEL", ""),
  "incremental_value_current_sofa_rowlevel_selected_days_with_currentday_SOFA_organsupport.csv"
)
ROWLEVEL_CANDIDATES <- ROWLEVEL_CANDIDATES[nzchar(ROWLEVEL_CANDIDATES)]

SUBPHENO_LEVELS_FULL <- c(
  "A1 low burden",
  "A2 mild burden",
  "B1 renal mild",
  "B2 renal severe",
  "G respiratory",
  "D1 lactate/shock",
  "D2 hepatic",
  "D3 shock/coagulation-overlap"
)
SUBPHENO_LEVELS_SHORT <- c(
  "A1 low burden",
  "A2 mild burden",
  "B1 renal mild",
  "B2 renal severe",
  "G respiratory",
  "D1 lactate/shock",
  "D2 hepatic",
  "D3 shock/coag overlap"
)

.stopif <- function(cond, msg) {
  if (cond) stop(msg, call. = FALSE)
}

first_existing <- function(paths) {
  paths <- paths[nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

as_bool <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x != 0)
  y <- tolower(trimws(as.character(x)))
  y %in% c("true", "t", "1", "yes", "y")
}

num <- function(x) suppressWarnings(as.numeric(x))

q25 <- function(x) {
  x <- num(x)
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE))
}
q75 <- function(x) {
  x <- num(x)
  if (all(is.na(x))) return(NA_real_)
  as.numeric(stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE))
}
med <- function(x) {
  x <- num(x)
  if (all(is.na(x))) return(NA_real_)
  stats::median(x, na.rm = TRUE)
}
mean_pct_bool <- function(x) 100 * mean(as_bool(x), na.rm = TRUE)

write_both <- function(df, stem) {
  readr::write_csv(df, paste0(stem, ".csv"), na = "")
  if (requireNamespace("arrow", quietly = TRUE)) {
    try(arrow::write_parquet(df, paste0(stem, ".parquet")), silent = TRUE)
  }
  invisible(df)
}

std_z <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - m) / s
}

harmonize_subphenotype <- function(x) {
  y <- as.character(x)
  map <- c(
    "A1 low burden" = "A1 low burden",
    "A2 mild burden" = "A2 mild burden",
    "B1 renal mild" = "B1 renal mild",
    "B2 renal severe" = "B2 renal severe",
    "G respiratory" = "G respiratory",
    "D1 lactate/shock" = "D1 lactate/shock",
    "D2 hepatic" = "D2 hepatic",
    "D3 shock/coagulation-overlap" = "D3 shock/coag overlap",
    "D3 shock/coag overlap" = "D3 shock/coag overlap"
  )
  out <- unname(map[y])
  out[is.na(out)] <- y[is.na(out)]
  factor(out, levels = SUBPHENO_LEVELS_SHORT)
}

row_fp <- first_existing(ROWLEVEL_CANDIDATES)
.stopif(is.na(row_fp), paste("Could not find row-level current-day SOFA/organsupport input. Tried:", paste(ROWLEVEL_CANDIDATES, collapse = "; ")))
message("[organ-support summary] Reading row-level input: ", row_fp)
row_df <- readr::read_csv(row_fp, show_col_types = FALSE)

required_cols <- c("stay_id", "day", "current_subphenotype")
missing_required <- setdiff(required_cols, names(row_df))
.stopif(length(missing_required) > 0, paste("Missing required columns:", paste(missing_required, collapse = ", ")))

row_df <- row_df %>%
  mutate(
    day = as.integer(day),
    current_subphenotype_short = harmonize_subphenotype(current_subphenotype),
    current_subphenotype = as.character(current_subphenotype_short)
  ) %>%
  filter(day %in% SELECTED_DAYS, !is.na(current_subphenotype_short))

support_candidate_cols <- c(
  "ventilation_active",
  "peep_max", "peep_median",
  "fio2_max", "fio2_median",
  "vasopressor_any", "vasopressor_sofa_any",
  "rrt_any", "crrt_any", "hemodialysis_any",
  "ecmo_any", "ecmo_flow_max",
  "uo_24h",
  "sofa_total_available_component", "sofa_total_complete6",
  "sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns", "sofa_renal"
)
missing_support_cols <- setdiff(support_candidate_cols, names(row_df))

for (cc in missing_support_cols) row_df[[cc]] <- NA

qc_overall <- tibble::tibble(
  metric = c(
    "rowlevel_input_file",
    "rows_selected_days",
    "unique_stays",
    "selected_days",
    "missing_optional_support_columns",
    "ventilation_active_rows_pct",
    "vasopressor_any_rows_pct",
    "rrt_any_rows_pct",
    "ecmo_any_rows_pct"
  ),
  value = c(
    row_fp,
    as.character(nrow(row_df)),
    as.character(dplyr::n_distinct(row_df$stay_id)),
    paste(sort(unique(row_df$day)), collapse = ","),
    ifelse(length(missing_support_cols) == 0, "none", paste(missing_support_cols, collapse = ",")),
    sprintf("%.3f", mean_pct_bool(row_df$ventilation_active)),
    sprintf("%.3f", mean_pct_bool(row_df$vasopressor_any)),
    sprintf("%.3f", mean_pct_bool(row_df$rrt_any)),
    sprintf("%.3f", mean_pct_bool(row_df$ecmo_any))
  )
)
write_both(qc_overall, paste0(OUT_PREFIX, "_QC_overall"))

qc_by_day <- row_df %>%
  group_by(day) %>%
  summarise(
    n = n(),
    unique_stays = n_distinct(stay_id),
    ventilation_active_pct = mean_pct_bool(ventilation_active),
    vasopressor_any_pct = mean_pct_bool(vasopressor_any),
    rrt_any_pct = mean_pct_bool(rrt_any),
    ecmo_any_pct = mean_pct_bool(ecmo_any),
    peep_max_nonmissing_pct = 100 * mean(is.finite(num(peep_max)), na.rm = TRUE),
    fio2_max_nonmissing_pct = 100 * mean(is.finite(num(fio2_max)), na.rm = TRUE),
    uo_24h_nonmissing_pct = 100 * mean(is.finite(num(uo_24h)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(day)
write_both(qc_by_day, paste0(OUT_PREFIX, "_QC_by_day"))

summary_by_day_subphenotype <- row_df %>%
  group_by(day, current_subphenotype_short) %>%
  summarise(
    n = n(),
    unique_stays = n_distinct(stay_id),
    ventilation_active_pct = mean_pct_bool(ventilation_active),
    peep_max_median = med(peep_max),
    peep_max_q25 = q25(peep_max),
    peep_max_q75 = q75(peep_max),
    fio2_max_median = med(fio2_max),
    fio2_max_q25 = q25(fio2_max),
    fio2_max_q75 = q75(fio2_max),
    vasopressor_any_pct = mean_pct_bool(vasopressor_any),
    vasopressor_sofa_any_pct = mean_pct_bool(vasopressor_sofa_any),
    rrt_any_pct = mean_pct_bool(rrt_any),
    crrt_any_pct = mean_pct_bool(crrt_any),
    hemodialysis_any_pct = mean_pct_bool(hemodialysis_any),
    ecmo_any_pct = mean_pct_bool(ecmo_any),
    ecmo_flow_max_median = med(ecmo_flow_max),
    urine_output_24h_median = med(uo_24h),
    urine_output_24h_q25 = q25(uo_24h),
    urine_output_24h_q75 = q75(uo_24h),
    sofa_total_available_component_median = med(sofa_total_available_component),
    sofa_total_available_component_q25 = q25(sofa_total_available_component),
    sofa_total_available_component_q75 = q75(sofa_total_available_component),
    sofa_total_complete6_median = med(sofa_total_complete6),
    sofa_total_complete6_q25 = q25(sofa_total_complete6),
    sofa_total_complete6_q75 = q75(sofa_total_complete6),
    .groups = "drop"
  ) %>%
  mutate(current_subphenotype = as.character(current_subphenotype_short)) %>%
  select(day, current_subphenotype, everything(), -current_subphenotype_short) %>%
  arrange(day, factor(current_subphenotype, levels = SUBPHENO_LEVELS_SHORT))
write_both(summary_by_day_subphenotype, paste0(OUT_PREFIX, "_by_day_subphenotype"))

summary_overall_subphenotype <- row_df %>%
  group_by(current_subphenotype_short) %>%
  summarise(
    n = n(),
    unique_stays = n_distinct(stay_id),
    ventilation_active_pct = mean_pct_bool(ventilation_active),
    peep_max_median = med(peep_max),
    peep_max_q25 = q25(peep_max),
    peep_max_q75 = q75(peep_max),
    fio2_max_median = med(fio2_max),
    fio2_max_q25 = q25(fio2_max),
    fio2_max_q75 = q75(fio2_max),
    vasopressor_any_pct = mean_pct_bool(vasopressor_any),
    vasopressor_sofa_any_pct = mean_pct_bool(vasopressor_sofa_any),
    rrt_any_pct = mean_pct_bool(rrt_any),
    crrt_any_pct = mean_pct_bool(crrt_any),
    hemodialysis_any_pct = mean_pct_bool(hemodialysis_any),
    ecmo_any_pct = mean_pct_bool(ecmo_any),
    urine_output_24h_median = med(uo_24h),
    urine_output_24h_q25 = q25(uo_24h),
    urine_output_24h_q75 = q75(uo_24h),
    sofa_total_available_component_median = med(sofa_total_available_component),
    sofa_total_available_component_q25 = q25(sofa_total_available_component),
    sofa_total_available_component_q75 = q75(sofa_total_available_component),
    sofa_total_complete6_median = med(sofa_total_complete6),
    sofa_total_complete6_q25 = q25(sofa_total_complete6),
    sofa_total_complete6_q75 = q75(sofa_total_complete6),
    .groups = "drop"
  ) %>%
  mutate(current_subphenotype = as.character(current_subphenotype_short)) %>%
  select(current_subphenotype, everything(), -current_subphenotype_short) %>%
  arrange(factor(current_subphenotype, levels = SUBPHENO_LEVELS_SHORT))
write_both(summary_overall_subphenotype, paste0(OUT_PREFIX, "_overall_subphenotype"))

reporting_summary <- summary_by_day_subphenotype %>%
  filter(day %in% PRIMARY_DAYS_FOR_PLOTS) %>%
  mutate(
    support_interpretation = "descriptive organ-support and current-day SOFA coherence summary"
  ) %>%
  arrange(day, factor(current_subphenotype, levels = SUBPHENO_LEVELS_SHORT))
write_both(reporting_summary, paste0(OUT_PREFIX, "_reporting_summary_selected_days"))

heatmap_wide <- summary_overall_subphenotype %>%
  transmute(
    current_subphenotype,
    `Ventilation active (%)` = ventilation_active_pct,
    `PEEP max median` = peep_max_median,
    `FiO2 max median` = fio2_max_median,
    `Vasopressor any (%)` = vasopressor_any_pct,
    `RRT any (%)` = rrt_any_pct,
    `ECMO any (%)` = ecmo_any_pct,
    `Lower urine output` = -urine_output_24h_median,
    `Available-component SOFA median` = sofa_total_available_component_median
  )

heatmap_long <- heatmap_wide %>%
  pivot_longer(-current_subphenotype, names_to = "metric", values_to = "value") %>%
  group_by(metric) %>%
  mutate(value_z = std_z(value)) %>%
  ungroup() %>%
  mutate(
    current_subphenotype = factor(current_subphenotype, levels = rev(SUBPHENO_LEVELS_SHORT)),
    metric = factor(metric, levels = c(
      "Ventilation active (%)",
      "PEEP max median",
      "FiO2 max median",
      "Vasopressor any (%)",
      "RRT any (%)",
      "ECMO any (%)",
      "Lower urine output",
      "Available-component SOFA median"
    ))
  )
write_both(heatmap_long, paste0(OUT_PREFIX, "_overall_support_burden_heatmap_values"))

selected_day_long <- summary_by_day_subphenotype %>%
  filter(day %in% PRIMARY_DAYS_FOR_PLOTS) %>%
  transmute(
    day,
    current_subphenotype,
    `Vasopressor any (%)` = vasopressor_any_pct,
    `RRT any (%)` = rrt_any_pct,
    `PEEP max median` = peep_max_median,
    `FiO2 max median` = fio2_max_median,
    `Ventilation active (%)` = ventilation_active_pct
  ) %>%
  pivot_longer(cols = -c(day, current_subphenotype), names_to = "metric", values_to = "value") %>%
  mutate(
    current_subphenotype = factor(current_subphenotype, levels = SUBPHENO_LEVELS_SHORT),
    metric = factor(metric, levels = c("Ventilation active (%)", "PEEP max median", "FiO2 max median", "Vasopressor any (%)", "RRT any (%)"))
  )
write_both(selected_day_long, paste0(OUT_PREFIX, "_selected_day_support_plot_values"))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))

  p_heat <- ggplot(heatmap_long, aes(x = metric, y = current_subphenotype, fill = value_z)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is.finite(value), sprintf("%.1f", value), "")), size = 2.4) +
    scale_fill_gradient2(name = "Within-metric\nstandardized\nburden", low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0, na.value = "grey90") +
    labs(
      title = "Organ-support and severity burden by current HMM subphenotype",
      subtitle = "Selected ICU days pooled; values are raw summaries, color is standardized within metric",
      x = NULL,
      y = "Current HMM subphenotype"
    ) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right"
    )
  ggsave(paste0(OUT_PREFIX, "_overall_support_burden_heatmap.pdf"), p_heat, width = 9.5, height = 4.8, bg = "white")
  ggsave(paste0(OUT_PREFIX, "_overall_support_burden_heatmap.png"), p_heat, width = 9.5, height = 4.8, dpi = 300, bg = "white")

  p_day <- ggplot(selected_day_long, aes(x = current_subphenotype, y = value, group = day)) +
    geom_line(aes(linetype = factor(day)), linewidth = 0.35) +
    geom_point(aes(shape = factor(day)), size = 1.3) +
    facet_wrap(~ metric, scales = "free_y", ncol = 2) +
    labs(
      title = "Selected-day organ-support descriptors by HMM subphenotype",
      subtitle = "Descriptive organ-support and current-day SOFA burden across HMM subphenotypes",
      x = NULL,
      y = "Value",
      linetype = "ICU day",
      shape = "ICU day"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  ggsave(paste0(OUT_PREFIX, "_selected_day_support_profiles.pdf"), p_day, width = 11, height = 7.2, bg = "white")
  ggsave(paste0(OUT_PREFIX, "_selected_day_support_profiles.png"), p_day, width = 11, height = 7.2, dpi = 300, bg = "white")
}

manifest <- tibble::tibble(
  item = c(
    "created_at",
    "work_dir",
    "rowlevel_input_file",
    "out_prefix",
    "selected_days",
    "primary_days_for_plots",
    "interpretation"
  ),
  value = c(
    as.character(Sys.time()),
    getwd(),
    row_fp,
    OUT_PREFIX,
    paste(SELECTED_DAYS, collapse = ","),
    paste(PRIMARY_DAYS_FOR_PLOTS, collapse = ","),
    "descriptive organ-support and current-day SOFA coherence summary"
  )
)
write_both(manifest, paste0(OUT_PREFIX, "_manifest"))

message("[organ-support summary] Done.")
print(qc_overall)
