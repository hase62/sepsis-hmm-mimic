
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
if (!exists("MIMIC_ACTIVE_ANNOT_MASTER", inherits = TRUE)) {
  source(code_path("17_analysis_config.R"))
}

if (!file.exists(code_path("17_analysis_config.R"))) stop("Missing shared analysis configuration: 17_analysis_config.R", call. = FALSE)
source(code_path("17_analysis_config.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(purrr)
})

select     <- dplyr::select
filter     <- dplyr::filter
rename     <- dplyr::rename
mutate     <- dplyr::mutate
summarise  <- dplyr::summarise
summarize  <- dplyr::summarise
arrange    <- dplyr::arrange
group_by   <- dplyr::group_by
ungroup    <- dplyr::ungroup
left_join  <- dplyr::left_join
inner_join <- dplyr::inner_join
bind_rows  <- dplyr::bind_rows
all_of     <- tidyselect::all_of
any_of     <- tidyselect::any_of


OUT_PREFIX <- "next_day_organ_burden"
ANNO_CANDIDATES <- c(
  MIMIC_ACTIVE_ANNOT_MASTER,
  MIMIC_ACTIVE_ANNOT_MASTER
)

MIN_NONMISSING_PER_CELL <- 20L

SUBPH_LEVELS <- c(
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

VAR_SPECS <- tibble::tribble(
  ~variable_id,       ~display_variable,  ~domain,          ~burden_direction, ~candidate_names,
  "PF",              "PaO2/FiO2",        "Respiratory",    -1,               list(c("PF", "PaO2_FiO2", "PFratio", "PF_ratio")),
  "FiO2",            "FiO2",             "Respiratory",     1,               list(c("FiO2", "FIO2", "fio2")),
  "MAP",             "MAP",              "Hemodynamic",    -1,               list(c("MAP", "MEAN", "MeanBP", "mean_bp")),
  "Lactate",         "Lactate",          "Acid-base",       1,               list(c("Lactate", "lactate", "LACTATE")),
  "pH",              "pH",               "Acid-base",      -1,               list(c("pH", "PH", "ph")),
  "HCO3",            "HCO3",             "Acid-base",      -1,               list(c("HCO3", "Bicarbonate", "bicarbonate")),
  "Creatinine",      "Creatinine",       "Renal",           1,               list(c("Creatinine", "creatinine", "CREATININE")),
  "UreaNitrogen",    "BUN",              "Renal",           1,               list(c("UreaNitrogen", "BUN", "BloodUreaNitrogen", "bun")),
  "TotalBilirubin",  "Total bilirubin",  "Hepatic",         1,               list(c("TotalBilirubin", "Bilirubin", "bilirubin", "total_bilirubin")),
  "INR",             "INR",              "Coagulation",     1,               list(c("INR", "INR.PT.", "INR_PT", "INR(PT)")),
  "PlateletCount",   "Platelets",        "Coagulation",    -1,               list(c("PlateletCount", "Platelet", "Platelets", "platelet", "platelet_count"))
)

DOMAIN_LEVELS <- c("Respiratory", "Hemodynamic", "Acid-base", "Renal", "Hepatic", "Coagulation")

TARGET_MAP <- tibble::tribble(
  ~current_subphenotype,             ~variable_id,      ~expected_pattern,
  "D1 lactate/shock",           "Lactate",        "higher burden",
  "D1 lactate/shock",           "pH",             "lower reserve",
  "D1 lactate/shock",           "HCO3",           "lower reserve",
  "D1 lactate/shock",           "MAP",            "lower reserve",
  "G respiratory",             "PF",             "lower reserve",
  "G respiratory",             "FiO2",           "higher burden",
  "B2 renal severe",            "Creatinine",     "higher burden",
  "B2 renal severe",            "UreaNitrogen",   "higher burden",
  "D2 hepatic",                "TotalBilirubin", "higher burden",
  "D2 hepatic",                "INR",            "higher burden",
  "D2 hepatic",                "PlateletCount",  "lower reserve",
  "D3 shock/coagulation-overlap",             "MAP",            "lower reserve",
  "D3 shock/coagulation-overlap",             "HCO3",           "lower reserve",
  "D3 shock/coagulation-overlap",             "pH",             "lower reserve",
  "D3 shock/coagulation-overlap",             "PlateletCount",  "lower reserve",
  "D3 shock/coagulation-overlap",             "INR",            "higher burden",
  "D3 shock/coagulation-overlap",             "PF",             "lower reserve",
  "A1 low burden",       "Lactate",        "lower burden",
  "A1 low burden",       "PF",             "higher reserve",
  "A1 low burden",       "Creatinine",     "lower burden",
  "A2 mild burden",        "Lactate",        "mild burden",
  "A2 mild burden",        "PF",             "mild burden",
  "A2 mild burden",        "Creatinine",     "mild burden"
)

.stopif <- function(ok, msg) if (!ok) stop(msg, call. = FALSE)
num <- function(x) suppressWarnings(as.numeric(x))

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)][1]
  if (length(hit) == 0 || is.na(hit) || !nzchar(hit)) return(NA_character_)
  hit
}

state_id_from_name <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove(x, "^S")
  suppressWarnings(as.integer(x))
}

make_safe_level <- function(x) {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "Other/NA"
  x
}

preferred_subphenotype_levels <- function(x = character(0)) {
  x <- x[!is.na(x) & nzchar(x)]
  c(intersect(SUBPH_LEVELS, unique(x)), setdiff(unique(x), SUBPH_LEVELS))
}

normalize_annotation_table <- function(tbl) {
  nm <- names(tbl)

  state_col <- NA_character_
  for (cand in c("state_id", "state", "State", "STATE", "state_name", "state_label")) {
    if (cand %in% nm) {
      state_col <- cand
      break
    }
  }
  .stopif(!is.na(state_col),
          paste("Annotation file must contain a state identifier column such as state_id or state. Available columns:", paste(nm, collapse = ", ")))

  subph_col <- NA_character_
  for (cand in c("subphenotype", "Subphenotype", "assigned_subphenotype", "hmm_subphenotype", "subtype_label")) {
    if (cand %in% nm) {
      subph_col <- cand
      break
    }
  }
  .stopif(!is.na(subph_col),
          paste("Annotation file must contain a subphenotype column. Available columns:", paste(nm, collapse = ", ")))

  out <- tibble::tibble(
    state_id = state_id_from_name(tbl[[state_col]]),
    subphenotype = make_safe_level(tbl[[subph_col]])
  ) %>%
    dplyr::filter(!is.na(state_id)) %>%
    dplyr::distinct(state_id, .keep_all = TRUE)

  .stopif(nrow(out) > 0, "No usable state_id/subphenotype mapping could be extracted from annotation file.")
  out
}

pick_col_name <- function(df, candidates) {
  candidates <- unlist(candidates, use.names = FALSE)
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

safe_iqr <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  out <- as.numeric(stats::IQR(x, na.rm = TRUE, type = 7))
  if (!is.finite(out) || out <= 0) return(NA_real_)
  out
}

write_both <- function(df, stem) {
  readr::write_csv(df, paste0(stem, ".csv"))
  readr::write_tsv(df, paste0(stem, ".tsv"))
  invisible(df)
}

make_png_if_possible <- function(plot, filename, width, height, dpi = 300) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))
  ggplot2::ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
  invisible(TRUE)
}

needed <- c("post_smooth_full", "combined", "ids")
missing_needed <- needed[!vapply(needed, exists, logical(1), envir = .GlobalEnv)]
.stopif(length(missing_needed) == 0,
        paste("Missing required objects:", paste(missing_needed, collapse = ", ")))

post_smooth_full <- as.matrix(get("post_smooth_full", envir = .GlobalEnv))
combined <- as.data.frame(get("combined", envir = .GlobalEnv))
ids <- as.character(get("ids", envir = .GlobalEnv))

.stopif(nrow(post_smooth_full) == length(ids), "nrow(post_smooth_full) must equal length(ids).")
.stopif(nrow(post_smooth_full) == nrow(combined), "nrow(post_smooth_full) must equal nrow(combined).")
.stopif("lab_time" %in% names(combined), "combined must contain lab_time.")

if (!"stay_id" %in% names(combined)) combined$stay_id <- ids
combined$stay_id <- as.character(combined$stay_id)
combined$lab_time <- suppressWarnings(as.integer(combined$lab_time))
message("Using combined columns: stay_id + lab_time; n rows = ", nrow(combined))

anno_path <- first_existing(ANNO_CANDIDATES)
if (is.na(anno_path)) {
  anno_hits <- list.files(pattern = "^state_annotation\\.csv$", full.names = FALSE)
  if (length(anno_hits) > 0) anno_path <- anno_hits[1]
}
.stopif(!is.na(anno_path), paste("Could not find annotation CSV. Tried:", paste(ANNO_CANDIDATES, collapse = ", ")))
state_master_raw <- readr::read_csv(anno_path, show_col_types = FALSE)
state_master <- normalize_annotation_table(state_master_raw)
message("Using annotation file: ", anno_path)
message("Using annotation columns resolved from file; standardized internally to: state_id + subphenotype")

state_names <- colnames(post_smooth_full)
if (is.null(state_names) || any(is.na(state_names)) || any(!nzchar(state_names))) {
  state_names <- paste0("S", seq_len(ncol(post_smooth_full)))
  colnames(post_smooth_full) <- state_names
}
state_ids <- state_id_from_name(state_names)
state_to_subph <- tibble::tibble(state_name = state_names, state_id = state_ids) %>%
  dplyr::left_join(state_master %>% dplyr::select(state_id, subphenotype), by = "state_id") %>%
  dplyr::mutate(subphenotype = make_safe_level(subphenotype))

if (any(state_to_subph$subphenotype == "Other/NA")) {
  warning("Some posterior states did not match annotation state_id values: ",
          paste(state_to_subph$state_name[state_to_subph$subphenotype == "Other/NA"], collapse = ", "))
}

subph_levels <- preferred_subphenotype_levels(state_to_subph$subphenotype)
.stopif(length(subph_levels) > 0, "No subphenotype levels found after matching annotation to posterior states.")

subph_prob <- matrix(0, nrow = nrow(post_smooth_full), ncol = length(subph_levels),
                     dimnames = list(NULL, subph_levels))
for (sp in subph_levels) {
  st <- state_to_subph %>% dplyr::filter(subphenotype == sp) %>% dplyr::pull(state_name)
  st <- intersect(st, colnames(post_smooth_full))
  if (length(st) > 0) subph_prob[, sp] <- rowSums(post_smooth_full[, st, drop = FALSE], na.rm = TRUE)
}

row_tot <- rowSums(subph_prob, na.rm = TRUE)
dom_idx <- max.col(subph_prob, ties.method = "first")
dominant_subphenotype <- colnames(subph_prob)[dom_idx]
dominant_subphenotype[!is.finite(row_tot) | row_tot <= 0] <- "Other/NA"

row_info <- tibble(
  row_id = seq_len(nrow(combined)),
  stay_id = combined$stay_id,
  day = combined$lab_time,
  current_subphenotype = dominant_subphenotype
) %>%
  filter(!is.na(day), !is.na(stay_id), nzchar(stay_id))

source_name <- "combined_original_scale_matrix"
phys_source <- combined

use_imputed <- identical(Sys.getenv("NEXTDAY_BURDEN_USE_IMPUTED_DATA"), "1")
if (use_imputed && exists("imputed_data", envir = .GlobalEnv)) {
  tmp <- as.data.frame(get("imputed_data", envir = .GlobalEnv))
  if (nrow(tmp) == nrow(combined)) {
    phys_source <- tmp
    source_name <- "imputed_data_HMM_observation_matrix_requested_by_env"
  } else {
    warning("NEXTDAY_BURDEN_USE_IMPUTED_DATA=1 but nrow(imputed_data) != nrow(combined); using combined instead.")
  }
}

if (!"PF" %in% names(phys_source)) {
  pao2_col <- pick_col_name(phys_source, c("PaO2", "pO2", "PO2", "Arterial_PO2"))
  fio2_col <- pick_col_name(phys_source, c("FiO2", "FIO2", "fio2"))
  if (!is.na(pao2_col) && !is.na(fio2_col)) {
    fio2 <- num(phys_source[[fio2_col]])
    fio2_frac <- ifelse(fio2 > 1.5, fio2 / 100, fio2)
    fio2_frac <- pmin(pmax(fio2_frac, 0.21), 1.0)
    phys_source$PF <- num(phys_source[[pao2_col]]) / fio2_frac
    message("Computed PF from ", pao2_col, " and ", fio2_col, ".")
  }
}

resolved_vars <- VAR_SPECS %>%
  rowwise() %>%
  mutate(source_column = pick_col_name(phys_source, candidate_names)) %>%
  ungroup() %>%
  mutate(
    available = !is.na(source_column),
    domain = factor(domain, levels = DOMAIN_LEVELS),
    variable_order = row_number()
  )

missing_vars <- resolved_vars %>% filter(!available)
if (nrow(missing_vars) > 0) {
  warning("Some requested physiological variables were not available and will be omitted: ",
          paste(missing_vars$variable_id, collapse = ", "))
}

available_specs <- resolved_vars %>% filter(available)
.stopif(nrow(available_specs) > 0, "No requested physiological variables were available.")
message("Physiology source: ", source_name)
message("Resolved next-day variables: ", paste(available_specs$variable_id, "=", available_specs$source_column, collapse = "; "))

phys_next <- tibble(
  stay_id = combined$stay_id,
  day_next = combined$lab_time,
  row_id_next = seq_len(nrow(combined))
)
for (i in seq_len(nrow(available_specs))) {
  phys_next[[available_specs$variable_id[i]]] <- num(phys_source[[available_specs$source_column[i]]])
}

pairs <- row_info %>%
  mutate(day_next = day + 1L) %>%
  inner_join(phys_next, by = c("stay_id", "day_next")) %>%
  mutate(
    current_subphenotype = factor(current_subphenotype,
                                  levels = c(intersect(SUBPH_LEVELS, unique(current_subphenotype)),
                                             setdiff(unique(current_subphenotype), SUBPH_LEVELS)))
  )

.stopif(nrow(pairs) > 0, "No adjacent active ICU-day pairs were found.")

var_ids <- available_specs$variable_id
pair_long <- pairs %>%
  select(row_id_current = row_id, row_id_next, stay_id, day, day_next, current_subphenotype, all_of(var_ids)) %>%
  pivot_longer(cols = all_of(var_ids), names_to = "variable_id", values_to = "next_day_value") %>%
  left_join(available_specs %>% select(variable_id, display_variable, domain, burden_direction, source_column, variable_order),
            by = "variable_id") %>%
  mutate(
    current_subphenotype = as.character(current_subphenotype),
    domain = as.character(domain)
  )

global_stats <- pair_long %>%
  group_by(variable_id, display_variable, domain, burden_direction, source_column, variable_order) %>%
  summarise(
    n_pairs_total = dplyr::n(),
    n_nonmissing_total = sum(!is.na(next_day_value)),
    global_median = median(next_day_value, na.rm = TRUE),
    global_q1 = as.numeric(stats::quantile(next_day_value, 0.25, na.rm = TRUE, names = FALSE)),
    global_q3 = as.numeric(stats::quantile(next_day_value, 0.75, na.rm = TRUE, names = FALSE)),
    global_iqr = safe_iqr(next_day_value),
    .groups = "drop"
  ) %>%
  mutate(global_median = ifelse(is.nan(global_median), NA_real_, global_median))

summary_median_iqr <- pair_long %>%
  group_by(current_subphenotype, variable_id, display_variable, domain, burden_direction, source_column, variable_order) %>%
  summarise(
    n_pairs = dplyr::n(),
    n_nonmissing = sum(!is.na(next_day_value)),
    median = median(next_day_value, na.rm = TRUE),
    q1 = as.numeric(stats::quantile(next_day_value, 0.25, na.rm = TRUE, names = FALSE)),
    q3 = as.numeric(stats::quantile(next_day_value, 0.75, na.rm = TRUE, names = FALSE)),
    iqr = safe_iqr(next_day_value),
    .groups = "drop"
  ) %>%
  mutate(
    median = ifelse(is.nan(median), NA_real_, median),
    q1 = ifelse(is.nan(q1), NA_real_, q1),
    q3 = ifelse(is.nan(q3), NA_real_, q3)
  ) %>%
  left_join(global_stats %>% select(variable_id, global_median, global_q1, global_q3, global_iqr, n_nonmissing_total),
            by = "variable_id") %>%
  mutate(
    scaled_median_difference = (median - global_median) / global_iqr,
    burden_scaled_median_difference = burden_direction * scaled_median_difference,
    cell_flag = dplyr::case_when(
      is.na(n_nonmissing) | n_nonmissing < MIN_NONMISSING_PER_CELL ~ "low_n",
      is.na(global_iqr) ~ "zero_or_missing_global_iqr",
      TRUE ~ "ok"
    )
  ) %>%
  arrange(factor(current_subphenotype, levels = SUBPH_LEVELS), variable_order)

heatmap_values <- summary_median_iqr %>%
  mutate(
    current_subphenotype = factor(current_subphenotype, levels = SUBPH_LEVELS),
    domain = factor(domain, levels = DOMAIN_LEVELS),
    display_variable = factor(display_variable, levels = available_specs$display_variable[order(available_specs$variable_order)])
  ) %>%
  arrange(current_subphenotype, domain, display_variable)

target_summary <- heatmap_values %>%
  mutate(current_subphenotype = as.character(current_subphenotype), display_variable = as.character(display_variable), domain = as.character(domain)) %>%
  inner_join(TARGET_MAP, by = c("current_subphenotype", "variable_id")) %>%
  arrange(factor(current_subphenotype, levels = SUBPH_LEVELS), variable_order)

domain_summary <- heatmap_values %>%
  filter(cell_flag == "ok") %>%
  group_by(current_subphenotype, domain) %>%
  summarise(
    n_variables = dplyr::n(),
    mean_burden_scaled_median_difference = mean(burden_scaled_median_difference, na.rm = TRUE),
    median_burden_scaled_median_difference = median(burden_scaled_median_difference, na.rm = TRUE),
    max_burden_scaled_median_difference = max(burden_scaled_median_difference, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(factor(current_subphenotype, levels = SUBPH_LEVELS), factor(domain, levels = DOMAIN_LEVELS))

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_heat <- heatmap_values %>%
    filter(cell_flag == "ok") %>%
    ggplot2::ggplot(ggplot2::aes(x = display_variable, y = current_subphenotype,
                                 fill = burden_scaled_median_difference)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", burden_scaled_median_difference)), size = 2.6) +
    ggplot2::facet_grid(. ~ domain, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0,
      name = "Burden-oriented\nmedian difference\n(per pooled IQR)"
    ) +
    ggplot2::labs(
      title = "Current HMM subphenotype and next-day organ-domain physiological burden",
      subtitle = "Day d dominant subphenotype versus day d+1 HMM observation variables; positive values indicate greater burden",
      x = "Next-day physiological variable",
      y = "Current HMM subphenotype"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
  ggplot2::ggsave(paste0(OUT_PREFIX, "_heatmap.pdf"), p_heat, width = 12.5, height = 5.8)
  ggplot2::ggsave(paste0(OUT_PREFIX, "_heatmap.png"), p_heat, width = 12.5, height = 5.8, dpi = 300)
}

write_both(pair_long, paste0(OUT_PREFIX, "_pairlevel_long"))
write_both(summary_median_iqr, paste0(OUT_PREFIX, "_summary_median_iqr"))
write_both(heatmap_values, paste0(OUT_PREFIX, "_heatmap_values"))
write_both(target_summary, paste0(OUT_PREFIX, "_target_domain_summary"))
write_both(domain_summary, paste0(OUT_PREFIX, "_domain_summary"))
write_both(resolved_vars %>% mutate(candidate_names = vapply(candidate_names, paste, character(1), collapse = ";")),
           paste0(OUT_PREFIX, "_variable_manifest"))

manifest <- tibble::tibble(
  key = c(
    "analysis", "out_prefix", "annotation_file", "physiology_source",
    "n_rows", "n_unique_stays", "n_adjacent_day_pairs", "n_pair_long_rows",
    "n_available_variables", "n_missing_requested_variables", "interpretation_note"
  ),
  value = c(
    "current_subphenotype_nextday_organ_domain_physiological_burden",
    OUT_PREFIX,
    anno_path,
    source_name,
    as.character(nrow(combined)),
    as.character(dplyr::n_distinct(combined$stay_id)),
    as.character(nrow(pairs)),
    as.character(nrow(pair_long)),
    as.character(nrow(available_specs)),
    as.character(nrow(missing_vars)),
    "Descriptive physiological-persistence check using HMM observation variables only; not treatment exposure or delivered organ-support analysis."
  )
)
write_both(manifest, paste0(OUT_PREFIX, "_manifest"))

message("Done: next-day organ-domain burden analysis written with prefix: ", OUT_PREFIX)
message("Physiology source: ", source_name)
message("Adjacent active ICU-day pairs: ", nrow(pairs))
message("Available variables: ", paste(available_specs$variable_id, collapse = ", "))
if (nrow(missing_vars) > 0) message("Missing variables: ", paste(missing_vars$variable_id, collapse = ", "))
