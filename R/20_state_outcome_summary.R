
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
  library(ggplot2)
  library(readr)
  library(stringr)
})

if (file.exists("18_plot_palette.R")) {
  source(code_path("18_plot_palette.R"))
}

ANNO_FILE <- MIMIC_ACTIVE_ANNOT_MASTER
MAX_DAY   <- 28L
LONGSTAY_DAY <- 21L
DAYS_KEEP <- c(1, 3, 5, 7, 10, 14, 21, 28)

.pick_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(df[[nm]])
  rep(NA, nrow(df))
}

safe_fname <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace("^_+", "") %>%
    stringr::str_replace("_+$", "")
}

shorten_tag <- function(tag) {
  if (is.na(tag) || tag == "") return("")
  tag2 <- gsub("\\([^)]*\\)", "", tag)
  tokens <- unlist(strsplit(tag2, "[^A-Za-z]+"))
  tokens <- tokens[nzchar(tokens)]
  tokens[tolower(tokens) == "hypoxemia"] <- "Hypoxemia"
  tokens_cap <- tokens[grepl("^[A-Z]", tokens)]
  if (!length(tokens_cap)) return("")
  parts <- substr(tokens_cap, 1, pmin(4, nchar(tokens_cap)))
  paste0(parts, collapse = "")
}

state_id_from_name <- function(x) as.integer(sub("^S", "", x))

sum_prob_by_group <- function(df, value_col, group_col) {
  df %>%
    group_by(.data[[group_col]]) %>%
    summarise(value = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop")
}

write_both <- function(df, stem) {
  write_tsv(df, paste0(stem, ".tsv"))
  write_csv(df, paste0(stem, ".csv"))
}

stopifnot(exists("post_smooth_full"), exists("combined"), exists("ids"), exists("icustays_intubated_first_patient_admission"))

icu <- icustays_intubated_first_patient_admission %>%
  mutate(stay_id = as.character(stay_id))

anno_tbl <- readr::read_csv(ANNO_FILE, show_col_types = FALSE) %>%
  rename(
    tags = tags_main,
    union_Orthogonal = orthogonal_types,
    union_ABGD = ABGD
  ) %>%
  mutate(state_id = as.integer(state_id))

required_anno_cols <- c("state_id", "subphenotype", "assigned_subtype", "tags", "union_Orthogonal")
missing_anno <- setdiff(required_anno_cols, names(anno_tbl))
if (length(missing_anno) > 0) {
  stop("Annotation file is missing required columns: ", paste(missing_anno, collapse = ", "))
}

N <- nrow(post_smooth_full)
stopifnot(length(ids) == N)

ids <- as.character(ids)
state_levels <- colnames(post_smooth_full)
K <- length(state_levels)
max_day <- MAX_DAY

survivor_ids <- icu %>% filter(is.na(dod)) %>% pull(stay_id)
death_ids    <- icu %>% filter(!is.na(dod)) %>% pull(stay_id)
is_death_stay <- ids %in% death_ids
death_flag    <- as.integer(is_death_stay)

state_base <- tibble(
  state = state_levels,
  state_id = state_id_from_name(state_levels)
) %>%
  left_join(
    anno_tbl %>% select(state_id, subphenotype, assigned_subtype, tags, union_Orthogonal),
    by = "state_id"
  )

subphenotype_base <- state_base %>%
  group_by(subphenotype) %>%
  summarise(n_states = n(), .groups = "drop")

N_state_any   <- colSums(post_smooth_full)
post_death    <- post_smooth_full * death_flag
N_state_death <- colSums(post_death)

state_prognosis <- tibble(
  state         = state_levels,
  state_id      = state_id_from_name(state_levels),
  N_state_any   = as.numeric(N_state_any),
  N_state_death = as.numeric(N_state_death)
) %>%
  mutate(
    p_death_eventual = ifelse(N_state_any > 0, N_state_death / N_state_any, NA_real_)
  ) %>%
  left_join(state_base %>% select(state_id, subphenotype, assigned_subtype, tags), by = "state_id")

subtype_prognosis <- state_prognosis %>%
  group_by(subphenotype) %>%
  summarise(
    n_states      = n(),
    N_state_any   = sum(N_state_any),
    N_state_death = sum(N_state_death),
    p_death_eventual = ifelse(N_state_any > 0, N_state_death / N_state_any, NA_real_),
    .groups = "drop"
  )

day_tables <- purrr::map(seq_len(max_day), function(d) {
  mask_d <- combined$lab_time == d
  if (!any(mask_d)) {
    return(tibble(subphenotype = subphenotype_base$subphenotype,
                  !!paste0("p_death_eventual_d", d) := NA_real_))
  }

  N_state_any_d   <- colSums(post_smooth_full[mask_d, , drop = FALSE])
  post_death_d    <- post_smooth_full[mask_d, , drop = FALSE] * death_flag[mask_d]
  N_state_death_d <- colSums(post_death_d)

  tibble(
    state = state_levels,
    state_id = state_id_from_name(state_levels),
    N_state_any   = as.numeric(N_state_any_d),
    N_state_death = as.numeric(N_state_death_d)
  ) %>%
    left_join(state_base %>% select(state_id, subphenotype), by = "state_id") %>%
    group_by(subphenotype) %>%
    summarise(
      N_state_any   = sum(N_state_any),
      N_state_death = sum(N_state_death),
      .groups = "drop"
    ) %>%
    mutate(
      !!paste0("p_death_eventual_d", d) := ifelse(N_state_any > 0, N_state_death / N_state_any, NA_real_)
    ) %>%
    select(subphenotype, starts_with("p_death_eventual_d"))
})

subphenotype_prognosis_by_day <- purrr::reduce(
  day_tables,
  .init = subphenotype_base,
  .f = ~ dplyr::left_join(.x, .y, by = "subphenotype")
)

plot_subphenotype_pdeath <- subphenotype_prognosis_by_day %>%
  pivot_longer(cols = starts_with("p_death_eventual_d"), names_to = "day", values_to = "p_death_eventual") %>%
  mutate(
    day = readr::parse_number(day),
    subphenotype = factor(subphenotype, levels = subphenotype_levels_order)
  )

p_subphenotype_pdeath <- ggplot(
  plot_subphenotype_pdeath,
  aes(x = day, y = p_death_eventual, color = subphenotype, linetype = subphenotype, shape = subphenotype, group = subphenotype)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  scale_color_manual(values = subphenotype_colors) +
  scale_linetype_manual(values = subphenotype_linetypes) +
  scale_shape_manual(values = subphenotype_shapes) +
  scale_x_continuous(breaks = seq_len(max_day), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "ICU day", y = "Prob. (eventual death)", title = "Day-wise eventual death probability by subphenotype") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA))

ggsave("state_outcomes_subphenotype_pdeath_timeseries.png", p_subphenotype_pdeath, width = 8, height = 5, dpi = 600, bg = "white")
ggsave("state_outcomes_subphenotype_pdeath_timeseries.pdf", p_subphenotype_pdeath, width = 8, height = 5, bg = "white", device = cairo_pdf)

icu_los <- icu %>%
  transmute(
    stay_id = as.character(stay_id),
    icu_intime = .pick_col(
      pick(everything()),
      c("intime", "icu_intime", "ICU_intime", "icu_in_time", "in_time", "admittime", "admit_time")
    ),
    icu_outtime = .pick_col(
      pick(everything()),
      c("outtime", "icu_outtime", "ICU_outtime", "out_time", "dischtime", "discharge_time")
    )
  ) %>%
  mutate(
    los_days = ifelse(
      !is.na(icu_intime) & !is.na(icu_outtime),
      as.numeric(difftime(icu_outtime, icu_intime, units = "days")),
      NA_real_
    ),
    los_days = ifelse(is.finite(los_days) & los_days >= 0, los_days, NA_real_)
  ) %>%
  distinct(stay_id, .keep_all = TRUE)

los_days_vec <- icu_los$los_days[match(ids, icu_los$stay_id)]
long_flag <- as.integer(!is.na(los_days_vec) & los_days_vec >= LONGSTAY_DAY)
stopifnot(length(long_flag) == nrow(post_smooth_full))

message(
  sprintf(
    "Long-stay definition: ICU LOS >= %d days from difftime(outtime, intime); missing LOS rows: %d/%d",
    LONGSTAY_DAY,
    sum(is.na(los_days_vec)),
    length(los_days_vec)
  )
)

N_state_long <- colSums(post_smooth_full * long_flag)

state_long_overall <- tibble(
  state = state_levels,
  state_id = state_id_from_name(state_levels),
  N_state_any = as.numeric(N_state_any),
  N_state_long = as.numeric(N_state_long),
  p_longstay_all = ifelse(N_state_any > 0, N_state_long / N_state_any, NA_real_)
) %>%
  left_join(state_base %>% select(state_id, subphenotype), by = "state_id")

subphenotype_long_overall <- state_long_overall %>%
  group_by(subphenotype) %>%
  summarise(
    n_states = n(),
    N_state_any = sum(N_state_any),
    N_state_long = sum(N_state_long),
    p_longstay_all = ifelse(N_state_any > 0, N_state_long / N_state_any, NA_real_),
    .groups = "drop"
  )

write_tsv(subphenotype_long_overall, sprintf("state_outcomes_subphenotype_longstay_overall_Y%d.tsv", LONGSTAY_DAY))
write_csv(subphenotype_long_overall, sprintf("state_outcomes_subphenotype_longstay_overall_Y%d.csv", LONGSTAY_DAY))

max_cond_day <- min(max_day, LONGSTAY_DAY - 1L)

long_day_tables <- purrr::map(seq_len(max_cond_day), function(d) {
  mask_d <- combined$lab_time == d
  if (!any(mask_d)) {
    return(tibble(subphenotype = subphenotype_base$subphenotype,
                  !!paste0("p_longstay_d", d) := NA_real_))
  }

  N_state_any_d <- colSums(post_smooth_full[mask_d, , drop = FALSE])
  post_long_d <- post_smooth_full[mask_d, , drop = FALSE] * long_flag[mask_d]
  N_state_long_d <- colSums(post_long_d)

  tibble(
    state = state_levels,
    state_id = state_id_from_name(state_levels),
    N_state_any = as.numeric(N_state_any_d),
    N_state_long = as.numeric(N_state_long_d)
  ) %>%
    left_join(state_base %>% select(state_id, subphenotype), by = "state_id") %>%
    group_by(subphenotype) %>%
    summarise(
      N_state_any = sum(N_state_any),
      N_state_long = sum(N_state_long),
      .groups = "drop"
    ) %>%
    mutate(
      !!paste0("p_longstay_d", d) := ifelse(N_state_any > 0, N_state_long / N_state_any, NA_real_)
    ) %>%
    select(subphenotype, starts_with("p_longstay_d"))
})

subphenotype_longstay_by_day <- purrr::reduce(
  long_day_tables,
  .init = subphenotype_base,
  .f = ~ dplyr::left_join(.x, .y, by = "subphenotype")
)

plot_subphenotype_long <- subphenotype_longstay_by_day %>%
  pivot_longer(cols = starts_with("p_longstay_d"), names_to = "day", values_to = "p_longstay") %>%
  mutate(
    day = readr::parse_number(day),
    subphenotype = factor(subphenotype, levels = subphenotype_levels_order)
  )

p_subphenotype_long <- ggplot(
  plot_subphenotype_long,
  aes(x = day, y = p_longstay, color = subphenotype, linetype = subphenotype, shape = subphenotype, group = subphenotype)
) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  scale_color_manual(values = subphenotype_colors) +
  scale_linetype_manual(values = subphenotype_linetypes) +
  scale_shape_manual(values = subphenotype_shapes) +
  scale_x_continuous(breaks = seq_len(max_cond_day), minor_breaks = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  labs(x = "ICU day", y = sprintf("P(LOS ≥ %d days)", LONGSTAY_DAY), title = sprintf("Day-wise long-stay probability (LOS ≥ %d days) by subphenotype", LONGSTAY_DAY)) +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA))

ggsave(sprintf("state_outcomes_subphenotype_longstay_toY%d_timeseries.png", LONGSTAY_DAY), p_subphenotype_long, width = 8, height = 5, dpi = 600, bg = "white")
ggsave(sprintf("state_outcomes_subphenotype_longstay_toY%d_timeseries.pdf", LONGSTAY_DAY), p_subphenotype_long, width = 8, height = 5, bg = "white", device = cairo_pdf)

state_day_list <- purrr::map(seq_len(max_day), function(d) {
  mask_d <- combined$lab_time == d
  if (!any(mask_d)) {
    return(tibble(state = state_levels, !!paste0("p_death_eventual_d", d) := NA_real_))
  }
  N_state_any_d <- colSums(post_smooth_full[mask_d, , drop = FALSE])
  post_death_d <- post_smooth_full[mask_d, , drop = FALSE] * death_flag[mask_d]
  N_state_death_d <- colSums(post_death_d)
  p_d <- ifelse(N_state_any_d > 0, N_state_death_d / N_state_any_d, NA_real_)
  tibble(state = state_levels, !!paste0("p_death_eventual_d", d) := as.numeric(p_d))
})

state_long_day_list <- purrr::map(seq_len(max_cond_day), function(d) {
  mask_d <- combined$lab_time == d
  if (!any(mask_d)) {
    return(tibble(state = state_levels, !!paste0("p_longstay_d", d) := NA_real_))
  }
  N_state_any_d <- colSums(post_smooth_full[mask_d, , drop = FALSE])
  post_long_d <- post_smooth_full[mask_d, , drop = FALSE] * long_flag[mask_d]
  N_state_long_d <- colSums(post_long_d)
  p_d <- ifelse(N_state_any_d > 0, N_state_long_d / N_state_any_d, NA_real_)
  tibble(state = state_levels, !!paste0("p_longstay_d", d) := as.numeric(p_d))
})

state_prognosis_by_day <- purrr::reduce(
  state_day_list,
  .init = tibble(state = state_levels),
  .f = ~ dplyr::left_join(.x, .y, by = "state")
) %>%
  mutate(state_id = state_id_from_name(state)) %>%
  left_join(state_base %>% select(state_id, subphenotype, tags), by = "state_id")

state_longstay_by_day <- purrr::reduce(
  state_long_day_list,
  .init = tibble(state = state_levels),
  .f = ~ dplyr::left_join(.x, .y, by = "state")
) %>%
  mutate(state_id = state_id_from_name(state)) %>%
  left_join(state_base %>% select(state_id, subphenotype), by = "state_id")

state_tags <- state_prognosis_by_day %>%
  select(state_id, subphenotype, tags) %>%
  filter(!is.na(subphenotype)) %>%
  mutate(tag_list = str_split(tags, "\\s*;\\s*"))

subphenotype_common_tags <- state_tags %>%
  group_by(subphenotype) %>%
  summarise(
    tags_common = list({
      lst <- tag_list[lengths(tag_list) > 0]
      if (length(lst) == 0) character(0) else Reduce(intersect, lst)
    }),
    .groups = "drop"
  )

state_tags_resid <- state_tags %>%
  left_join(subphenotype_common_tags, by = "subphenotype") %>%
  mutate(
    residual_tags = purrr::map2(tag_list, tags_common, setdiff),
    residual_short = purrr::map_chr(
      residual_tags,
      ~ {
        if (length(.x) == 0) return("")
        shorts <- vapply(.x, shorten_tag, character(1))
        shorts <- shorts[shorts != ""]
        if (!length(shorts)) "" else paste(shorts, collapse = "; ")
      }
    )
  )

state_labels <- state_tags_resid %>%
  mutate(
    state = paste0("S", state_id),
    state_label = if_else(
      residual_short == "" | is.na(residual_short),
      state,
      paste0(state, " (", residual_short, ")")
    )
  ) %>%
  select(state_id, subphenotype, state, state_label)

subphenotype_common_labels_short <- subphenotype_common_tags %>%
  mutate(
    common_short = purrr::map_chr(
      tags_common,
      ~ {
        if (length(.x) == 0) return("")
        shorts <- vapply(.x, shorten_tag, character(1))
        shorts <- shorts[shorts != ""]
        if (!length(shorts)) "" else paste0(shorts, collapse = ";")
      }
    )
  )

state_prognosis_labeled <- state_prognosis_by_day %>%
  left_join(state_labels, by = c("state_id", "subphenotype"))

state_longstay_labeled <- state_longstay_by_day %>%
  left_join(state_labels, by = c("state_id", "subphenotype"))

plot_df_all <- state_prognosis_labeled %>%
  pivot_longer(cols = starts_with("p_death_eventual_d"), names_to = "day", values_to = "p_death_eventual") %>%
  mutate(day = readr::parse_number(day))

plot_long_all <- state_longstay_labeled %>%
  pivot_longer(cols = starts_with("p_longstay_d"), names_to = "day", values_to = "p_longstay") %>%
  mutate(day = readr::parse_number(day))

subphenotype_levels <- unique(plot_df_all$subphenotype)
subphenotype_levels_long <- unique(plot_long_all$subphenotype)

for (mg in subphenotype_levels_long) {
  df_g <- plot_long_all %>% filter(subphenotype == mg)
  if (nrow(df_g) == 0) next

  lab_common_short_vec <- subphenotype_common_labels_short %>%
    filter(subphenotype == mg) %>%
    pull(common_short)
  lab_common_short <- if (length(lab_common_short_vec) == 0 || is.na(lab_common_short_vec[1])) "" else lab_common_short_vec[1]
  title_mg <- if (lab_common_short == "") mg else paste0(mg, " (", lab_common_short, ")")

  n_states_g <- n_distinct(df_g$state_label)
  ncol_legend <- if (mg == "D1 lactate/shock") 1L else dplyr::case_when(
    n_states_g <= 6 ~ 2L,
    n_states_g <= 12 ~ 3L,
    TRUE ~ 3L
  )

  p_g_long <- ggplot(df_g, aes(x = day, y = p_longstay, colour = state_label, group = state_label)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 10, 14, 21, 28)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "ICU day", y = sprintf("P(LOS ≥ %d days)", LONGSTAY_DAY), colour = "State (residual tags)", title = title_mg) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8), legend.title = element_text(size = 9), legend.key.width = grid::unit(0.8, "lines"), legend.key.height = grid::unit(0.8, "lines"), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA)) +
    guides(colour = guide_legend(ncol = ncol_legend, byrow = TRUE, title.position = "top"))

  fname_safe <- safe_fname(mg)
  ggsave(paste0("state_outcomes_state_longstay_timeseries_", fname_safe, ".png"), p_g_long, width = 8, height = 5, dpi = 600, bg = "white")
  ggsave(paste0("state_outcomes_state_longstay_timeseries_", fname_safe, ".pdf"), p_g_long, width = 8, height = 5, bg = "white", device = cairo_pdf)
}

for (mg in subphenotype_levels) {
  df_g <- plot_df_all %>% filter(subphenotype == mg)
  if (nrow(df_g) == 0) next

  lab_common_short_vec <- subphenotype_common_labels_short %>%
    filter(subphenotype == mg) %>%
    pull(common_short)
  lab_common_short <- if (length(lab_common_short_vec) == 0 || is.na(lab_common_short_vec[1])) "" else lab_common_short_vec[1]
  title_mg <- if (lab_common_short == "") mg else paste0(mg, " (", lab_common_short, ")")

  n_states_g <- n_distinct(df_g$state_label)
  ncol_legend <- if (mg == "D1 lactate/shock") 1L else dplyr::case_when(
    n_states_g <= 6 ~ 2L,
    n_states_g <= 12 ~ 3L,
    TRUE ~ 3L
  )

  p_g <- ggplot(df_g, aes(x = day, y = p_death_eventual, colour = state_label, group = state_label)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 10, 14, 21, 28)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "ICU day", y = "P(eventual death)", colour = "State (residual tags)", title = title_mg) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8), legend.title = element_text(size = 9), legend.key.width = grid::unit(0.8, "lines"), legend.key.height = grid::unit(0.8, "lines"), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA)) +
    guides(colour = guide_legend(ncol = ncol_legend, byrow = TRUE, title.position = "top"))

  fname_safe <- safe_fname(mg)
  ggsave(paste0("state_outcomes_state_pdeath_timeseries_", fname_safe, ".png"), p_g, width = 8, height = 5, dpi = 600, bg = "white")
  ggsave(paste0("state_outcomes_state_pdeath_timeseries_", fname_safe, ".pdf"), p_g, width = 8, height = 5, bg = "white", device = cairo_pdf)
}

state_day_tables <- purrr::map(seq_len(max_day), function(d) {
  mask_d <- combined$lab_time == d
  if (!any(mask_d)) {
    return(tibble(state = state_levels, !!paste0("p_state_d", d) := NA_real_))
  }
  n_d <- sum(mask_d)
  N_state_any_d <- colSums(post_smooth_full[mask_d, , drop = FALSE])
  p_state_d <- N_state_any_d / n_d
  tibble(state = state_levels, !!paste0("p_state_d", d) := as.numeric(p_state_d))
})

state_occ_by_day <- purrr::reduce(
  state_day_tables,
  .init = tibble(state = state_levels),
  .f = ~ dplyr::left_join(.x, .y, by = "state")
)

subphenotype_occ_by_day <- state_occ_by_day %>%
  mutate(state_id = state_id_from_name(state)) %>%
  left_join(state_base %>% select(state_id, subphenotype), by = "state_id") %>%
  group_by(subphenotype) %>%
  summarise(
    n_states = n(),
    across(starts_with("p_state_d"), ~ sum(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename_with(~ sub("^p_state_", "p_subphenotype_", .x), starts_with("p_state_d"))

plot_subphenotype_occ <- subphenotype_occ_by_day %>%
  pivot_longer(cols = starts_with("p_subphenotype_d"), names_to = "day", values_to = "p_subphenotype") %>%
  mutate(day = readr::parse_number(day), subphenotype = factor(subphenotype, levels = subphenotype_levels_order))

p_subphenotype_occ <- ggplot(plot_subphenotype_occ, aes(x = day, y = p_subphenotype, color = subphenotype, linetype = subphenotype, shape = subphenotype, group = subphenotype)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.8) +
  scale_color_manual(values = subphenotype_colors) +
  scale_linetype_manual(values = subphenotype_linetypes) +
  scale_shape_manual(values = subphenotype_shapes) +
  scale_x_continuous(breaks = seq_len(max_day)) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "ICU day", y = "P(being in this subphenotype)", title = "Day-wise occupancy probability by subphenotype") +
  theme_bw(base_size = 13) +
  theme(panel.grid.minor = element_blank(), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA))

ggsave("state_outcomes_subphenotype_occupancy_timeseries.png", p_subphenotype_occ, width = 8, height = 5, dpi = 600, bg = "white")
ggsave("state_outcomes_subphenotype_occupancy_timeseries.pdf", p_subphenotype_occ, width = 8, height = 5, bg = "white", device = cairo_pdf)

state_ortho <- state_prognosis_by_day %>%
  select(state, state_id, subphenotype) %>%
  left_join(state_base %>% select(state_id, union_Orthogonal), by = "state_id") %>%
  separate_rows(union_Orthogonal, sep = "\\s*;\\s*") %>%
  filter(!is.na(union_Orthogonal), union_Orthogonal != "") %>%
  rename(ortho = union_Orthogonal) %>%
  distinct()

ortho_levels <- sort(unique(state_ortho$ortho))

state_pdeath_long <- state_prognosis_by_day %>%
  pivot_longer(cols = starts_with("p_death_eventual_d"), names_to = "day", values_to = "p_death_eventual") %>%
  mutate(day = readr::parse_number(day))

for (o in ortho_levels) {
  df_o <- state_pdeath_long %>%
    inner_join(state_ortho %>% filter(ortho == o), by = c("state", "state_id", "subphenotype")) %>%
    mutate(state_subphenotype_label = paste0(state, " (", subphenotype, ")"))
  if (nrow(df_o) == 0) next

  n_states_o <- n_distinct(df_o$state_subphenotype_label)
  ncol_legend <- dplyr::case_when(n_states_o <= 6 ~ 1L, n_states_o <= 12 ~ 2L, TRUE ~ 3L)

  p_o <- ggplot(df_o, aes(x = day, y = p_death_eventual, colour = state_subphenotype_label, group = state_subphenotype_label)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = c(1, 3, 5, 7, 10, 14, 21, 28)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(x = "ICU day", y = "P(eventual death)", colour = "State (subphenotype)", title = paste0("Orthogonal subtype: ", o)) +
    theme_bw(base_size = 13) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom", legend.text = element_text(size = 8), legend.title = element_text(size = 9), legend.key.width = grid::unit(0.8, "lines"), legend.key.height = grid::unit(0.8, "lines"), plot.background = element_rect(fill = "white", colour = NA), panel.background = element_rect(fill = "white", colour = NA)) +
    guides(colour = guide_legend(ncol = ncol_legend, byrow = TRUE, title.position = "top"))

  fname_safe <- safe_fname(o)
  ggsave(paste0("state_outcomes_ortho_pdeath_timeseries_", fname_safe, ".png"), p_o, width = 8, height = 5, dpi = 600, bg = "white")
  ggsave(paste0("state_outcomes_ortho_pdeath_timeseries_", fname_safe, ".pdf"), p_o, width = 8, height = 5, bg = "white", device = cairo_pdf)
}

death_cols <- paste0("p_death_eventual_d", DAYS_KEEP)
occ_cols   <- paste0("p_subphenotype_d", DAYS_KEEP)

subphenotype_death_table <- subphenotype_prognosis_by_day %>%
  select(subphenotype, all_of(death_cols)) %>%
  mutate(subphenotype = factor(subphenotype, levels = subphenotype_levels_order)) %>%
  arrange(subphenotype) %>%
  mutate(across(starts_with("p_death_eventual_d"), ~ round(.x, 3)))

write_both(subphenotype_death_table, "state_outcomes_table_subphenotype_pdeath_timeseries_selected")

ortho_base <- tibble(ortho = ortho_levels)
ortho_day_tables <- purrr::map(seq_len(max_day), function(d) {
  mask_d <- combined$lab_time == d
  if (!any(mask_d)) {
    return(tibble(ortho = ortho_base$ortho, !!paste0("p_death_eventual_d", d) := NA_real_))
  }

  N_state_any_d   <- colSums(post_smooth_full[mask_d, , drop = FALSE])
  post_death_d    <- post_smooth_full[mask_d, , drop = FALSE] * death_flag[mask_d]
  N_state_death_d <- colSums(post_death_d)

  tibble(
    state = state_levels,
    state_id = state_id_from_name(state_levels),
    N_state_any   = as.numeric(N_state_any_d),
    N_state_death = as.numeric(N_state_death_d)
  ) %>%
    left_join(state_ortho %>% select(state_id, ortho) %>% distinct(), by = "state_id") %>%
    filter(!is.na(ortho), ortho != "") %>%
    group_by(ortho) %>%
    summarise(
      N_state_any   = sum(N_state_any),
      N_state_death = sum(N_state_death),
      .groups = "drop"
    ) %>%
    mutate(!!paste0("p_death_eventual_d", d) := ifelse(N_state_any > 0, N_state_death / N_state_any, NA_real_)) %>%
    select(ortho, starts_with("p_death_eventual_d")) %>%
    right_join(ortho_base, by = "ortho")
})

ortho_prognosis_by_day <- purrr::reduce(
  ortho_day_tables,
  .init = ortho_base,
  .f = ~ dplyr::left_join(.x, .y, by = "ortho")
)

ortho_death_table <- ortho_prognosis_by_day %>%
  select(ortho, all_of(death_cols)) %>%
  arrange(ortho) %>%
  mutate(across(starts_with("p_death_eventual_d"), ~ round(.x, 3)))

write_both(ortho_death_table, "state_outcomes_table_ortho_pdeath_timeseries_selected")

subphenotype_occ_table <- subphenotype_occ_by_day %>%
  select(subphenotype, all_of(occ_cols)) %>%
  mutate(subphenotype = factor(subphenotype, levels = subphenotype_levels_order)) %>%
  arrange(subphenotype) %>%
  mutate(across(starts_with("p_subphenotype_d"), ~ round(.x, 3)))

write_both(subphenotype_occ_table, "state_outcomes_table_subphenotype_occupancy_timeseries_selected")

icu_intime <- .pick_col(
  icu,
  c("intime", "icu_intime", "ICU_intime", "icu_in_time", "in_time", "admittime", "admit_time")
)

icu_for_death <- icu %>% mutate(icu_intime = icu_intime)

death_day_tbl <- icu_for_death %>%
  transmute(
    stay_id = as.character(stay_id),
    dod = dod,
    icu_intime = icu_intime
  ) %>%
  mutate(
    death_days = ifelse(!is.na(dod) & !is.na(icu_intime), as.numeric(difftime(dod, icu_intime, units = "days")), NA_real_),
    death_days = ifelse(!is.na(death_days) & death_days < 0, NA_real_, death_days)
  ) %>%
  select(stay_id, death_days)

death_days_vec <- death_day_tbl$death_days[match(ids, death_day_tbl$stay_id)]
death3_flag <- as.integer(!is.na(death_days_vec) & death_days_vec <= (combined$lab_time + 2L))
death7_flag <- as.integer(!is.na(death_days_vec) & death_days_vec <= (combined$lab_time + 6L))

stopifnot(length(death3_flag) == nrow(post_smooth_full))
stopifnot(length(death7_flag) == nrow(post_smooth_full))

post_death3 <- post_smooth_full * death3_flag
post_death7 <- post_smooth_full * death7_flag
N_state_death3 <- colSums(post_death3)
N_state_death7 <- colSums(post_death7)

state_prognosis_3d7d <- state_prognosis %>%
  mutate(
    N_state_death3 = as.numeric(N_state_death3),
    N_state_death7 = as.numeric(N_state_death7),
    p_death_3d = ifelse(N_state_any > 0, N_state_death3 / N_state_any, NA_real_),
    p_death_7d = ifelse(N_state_any > 0, N_state_death7 / N_state_any, NA_real_)
  )

subphenotype_prognosis_overall_3d7d <- state_prognosis_3d7d %>%
  group_by(subphenotype) %>%
  summarise(
    n_states = n(),
    N_state_any = sum(N_state_any),
    N_state_death3 = sum(N_state_death3),
    N_state_death7 = sum(N_state_death7),
    p_death_3d = ifelse(N_state_any > 0, N_state_death3 / N_state_any, NA_real_),
    p_death_7d = ifelse(N_state_any > 0, N_state_death7 / N_state_any, NA_real_),
    .groups = "drop"
  )

write_both(subphenotype_prognosis_overall_3d7d, "state_outcomes_table_subphenotype_pdeath_3d7d_overall")

day_tables_3d <- purrr::map(seq_len(max_day), function(d) {
  mask_d <- combined$lab_time == d
  if (!any(mask_d)) {
    return(tibble(subphenotype = subphenotype_base$subphenotype, !!paste0("p_death_3d_d", d) := NA_real_))
  }

  N_state_any_d <- colSums(post_smooth_full[mask_d, , drop = FALSE])
  post_death3_d <- post_smooth_full[mask_d, , drop = FALSE] * death3_flag[mask_d]
  N_state_death3_d <- colSums(post_death3_d)

  tibble(
    state = state_levels,
    state_id = state_id_from_name(state_levels),
    N_state_any = as.numeric(N_state_any_d),
    N_state_death3 = as.numeric(N_state_death3_d)
  ) %>%
    left_join(state_base %>% select(state_id, subphenotype), by = "state_id") %>%
    group_by(subphenotype) %>%
    summarise(
      N_state_any = sum(N_state_any),
      N_state_death3 = sum(N_state_death3),
      .groups = "drop"
    ) %>%
    mutate(!!paste0("p_death_3d_d", d) := ifelse(N_state_any > 0, N_state_death3 / N_state_any, NA_real_)) %>%
    select(subphenotype, starts_with("p_death_3d_d"))
})

day_tables_7d <- purrr::map(seq_len(max_day), function(d) {
  mask_d <- combined$lab_time == d
  if (!any(mask_d)) {
    return(tibble(subphenotype = subphenotype_base$subphenotype, !!paste0("p_death_7d_d", d) := NA_real_))
  }

  N_state_any_d <- colSums(post_smooth_full[mask_d, , drop = FALSE])
  post_death7_d <- post_smooth_full[mask_d, , drop = FALSE] * death7_flag[mask_d]
  N_state_death7_d <- colSums(post_death7_d)

  tibble(
    state = state_levels,
    state_id = state_id_from_name(state_levels),
    N_state_any = as.numeric(N_state_any_d),
    N_state_death7 = as.numeric(N_state_death7_d)
  ) %>%
    left_join(state_base %>% select(state_id, subphenotype), by = "state_id") %>%
    group_by(subphenotype) %>%
    summarise(
      N_state_any = sum(N_state_any),
      N_state_death7 = sum(N_state_death7),
      .groups = "drop"
    ) %>%
    mutate(!!paste0("p_death_7d_d", d) := ifelse(N_state_any > 0, N_state_death7 / N_state_any, NA_real_)) %>%
    select(subphenotype, starts_with("p_death_7d_d"))
})

subphenotype_pdeath3_by_day <- purrr::reduce(day_tables_3d, .init = subphenotype_base, .f = ~ dplyr::left_join(.x, .y, by = "subphenotype"))
subphenotype_pdeath7_by_day <- purrr::reduce(day_tables_7d, .init = subphenotype_base, .f = ~ dplyr::left_join(.x, .y, by = "subphenotype"))

death3_cols <- paste0("p_death_3d_d", DAYS_KEEP)
death7_cols <- paste0("p_death_7d_d", DAYS_KEEP)

subphenotype_death3_table <- subphenotype_pdeath3_by_day %>%
  select(subphenotype, all_of(death3_cols)) %>%
  mutate(subphenotype = factor(subphenotype, levels = subphenotype_levels_order)) %>%
  arrange(subphenotype) %>%
  mutate(across(starts_with("p_death_3d_d"), ~ round(.x, 3)))

subphenotype_death7_table <- subphenotype_pdeath7_by_day %>%
  select(subphenotype, all_of(death7_cols)) %>%
  mutate(subphenotype = factor(subphenotype, levels = subphenotype_levels_order)) %>%
  arrange(subphenotype) %>%
  mutate(across(starts_with("p_death_7d_d"), ~ round(.x, 3)))

write_both(subphenotype_death3_table, "state_outcomes_table_subphenotype_pdeath_3d_timeseries_selected")
write_both(subphenotype_death7_table, "state_outcomes_table_subphenotype_pdeath_7d_timeseries_selected")

.make_long <- function(df_wide, prefix) {
  df_wide %>%
    pivot_longer(cols = starts_with(prefix), names_to = "metric_day", values_to = "p_death") %>%
    mutate(day = as.integer(str_extract(metric_day, "(?<=_d)\\d+")), p_death = as.numeric(p_death)) %>%
    select(subphenotype, day, p_death) %>%
    arrange(subphenotype, day)
}

.plot_pdeath_ts <- function(df_long, out_prefix, title, ylab) {
  lv <- c("A1 low burden", "A2 mild burden", "B1 renal mild", "B2 renal severe", "G respiratory", "D1 lactate/shock", "D2 hepatic", "D3 shock/coagulation-overlap")
  df_long <- df_long %>% mutate(subphenotype = factor(subphenotype, levels = lv[lv %in% unique(subphenotype)]))

  p <- ggplot(df_long, aes(x = day, y = p_death, group = subphenotype, color = subphenotype)) +
    geom_line(linewidth = 1.0, alpha = 0.9) +
    geom_point(size = 2.0, alpha = 0.95) +
    scale_x_continuous(breaks = sort(unique(df_long$day))) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = title, x = "ICU day", y = ylab, color = "Subphenotype") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())

  if (exists("SUBPHENO_COLORS")) {
    cols <- get("SUBPHENO_COLORS")
    if (is.character(cols) && !is.null(names(cols))) {
      keep <- intersect(levels(df_long$subphenotype), names(cols))
      if (length(keep) >= 2) {
        p <- p + scale_color_manual(values = cols[keep], drop = FALSE)
      }
    }
  }

  ggsave(paste0(out_prefix, ".pdf"), p, width = 10, height = 6, device = cairo_pdf)
  ggsave(paste0(out_prefix, ".png"), p, width = 10, height = 6, dpi = 300)
  invisible(p)
}

l3 <- .make_long(subphenotype_death3_table, prefix = "p_death_3d_d")
.plot_pdeath_ts(l3, "state_outcomes_subphenotype_pdeath_3d_timeseries", "Subphenotype time series: 3-day death probability", "P(death within 3 days | at day d, weighted)")

l7 <- .make_long(subphenotype_death7_table, prefix = "p_death_7d_d")
.plot_pdeath_ts(l7, "state_outcomes_subphenotype_pdeath_7d_timeseries", "Subphenotype time series: 7-day death probability", "P(death within 7 days | at day d, weighted)")

message("State outcome summary completed.")
