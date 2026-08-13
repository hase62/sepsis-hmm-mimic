
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(lubridate)
})

files <- c(
  lab = "mimiciv.lab.wide.csv",
  vital = "mimiciv.vital.wide.csv",
  procedure = "mimiciv.proc.wide.csv"
)
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0L) {
  stop("Missing required wide files: ", paste(missing_files, collapse = ", "))
}

read_keys <- function(path, source_name) {
  read_csv(
    path,
    col_select = c(stay_id, hadm_id, charttime_year_day),
    show_col_types = FALSE
  ) %>%
    transmute(
      stay_id = as.character(stay_id),
      hadm_id = as.character(hadm_id),
      day = as.Date(charttime_year_day),
      source = source_name,
      present = 1L
    ) %>%
    distinct()
}

keys <- bind_rows(
  read_keys(files[["lab"]], "lab"),
  read_keys(files[["vital"]], "vital"),
  read_keys(files[["procedure"]], "procedure")
) %>%
  pivot_wider(names_from = source, values_from = present, values_fill = 0L)

for (nm in c("lab", "vital", "procedure")) {
  if (!nm %in% names(keys)) keys[[nm]] <- 0L
}

lab_range <- keys %>%
  filter(lab == 1L) %>%
  group_by(stay_id) %>%
  summarise(
    first_lab_day = min(day),
    last_lab_day = max(day),
    n_lab_days = n(),
    .groups = "drop"
  )

qc <- keys %>%
  left_join(lab_range, by = "stay_id") %>%
  mutate(
    included_primary = lab == 1L,
    source_pattern = paste0("L", lab, "_V", vital, "_P", procedure),
    excluded_day_position = case_when(
      lab == 1L ~ "included_lab_anchor",
      is.na(first_lab_day) ~ "no_lab_day_in_stay_without_any_lab_anchor",
      day < first_lab_day ~ "before_first_lab_anchor",
      day > last_lab_day ~ "after_last_lab_anchor",
      TRUE ~ "between_lab_anchor_days"
    ),
    days_from_first_lab = as.integer(day - first_lab_day)
  )

summary1 <- qc %>%
  count(source_pattern, included_primary, excluded_day_position, name = "n_stay_days") %>%
  arrange(desc(n_stay_days))

summary2 <- qc %>%
  group_by(stay_id) %>%
  summarise(
    n_union_days = n(),
    n_included_lab_days = sum(included_primary),
    n_excluded_vital_or_proc_days = sum(!included_primary),
    n_before_first_lab = sum(excluded_day_position == "before_first_lab_anchor"),
    n_between_lab_days = sum(excluded_day_position == "between_lab_anchor_days"),
    n_after_last_lab = sum(excluded_day_position == "after_last_lab_anchor"),
    .groups = "drop"
  )

write_csv(summary1, "analytical_day_detailed_summary.csv")
write_csv(summary2, "analytical_day_by_stay.csv")

cat("Union stay-days:", nrow(qc), "\n")
cat("Primary lab-anchored stay-days:", sum(qc$included_primary), "\n")
cat("Excluded vital/procedure-only stay-days:", sum(!qc$included_primary), "\n")
cat("Excluded before first lab:", sum(qc$excluded_day_position == "before_first_lab_anchor"), "\n")
cat("Excluded between lab days:", sum(qc$excluded_day_position == "between_lab_anchor_days"), "\n")
cat("Excluded after last lab:", sum(qc$excluded_day_position == "after_last_lab_anchor"), "\n")
print(summary1, n = Inf)
