
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
if (!exists("MIMIC_ACTIVE_ANNOT_MASTER", inherits = TRUE)) {
  source(code_path("17_analysis_config.R"))
}

if (!file.exists(code_path("17_analysis_config.R"))) stop("Missing shared analysis configuration: 17_analysis_config.R", call. = FALSE)
source(code_path("17_analysis_config.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(lubridate)
})

options(datatable.integer64 = "double")

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && !is.na(a) && nzchar(a)) a else b

WORK_DIR <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

MIMIC_DIR <- normalizePath(mimic_root, winslash = "/", mustWork = TRUE)

HMM_ROW_SOURCE <- Sys.getenv("HMM_ROW_SOURCE", unset = "mimiciv.combined.wide.HMM_raw.csv")
if (!file.exists(file.path(WORK_DIR, HMM_ROW_SOURCE))) {
  stop("Final analytical row source is missing: ", HMM_ROW_SOURCE, call. = FALSE)
}

COHORT_FILE <- file.path(WORK_DIR, "icustays_intubated_first_patient_admission.HMM_eligible.csv")
stopifnot(file.exists(COHORT_FILE))

OUT_PREFIX <- "patient_day_sofa_organsupport"

LAB_ITEMIDS <- list(
  pao2       = c(50821),
  platelets  = c(51265, 51266),
  bilirubin  = c(50885, 225690),
  creatinine = c(50912, 52546)
)
LAB_ALL <- sort(unique(unlist(LAB_ITEMIDS)))

CH_ITEMIDS <- list(
  fio2_chart      = c(223835),
  map             = c(220052, 220181, 225312),
  gcs_motor       = c(223901),
  gcs_verbal      = c(223900),
  gcs_eyes        = c(220739),
  peep_set        = c(220339),
  peep_backup     = c(224700, 224699),
  vent_mode       = c(223849, 229314),
  vent_type       = c(223848),
  o2_device       = c(226732),
  ecmo_chart      = c(224660, 229270, 229280, 229841, 229842)
)
CH_ALL <- sort(unique(unlist(CH_ITEMIDS)))

VENT_PROC_PRIMARY <- c(225792)
INTUB_PROC_EVENT  <- c(224385)
VENT_PROC_ALL     <- c(VENT_PROC_PRIMARY, INTUB_PROC_EVENT)

VASO_ITEMIDS <- list(
  norepi = c(221906),
  epi    = c(221289, 229617),
  dopa   = c(221662),
  dobu   = c(221653),
  phenyl = c(221749, 229630, 229631, 229632, 229789),
  vaso   = c(222315),
  angii  = c(229709, 229764)
)
VASO_SOFA <- sort(unique(unlist(VASO_ITEMIDS[c("norepi", "epi", "dopa", "dobu")])))
VASO_ALL  <- sort(unique(unlist(VASO_ITEMIDS)))

URINE_IDS <- c(226559, 226560, 226561, 226584, 226563, 226564, 226565, 226567,
               226557, 226558, 227488, 227489)
IRRIGANT_IDS <- c(227488)

RRT_ITEMIDS <- c(
  225441,
  225802,
  225803,
  225809,
  225955,
  225805,
  227290,
  226457,
  226499
)
RRT_CRRT_ITEMIDS <- c(225802, 225803, 225809, 225955, 227290, 226457)
RRT_HD_ITEMIDS   <- c(225441, 226499)

ECMO_ITEMIDS <- c(224660, 229270, 229280, 229841, 229842)

clean_fio2_fraction <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[!is.na(x) & x == 0] <- NA_real_

  x[!is.na(x) & x > 1 & x <= 100] <- x[!is.na(x) & x > 1 & x <= 100] / 100
  x[!is.na(x) & x < 0.21] <- NA_real_
  x[!is.na(x) & x > 1.00] <- NA_real_
  x
}

score_resp <- function(pfr_min, vent_active) {
  fcase(
    !is.na(pfr_min) & pfr_min < 100 & vent_active, 4,
    !is.na(pfr_min) & pfr_min < 200 & vent_active, 3,
    !is.na(pfr_min) & pfr_min < 300, 2,
    !is.na(pfr_min) & pfr_min < 400, 1,
    default = 0
  )
}

score_coag <- function(plt_min) {
  fcase(
    !is.na(plt_min) & plt_min < 20, 4,
    !is.na(plt_min) & plt_min < 50, 3,
    !is.na(plt_min) & plt_min < 100, 2,
    !is.na(plt_min) & plt_min < 150, 1,
    default = 0
  )
}

score_liver <- function(bili_max) {
  fcase(
    !is.na(bili_max) & bili_max >= 12.0, 4,
    !is.na(bili_max) & bili_max >= 6.0, 3,
    !is.na(bili_max) & bili_max >= 2.0, 2,
    !is.na(bili_max) & bili_max >= 1.2, 1,
    default = 0
  )
}

score_cv <- function(map_min, dopa_max, epi_max, norepi_max, dobu_max) {
  out <- rep(0L, length(map_min))
  out[!is.na(map_min) & map_min < 70] <- 1L
  out[(!is.na(dopa_max) & dopa_max > 0) | (!is.na(dobu_max) & dobu_max > 0)] <- 2L
  out[(!is.na(dopa_max) & dopa_max > 5) |
        (!is.na(epi_max) & epi_max > 0 & epi_max <= 0.1) |
        (!is.na(norepi_max) & norepi_max > 0 & norepi_max <= 0.1)] <- 3L
  out[(!is.na(dopa_max) & dopa_max > 15) |
        (!is.na(epi_max) & epi_max > 0.1) |
        (!is.na(norepi_max) & norepi_max > 0.1)] <- 4L
  out
}

score_cns <- function(gcs_min) {
  fcase(
    !is.na(gcs_min) & gcs_min < 6, 4,
    !is.na(gcs_min) & gcs_min < 10, 3,
    !is.na(gcs_min) & gcs_min < 13, 2,
    !is.na(gcs_min) & gcs_min < 15, 1,
    default = 0
  )
}

score_renal <- function(crea_max, uo_24h) {
  out <- rep(0L, length(crea_max))
  out[!is.na(crea_max) & crea_max >= 1.2] <- 1L
  out[!is.na(crea_max) & crea_max >= 2.0] <- 2L
  out[!is.na(crea_max) & crea_max >= 3.5] <- 3L
  out[!is.na(crea_max) & crea_max >= 5.0] <- 4L
  out[!is.na(uo_24h) & uo_24h < 500 & out < 3] <- 3L
  out[!is.na(uo_24h) & uo_24h < 200] <- 4L
  out
}

as_posix_utc <- function(x) as.POSIXct(x, tz = "UTC")

read_chartevents_filtered_duckdb <- function(csv_gz,
                                              target_stays,
                                              itemids,
                                              cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file)) {
    cat("[", OUT_PREFIX, "] Reading cached filtered chartevents: ", cache_file, "\n", sep = "")
    return(as.data.table(readRDS(cache_file)))
  }

  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("duckdb", quietly = TRUE)) {
    stop(
      "Packages 'DBI' and 'duckdb' are required to stream-filter chartevents.csv.gz. ",
      "Install them with install.packages(c('DBI', 'duckdb')).",
      call. = FALSE
    )
  }

  csv_gz <- normalizePath(csv_gz, winslash = "/", mustWork = TRUE)
  sql_path <- gsub("'", "''", csv_gz, fixed = TRUE)

  duckdb_tmp <- file.path(WORK_DIR, paste0(OUT_PREFIX, "_duckdb_tmp"))
  dir.create(duckdb_tmp, recursive = TRUE, showWarnings = FALSE)
  sql_tmp <- gsub("'", "''", normalizePath(duckdb_tmp, winslash = "/", mustWork = TRUE), fixed = TRUE)

  memory_limit <- Sys.getenv("DUCKDB_MEMORY_LIMIT", unset = "8GB")
  n_threads <- suppressWarnings(as.integer(Sys.getenv("DUCKDB_THREADS", unset = "4")))
  if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 4L

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit({
    try(duckdb::duckdb_unregister(con, "target_stays_filter"), silent = TRUE)
    try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }, add = TRUE)

  DBI::dbExecute(con, sprintf("SET memory_limit = '%s'", memory_limit))
  DBI::dbExecute(con, sprintf("SET threads = %d", n_threads))
  DBI::dbExecute(con, "SET preserve_insertion_order = false")
  DBI::dbExecute(con, sprintf("SET temp_directory = '%s'", sql_tmp))

  stay_filter <- data.frame(stay_id = as.numeric(target_stays))
  duckdb::duckdb_register(con, "target_stays_filter", stay_filter)

  item_sql <- paste(as.integer(itemids), collapse = ",")
  query <- sprintf(
    paste0(
      "SELECT ",
      "TRY_CAST(c.subject_id AS DOUBLE) AS subject_id, ",
      "TRY_CAST(c.hadm_id AS DOUBLE) AS hadm_id, ",
      "TRY_CAST(c.stay_id AS DOUBLE) AS stay_id, ",
      "c.charttime AS charttime, ",
      "TRY_CAST(c.itemid AS INTEGER) AS itemid, ",
      "TRY_CAST(c.valuenum AS DOUBLE) AS valuenum, ",
      "c.value AS value ",
      "FROM read_csv_auto('%s', header = true, all_varchar = true, compression = 'gzip') c ",
      "INNER JOIN target_stays_filter t ",
      "ON TRY_CAST(c.stay_id AS DOUBLE) = t.stay_id ",
      "WHERE TRY_CAST(c.itemid AS INTEGER) IN (%s)"
    ),
    sql_path,
    item_sql
  )

  cat(
    "[", OUT_PREFIX, "] Streaming chartevents with DuckDB; ",
    "memory_limit=", memory_limit, ", threads=", n_threads, "\n",
    sep = ""
  )

  out <- as.data.table(DBI::dbGetQuery(con, query))
  out <- out[!is.na(stay_id) & !is.na(itemid)]

  if (!is.null(cache_file)) {
    saveRDS(out, cache_file, compress = FALSE)
    cat("[", OUT_PREFIX, "] Cached filtered chartevents: ", cache_file, "\n", sep = "")
  }

  out
}

read_labevents_filtered_duckdb <- function(csv_gz,
                                            target_pairs,
                                            itemids,
                                            cache_file = NULL) {
  if (!is.null(cache_file) && file.exists(cache_file)) {
    cat("[", OUT_PREFIX, "] Reading cached filtered labevents: ", cache_file, "\n", sep = "")
    out <- as.data.table(readRDS(cache_file))
    required_cols <- c("subject_id", "hadm_id", "charttime", "itemid", "valuenum")
    if (all(required_cols %in% names(out))) return(out)
    warning("Ignoring invalid labevents cache: ", cache_file, call. = FALSE)
  }

  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("duckdb", quietly = TRUE)) {
    stop(
      "Packages 'DBI' and 'duckdb' are required to stream-filter labevents.csv.gz. ",
      "Install them with install.packages(c('DBI', 'duckdb')).",
      call. = FALSE
    )
  }

  csv_gz <- normalizePath(csv_gz, winslash = "/", mustWork = TRUE)
  sql_path <- gsub("'", "''", csv_gz, fixed = TRUE)

  duckdb_tmp_root <- Sys.getenv("DUCKDB_TEMP_DIR", unset = tempdir())
  duckdb_tmp <- file.path(duckdb_tmp_root, paste0(OUT_PREFIX, "_duckdb_labevents_tmp"))
  dir.create(duckdb_tmp, recursive = TRUE, showWarnings = FALSE)
  sql_tmp <- gsub("'", "''", normalizePath(duckdb_tmp, winslash = "/", mustWork = TRUE), fixed = TRUE)

  memory_limit <- Sys.getenv("DUCKDB_MEMORY_LIMIT", unset = "4GB")
  n_threads <- suppressWarnings(as.integer(Sys.getenv("DUCKDB_THREADS", unset = "4")))
  if (!is.finite(n_threads) || n_threads < 1L) n_threads <- 4L

  target_pairs <- unique(as.data.frame(target_pairs[, c("subject_id", "hadm_id"), drop = FALSE]))
  target_pairs$subject_id <- suppressWarnings(as.numeric(target_pairs$subject_id))
  target_pairs$hadm_id <- suppressWarnings(as.numeric(target_pairs$hadm_id))
  target_pairs <- target_pairs[is.finite(target_pairs$subject_id) & is.finite(target_pairs$hadm_id), , drop = FALSE]
  if (nrow(target_pairs) == 0L) stop("No valid subject_id/hadm_id pairs were available for labevents filtering.")

  itemids <- unique(suppressWarnings(as.integer(itemids)))
  itemids <- itemids[is.finite(itemids)]
  if (length(itemids) == 0L) stop("No valid labevents itemids were supplied.")

  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit({
    try(duckdb::duckdb_unregister(con, "target_lab_pairs"), silent = TRUE)
    try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
  }, add = TRUE)

  DBI::dbExecute(con, sprintf("SET memory_limit = '%s'", memory_limit))
  DBI::dbExecute(con, sprintf("SET threads = %d", n_threads))
  DBI::dbExecute(con, "SET preserve_insertion_order = false")
  DBI::dbExecute(con, sprintf("SET temp_directory = '%s'", sql_tmp))

  duckdb::duckdb_register(con, "target_lab_pairs", target_pairs)
  item_sql <- paste(itemids, collapse = ",")

  query <- sprintf(
    paste0(
      "SELECT ",
      "TRY_CAST(c.subject_id AS DOUBLE) AS subject_id, ",
      "TRY_CAST(c.hadm_id AS DOUBLE) AS hadm_id, ",
      "c.charttime AS charttime, ",
      "TRY_CAST(c.itemid AS INTEGER) AS itemid, ",
      "TRY_CAST(c.valuenum AS DOUBLE) AS valuenum ",
      "FROM read_csv_auto('%s', header = true, all_varchar = true, compression = 'gzip') c ",
      "INNER JOIN target_lab_pairs t ",
      "ON TRY_CAST(c.subject_id AS DOUBLE) = t.subject_id ",
      "AND TRY_CAST(c.hadm_id AS DOUBLE) = t.hadm_id ",
      "WHERE TRY_CAST(c.itemid AS INTEGER) IN (%s)"
    ),
    sql_path,
    item_sql
  )

  cat(
    "[", OUT_PREFIX, "] Streaming labevents with DuckDB; ",
    "target admissions=", nrow(target_pairs),
    ", memory_limit=", memory_limit,
    ", threads=", n_threads, "\n",
    sep = ""
  )

  out <- as.data.table(DBI::dbGetQuery(con, query))
  out <- out[
    !is.na(subject_id) & !is.na(hadm_id) &
      !is.na(charttime) & !is.na(itemid) & !is.na(valuenum)
  ]

  if (!is.null(cache_file)) {
    saveRDS(out, cache_file, compress = FALSE)
    cat("[", OUT_PREFIX, "] Cached filtered labevents: ", cache_file, "\n", sep = "")
  }

  out
}

cat("[", OUT_PREFIX, "] Reading HMM row source: ", HMM_ROW_SOURCE, "\n", sep = "")
hmm_row_source_dt <- fread(file.path(WORK_DIR, HMM_ROW_SOURCE))
if (!"stay_id" %in% names(hmm_row_source_dt)) stop("HMM row source has no stay_id")
if (!"charttime_year_day" %in% names(hmm_row_source_dt)) stop("HMM row source has no charttime_year_day")
if (!"hadm_id" %in% names(hmm_row_source_dt)) stop("HMM row source has no hadm_id")

hmm_row_source_dt[, stay_id := as.character(stay_id)]
hmm_row_source_dt[, hadm_id := as.numeric(hadm_id)]
hmm_row_source_dt[, charttime_year_day := as.IDate(charttime_year_day)]

count_dt <- hmm_row_source_dt[, .N, by = stay_id]
morethan2 <- count_dt[N >= 2, stay_id]

cohort <- fread(COHORT_FILE)
cohort[, stay_id := as.character(stay_id)]
cohort[, hadm_id := as.numeric(hadm_id)]
cohort[, subject_id := as.numeric(subject_id)]
cohort[, `:=`(
  intime  = as_posix_utc(intime),
  outtime = as_posix_utc(outtime)
)]
cohort <- cohort[age >= 18 & age < 80]

row_key <- hmm_row_source_dt[stay_id %in% morethan2 & stay_id %in% cohort$stay_id]
row_key <- merge(
  row_key[, .(stay_id, hadm_id, charttime_year_day,
              dialysis_indicator_input = if ("Dialysis" %in% names(row_key)) Dialysis else NA_integer_,
              iabp_indicator_input     = if ("IABP" %in% names(row_key)) IABP else NA_integer_,
              ecmo_indicator_input     = if ("ECMO" %in% names(row_key)) ECMO else NA_integer_)],
  cohort[, .(stay_id, subject_id, hadm_id, intime, outtime, age, sex)],
  by = c("stay_id", "hadm_id"),
  all.x = TRUE,
  sort = FALSE
)

if (any(is.na(row_key$subject_id))) {
  warning("Some HMM rows did not match cohort metadata; check stay_id/hadm_id alignment.", call. = FALSE)
}

row_key[, `:=`(
  day_start_raw = as.POSIXct(charttime_year_day, tz = "UTC"),
  day_end_raw   = as.POSIXct(charttime_year_day + 1L, tz = "UTC")
)]
row_key[, `:=`(
  day_start = pmax(day_start_raw, intime, na.rm = TRUE),
  day_end   = pmin(day_end_raw, outtime, na.rm = TRUE)
)]
row_key <- row_key[!is.na(day_start) & !is.na(day_end) & day_end > day_start]

setorder(row_key, stay_id, charttime_year_day)
row_key[, lab_time := as.numeric(charttime_year_day - min(charttime_year_day)), by = stay_id]
row_key[, lab_time := lab_time + 1]
row_key[, hmm_row_index := .I]
setcolorder(row_key, c("hmm_row_index", "stay_id", "subject_id", "hadm_id", "charttime_year_day", "lab_time",
                       "day_start", "day_end", "intime", "outtime", "age", "sex",
                       "dialysis_indicator_input", "iabp_indicator_input", "ecmo_indicator_input"))

fwrite(row_key, "patient_day_row_key.csv")
cat("[", OUT_PREFIX, "] Row key rows: ", nrow(row_key), "; stays: ", uniqueN(row_key$stay_id), "\n", sep = "")

key_days <- row_key[, .(hmm_row_index, stay_id, subject_id, hadm_id, charttime_year_day, day_start, day_end)]
setkey(key_days, stay_id, hmm_row_index)

target_stays <- unique(row_key$stay_id)
target_hadms <- unique(row_key$hadm_id)
target_subjects <- unique(row_key$subject_id)

cat("[", OUT_PREFIX, "] Loading raw MIMIC tables from: ", MIMIC_DIR, "\n", sep = "")

labevents_cache <- file.path(WORK_DIR, paste0(OUT_PREFIX, "_labevents_filtered.rds"))
labevents <- read_labevents_filtered_duckdb(
  csv_gz = file.path(MIMIC_DIR, "hosp/labevents.csv.gz"),
  target_pairs = unique(row_key[, .(subject_id, hadm_id)]),
  itemids = LAB_ALL,
  cache_file = labevents_cache
)
labevents[, charttime := as_posix_utc(charttime)]
labevents[, charttime_year_day := as.IDate(charttime)]

chartevents_cache <- file.path(WORK_DIR, paste0(OUT_PREFIX, "_chartevents_filtered.rds"))
chartevents <- read_chartevents_filtered_duckdb(
  csv_gz = file.path(MIMIC_DIR, "icu/chartevents.csv.gz"),
  target_stays = target_stays,
  itemids = CH_ALL,
  cache_file = chartevents_cache
)
chartevents[, stay_id := as.character(stay_id)]
chartevents[, charttime := as_posix_utc(charttime)]
chartevents[, charttime_year_day := as.IDate(charttime)]

inputevents <- fread(
  file.path(MIMIC_DIR, "icu/inputevents.csv.gz"),
  select = c("subject_id", "hadm_id", "stay_id", "starttime", "endtime", "itemid", "rate", "rateuom")
)
inputevents[, stay_id := as.character(stay_id)]
inputevents <- inputevents[stay_id %in% target_stays & itemid %in% VASO_ALL]
inputevents[, `:=`(starttime = as_posix_utc(starttime), endtime = as_posix_utc(endtime))]
inputevents[is.na(endtime), endtime := starttime]

procedureevents <- fread(
  file.path(MIMIC_DIR, "icu/procedureevents.csv.gz"),
  select = c("subject_id", "hadm_id", "stay_id", "starttime", "endtime", "itemid")
)
procedureevents[, stay_id := as.character(stay_id)]
procedureevents <- procedureevents[stay_id %in% target_stays & itemid %in% unique(c(VENT_PROC_ALL, RRT_ITEMIDS, ECMO_ITEMIDS))]
procedureevents[, `:=`(starttime = as_posix_utc(starttime), endtime = as_posix_utc(endtime))]
procedureevents[is.na(endtime), endtime := starttime]

outputevents <- fread(
  file.path(MIMIC_DIR, "icu/outputevents.csv.gz"),
  select = c("subject_id", "hadm_id", "stay_id", "charttime", "itemid", "value")
)
outputevents[, stay_id := as.character(stay_id)]
outputevents <- outputevents[stay_id %in% target_stays & itemid %in% unique(c(URINE_IDS, RRT_ITEMIDS))]
outputevents[, charttime := as_posix_utc(charttime)]
outputevents[, charttime_year_day := as.IDate(charttime)]

cat("[", OUT_PREFIX, "] Aggregating daily SOFA components...\n", sep = "")

lab_day <- merge(
  labevents,
  key_days[, .(hmm_row_index, stay_id, subject_id, hadm_id, charttime_year_day, day_start, day_end)],
  by = c("subject_id", "hadm_id", "charttime_year_day"),
  allow.cartesian = TRUE
)
lab_day <- lab_day[charttime >= day_start & charttime < day_end]

sofa_labs <- lab_day[, .(
  plt_min  = suppressWarnings(min(valuenum[itemid %in% LAB_ITEMIDS$platelets], na.rm = TRUE)),
  bili_max = suppressWarnings(max(valuenum[itemid %in% LAB_ITEMIDS$bilirubin], na.rm = TRUE)),
  crea_max = suppressWarnings(max(valuenum[itemid %in% LAB_ITEMIDS$creatinine], na.rm = TRUE))
), by = hmm_row_index]
for (cc in c("plt_min", "bili_max", "crea_max")) sofa_labs[is.infinite(get(cc)), (cc) := NA_real_]

chart_day <- merge(
  chartevents,
  key_days[, .(hmm_row_index, stay_id, charttime_year_day, day_start, day_end)],
  by = c("stay_id", "charttime_year_day"),
  allow.cartesian = TRUE
)
chart_day <- chart_day[charttime >= day_start & charttime < day_end]

sofa_map <- chart_day[itemid %in% CH_ITEMIDS$map,
                      .(map_min = suppressWarnings(min(valuenum, na.rm = TRUE))),
                      by = hmm_row_index]
sofa_map[is.infinite(map_min), map_min := NA_real_]

pao2_dt <- lab_day[itemid %in% LAB_ITEMIDS$pao2 & !is.na(valuenum),
                   .(hmm_row_index, stay_id, t_lab = charttime, pao2 = suppressWarnings(as.numeric(valuenum)))]
pao2_dt <- pao2_dt[!is.na(pao2)]
pao2_dt <- pao2_dt[, .(pao2 = min(pao2, na.rm = TRUE)), by = .(hmm_row_index, stay_id, t_lab)]
pao2_dt[is.infinite(pao2), pao2 := NA_real_]
pao2_dt <- pao2_dt[!is.na(pao2)]

fio2_dt <- chart_day[itemid %in% CH_ITEMIDS$fio2_chart,
                     .(stay_id, t_chart = charttime, fio2 = clean_fio2_fraction(valuenum))]
fio2_dt <- fio2_dt[!is.na(fio2)]
fio2_dt <- fio2_dt[, .(fio2 = max(fio2, na.rm = TRUE)), by = .(stay_id, t_chart)]
fio2_dt[is.infinite(fio2), fio2 := NA_real_]
fio2_dt <- fio2_dt[!is.na(fio2)]

if (nrow(pao2_dt) > 0 && nrow(fio2_dt) > 0) {
  setkey(fio2_dt, stay_id, t_chart)
  setkey(pao2_dt, stay_id, t_lab)
  pf_merged <- fio2_dt[pao2_dt, on = .(stay_id, t_chart = t_lab), roll = "nearest", mult = "first",
                        .(hmm_row_index = i.hmm_row_index,
                          stay_id = i.stay_id,
                          t_lab = i.t_lab,
                          t_fio2 = x.t_chart,
                          pao2 = i.pao2,
                          fio2 = x.fio2)]
  pf_merged <- pf_merged[!is.na(t_fio2) & abs(as.numeric(difftime(t_lab, t_fio2, units = "hours"))) <= 2]
  pf_merged[, pf_ratio := pao2 / fio2]
  sofa_resp <- pf_merged[, .(pfr_min = suppressWarnings(min(pf_ratio, na.rm = TRUE))), by = hmm_row_index]
  sofa_resp[is.infinite(pfr_min), pfr_min := NA_real_]
} else {
  sofa_resp <- data.table(hmm_row_index = integer(), pfr_min = numeric())
}

gcs_raw <- chart_day[itemid %in% c(CH_ITEMIDS$gcs_motor, CH_ITEMIDS$gcs_verbal, CH_ITEMIDS$gcs_eyes),
                     .(hmm_row_index, stay_id, charttime, itemid, valuenum, value)]

gcs_m <- gcs_raw[itemid %in% CH_ITEMIDS$gcs_motor,
                 .(stay_id, hmm_row_index, charttime, gcs_m = suppressWarnings(as.numeric(valuenum)))]
gcs_m <- gcs_m[!is.na(gcs_m)]
gcs_m <- gcs_m[, .(gcs_m = min(gcs_m, na.rm = TRUE)), by = .(stay_id, hmm_row_index, charttime)]
gcs_m[is.infinite(gcs_m), gcs_m := NA_real_]
gcs_m <- gcs_m[!is.na(gcs_m)]

gcs_v <- gcs_raw[itemid %in% CH_ITEMIDS$gcs_verbal,
                 .(stay_id, charttime,
                   gcs_v = suppressWarnings(as.numeric(valuenum)),
                   v_txt = tolower(trimws(fifelse(is.na(value), "", as.character(value)))))]
gcs_v[is.na(gcs_v) & nzchar(v_txt) & grepl("ett|intubat|tube", v_txt), gcs_v := 1]
gcs_v[, v_txt := NULL]
gcs_v <- gcs_v[!is.na(gcs_v)]
gcs_v <- gcs_v[, .(gcs_v = min(gcs_v, na.rm = TRUE)), by = .(stay_id, charttime)]
gcs_v[is.infinite(gcs_v), gcs_v := NA_real_]
gcs_v <- gcs_v[!is.na(gcs_v)]

gcs_e <- gcs_raw[itemid %in% CH_ITEMIDS$gcs_eyes,
                 .(stay_id, charttime, gcs_e = suppressWarnings(as.numeric(valuenum)))]
gcs_e <- gcs_e[!is.na(gcs_e)]
gcs_e <- gcs_e[, .(gcs_e = min(gcs_e, na.rm = TRUE)), by = .(stay_id, charttime)]
gcs_e[is.infinite(gcs_e), gcs_e := NA_real_]
gcs_e <- gcs_e[!is.na(gcs_e)]

if (nrow(gcs_m) > 0 && nrow(gcs_v) > 0 && nrow(gcs_e) > 0) {
  setkey(gcs_m, stay_id, charttime)
  setkey(gcs_v, stay_id, charttime)
  setkey(gcs_e, stay_id, charttime)
  gcs_mv <- gcs_v[gcs_m, on = .(stay_id, charttime), roll = "nearest", mult = "first",
                  .(stay_id = i.stay_id,
                    hmm_row_index = i.hmm_row_index,
                    t_m = i.charttime,
                    t_v = x.charttime,
                    m = i.gcs_m,
                    v = x.gcs_v)]
  gcs_mv <- gcs_mv[!is.na(t_v) & abs(as.numeric(difftime(t_m, t_v, units = "hours"))) <= 1]
  gcs_mve <- gcs_e[gcs_mv, on = .(stay_id, charttime = t_m), roll = "nearest", mult = "first",
                   .(stay_id = i.stay_id,
                     hmm_row_index = i.hmm_row_index,
                     t_m = i.t_m,
                     t_v = i.t_v,
                     t_e = x.charttime,
                     m = i.m,
                     v = i.v,
                     e = x.gcs_e)]
  gcs_mve <- gcs_mve[!is.na(t_e) & abs(as.numeric(difftime(t_m, t_e, units = "hours"))) <= 1]
  gcs_mve[, gcs_sum := m + v + e]
  sofa_cns <- gcs_mve[, .(gcs_min = suppressWarnings(min(gcs_sum, na.rm = TRUE))), by = hmm_row_index]
  sofa_cns[is.infinite(gcs_min), gcs_min := NA_real_]
} else {
  sofa_cns <- data.table(hmm_row_index = integer(), gcs_min = numeric())
}

cat("[", OUT_PREFIX, "] Aggregating interval-based organ support...\n", sep = "")

interval_days <- key_days[, .(hmm_row_index, stay_id, day_start, day_end)]

proc_day <- merge(
  procedureevents,
  interval_days,
  by = "stay_id",
  allow.cartesian = TRUE
)
proc_day <- proc_day[starttime < day_end & endtime >= day_start]

vent_day <- proc_day[itemid %in% VENT_PROC_ALL,
                     .(ventilation_active = any(itemid %in% VENT_PROC_PRIMARY),
                       intubation_event_day = any(itemid %in% INTUB_PROC_EVENT),
                       vent_or_intub_proc_any = TRUE),
                     by = hmm_row_index]

rrt_proc_day <- proc_day[itemid %in% RRT_ITEMIDS,
                         .(rrt_proc_any = TRUE,
                           crrt_proc_any = any(itemid %in% RRT_CRRT_ITEMIDS),
                           hemodialysis_proc_any = any(itemid %in% RRT_HD_ITEMIDS)),
                         by = hmm_row_index]

ecmo_proc_day <- proc_day[itemid %in% ECMO_ITEMIDS,
                          .(ecmo_proc_any = TRUE),
                          by = hmm_row_index]

input_day <- merge(inputevents, interval_days, by = "stay_id", allow.cartesian = TRUE)
input_day <- input_day[starttime < day_end & endtime >= day_start]

rateuom_summary <- input_day[, .N, by = .(itemid, rateuom)][order(itemid, -N)]
fwrite(rateuom_summary, "patient_day_vasopressor_rate_unit_summary.csv")

sofa_vaso <- input_day[, .(
  dopa_max   = suppressWarnings(max(rate[itemid %in% VASO_ITEMIDS$dopa], na.rm = TRUE)),
  epi_max    = suppressWarnings(max(rate[itemid %in% VASO_ITEMIDS$epi], na.rm = TRUE)),
  norepi_max = suppressWarnings(max(rate[itemid %in% VASO_ITEMIDS$norepi], na.rm = TRUE)),
  dobu_max   = suppressWarnings(max(rate[itemid %in% VASO_ITEMIDS$dobu], na.rm = TRUE)),
  phenyl_max = suppressWarnings(max(rate[itemid %in% VASO_ITEMIDS$phenyl], na.rm = TRUE)),
  vaso_max   = suppressWarnings(max(rate[itemid %in% VASO_ITEMIDS$vaso], na.rm = TRUE)),
  angii_max  = suppressWarnings(max(rate[itemid %in% VASO_ITEMIDS$angii], na.rm = TRUE))
), by = hmm_row_index]
for (cc in setdiff(names(sofa_vaso), "hmm_row_index")) sofa_vaso[is.infinite(get(cc)), (cc) := NA_real_]

out_day <- merge(
  outputevents,
  key_days[, .(hmm_row_index, stay_id, charttime_year_day, day_start, day_end)],
  by = c("stay_id", "charttime_year_day"),
  allow.cartesian = TRUE
)
out_day <- out_day[charttime >= day_start & charttime < day_end]
out_day[, value_num := suppressWarnings(as.numeric(value))]

uo_day <- out_day[itemid %in% URINE_IDS & !is.na(value_num)]
uo_day[, val_signed := fifelse(itemid %in% IRRIGANT_IDS & value_num > 0, -value_num, value_num)]
uo_agg <- uo_day[, .(uo_24h = sum(val_signed, na.rm = TRUE)), by = hmm_row_index]
uo_agg[uo_24h < 0, uo_24h := 0]

rrt_out_day <- out_day[itemid %in% RRT_ITEMIDS & !is.na(value_num) & value_num != 0,
                       .(rrt_output_any = TRUE,
                         ultrafiltrate_output_ml = sum(value_num[itemid == 226457], na.rm = TRUE),
                         hemodialysis_output_ml = sum(value_num[itemid == 226499], na.rm = TRUE)),
                       by = hmm_row_index]

peep_day <- chart_day[itemid %in% c(CH_ITEMIDS$peep_set, CH_ITEMIDS$peep_backup),
                      .(peep_max = suppressWarnings(max(as.numeric(valuenum), na.rm = TRUE)),
                        peep_median = suppressWarnings(median(as.numeric(valuenum), na.rm = TRUE)),
                        peep_set_max = suppressWarnings(max(as.numeric(valuenum[itemid %in% CH_ITEMIDS$peep_set]), na.rm = TRUE)),
                        peep_backup_max = suppressWarnings(max(as.numeric(valuenum[itemid %in% CH_ITEMIDS$peep_backup]), na.rm = TRUE))),
                      by = hmm_row_index]
for (cc in setdiff(names(peep_day), "hmm_row_index")) peep_day[is.infinite(get(cc)), (cc) := NA_real_]

fio2_day <- chart_day[itemid %in% CH_ITEMIDS$fio2_chart,
                      .(fio2_max = suppressWarnings(max(clean_fio2_fraction(valuenum), na.rm = TRUE)),
                        fio2_median = suppressWarnings(median(clean_fio2_fraction(valuenum), na.rm = TRUE))),
                      by = hmm_row_index]
for (cc in setdiff(names(fio2_day), "hmm_row_index")) fio2_day[is.infinite(get(cc)), (cc) := NA_real_]

vent_mode_day <- chart_day[itemid %in% c(CH_ITEMIDS$vent_mode, CH_ITEMIDS$vent_type, CH_ITEMIDS$o2_device),
                           .(vent_mode_chart_any = any(itemid %in% CH_ITEMIDS$vent_mode & !is.na(value) & nzchar(trimws(as.character(value)))),
                             vent_type_chart_any = any(itemid %in% CH_ITEMIDS$vent_type & !is.na(value) & nzchar(trimws(as.character(value)))),
                             o2_device_chart_any = any(itemid %in% CH_ITEMIDS$o2_device & !is.na(value) & nzchar(trimws(as.character(value))))),
                           by = hmm_row_index]

rrt_chart_day <- chart_day[itemid %in% RRT_ITEMIDS,
                           .(rrt_chart_any = any(!is.na(value) | !is.na(valuenum)),
                             crrt_chart_any = any(itemid %in% RRT_CRRT_ITEMIDS & (!is.na(value) | !is.na(valuenum))),
                             hemodialysis_chart_any = any(itemid %in% RRT_HD_ITEMIDS & (!is.na(value) | !is.na(valuenum)))),
                           by = hmm_row_index]

ecmo_chart_day <- chart_day[itemid %in% ECMO_ITEMIDS,
                            .(ecmo_chart_any = any(!is.na(value) | !is.na(valuenum)),
                              ecmo_flow_max = suppressWarnings(max(as.numeric(valuenum[itemid %in% c(224660, 229270, 229842)]), na.rm = TRUE))),
                            by = hmm_row_index]
ecmo_chart_day[is.infinite(ecmo_flow_max), ecmo_flow_max := NA_real_]

cat("[", OUT_PREFIX, "] Combining outputs and calculating scores...\n", sep = "")

dt <- row_key[, .(hmm_row_index, stay_id, subject_id, hadm_id, charttime_year_day, lab_time, age, sex,
                  day_start, day_end, dialysis_indicator_input, iabp_indicator_input, ecmo_indicator_input)]
for (x in list(sofa_labs, sofa_resp, sofa_map, sofa_cns, sofa_vaso, uo_agg, vent_day,
               rrt_proc_day, rrt_out_day, rrt_chart_day, ecmo_proc_day, ecmo_chart_day,
               peep_day, fio2_day, vent_mode_day)) {
  if (nrow(x) > 0) dt <- merge(dt, x, by = "hmm_row_index", all.x = TRUE, sort = FALSE)
}

logical_cols <- c("ventilation_active", "intubation_event_day", "vent_or_intub_proc_any",
                  "rrt_proc_any", "crrt_proc_any", "hemodialysis_proc_any",
                  "rrt_output_any", "rrt_chart_any", "crrt_chart_any", "hemodialysis_chart_any",
                  "ecmo_proc_any", "ecmo_chart_any",
                  "vent_mode_chart_any", "vent_type_chart_any", "o2_device_chart_any")
for (cc in intersect(logical_cols, names(dt))) dt[is.na(get(cc)), (cc) := FALSE]

expected_numeric <- c("pfr_min", "plt_min", "bili_max", "crea_max", "map_min", "gcs_min",
                      "dopa_max", "epi_max", "norepi_max", "dobu_max", "phenyl_max", "vaso_max", "angii_max",
                      "uo_24h", "peep_max", "peep_median", "peep_set_max", "peep_backup_max",
                      "fio2_max", "fio2_median", "ultrafiltrate_output_ml", "hemodialysis_output_ml", "ecmo_flow_max")
for (cc in setdiff(expected_numeric, names(dt))) dt[, (cc) := NA_real_]
for (cc in setdiff(logical_cols, names(dt))) dt[, (cc) := FALSE]

dt[, `:=`(
  vasopressor_sofa_any = (!is.na(dopa_max) & dopa_max > 0) |
    (!is.na(dobu_max) & dobu_max > 0) |
    (!is.na(epi_max) & epi_max > 0) |
    (!is.na(norepi_max) & norepi_max > 0),
  vasopressor_any = (!is.na(dopa_max) & dopa_max > 0) |
    (!is.na(dobu_max) & dobu_max > 0) |
    (!is.na(epi_max) & epi_max > 0) |
    (!is.na(norepi_max) & norepi_max > 0) |
    (!is.na(phenyl_max) & phenyl_max > 0) |
    (!is.na(vaso_max) & vaso_max > 0) |
    (!is.na(angii_max) & angii_max > 0),
  rrt_any = rrt_proc_any | rrt_output_any | rrt_chart_any | (dialysis_indicator_input %in% 1),
  crrt_any = crrt_proc_any | crrt_chart_any,
  hemodialysis_any = hemodialysis_proc_any | hemodialysis_chart_any | (!is.na(hemodialysis_output_ml) & hemodialysis_output_ml > 0),
  ecmo_any = ecmo_proc_any | ecmo_chart_any | (ecmo_indicator_input %in% 1)
)]

dt[, `:=`(
  resp_available  = !is.na(pfr_min),
  coag_available  = !is.na(plt_min),
  liver_available = !is.na(bili_max),
  cv_available    = !is.na(map_min) | vasopressor_sofa_any,
  cns_available   = !is.na(gcs_min),
  renal_available = !is.na(crea_max) | !is.na(uo_24h)
)]
dt[, sofa_component_available_n := rowSums(.SD), .SDcols = c("resp_available", "coag_available", "liver_available", "cv_available", "cns_available", "renal_available")]

dt[, `:=`(
  sofa_resp  = score_resp(pfr_min, ventilation_active),
  sofa_coag  = score_coag(plt_min),
  sofa_liver = score_liver(bili_max),
  sofa_cv    = score_cv(map_min, dopa_max, epi_max, norepi_max, dobu_max),
  sofa_cns   = score_cns(gcs_min),
  sofa_renal = score_renal(crea_max, uo_24h)
)]
dt[, sofa_total_available_component := rowSums(.SD, na.rm = TRUE),
   .SDcols = c("sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns", "sofa_renal")]
dt[, sofa_total_complete6 := fifelse(sofa_component_available_n == 6L, sofa_total_available_component, NA_real_)]

sofa_cols <- c("hmm_row_index", "stay_id", "subject_id", "hadm_id", "charttime_year_day", "lab_time",
               "sofa_resp", "sofa_coag", "sofa_liver", "sofa_cv", "sofa_cns", "sofa_renal",
               "sofa_total_available_component", "sofa_total_complete6", "sofa_component_available_n",
               "resp_available", "coag_available", "liver_available", "cv_available", "cns_available", "renal_available",
               "pfr_min", "plt_min", "bili_max", "crea_max", "map_min", "gcs_min",
               "dopa_max", "epi_max", "norepi_max", "dobu_max", "uo_24h", "ventilation_active")

support_cols <- c("hmm_row_index", "stay_id", "subject_id", "hadm_id", "charttime_year_day", "lab_time",
                  "ventilation_active", "intubation_event_day", "vent_or_intub_proc_any",
                  "peep_max", "peep_median", "peep_set_max", "peep_backup_max", "fio2_max", "fio2_median",
                  "vasopressor_sofa_any", "vasopressor_any", "dopa_max", "dobu_max", "epi_max", "norepi_max",
                  "phenyl_max", "vaso_max", "angii_max",
                  "rrt_any", "crrt_any", "hemodialysis_any", "rrt_proc_any", "rrt_output_any", "rrt_chart_any",
                  "ultrafiltrate_output_ml", "hemodialysis_output_ml",
                  "ecmo_any", "ecmo_proc_any", "ecmo_chart_any", "ecmo_flow_max",
                  "uo_24h", "dialysis_indicator_input", "iabp_indicator_input", "ecmo_indicator_input",
                  "vent_mode_chart_any", "vent_type_chart_any", "o2_device_chart_any")

fwrite(dt[, ..sofa_cols], "patient_day_sofa_components.csv")
fwrite(dt[, ..support_cols], "patient_day_organsupport.csv")
fwrite(dt, "patient_day_sofa_organsupport_joined.csv")

qc_overall <- data.table(
  metric = c(
    "hmm_patientday_rows",
    "unique_stays",
    "rows_resp_available_pct",
    "rows_coag_available_pct",
    "rows_liver_available_pct",
    "rows_cv_available_pct",
    "rows_cns_available_pct",
    "rows_renal_available_pct",
    "rows_complete6_sofa_pct",
    "rows_ventilation_active_pct",
    "rows_vasopressor_any_pct",
    "rows_rrt_any_pct",
    "rows_ecmo_any_pct"
  ),
  value = c(
    nrow(dt),
    uniqueN(dt$stay_id),
    mean(dt$resp_available) * 100,
    mean(dt$coag_available) * 100,
    mean(dt$liver_available) * 100,
    mean(dt$cv_available) * 100,
    mean(dt$cns_available) * 100,
    mean(dt$renal_available) * 100,
    mean(dt$sofa_component_available_n == 6L) * 100,
    mean(dt$ventilation_active) * 100,
    mean(dt$vasopressor_any) * 100,
    mean(dt$rrt_any) * 100,
    mean(dt$ecmo_any) * 100
  )
)

qc_by_day <- dt[, .(
  rows = .N,
  unique_stays = uniqueN(stay_id),
  resp_available_pct = mean(resp_available) * 100,
  complete6_sofa_pct = mean(sofa_component_available_n == 6L) * 100,
  ventilation_active_pct = mean(ventilation_active) * 100,
  vasopressor_any_pct = mean(vasopressor_any) * 100,
  rrt_any_pct = mean(rrt_any) * 100,
  ecmo_any_pct = mean(ecmo_any) * 100,
  sofa_available_median = median(sofa_total_available_component, na.rm = TRUE),
  sofa_complete6_median = median(sofa_total_complete6, na.rm = TRUE)
), by = lab_time][order(lab_time)]

qc_by_component_count <- dt[, .N, by = sofa_component_available_n][order(sofa_component_available_n)]

qc <- rbindlist(list(
  cbind(section = "overall", qc_overall),
  data.table(section = "by_component_available_n", metric = paste0("n_components_", qc_by_component_count$sofa_component_available_n), value = qc_by_component_count$N)
), fill = TRUE)

fwrite(qc, "patient_day_sofa_organsupport_qc.csv")
fwrite(qc_by_day, "patient_day_sofa_organsupport_qc_by_day.csv")

manifest <- data.table(
  item = c("created_at", "work_dir", "mimic_dir", "hmm_row_source", "cohort_file",
           "row_key_csv", "sofa_csv", "support_csv", "joined_csv", "qc_csv", "qc_by_day_csv",
           "note_resp_support", "note_cv_sofa"),
  value = c(as.character(Sys.time()), WORK_DIR, MIMIC_DIR, HMM_ROW_SOURCE, COHORT_FILE,
            "patient_day_row_key.csv",
            "patient_day_sofa_components.csv",
            "patient_day_organsupport.csv",
            "patient_day_sofa_organsupport_joined.csv",
            "patient_day_sofa_organsupport_qc.csv",
            "patient_day_sofa_organsupport_qc_by_day.csv",
            "Respiratory SOFA support uses invasive-ventilation procedure intervals overlapping each analytical calendar day.",
            "Cardiovascular SOFA uses dopamine, dobutamine, epinephrine, and norepinephrine infusion rates; other vasoactive agents are reported descriptively.")
)
fwrite(manifest, paste0(OUT_PREFIX, "_manifest.csv"))

cat("[", OUT_PREFIX, "] Done.\n", sep = "")
print(qc_overall)
