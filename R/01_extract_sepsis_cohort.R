
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
library(readr)
library(data.table)
library(lubridate)
options(datatable.integer64 = "double")

input_dir <- mimic_root
abx_csv <- resource_path("antibiotics_list.txt")
sepsis_out <- "mimic_iv_sepsis3_cohort_final.csv.gz"

LAB_ITEMIDS <- list(
  pao2       = c(50821),
  platelets  = c(51265, 51266),
  bilirubin  = c(50885, 225690),
  creatinine = c(50912, 52546)
)

CH_ITEMIDS <- list(
  fio2_chart = c(223835),
  map        = c(220052, 220181, 225312),
  gcs_motor  = c(223901),
  gcs_verbal = c(223900),
  gcs_eyes   = c(220739)
)

INTUB_POINT_ITEMIDS <- c(224385)
VENT_DURATION_ITEMIDS <- c(225792)

VASO_ITEMIDS <- list(
  norepi = c(221906),
  epi    = c(221289),
  dopa   = c(221662),
  dobu   = c(221653)
)
VASO_ALL <- sort(unique(unlist(VASO_ITEMIDS)))

URINE_IDS <- c(
  226559, 226560, 226561, 226584, 226563, 226564, 226565, 226567,
  226557, 226558, 227488, 227489
)

IRRIGANT_IDS <- c(227488)

cat("Loading tables...\n")

if (!file.exists(abx_csv)) stop("abx_csv missing")
escape_regex <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x, perl = TRUE)

abx_terms <- tolower(trimws(readLines(abx_csv, warn = FALSE)))
abx_terms <- abx_terms[nzchar(abx_terms)]
if (length(abx_terms) == 0) stop("STOP: antibiotics_list.txt has no usable terms (empty after trimming).", call. = FALSE)

regex_terms <- vapply(
  abx_terms,
  function(x) {
    x2 <- escape_regex(x)
    x2 <- gsub("[[:space:]-]+", "[\\\\s\\\\-]*", x2)
    paste0("\\b", x2, "\\b")
  },
  character(1)
)
abx_pattern <- paste(regex_terms, collapse = "|")

prescriptions  <- fread(file.path(input_dir, "hosp/prescriptions.csv.gz"), select = c("subject_id","hadm_id","starttime","route","drug"))
microbioevents <- fread(file.path(input_dir, "hosp/microbiologyevents.csv.gz"), select = c("subject_id","hadm_id","charttime","chartdate","spec_type_desc"))
admissions     <- fread(file.path(input_dir, "hosp/admissions.csv.gz"), select = c("subject_id","hadm_id","admittime","dischtime"))
icustays       <- fread(file.path(input_dir, "icu/icustays.csv.gz"), select = c("subject_id","hadm_id","stay_id","intime","outtime"))

labevents      <- fread(file.path(input_dir, "hosp/labevents.csv.gz"), select = c("subject_id","hadm_id","charttime","itemid","valuenum"))
chartevents    <- fread(file.path(input_dir, "icu/chartevents.csv.gz"), select = c("subject_id","hadm_id","stay_id","charttime","itemid","valuenum","value"))
inputevents    <- fread(file.path(input_dir, "icu/inputevents.csv.gz"), select = c("subject_id","hadm_id","stay_id","starttime","endtime","itemid","rate","rateuom","patientweight"))
procedureevents<- fread(file.path(input_dir, "icu/procedureevents.csv.gz"), select = c("subject_id","hadm_id","stay_id","starttime","endtime","itemid"))
outputevents   <- fread(file.path(input_dir, "icu/outputevents.csv.gz"), select = c("subject_id","hadm_id","stay_id","charttime","itemid","value"))

cat("Normalizing times...\n")
admissions[, `:=`(admittime = as.POSIXct(admittime, tz="UTC"), dischtime = as.POSIXct(dischtime, tz="UTC"))]
admissions[, dischtime_fix := fifelse(is.na(dischtime), admittime, dischtime)]

icustays[, `:=`(intime = as.POSIXct(intime, tz="UTC"), outtime = as.POSIXct(outtime, tz="UTC"))]
icustays[, outtime_fix := fifelse(is.na(outtime), intime, outtime)]

prescriptions[, starttime := as.POSIXct(starttime, tz="UTC")]
microbioevents[, charttime := as.POSIXct(charttime, tz="UTC")]
labevents[, charttime := as.POSIXct(charttime, tz="UTC")]
chartevents[, charttime := as.POSIXct(charttime, tz="UTC")]
outputevents[, charttime := as.POSIXct(charttime, tz="UTC")]

inputevents[, `:=`(starttime = as.POSIXct(starttime, tz="UTC"), endtime = as.POSIXct(endtime, tz="UTC"))]
inputevents[, endtime_fix := fifelse(is.na(endtime), starttime, endtime)]

procedureevents[, `:=`(starttime = as.POSIXct(starttime, tz="UTC"), endtime = as.POSIXct(endtime, tz="UTC"))]
procedureevents[, endtime_fix := fifelse(is.na(endtime), starttime, endtime)]

cat("Identifying suspected infection...\n")

prescriptions[, route_clean := tolower(trimws(fifelse(is.na(route), "", route)))]
abx_any <- prescriptions[
  !is.na(starttime) &
    (route_clean == "" | grepl("iv|po|im|intravenous|oral", route_clean)) &
    grepl(abx_pattern, tolower(drug), perl = TRUE),
  .(subject_id, hadm_id, abx_time = starttime)
]
abx_any[, abx_end24 := abx_time + hours(24)]

microbioevents[, chartdate_fix := as.POSIXct(as.IDate(chartdate), tz="UTC") + hours(12)]
culture_any <- microbioevents[
  (is.na(spec_type_desc) | !grepl("screen|surveillance|mrsa", tolower(spec_type_desc))),
  .(subject_id, hadm_id, cult_time = fcoalesce(charttime, chartdate_fix))
]

culture_ok <- culture_any[!is.na(hadm_id)]
culture_na <- culture_any[is.na(hadm_id)]

if(nrow(culture_na) > 0) {
  culture_na[, cult_limit := cult_time + hours(24)]
  culture_imputed <- admissions[
    culture_na,
    on = .(subject_id, admittime <= cult_limit, dischtime_fix >= cult_time),
    nomatch = 0,
    .(subject_id, hadm_id, cult_time = i.cult_time, admittime)
  ]
  setorder(culture_imputed, subject_id, cult_time, -admittime)
  culture_imputed <- culture_imputed[, .SD[1], by = .(subject_id, cult_time)][, admittime := NULL]
  culture_any <- rbind(culture_ok, culture_imputed)
} else {
  culture_any <- culture_ok
}
culture_any <- unique(culture_any)
culture_any[, cult_end72 := cult_time + hours(72)]

abx_first <- culture_any[
  abx_any,
  on = .(subject_id, hadm_id, cult_time >= abx_time, cult_time <= abx_end24),
  nomatch = 0,
  allow.cartesian = TRUE,
  .(subject_id, hadm_id, abx_time = i.abx_time, cult_time = x.cult_time)
]
cult_first <- abx_any[
  culture_any,
  on = .(subject_id, hadm_id, abx_time >= cult_time, abx_time <= cult_end72),
  nomatch = 0,
  allow.cartesian = TRUE,
  .(subject_id, hadm_id, abx_time = x.abx_time, cult_time = i.cult_time)
]

sus_inf_all <- unique(rbind(abx_first, cult_first))
sus_inf_all[, T0 := pmin(abx_time, cult_time)]
setorder(sus_inf_all, subject_id, hadm_id, T0)

sus_inf_all[, gap_hours := as.numeric(difftime(T0, shift(T0), units = "hours")),
            by = .(subject_id, hadm_id)]
sus_inf_all[, infection_cluster := cumsum(is.na(gap_hours) | gap_hours > 72),
            by = .(subject_id, hadm_id)]
fwrite(
  sus_inf_all,
  "mimic_iv_suspected_infection_candidates_qc.csv.gz"
)

sus_inf <- sus_inf_all[, .SD[1], by = .(subject_id, hadm_id)]
sus_inf[, `:=`(start_win = T0 - hours(48), end_win = T0 + hours(24))]

cat("Assigning ICU stays...\n")
icustays[, `:=`(match_start = intime - hours(24), match_end = outtime_fix)]

si_stay_cand <- icustays[
  sus_inf,
  on = .(subject_id, hadm_id, match_start <= T0, match_end >= T0),
  nomatch = 0,
  allow.cartesian = TRUE,
  .(
    subject_id, hadm_id, stay_id,
    intime, outtime_fix,
    T0        = i.T0,
    start_win = i.start_win,
    end_win   = i.end_win
  )
]

si_stay_cand[, overlap := pmin(outtime_fix, end_win) - pmax(intime, start_win)]
setorder(si_stay_cand, subject_id, hadm_id, -overlap)
cohort <- si_stay_cand[, .SD[1], by = .(subject_id, hadm_id)]
cohort[, overlap := NULL]

cat(sprintf("Cohort size: %d admissions\n", nrow(cohort)))

target_hadms <- unique(cohort$hadm_id)
target_stays <- unique(cohort$stay_id)

proc_vent <- procedureevents[
  stay_id %in% target_stays & itemid %in% VENT_DURATION_ITEMIDS,
  .(stay_id, starttime, endtime_fix)
]
vent_flags <- proc_vent[
  cohort,
  on = .(stay_id, endtime_fix >= start_win, starttime <= end_win),
  nomatch = 0,
  .(stay_id, is_vent = TRUE)
]
vent_flags <- unique(vent_flags)

lab_dt <- labevents[
  cohort,
  on = .(subject_id, hadm_id, charttime >= start_win, charttime <= end_win),
  nomatch = 0,
  .(hadm_id, itemid, valuenum, charttime)
]
sofa_labs <- lab_dt[, .(
  plt_min  = suppressWarnings(min(valuenum[itemid %in% LAB_ITEMIDS$platelets], na.rm=TRUE)),
  bili_max = suppressWarnings(max(valuenum[itemid %in% LAB_ITEMIDS$bilirubin], na.rm=TRUE)),
  crea_max = suppressWarnings(max(valuenum[itemid %in% LAB_ITEMIDS$creatinine], na.rm=TRUE))
), by = hadm_id]

pao2_dt <- lab_dt[itemid %in% LAB_ITEMIDS$pao2,
                  .(hadm_id, t_lab = charttime, pao2 = valuenum)]

fio2_dt <- chartevents[
  cohort,
  on = .(stay_id, charttime >= start_win, charttime <= end_win),
  nomatch = 0,
  .(hadm_id, t_chart = charttime, itemid, valuenum)
][itemid %in% CH_ITEMIDS$fio2_chart]

fio2_dt[, fio2 := ifelse(valuenum > 1.0, valuenum / 100, valuenum)]
fio2_dt <- fio2_dt[fio2 >= 0.21 & fio2 <= 1.0]

pao2_dt <- pao2_dt[!is.na(pao2)]
pao2_dt <- pao2_dt[, .(pao2 = suppressWarnings(min(pao2, na.rm = TRUE))),
                   by = .(hadm_id, t_lab)]
pao2_dt[is.infinite(pao2), pao2 := NA_real_]
pao2_dt <- pao2_dt[!is.na(pao2)]

fio2_dt <- fio2_dt[!is.na(fio2)]
fio2_dt <- fio2_dt[, .(fio2 = suppressWarnings(max(fio2, na.rm = TRUE))),
                   by = .(hadm_id, t_chart)]
fio2_dt[is.infinite(fio2), fio2 := NA_real_]
fio2_dt <- fio2_dt[!is.na(fio2)]

setkey(pao2_dt, hadm_id, t_lab)
setkey(fio2_dt, hadm_id, t_chart)

pf_merged <- fio2_dt[
  pao2_dt,
  on   = .(hadm_id, t_chart = t_lab),
  roll = "nearest",
  mult = "first",
  .(hadm_id,
    t_lab  = i.t_lab,
    t_fio2 = x.t_chart,
    pao2   = i.pao2,
    fio2   = x.fio2)
]

pf_merged <- pf_merged[
  !is.na(t_fio2) &
    abs(as.numeric(difftime(t_lab, t_fio2, units = "hours"))) <= 2
]

pf_merged[, pf_ratio := pao2 / fio2]
sofa_resp <- pf_merged[, .(pfr_min = suppressWarnings(min(pf_ratio, na.rm = TRUE))), by = hadm_id]
sofa_resp[is.infinite(pfr_min), pfr_min := NA_real_]

map_dt <- chartevents[
  cohort,
  on = .(stay_id, charttime >= start_win, charttime <= end_win),
  nomatch = 0,
  .(stay_id, itemid, valuenum)
][itemid %in% CH_ITEMIDS$map]
sofa_map <- map_dt[, .(map_min = suppressWarnings(min(valuenum, na.rm=TRUE))), by = stay_id]

gcs_raw <- chartevents[
  cohort,
  on = .(stay_id, charttime >= start_win, charttime <= end_win),
  nomatch = 0,
  .(stay_id, itemid, charttime, valuenum, value)
][itemid %in% c(CH_ITEMIDS$gcs_motor, CH_ITEMIDS$gcs_verbal, CH_ITEMIDS$gcs_eyes)]

dedupe_gcs <- function(dt, value_col) {
  empty <- data.table(
    stay_id = numeric(),
    charttime = as.POSIXct(character()),
    value = numeric()
  )
  if (nrow(dt) == 0L) return(empty)
  dt <- dt[!is.na(get(value_col))]
  if (nrow(dt) == 0L) return(empty)
  dt <- dt[, .(value = suppressWarnings(min(get(value_col), na.rm = TRUE))),
           by = .(stay_id, charttime)]
  dt[is.infinite(value), value := NA_real_]
  dt[!is.na(value)]
}

gcs_m <- gcs_raw[itemid %in% CH_ITEMIDS$gcs_motor,
                 .(stay_id, charttime, score = suppressWarnings(as.numeric(valuenum)))]
gcs_e <- gcs_raw[itemid %in% CH_ITEMIDS$gcs_eyes,
                 .(stay_id, charttime, score = suppressWarnings(as.numeric(valuenum)))]

gcs_v0 <- gcs_raw[itemid %in% CH_ITEMIDS$gcs_verbal,
                  .(stay_id, charttime,
                    score = suppressWarnings(as.numeric(valuenum)),
                    value_text = tolower(trimws(fifelse(is.na(value), "", as.character(value)))))]
gcs_v0[, is_ett := nzchar(value_text) & grepl("ett|intubat|endotracheal|tube", value_text, perl = TRUE)]

gcs_v_primary <- gcs_v0[!is_ett, .(stay_id, charttime, score)]

gcs_v_ett1 <- copy(gcs_v0)
gcs_v_ett1[is_ett, score := 1]
gcs_v_ett1 <- gcs_v_ett1[, .(stay_id, charttime, score)]

gcs_m <- dedupe_gcs(gcs_m, "score")
gcs_e <- dedupe_gcs(gcs_e, "score")
gcs_v_primary <- dedupe_gcs(gcs_v_primary, "score")
gcs_v_ett1 <- dedupe_gcs(gcs_v_ett1, "score")
setnames(gcs_m, "value", "gcs_m")
setnames(gcs_e, "value", "gcs_e")
setnames(gcs_v_primary, "value", "gcs_v")
setnames(gcs_v_ett1, "value", "gcs_v")

calculate_gcs_min <- function(gcs_m, gcs_v, gcs_e, output_name) {
  empty <- data.table(stay_id = numeric(), value = numeric())
  setnames(empty, "value", output_name)
  if (nrow(gcs_m) == 0L || nrow(gcs_v) == 0L || nrow(gcs_e) == 0L) return(empty)

  setkey(gcs_m, stay_id, charttime)
  setkey(gcs_v, stay_id, charttime)
  setkey(gcs_e, stay_id, charttime)

  mv <- gcs_v[
    gcs_m,
    on = .(stay_id, charttime),
    roll = "nearest",
    mult = "first",
    .(stay_id, t_m = i.charttime, t_v = x.charttime,
      m = i.gcs_m, v = x.gcs_v)
  ]
  mv <- mv[!is.na(t_v) & abs(as.numeric(difftime(t_m, t_v, units = "hours"))) <= 1]
  if (nrow(mv) == 0L) return(empty)

  mve <- gcs_e[
    mv,
    on = .(stay_id, charttime = t_m),
    roll = "nearest",
    mult = "first",
    .(stay_id, t_m = i.t_m, t_e = x.charttime,
      m = i.m, v = i.v, e = x.gcs_e)
  ]
  mve <- mve[!is.na(t_e) & abs(as.numeric(difftime(t_m, t_e, units = "hours"))) <= 1]
  if (nrow(mve) == 0L) return(empty)

  mve[, gcs_sum := m + v + e]
  ans <- mve[, .(value = suppressWarnings(min(gcs_sum, na.rm = TRUE))), by = stay_id]
  ans[is.infinite(value), value := NA_real_]
  setnames(ans, "value", output_name)
  ans
}

sofa_cns_primary <- calculate_gcs_min(gcs_m, gcs_v_primary, gcs_e, "gcs_min")
sofa_cns_ett1 <- calculate_gcs_min(gcs_m, gcs_v_ett1, gcs_e, "gcs_min_ett1")
sofa_cns <- merge(sofa_cns_primary, sofa_cns_ett1, by = "stay_id", all = TRUE)

fwrite(
  gcs_v0[, .(
    n_verbal_rows = .N,
    n_ett_rows = sum(is_ett, na.rm = TRUE),
    n_numeric_nonett = sum(!is_ett & !is.na(score))
  )],
  "gcs_ett_handling_qc.csv"
)

normalize_vaso_rate <- function(rate, rateuom, patientweight) {
  r <- suppressWarnings(as.numeric(rate))
  w <- suppressWarnings(as.numeric(patientweight))
  u <- tolower(trimws(fifelse(is.na(rateuom), "", as.character(rateuom))))
  u <- gsub("[[:space:]]+", "", u)
  u <- gsub("[µμ]", "u", u)
  u <- gsub("micrograms?", "mcg", u)
  u <- gsub("ug", "mcg", u, fixed = TRUE)
  u <- gsub("minutes?", "min", u)
  u <- gsub("kilograms?", "kg", u)

  out <- rep(NA_real_, length(r))
  unit_class <- rep("unsupported_or_missing", length(r))

  hit <- grepl("^mcg/kg/min$", u)
  out[hit] <- r[hit]
  unit_class[hit] <- "mcg/kg/min"

  hit <- grepl("^mg/kg/min$", u)
  out[hit] <- r[hit] * 1000
  unit_class[hit] <- "mg/kg/min_to_mcg/kg/min"

  hit <- grepl("^mcg/min$", u) & is.finite(w) & w > 0
  out[hit] <- r[hit] / w[hit]
  unit_class[hit] <- "mcg/min_weight_converted"

  hit <- grepl("^mg/min$", u) & is.finite(w) & w > 0
  out[hit] <- r[hit] * 1000 / w[hit]
  unit_class[hit] <- "mg/min_weight_converted"

  out[!is.finite(out) | out < 0] <- NA_real_
  data.table(rate_mcg_kg_min = out, vaso_unit_class = unit_class)
}

input_dt <- inputevents[
  cohort,
  on = .(stay_id, endtime_fix >= start_win, starttime <= end_win),
  nomatch = 0,
  .(stay_id, itemid, rate, rateuom, patientweight)
][itemid %in% VASO_ALL]

if (nrow(input_dt) > 0L) {
  normalized <- normalize_vaso_rate(input_dt$rate, input_dt$rateuom, input_dt$patientweight)
  input_dt[, `:=`(
    rate_mcg_kg_min = normalized$rate_mcg_kg_min,
    vaso_unit_class = normalized$vaso_unit_class
  )]
} else {
  input_dt[, `:=`(rate_mcg_kg_min = numeric(), vaso_unit_class = character())]
}

vaso_unit_qc <- input_dt[, .(
  n_rows = .N,
  n_rate_nonmissing = sum(!is.na(rate)),
  n_standardized = sum(!is.na(rate_mcg_kg_min))
), by = .(itemid, rateuom, vaso_unit_class)]
fwrite(vaso_unit_qc, "vasopressor_rate_unit_qc.csv")

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  max(x)
}

sofa_vaso <- input_dt[, .(
  dopa_max   = safe_max(rate_mcg_kg_min[itemid %in% VASO_ITEMIDS$dopa]),
  epi_max    = safe_max(rate_mcg_kg_min[itemid %in% VASO_ITEMIDS$epi]),
  norepi_max = safe_max(rate_mcg_kg_min[itemid %in% VASO_ITEMIDS$norepi]),
  dobu_max   = safe_max(rate_mcg_kg_min[itemid %in% VASO_ITEMIDS$dobu]),
  any_vaso_record = .N > 0L
), by = stay_id]

MIN_UO_OBSERVABLE_HOURS <- 23.5

uo_window <- cohort[, .(stay_id, T0, intime, outtime_fix, end_win)]
uo_window[, uo_observation_start := fifelse(T0 > intime, T0, intime)]
uo_window[, uo_observation_end := fifelse(end_win < outtime_fix, end_win, outtime_fix)]
uo_window[, uo_observable_hours := pmax(
  0,
  as.numeric(difftime(uo_observation_end, uo_observation_start, units = "hours"))
)]
uo_window[, uo_full_24h_window := uo_observable_hours >= MIN_UO_OBSERVABLE_HOURS]

uo_dt <- outputevents[
  cohort,
  on = .(stay_id, charttime >= T0, charttime <= end_win),
  nomatch = 0,
  .(
    stay_id = i.stay_id,
    charttime = x.charttime,
    itemid = x.itemid,
    value = x.value
  )
][itemid %in% URINE_IDS]

uo_dt[, value_num := suppressWarnings(as.numeric(value))]
bad_n <- uo_dt[is.na(value_num) & !is.na(value), .N]
if (bad_n > 0) {
  warning(sprintf("Urine rows with non-numeric values were dropped: %d", bad_n), call. = FALSE)
}
uo_dt <- uo_dt[!is.na(value_num)]
uo_dt[, val_signed := fifelse(itemid %in% IRRIGANT_IDS & value_num > 0, -value_num, value_num)]

uo_sum <- uo_dt[, .(
  uo_24h_raw = sum(val_signed, na.rm = TRUE),
  n_uo_events = .N,
  first_uo_time = min(charttime),
  last_uo_time = max(charttime)
), by = stay_id]
uo_sum[uo_24h_raw < 0, uo_24h_raw := 0]

uo_agg <- merge(uo_window, uo_sum, by = "stay_id", all.x = TRUE)
uo_agg[, uo_24h_valid := fifelse(
  uo_full_24h_window & !is.na(n_uo_events) & n_uo_events > 0,
  uo_24h_raw,
  NA_real_
)]

fwrite(
  uo_agg[, .(
    stay_id,
    uo_observable_hours,
    uo_full_24h_window,
    n_uo_events,
    uo_24h_raw,
    uo_24h_valid,
    first_uo_time,
    last_uo_time
  )],
  "urine_output_window_qc.csv.gz"
)

cat("Calculating SOFA...\n")

dt <- merge(cohort, sofa_labs, by="hadm_id", all.x=TRUE)
dt <- merge(dt, sofa_resp, by="hadm_id", all.x=TRUE)
dt <- merge(dt, sofa_map, by="stay_id", all.x=TRUE)
dt <- merge(dt, sofa_cns, by="stay_id", all.x=TRUE)
dt <- merge(dt, sofa_vaso, by="stay_id", all.x=TRUE)
dt <- merge(dt, vent_flags, by="stay_id", all.x=TRUE)
dt <- merge(dt, uo_agg, by="stay_id", all.x=TRUE)

num_cols <- c("pfr_min", "plt_min", "bili_max", "crea_max",
              "map_min", "gcs_min", "gcs_min_ett1",
              "dopa_max", "epi_max", "norepi_max", "dobu_max")
for(col in num_cols) {
  if(col %in% names(dt)) set(dt, i=which(is.infinite(dt[[col]])), j=col, value=NA)
}
dt[is.na(is_vent), is_vent := FALSE]

dt[, resp_score := fcase(
  is.na(pfr_min), NA_real_,
  pfr_min < 100 & is_vent, 4,
  pfr_min < 200 & is_vent, 3,
  pfr_min < 300, 2,
  pfr_min < 400, 1,
  default = 0
)]

dt[, coag_score := fcase(
  is.na(plt_min), NA_real_,
  plt_min < 20, 4,
  plt_min < 50, 3,
  plt_min < 100, 2,
  plt_min < 150, 1,
  default = 0
)]

dt[, liver_score := fcase(
  is.na(bili_max), NA_real_,
  bili_max >= 12.0, 4,
  bili_max >= 6.0, 3,
  bili_max >= 2.0, 2,
  bili_max >= 1.2, 1,
  default = 0
)]

dt[, cv_score := NA_real_]
dt[!is.na(map_min), cv_score := fifelse(map_min < 70, 1, 0)]
dt[(!is.na(dopa_max) & dopa_max > 0) |
   (!is.na(dobu_max) & dobu_max > 0), cv_score := 2]
dt[(!is.na(dopa_max) & dopa_max > 5) |
   (!is.na(epi_max) & epi_max > 0 & epi_max <= 0.1) |
   (!is.na(norepi_max) & norepi_max > 0 & norepi_max <= 0.1), cv_score := 3]
dt[(!is.na(dopa_max) & dopa_max > 15) |
   (!is.na(epi_max) & epi_max > 0.1) |
   (!is.na(norepi_max) & norepi_max > 0.1), cv_score := 4]

dt[, cns_score := fcase(
  is.na(gcs_min), NA_real_,
  gcs_min < 6, 4,
  gcs_min < 10, 3,
  gcs_min < 13, 2,
  gcs_min < 15, 1,
  default = 0
)]
dt[, cns_score_ett1 := fcase(
  is.na(gcs_min_ett1), NA_real_,
  gcs_min_ett1 < 6, 4,
  gcs_min_ett1 < 10, 3,
  gcs_min_ett1 < 13, 2,
  gcs_min_ett1 < 15, 1,
  default = 0
)]

dt[, renal_score := NA_real_]
dt[!is.na(crea_max), renal_score := 0]
dt[!is.na(crea_max) & crea_max >= 1.2, renal_score := 1]
dt[!is.na(crea_max) & crea_max >= 2.0, renal_score := 2]
dt[!is.na(crea_max) & crea_max >= 3.5, renal_score := 3]
dt[!is.na(crea_max) & crea_max >= 5.0, renal_score := 4]
dt[!is.na(uo_24h_valid) & uo_24h_valid < 500,
   renal_score := pmax(fcoalesce(renal_score, 0.0), 3)]
dt[!is.na(uo_24h_valid) & uo_24h_valid < 200,
   renal_score := pmax(fcoalesce(renal_score, 0.0), 4)]

score_cols <- c(
  "resp_score", "coag_score", "liver_score",
  "cv_score", "cns_score", "renal_score"
)

dt[, n_sofa_components_available := rowSums(!is.na(.SD)), .SDcols = score_cols]
dt[, total_sofa_lower_bound := rowSums(.SD, na.rm = TRUE), .SDcols = score_cols]
dt[, total_sofa_complete := fifelse(
  n_sofa_components_available == length(score_cols),
  total_sofa_lower_bound,
  NA_real_
)]

dt[, total_sofa_lower_bound_ett1 := rowSums(
  cbind(resp_score, coag_score, liver_score, cv_score, cns_score_ett1, renal_score),
  na.rm = TRUE
)]

dt[, total_sofa := total_sofa_lower_bound]
dt[, sofa_below2_with_missing := total_sofa_lower_bound < 2 & n_sofa_components_available < 6]

sepsis3_cohort <- dt[total_sofa_lower_bound >= 2]

fwrite(
  dt[, .N, by = .(
    n_sofa_components_available,
    total_sofa_ge2 = total_sofa_lower_bound >= 2,
    sofa_below2_with_missing
  )][order(n_sofa_components_available, -total_sofa_ge2)],
  "sofa_component_availability_qc.csv"
)

cat(sprintf("Sepsis-3 Cohort extracted: %d stays (from %d candidates)\n", nrow(sepsis3_cohort), nrow(dt)))

fwrite(sepsis3_cohort, sepsis_out)
cat("Done.\n")
