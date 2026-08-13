
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
gc()

library(lubridate)

input_dir <- mimic_root

admissions <- read_csv(
  file.path(input_dir, "hosp/admissions.csv.gz"),
  col_select = c(subject_id, hadm_id, admittime, dischtime)
)

patients <- read_csv(
  file.path(input_dir, "hosp/patients.csv.gz"),
  col_select = c(subject_id, anchor_year, anchor_age, gender, dod)
)

labevents <- read_csv(
  file.path(input_dir, "hosp/labevents.csv.gz"),
  col_select = c(subject_id, hadm_id, charttime, itemid, value, valueuom)
)

omr <- read_csv(
  file.path(input_dir, "hosp/omr.csv.gz"),
  col_select = c(subject_id, result_name, result_value)
)

icustays <- read_csv(
  file.path(input_dir, "icu/icustays.csv.gz"),
  col_select = c(subject_id, hadm_id, stay_id, intime, outtime)
)

chartevents <- read_csv(
  file.path(input_dir, "icu/chartevents.csv.gz"),
  col_select = c(subject_id, hadm_id, stay_id, charttime, itemid, value, valueuom)
) %>%
  mutate(value = suppressWarnings(as.numeric(value)))

d_items <- read_csv(
  file.path(input_dir, "icu/d_items.csv.gz"),
  col_select = c(itemid, label)
)

core_ids <- c(223848, 223849, 229314, 225792, 226260, 225411, 228719)

intubation_itemids <- d_items %>%
  filter(itemid %in% core_ids) %>%
  mutate(label_lower = str_to_lower(label)) %>%
  dplyr::select(itemid, label, label_lower)

rm(d_items)

bmi <- omr %>%
  filter(str_detect(result_name, regex("BMI", ignore_case = TRUE))) %>%
  group_by(subject_id) %>%
  summarise(
    avg_value = median(as.numeric(result_value), na.rm = TRUE),
    .groups = "drop"
  )
rm(omr)

PRE_ICU_INTUB_HOURS <- 6

SEPSIS_COHORT_FILE <- Sys.getenv(
  "MIMIC_SEPSIS_COHORT_FILE",
  unset = "mimic_iv_sepsis3_cohort_final.csv.gz"
)
if (!file.exists(SEPSIS_COHORT_FILE)) {
  stop("Sepsis cohort file not found: ", SEPSIS_COHORT_FILE)
}
sepsis3_raw <- data.table::fread(
  SEPSIS_COHORT_FILE,
  select = c("subject_id", "hadm_id", "stay_id", "T0", "total_sofa")
)

cat("sepsis3_raw subjects:", uniqueN(sepsis3_raw$subject_id), "\n")
cat("sepsis3_raw stays   :", uniqueN(sepsis3_raw$stay_id), "\n")

sepsis3_raw[, T0 := as.POSIXct(T0)]
setorder(sepsis3_raw, subject_id, hadm_id, stay_id, T0)
sepsis_stays <- sepsis3_raw[
  , .SD[1],
  by = .(subject_id, hadm_id, stay_id)
][, .(subject_id, hadm_id, stay_id, T0, total_sofa)]

cat("sepsis_stays rows   :", nrow(sepsis_stays), "\n")

stays_dt <- data.table::as.data.table(icustays)[
  data.table::as.data.table(sepsis_stays),
  on = .(subject_id, hadm_id, stay_id),
  nomatch = 0
]

cat("stays_dt rows (after join icustays):", nrow(stays_dt), "\n")
cat("stays_dt unique stays:", uniqueN(stays_dt$stay_id), "\n")
cat("stays_dt unique subjects:", uniqueN(stays_dt$subject_id), "\n")

candidate_stays <- unique(stays_dt$stay_id)

PRE_H  <- 6
POST_H <- 6

proc_ev <- data.table::fread(
  file.path(input_dir, "icu/procedureevents.csv.gz"),
  select = c("stay_id", "itemid", "starttime", "endtime")
)[stay_id %in% candidate_stays]

proc_ev[, starttime := as.POSIXct(starttime)]
proc_ev[, endtime   := as.POSIXct(endtime)]

B1_vent_window <- proc_ev[itemid == 225792][
  stays_dt[, .(stay_id, intime)],
  on = "stay_id",
  nomatch = 0
][
  starttime <= (intime + lubridate::hours(POST_H)) &
    (is.na(endtime) | endtime >= (intime - lubridate::hours(PRE_H))),
  .(stay_id)
]
B1_vent_window <- unique(B1_vent_window)

mode_itemids <- c(223848, 223849, 226732)
inv_mode_pat <- "\\b(assist[-\\s]*control|assist/control|a/c|cmv|simv|prvc|aprv|bilevel|bi-level)\\b"

chart_mode <- data.table::fread(
  file.path(input_dir, "icu/chartevents.csv.gz"),
  select = c("stay_id", "charttime", "itemid", "value")
)
chart_mode <- chart_mode[stay_id %in% candidate_stays & itemid %in% mode_itemids]
chart_mode[, charttime := as.POSIXct(charttime)]

mode_on_admit <- chart_mode[
  stays_dt[, .(stay_id, intime)],
  on = "stay_id",
  nomatch = 0
][
  charttime >= (intime - lubridate::hours(1)) &
    charttime <= (intime + lubridate::hours(1))
][
  !is.na(value) & grepl(inv_mode_pat, tolower(as.character(value)), perl = TRUE),
  .(stay_id)
]
mode_on_admit <- unique(mode_on_admit)

tube_near <- proc_ev[itemid == 224385][
  stays_dt[, .(stay_id, intime)],
  on = "stay_id",
  nomatch = 0
][
  starttime <= intime &
    starttime >= (intime - lubridate::hours(PRE_ICU_INTUB_HOURS)),
  .(stay_id)
]
tube_near <- unique(tube_near)

cat("Peri-admission invasive-ventilation stays:", uniqueN(B1_vent_window$stay_id), "\n")
cat("Mode-supported stays:", uniqueN(mode_on_admit$stay_id), "\n")
cat("Intubation-event-supported stays:", uniqueN(tube_near$stay_id), "\n")

intubated_on_admit_stays <- unique(B1_vent_window)
cat("Primary peri-admission ventilation stays:", uniqueN(intubated_on_admit_stays$stay_id), "\n")
cat("Supportive mode-on-admission stays:", uniqueN(mode_on_admit$stay_id), "\n")
cat("Supportive near-intubation stays:", uniqueN(tube_near$stay_id), "\n")

stays_intub <- stays_dt[stay_id %in% intubated_on_admit_stays$stay_id]

icustays_intubated_first <- stays_intub[
  order(subject_id, T0)
][
  , .SD[1], by = subject_id
][
  , .(stay_id, subject_id, hadm_id, intime, outtime)
]

cat("final (per-subject) subjects:", uniqueN(icustays_intubated_first$subject_id), "\n")
cat("final (per-subject) stays:", uniqueN(icustays_intubated_first$stay_id), "\n")

rm(sepsis3_raw, sepsis_stays, stays_dt, candidate_stays,
   proc_ev, chart_mode, mode_on_admit, tube_near,
   intubated_on_admit_stays, stays_intub)
gc()

rm(intubation_itemids)

chartevents_sepsis <- chartevents %>%
  filter(subject_id %in% unlist(icustays_intubated_first$subject_id))
rm(chartevents)

labevents_sepsis <- labevents %>%
  filter(subject_id %in% unlist(icustays_intubated_first$subject_id))
rm(labevents)

patients_sepsis <- patients %>%
  filter(subject_id %in% unlist(icustays_intubated_first$subject_id))
rm(patients)

icustays_intubated_first_patient_admission <- icustays_intubated_first %>%
  left_join(patients_sepsis, by = "subject_id") %>%
  left_join(admissions, by = c("subject_id", "hadm_id")) %>%
  left_join(bmi, by = "subject_id") %>%
  mutate(
    age = as.numeric(format(intime, "%Y")) - anchor_year + anchor_age,
    sex = gender
  ) %>%
  filter(age >= 18)

rm(icustays_intubated_first)
rm(patients_sepsis)
rm(admissions)
rm(bmi)

write_csv(icustays_intubated_first_patient_admission, "icustays_intubated_first_patient_admission.csv")
icustays_intubated_first_patient_admission <- read_csv("icustays_intubated_first_patient_admission.csv")

procedureevents <- read_csv(
  file.path(input_dir, "icu/procedureevents.csv.gz"),
  col_select = c(subject_id, hadm_id, stay_id, starttime, itemid, value, valueuom)
)
procedureevents_sepsis <- procedureevents %>%
  filter(subject_id %in% unlist(icustays_intubated_first_patient_admission$subject_id))
rm(procedureevents)

if (!file.exists(resource_path("mimic_item_map.csv"))) stop("resources/mimic_item_map.csv is missing.", call. = FALSE)
itemid_name_table <- read_csv(file = resource_path("mimic_item_map.csv"))
lab_items <- unlist(itemid_name_table$itemid)

lab_sub <- labevents_sepsis %>%
  filter(itemid %in% lab_items) %>%
  dplyr::select(subject_id, hadm_id, charttime, itemid, value, valueuom)
rm(labevents_sepsis)

vital_sub <- chartevents_sepsis %>%
  filter(itemid %in% lab_items) %>%
  dplyr::select(stay_id, subject_id, hadm_id, charttime, itemid, value, valueuom)
rm(chartevents_sepsis)

proc_sub <- procedureevents_sepsis %>%
  filter(itemid %in% lab_items) %>%
  dplyr::select(stay_id, subject_id, hadm_id, starttime, itemid, value, valueuom)
rm(procedureevents_sepsis)

rm(lab_items)

lab_sub <- lab_sub %>% filter(!is.na(hadm_id))
lab_icustay <- lab_sub %>%
  inner_join(
    icustays_intubated_first_patient_admission,
    by = c("subject_id", "hadm_id")
  ) %>%
  filter(charttime >= intime & charttime <= outtime) %>%
  mutate(stay_id = stay_id) %>%
  dplyr::select(stay_id, subject_id, hadm_id, charttime, itemid, value, valueuom, intime, outtime)
rm(lab_sub)

vital_sub <- vital_sub %>% filter(!is.na(hadm_id))
vital_icustay <- vital_sub %>%
  inner_join(
    icustays_intubated_first_patient_admission,
    by = c("stay_id", "subject_id", "hadm_id")
  ) %>%
  filter(charttime >= intime & charttime <= outtime) %>%
  dplyr::select(stay_id, subject_id, hadm_id, charttime, itemid, value, valueuom, intime, outtime)
rm(vital_sub)

proc_icustay <- proc_sub %>%
  inner_join(
    icustays_intubated_first_patient_admission,
    by = c("stay_id", "subject_id", "hadm_id")
  ) %>%
  filter(starttime >= intime & starttime <= outtime) %>%
  dplyr::select(stay_id, subject_id, hadm_id, starttime, itemid, value, valueuom, intime, outtime) %>%
  rename(charttime = starttime)
rm(proc_sub)

lab_icustay %>% filter(itemid == 50825) %>% dplyr::select(value)
vital_icustay <- vital_icustay %>%
  mutate(value = ifelse(itemid == 223761, 5 / 9 * (as.numeric(value) - 32), value))

lab_data_wide <- lab_icustay %>%
  filter(charttime >= intime & charttime <= outtime) %>%
  mutate(charttime_year_day = format(charttime, "%Y-%m-%d")) %>%
  dplyr::select(stay_id, charttime_year_day, itemid, hadm_id, itemid, value) %>%
  group_by(stay_id, charttime_year_day, itemid, hadm_id) %>%
  summarise(value = median(as.numeric(value), na.rm = TRUE), .groups="drop") %>%
  pivot_wider(names_from = itemid, values_from = value)

write_csv(lab_data_wide, "mimiciv.lab.wide.csv")

vital_data_wide <- vital_icustay %>%
  filter(charttime >= intime & charttime <= outtime) %>%
  mutate(charttime_year_day = format(charttime, "%Y-%m-%d")) %>%
  dplyr::select(stay_id, charttime_year_day, itemid, hadm_id, itemid, value) %>%
  group_by(stay_id, charttime_year_day, itemid, hadm_id) %>%
  summarise(value = median(as.numeric(value), na.rm = TRUE), .groups="drop") %>%
  pivot_wider(names_from = itemid, values_from = value)

write_csv(vital_data_wide, "mimiciv.vital.wide.csv")

proc_data_wide <- proc_icustay %>%
  filter(charttime >= intime & charttime <= outtime) %>%
  mutate(charttime_year_day = format(charttime, "%Y-%m-%d")) %>%
  dplyr::select(stay_id, charttime_year_day, itemid, hadm_id, itemid, value) %>%
  group_by(stay_id, charttime_year_day, itemid, hadm_id) %>%
  summarise(value = median(as.numeric(value), na.rm = TRUE), .groups="drop") %>%
  pivot_wider(names_from = itemid, values_from = value)

write_csv(proc_data_wide, "mimiciv.proc.wide.csv")

lab_data_wide <- read_csv("mimiciv.lab.wide.csv")
vital_data_wide <- read_csv("mimiciv.vital.wide.csv")
proc_data_wide <- read_csv("mimiciv.proc.wide.csv")

all_day_keys <- bind_rows(
  lab_data_wide %>% transmute(stay_id, hadm_id, charttime_year_day, source = "lab"),
  vital_data_wide %>% transmute(stay_id, hadm_id, charttime_year_day, source = "vital"),
  proc_data_wide %>% transmute(stay_id, hadm_id, charttime_year_day, source = "procedure")
) %>%
  distinct() %>%
  mutate(present = 1L) %>%
  tidyr::pivot_wider(names_from = source, values_from = present, values_fill = 0L)
for (nm in c("lab", "vital", "procedure")) {
  if (!nm %in% names(all_day_keys)) all_day_keys[[nm]] <- 0L
}
analytical_day_qc <- all_day_keys %>%
  mutate(
    included_primary = lab == 1L,
    source_pattern = paste0("L", lab, "_V", vital, "_P", procedure)
  ) %>%
  count(source_pattern, included_primary, name = "n_stay_days")
write_csv(analytical_day_qc, "analytical_day_source_qc.csv")

combined <- lab_data_wide %>%
  left_join(vital_data_wide, by = c("stay_id", "hadm_id", "charttime_year_day")) %>%
  left_join(proc_data_wide, by = c("stay_id", "hadm_id", "charttime_year_day"))
combined_keep <- combined

fio2_itemids <- intersect(c("223835","229280","229841"), colnames(combined_keep))
clean_fio2 <- function(x){
  x <- suppressWarnings(as.numeric(x))

  x[!is.na(x) & x == 0] <- NA_real_

  hit_frac <- !is.na(x) & x > 0 & x <= 1
  x[hit_frac] <- x[hit_frac] * 100

  x[!is.na(x) & x < 21]  <- NA_real_
  x[!is.na(x) & x > 100] <- NA_real_

  x
}

for(cc in fio2_itemids){
  combined_keep[[cc]] <- clean_fio2(combined_keep[[cc]])
}

if (!file.exists(resource_path("mimic_item_map.csv"))) stop("resources/mimic_item_map.csv is missing.", call. = FALSE)
itemid_name_table <- read_csv(file = resource_path("mimic_item_map.csv"))
key_cols <- intersect(c("stay_id", "charttime_year_day", "hadm_id"), names(combined_keep))
hmm_item_cols <- intersect(as.character(itemid_name_table$itemid), names(combined_keep))
combined <- combined_keep[, unique(c(key_cols, hmm_item_cols)), drop = FALSE]

mad_outliers <- function(vec, threshold = 100) {
  median_x <- median(vec, na.rm = TRUE)
  mad_x <- mad(vec, constant = 1.4826, na.rm = TRUE)
  if(mad_x == 0) return(vec)
  z_score <- abs(vec - median_x) / mad_x
  vec[which(z_score > threshold)] <- NA
  return(vec)
}

replace_outliers_with_na_grubbs <- function(vec, alpha = 0.05) {

  if (sum(!is.na(vec)) < 3) return(vec)

  repeat {
    med <- median(vec, na.rm = TRUE)
    if (is.na(med)) break

    one <- sort(vec, decreasing = TRUE, na.last = NA)[1]
    if (is.na(one)) break

    two_vec <- vec[vec != one]
    two <- sort(two_vec, decreasing = TRUE, na.last = NA)[1]
    if (is.na(two)) break

    if (one < med | two < med) break

    if (one - two > two - med) {
      vec[which(vec == one)] <- NA
    } else {
      break
    }
  }

  repeat {
    med <- median(vec, na.rm = TRUE)
    if (is.na(med)) break

    one <- sort(vec, decreasing = TRUE, na.last = NA)[1]
    if (is.na(one)) break

    two_vec <- vec[vec != one]
    two <- sort(two_vec, decreasing = FALSE, na.last = NA)[1]
    if (is.na(two)) break

    if (one < med | two < med) break

    if (two - one > med - two) {
      vec[which(vec == one)] <- NA
    } else {
      break
    }
  }

  if (length(which(!is.na(vec))) > 50) {
    return(mad_outliers(vec, threshold = 50))
  }

  vec_keep <- vec
  count <- 0
  repeat {
    if (sum(!is.na(vec)) < 3) break

    test_result <- tryCatch(grubbs.test(as.numeric(vec)), error = function(e) NULL)
    if (is.null(test_result)) break

    p_value <- test_result$p.value
    if (p_value < alpha) {
      outlier_value <- as.numeric(strsplit(test_result$alternative, " ")[[1]][3])
      count <- count + 1
      vec[which(vec == outlier_value)] <- NA
    } else {
      break
    }
    if (count > 10) return(vec_keep)
  }
  return(vec)
}

combined[combined == 999999] <- NA
combined <- combined[, colSums(is.na(combined)) < (nrow(combined) - 2)]
combined <- combined[, unlist(lapply(combined, class)) != "logical"]

for(j in 4:ncol(combined)){
  vec <- unlist(combined[, j])
  long_tail <- FALSE
  if(colnames(combined)[j] == "50806") vec[vec < 20  & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "50810") vec[vec > 100 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "50811") vec[vec > 80  & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "50813") long_tail <- TRUE
  if(colnames(combined)[j] == "50824") vec[vec < 50 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "50910") long_tail <- TRUE
  if(colnames(combined)[j] == "50912") long_tail <- TRUE
  if(colnames(combined)[j] == "51237") long_tail <- TRUE
  if(colnames(combined)[j] == "51274") long_tail <- TRUE
  if(colnames(combined)[j] == "51301") long_tail <- TRUE
  if(colnames(combined)[j] == "50861") long_tail <- TRUE
  if(colnames(combined)[j] == "50861") vec[vec > 10000 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "50878") long_tail <- TRUE
  if(colnames(combined)[j] == "50863") long_tail <- TRUE
  if(colnames(combined)[j] == "50885") long_tail <- TRUE
  if(colnames(combined)[j] == "51133") long_tail <- TRUE
  if(colnames(combined)[j] == "51200") long_tail <- TRUE
  if(colnames(combined)[j] == "52069") long_tail <- TRUE
  if(colnames(combined)[j] == "52073") long_tail <- TRUE
  if(colnames(combined)[j] == "52074") long_tail <- TRUE
  if(colnames(combined)[j] == "50954") long_tail <- TRUE
  if(colnames(combined)[j] == "50924") long_tail <- TRUE
  if(colnames(combined)[j] == "50825") vec[vec  < 25 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "223761") vec[vec < 25 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "223762") vec[vec < 25 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "223762") vec[vec > 50 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "229280") vec[vec < 10 & !is.na(vec)] <- NA
  if(!long_tail) {
    vec <- replace_outliers_with_na_grubbs(vec, 1.0e-16)
  }
  combined[, j] <- vec
}

tmp <- itemid_name_table$matched_test[match(colnames(combined)[-c(1, 2, 3)], itemid_name_table$itemid)]
for(j in 4:ncol(combined)){
  vec <- unlist(combined[, j])
  long_tail <- FALSE
  if(colnames(combined)[j] == "50813") long_tail <- TRUE
  if(colnames(combined)[j] == "50910") long_tail <- TRUE
  if(colnames(combined)[j] == "50912") long_tail <- TRUE
  if(colnames(combined)[j] == "51237") long_tail <- TRUE
  if(colnames(combined)[j] == "51274") long_tail <- TRUE
  if(colnames(combined)[j] == "51301") long_tail <- TRUE
  if(colnames(combined)[j] == "50861") long_tail <- TRUE
  if(colnames(combined)[j] == "50878") long_tail <- TRUE
  if(colnames(combined)[j] == "50863") long_tail <- TRUE
  if(colnames(combined)[j] == "50885") long_tail <- TRUE
  if(colnames(combined)[j] == "51133") long_tail <- TRUE
  if(colnames(combined)[j] == "51200") long_tail <- TRUE
  if(colnames(combined)[j] == "52069") long_tail <- TRUE
  if(colnames(combined)[j] == "52073") long_tail <- TRUE
  if(colnames(combined)[j] == "52074") long_tail <- TRUE
  if(colnames(combined)[j] == "50954") long_tail <- TRUE
  if(colnames(combined)[j] == "50924") long_tail <- TRUE
  if(colnames(combined)[j] == "51249")  vec[vec <  20 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "220602") vec[vec <   30 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "223830") vec[vec <   5 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "226534") vec[vec <  50 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "226536") vec[vec <  50 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "226540") vec[vec > 100 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "227443") vec[vec >  80 & !is.na(vec)] <- NA
  if(colnames(combined)[j] == "227464") vec[vec >  15 & !is.na(vec)] <- NA
  if(!long_tail) {
    vec <- replace_outliers_with_na_grubbs(vec, 1.0e-16)
  }
  combined[, j] <- vec
}

for(id in unique(tmp)){
  hit <- which(id == tmp)
  if(length(hit) > 1){
  }
}

tmp <- itemid_name_table$matched_test[match(colnames(combined)[-c(1, 2, 3)], itemid_name_table$itemid)]
for(id in unique(tmp)){
  hit <- which(id == tmp)
  if(length(hit) > 1){
    combined <- combined %>% mutate(!!id := apply(combined[, hit + 3], 1, function(x) median(x, na.rm = TRUE)))
  } else {
    colnames(combined)[hit + 3] <- id
  }
}
combined <- combined[, -(which(!is.na(match(colnames(combined)[-c(1, 2, 3)], itemid_name_table$itemid))) + 3)]
combined$ART.D. <- ifelse(unlist(combined$ART.D.) < 25, NA, unlist(combined$ART.D.))
combined$ART.S. <- ifelse(unlist(combined$ART.S.) < 40, NA, unlist(combined$ART.S.))
combined$ART.M. <- ifelse(unlist(combined$ART.M.) < 30, NA, unlist(combined$ART.M.))
combined$MEAN <- ifelse(unlist(combined$MEAN) < 20, NA, unlist(combined$MEAN))
combined$SpO2 <- ifelse(unlist(combined$SpO2) < 50, NA, unlist(combined$SpO2))
combined$RR <- ifelse(unlist(combined$RR) < 1, NA, unlist(combined$RR))
combined$AnionGap <- ifelse(unlist(combined$AnionGap) < 0, NA, unlist(combined$AnionGap))

combined <- combined %>% mutate(SYS. = ifelse(is.na(unlist(SYS.)), ART.S., SYS.))
combined <- combined %>% mutate(DIAS. = ifelse(is.na(unlist(DIAS.)), ART.D., DIAS.))
combined <- combined %>% mutate(MEAN = ifelse(is.na(unlist(MEAN)), ART.M., MEAN))
combined <- combined[, order(colSums(!is.na(combined)), decreasing = T)]

source(code_path("02_fio2_imputation.R"))
combined <- impute_fio2_from_o2flow(combined, input_dir = input_dir, fio2_col = "FiO2")

write_csv(combined, "mimiciv.combined.wide.csv")

