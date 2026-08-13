
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(zoo)
})

impute_fio2_from_o2flow <- function(combined,
                                    input_dir,
                                    fio2_col = "FiO2",
                                    stay_col = "stay_id",
                                    day_col  = "charttime_year_day",

                                    o2flow_itemids = c(223834, 227287, 227582),

                                    o2device_itemids = c(226732),
                                    room_air_device_regex = "(room\\s*air|none|no\\s*o2|ra\\b)",

                                    vent_itemids = c(225792),

                                    locf_maxgap_days = 1,

                                    fill_room_air = TRUE,
                                    room_air_fio2 = 0.21,

                                    require_no_oxygen_evidence_for_room_air = TRUE,

                                    force_percent_output = NA,

                                    keep_flags = FALSE) {

  dt <- as.data.table(combined)
  if(!all(c(stay_col, day_col, fio2_col) %in% names(dt))) {
    stop("required columns not found: stay_id/charttime_year_day/FiO2")
  }

  dayvec <- dt[[day_col]]
  if(inherits(dayvec, "Date")) {
    dt[, .day := as.Date(dayvec)]
  } else if(inherits(dayvec, "POSIXct")) {
    dt[, .day := as.Date(dayvec)]
  } else if(is.numeric(dayvec)) {
    dt[, .day := as.Date(dayvec, origin = "1970-01-01")]
  } else {
    dt[, .day := as.Date(as.character(dayvec))]
  }
  dt[, .stay := get(stay_col)]

  fio2_to_frac <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    y <- ifelse(is.na(x), NA_real_, ifelse(x > 1.5, x/100, x))
    y[y == 0] <- NA_real_
    y[y < 0.05 | y > 1.5] <- NA_real_
    y
  }

  flow_to_fio2 <- function(flow) {
    f <- suppressWarnings(as.numeric(flow))
    y <- 0.21 + 0.04 * pmax(f, 0)
    y <- pmin(pmax(y, 0.21), 0.80)
    y[!is.finite(y)] <- NA_real_
    y
  }

  fio2_max <- suppressWarnings(max(as.numeric(dt[[fio2_col]]), na.rm = TRUE))
  auto_is_percent <- is.finite(fio2_max) && fio2_max > 1.5
  is_percent <- if (is.na(force_percent_output)) auto_is_percent else isTRUE(force_percent_output)

  target_stays <- unique(dt[[stay_col]])

  ce_flow <- fread(file.path(input_dir, "icu/chartevents.csv.gz"),
                   select = c("stay_id", "charttime", "itemid", "value"))
  ce_flow <- ce_flow[stay_id %in% target_stays & itemid %in% o2flow_itemids]
  if(nrow(ce_flow) > 0) {
    ce_flow[, charttime := as.POSIXct(charttime)]
    ce_flow[, .day := as.Date(charttime)]
    ce_flow[, .stay := stay_id]

    flow_daily <- ce_flow[, .(flow = median(suppressWarnings(as.numeric(value)), na.rm = TRUE)),
                          by = .(.stay, .day, itemid)]
    flow_daily <- flow_daily[is.finite(flow)]
    flow_daily <- flow_daily[, .(O2Flow = max(flow, na.rm = TRUE)), by = .(.stay, .day)]
    flow_daily[!is.finite(O2Flow), O2Flow := NA_real_]
  } else {
    flow_daily <- data.table(.stay = integer(), .day = as.Date(character()), O2Flow = numeric())
  }

  ce_dev <- fread(file.path(input_dir, "icu/chartevents.csv.gz"),
                  select = c("stay_id", "charttime", "itemid", "value"))
  ce_dev <- ce_dev[stay_id %in% target_stays & itemid %in% o2device_itemids]
  if(nrow(ce_dev) > 0) {
    ce_dev[, charttime := as.POSIXct(charttime)]
    ce_dev[, .day := as.Date(charttime)]
    ce_dev[, .stay := stay_id]
    ce_dev[, dev := tolower(as.character(value))]

    dev_daily <- ce_dev[
      !is.na(dev) & dev != "",
      .(
        O2Device_any = 1L,
        O2Device_nonroomair = as.integer(any(!grepl(room_air_device_regex, dev, perl = TRUE)))
      ),
      by = .(.stay, .day)
    ]
  } else {
    dev_daily <- data.table(.stay = integer(), .day = as.Date(character()),
                            O2Device_any = integer(), O2Device_nonroomair = integer())
  }

  pe <- fread(file.path(input_dir, "icu/procedureevents.csv.gz"),
              select = c("stay_id", "itemid", "starttime", "endtime"))
  pe <- pe[stay_id %in% target_stays & itemid %in% vent_itemids]
  if(nrow(pe) > 0) {

    pe[, starttime := as.POSIXct(starttime)]
    pe[, endtime   := as.POSIXct(endtime)]
    pe[, endtime   := data.table::fcoalesce(endtime, starttime)]

    pe[, start_day := as.Date(starttime)]
    pe[, end_day   := as.Date(endtime)]

    pe[, start_i := as.integer(start_day)]
    pe[, end_i   := as.integer(end_day)]
    pe <- pe[!is.na(start_i) & !is.na(end_i) & end_i >= start_i]

    vent_daily <- pe[
      ,
      .(.day = as.Date(seq.int(from = start_i, to = end_i, by = 1L), origin = "1970-01-01")),
      by = .(stay_id, start_i, end_i)
    ]

    vent_daily[, .stay := stay_id]
    vent_daily[, Vent_any := 1L]
    vent_daily <- unique(vent_daily[, .(.stay, .day, Vent_any)])

    pe[, c("start_i", "end_i") := NULL]

  } else {
    vent_daily <- data.table(.stay = integer(), .day = as.Date(character()), Vent_any = integer())
  }

  dt <- flow_daily[dt, on = c(".stay", ".day")]

  dt <- dev_daily[dt, on = c(".stay", ".day")]
  dt[is.na(O2Device_any), O2Device_any := 0L]
  dt[is.na(O2Device_nonroomair), O2Device_nonroomair := 0L]

  dt <- vent_daily[dt, on = c(".stay", ".day")]
  dt[is.na(Vent_any), Vent_any := 0L]

  dt[, FiO2_frac := fio2_to_frac(get(fio2_col))]

  if(locf_maxgap_days > 0) {
    setorder(dt, .stay, .day)
    dt[, FiO2_frac := zoo::na.locf(FiO2_frac, na.rm = FALSE, maxgap = locf_maxgap_days), by = .stay]
    dt[, FiO2_frac := zoo::na.locf(FiO2_frac, na.rm = FALSE, maxgap = locf_maxgap_days, fromLast = TRUE), by = .stay]
  }

  dt[, FiO2_from_flow := flow_to_fio2(O2Flow)]
  dt[is.na(FiO2_frac) & !is.na(FiO2_from_flow), FiO2_frac := FiO2_from_flow]

  if(fill_room_air) {
    if(require_no_oxygen_evidence_for_room_air) {

      dt[, oxygen_evidence := (Vent_any == 1L) | (O2Device_nonroomair == 1L) | (!is.na(O2Flow) & O2Flow > 0)]
      dt[is.na(FiO2_frac) & !oxygen_evidence, FiO2_frac := room_air_fio2]
    } else {
      dt[is.na(FiO2_frac), FiO2_frac := room_air_fio2]
    }
  }

  if(is_percent) {
    dt[, (fio2_col) := FiO2_frac * 100]
  } else {
    dt[, (fio2_col) := FiO2_frac]
  }

  if(!keep_flags) {
    rm_cols <- c("O2Flow", "O2Device_any", "O2Device_nonroomair", "Vent_any",
                 "FiO2_frac", "FiO2_from_flow", "oxygen_evidence", ".stay", ".day")
    rm_cols <- intersect(rm_cols, names(dt))
    dt[, (rm_cols) := NULL]
  } else {
    rm_cols <- intersect(c("FiO2_frac", "FiO2_from_flow", ".stay", ".day"), names(dt))
    dt[, (rm_cols) := NULL]
  }

  return(as.data.frame(dt))
}
