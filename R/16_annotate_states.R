source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

if (!exists("mimic_state_profile_bundle", inherits = TRUE)) {
  if (!file.exists("state_profile_bundle.rds")) stop("state_profile_bundle.rds not found. Run 15_reconstruct_state_profiles.R first.", call. = FALSE)
  mimic_state_profile_bundle <- readRDS("state_profile_bundle.rds")
}
if (!identical(mimic_state_profile_bundle$version, "state_profile_v1")) stop("Unexpected state-profile bundle version.", call. = FALSE)
AGES <- as.numeric(mimic_state_profile_bundle$ages)
SEXS <- as.numeric(mimic_state_profile_bundle$sexes)
if (!identical(AGES, c(55, 65, 75)) || !identical(SEXS, c(0, 1))) {
  stop("Unexpected annotation reference strata.", call. = FALSE)
}
num <- function(x) suppressWarnings(as.numeric(x))
get_num <- function(r, nm) {
  if (!nm %in% names(r)) return(NA_real_)
  num(r[[nm]])
}

empty_a_split_thresholds <- function() {
  list(
    MAP_A2 = NA_real_, PF_A2 = NA_real_, HCO3_A2 = NA_real_,
    HGB_A2 = NA_real_, FIO2_A2 = NA_real_, SCORE_CUT = Inf, N_A_REF = 0L
  )
}

REL_THRESHOLD_NAMES <- c(
  "MAP_SOFT", "MAP_STR",
  "LAC_SOFT", "LAC_STR",
  "CR_SOFT", "CR_STR", "BUN_SOFT", "BUN_STR",
  "TB_SOFT", "TB_STR", "TB_MIX",
  "AST_SOFT", "AST_STR", "AST_MIX",
  "ALT_SOFT", "ALT_STR", "ALT_MIX",
  "INR_SOFT", "INR_STR",
  "PLT_SOFT", "PLT_STR",
  "HCO3_SOFT", "HCO3_STR",
  "CORE_AST_SOFT", "CORE_AST_STR",
  "PF_MILD", "PF_MOD", "PF_STR",
  "FIO2_MILD", "FIO2_MOD", "FIO2_HIGH", "FIO2_VHIGH",
  "RR_TACHY", "FIO2_ROOM"
)

ABS_THRESHOLD_NAMES <- c(
  "PF_REFRACTORY",
  "PH_CRISIS_LOW",
  "PH_CRISIS_HIGH",
  "TB_FAILURE",
  "AST_FAILURE",
  "ALT_FAILURE",
  "CR_FAILURE",
  "PLT_FAILURE",
  "HGB_FAILURE",
  "INR_FAILURE",
  "NA_CRISIS_LOW",
  "NA_CRISIS_HIGH",
  "K_CRISIS_LOW",
  "K_CRISIS_HIGH",
  "CL_CRISIS_LOW",
  "CL_CRISIS_HIGH",
  "GLU_SEV",
  "GLU_HIGH",
  "RR_RESP",
  "PH_RESP_LOW",
  "PH_RESP_HIGH",
  "PCO2_HIGH",
  "PCO2_LOW",
  "HGB_SOFT"
)
THR_ABS <- NULL

thresholds_to_table <- function(THR_REL, THR_ABS, THR_A_SPLIT = NULL) {
  rel_tbl <- tibble(
    threshold_name = names(THR_REL),
    value = unname(unlist(THR_REL, use.names = FALSE)),
    layer = "relative_main"
  )
  abs_tbl <- tibble(
    threshold_name = names(THR_ABS),
    value = unname(unlist(THR_ABS, use.names = FALSE)),
    layer = "absolute_orthogonal"
  )
  if (is.null(THR_A_SPLIT) || !length(THR_A_SPLIT)) return(bind_rows(rel_tbl, abs_tbl))
  a_split_tbl <- tibble(
    threshold_name = names(THR_A_SPLIT),
    value = unname(unlist(THR_A_SPLIT, use.names = FALSE)),
    layer = "a1_a2_split"
  )
  bind_rows(rel_tbl, abs_tbl, a_split_tbl)
}

is_gamma_rel <- function(tags_rel) {
  resp_sevmod <- any(grepl("^Severe hypoxemia|^Moderate hypoxemia", tags_rel))
  resp_mild <- any(grepl("^Mild hypoxemia", tags_rel))
  any_fio2_need <- any(grepl("FiO2 need", tags_rel))
  support_minor <- ("Respiratory acidosis" %in% tags_rel) ||
    ("Respiratory alkalosis" %in% tags_rel) ||
    any(grepl("^Tachypnea", tags_rel))
  resp_sevmod || (resp_mild && (any_fio2_need || support_minor))
}

beta_subtype_from_flags <- function(f) {
  renal_soft <- isTRUE(f$renal_soft)
  renal_str  <- isTRUE(f$renal_str)
  azot_str   <- isTRUE(f$azot_str)

  if (renal_str || (renal_soft && azot_str)) return("B2 renal severe")
  if (renal_soft) return("B1 renal mild")
  NA_character_
}

delta_subtype_from_flags <- function(f) {
  core_soft  <- isTRUE(f$core_soft)
  core_str   <- isTRUE(f$core_str)
  shock_soft <- isTRUE(f$shock_soft)
  shock_str  <- isTRUE(f$shock_str)
  lac_soft   <- isTRUE(f$lac_soft)
  lac_str    <- isTRUE(f$lac_str)
  hep_soft   <- isTRUE(f$hep_soft)
  hep_str    <- isTRUE(f$hep_str)
  coag_soft  <- isTRUE(f$coag_soft)
  coag_str   <- isTRUE(f$coag_str)

  d1_cond <- core_str || core_soft
  d2_cond <- hep_str || (hep_soft && (coag_soft || lac_str))

  d3_cond <- shock_str || coag_str ||
    (shock_soft && (lac_soft || coag_soft || hep_soft)) ||
    (coag_soft && (shock_soft || lac_str))

  is_delta_entry <- d1_cond || d2_cond || d3_cond
  if (!is_delta_entry) return(NA_character_)
  if (d1_cond) return("D1 lactate/shock")
  if (d2_cond) return("D2 hepatic")
  if (d3_cond) return("D3 shock/coagulation-overlap")
  stop("unreachable: is_delta_entry TRUE but no subtype matched")
}

domains_from_flags <- function(f, gam) {
  dom <- character(0)
  if (isTRUE(f$shock_soft) || isTRUE(f$shock_str) ||
      isTRUE(f$lac_soft)   || isTRUE(f$lac_str)   ||
      isTRUE(f$core_soft)  || isTRUE(f$core_str)) dom <- c(dom, "shock")
  if (isTRUE(gam)) dom <- c(dom, "resp")
  if (isTRUE(f$renal_soft) || isTRUE(f$renal_str)) dom <- c(dom, "renal")
  if (isTRUE(f$hep_soft)   || isTRUE(f$hep_str))   dom <- c(dom, "hepatic")
  if (isTRUE(f$coag_soft)  || isTRUE(f$coag_str))  dom <- c(dom, "coagulation")
  dom <- unique(dom)
  list(
    primary   = if (length(dom) >= 1) dom[1] else "",
    secondary = if (length(dom) >= 2) dom[2] else ""
  )
}


tag_one <- function(r, THR_REL, THR_ABS) {
  if (is.data.frame(r)) r <- as.list(r[1, , drop = FALSE])

  MAP  <- get_num(r, "MAP")
  Lac  <- get_num(r, "Lactate")
  Cr   <- get_num(r, "Creatinine")
  BUN  <- get_num(r, "UreaNitrogen")
  Tb   <- get_num(r, "TotalBilirubin")
  AST  <- get_num(r, "AST")
  ALT  <- get_num(r, "ALT")
  INR  <- get_num(r, "INR")
  Plt  <- get_num(r, "PlateletCount")
  HCO3 <- get_num(r, "HCO3")
  FiO2 <- get_num(r, "FiO2")
  RR   <- get_num(r, "RR")
  pH   <- get_num(r, "pH")
  pCO2 <- get_num(r, "pCO2")
  Hgb  <- get_num(r, "Hemoglobin")
  Na   <- get_num(r, "Sodium")
  K    <- get_num(r, "Potassium")
  Cl   <- get_num(r, "Chloride")
  Glu  <- get_num(r, "Glucose")
  WBC  <- get_num(r, "WBC")
  PF   <- get_num(r, "PF")
  SpO2 <- get_num(r, "SpO2")

  tags_rel <- character(0)
  orth_abs <- character(0)
  rules_rel <- character(0)

  if (!is.na(MAP) && MAP <= THR_REL$MAP_STR) {
    tags_rel <- c(tags_rel, "Shock-like:Severe")
    rules_rel <- c(rules_rel, sprintf("MAP=%.1f<=%.1f", MAP, THR_REL$MAP_STR))
  } else if (!is.na(MAP) && MAP <= THR_REL$MAP_SOFT) {
    tags_rel <- c(tags_rel, "Hypotension")
    rules_rel <- c(rules_rel, sprintf("MAP=%.1f<=%.1f", MAP, THR_REL$MAP_SOFT))
  }

  if (!is.na(Lac) && Lac >= THR_REL$LAC_STR) {
    tags_rel <- c(tags_rel, "Hyperlactatemia:Severe")
    rules_rel <- c(rules_rel, sprintf("Lactate=%.2f>=%.2f", Lac, THR_REL$LAC_STR))
  } else if (!is.na(Lac) && Lac >= THR_REL$LAC_SOFT) {
    tags_rel <- c(tags_rel, "Hyperlactatemia")
    rules_rel <- c(rules_rel, sprintf("Lactate=%.2f>=%.2f", Lac, THR_REL$LAC_SOFT))
  }

  if (!is.na(PF)) {
    if (PF < THR_REL$PF_STR) {
      tags_rel <- c(tags_rel, sprintf("Severe hypoxemia (PF<%.0f)", THR_REL$PF_STR))
      rules_rel <- c(rules_rel, sprintf("PF=%.1f<%.1f", PF, THR_REL$PF_STR))
    } else if (PF < THR_REL$PF_MOD) {
      tags_rel <- c(tags_rel, sprintf("Moderate hypoxemia (PF<%.0f)", THR_REL$PF_MOD))
      rules_rel <- c(rules_rel, sprintf("PF=%.1f<%.1f", PF, THR_REL$PF_MOD))
    } else if (PF < THR_REL$PF_MILD) {
      tags_rel <- c(tags_rel, sprintf("Mild hypoxemia (PF<%.0f)", THR_REL$PF_MILD))
      rules_rel <- c(rules_rel, sprintf("PF=%.1f<%.1f", PF, THR_REL$PF_MILD))
    } else if (!is.na(FiO2) && FiO2 <= THR_REL$FIO2_ROOM) {
      tags_rel <- c(tags_rel, "Room-air oxygenation")
      rules_rel <- c(rules_rel, sprintf("PF=%.1f high & FiO2=%.1f<=%.1f", PF, FiO2, THR_REL$FIO2_ROOM))
    }
  }

  if (!is.na(FiO2)) {
    if (FiO2 >= THR_REL$FIO2_VHIGH) {
      tags_rel <- c(tags_rel, sprintf("Very high FiO2 need (>=%.0f%%)", THR_REL$FIO2_VHIGH))
      rules_rel <- c(rules_rel, sprintf("FiO2=%.1f>=%.1f", FiO2, THR_REL$FIO2_VHIGH))
    } else if (FiO2 >= THR_REL$FIO2_HIGH) {
      tags_rel <- c(tags_rel, sprintf("High FiO2 need (>=%.0f%%)", THR_REL$FIO2_HIGH))
      rules_rel <- c(rules_rel, sprintf("FiO2=%.1f>=%.1f", FiO2, THR_REL$FIO2_HIGH))
    } else if (FiO2 >= THR_REL$FIO2_MOD) {
      tags_rel <- c(tags_rel, sprintf("Moderate FiO2 need (>=%.0f%%)", THR_REL$FIO2_MOD))
      rules_rel <- c(rules_rel, sprintf("FiO2=%.1f>=%.1f", FiO2, THR_REL$FIO2_MOD))
    } else if (FiO2 >= THR_REL$FIO2_MILD) {
      tags_rel <- c(tags_rel, sprintf("Mild FiO2 need (>=%.0f%%)", THR_REL$FIO2_MILD))
      rules_rel <- c(rules_rel, sprintf("FiO2=%.1f>=%.1f", FiO2, THR_REL$FIO2_MILD))
    }
  }

  if (!is.na(RR) && RR >= THR_REL$RR_TACHY) {
    tags_rel <- c(tags_rel, sprintf("Tachypnea (RR>=%.0f)", THR_REL$RR_TACHY))
    rules_rel <- c(rules_rel, sprintf("RR=%.1f>=%.1f", RR, THR_REL$RR_TACHY))
  }

  if (!is.na(pH) && !is.na(pCO2)) {
    if (pH < THR_ABS$PH_RESP_LOW && pCO2 > THR_ABS$PCO2_HIGH) tags_rel <- c(tags_rel, "Respiratory acidosis")
    if (pH > THR_ABS$PH_RESP_HIGH && pCO2 < THR_ABS$PCO2_LOW) tags_rel <- c(tags_rel, "Respiratory alkalosis")
  }

  if (!is.na(Cr)) {
    if (Cr >= THR_REL$CR_STR) tags_rel <- c(tags_rel, "Renal dysfunction:Severe")
    else if (Cr >= THR_REL$CR_SOFT) tags_rel <- c(tags_rel, "Renal dysfunction")
  }
  if (!is.na(BUN)) {
    if (BUN >= THR_REL$BUN_STR) tags_rel <- c(tags_rel, "Azotemia:Severe")
    else if (BUN >= THR_REL$BUN_SOFT) tags_rel <- c(tags_rel, "Azotemia")
  }

  hepatic_dys <- (
    (!is.na(Tb) && Tb >= THR_REL$TB_SOFT && (((!is.na(AST) && AST >= THR_REL$AST_MIX) || (!is.na(ALT) && ALT >= THR_REL$ALT_MIX)))) ||
    (((!is.na(AST) && AST >= THR_REL$AST_SOFT) || (!is.na(ALT) && ALT >= THR_REL$ALT_SOFT)) && (!is.na(Tb) && Tb >= THR_REL$TB_MIX))
  )
  hepatic_sev <- (
    (!is.na(Tb) && Tb >= THR_REL$TB_STR) ||
    (((!is.na(AST) && AST >= THR_REL$AST_STR) || (!is.na(ALT) && ALT >= THR_REL$ALT_STR)) && (!is.na(Tb) && Tb >= THR_REL$TB_MIX))
  )
  if (hepatic_sev) tags_rel <- c(tags_rel, "Hepatic dysfunction:Severe")
  else if (hepatic_dys) tags_rel <- c(tags_rel, "Hepatic dysfunction")

  coag_soft_any <- ((!is.na(Plt) && Plt <= THR_REL$PLT_SOFT) || (!is.na(INR) && INR >= THR_REL$INR_SOFT))
  coag_str_any  <- ((!is.na(Plt) && Plt <= THR_REL$PLT_STR) || (!is.na(INR) && INR >= THR_REL$INR_STR))
  if (coag_str_any) {
    tags_rel <- c(tags_rel, "Coagulation dysfunction:Severe")
    if (!is.na(INR) && INR >= THR_REL$INR_STR) tags_rel <- c(tags_rel, "Coagulation dysfunction:Factor:Severe")
  } else if (coag_soft_any) {
    tags_rel <- c(tags_rel, "Coagulation dysfunction")
    if (!is.na(INR) && INR >= THR_REL$INR_SOFT) tags_rel <- c(tags_rel, "Coagulation dysfunction:Factor")
  }

  delta_core <- (!is.na(AST) && AST >= THR_REL$CORE_AST_SOFT) &&
    (!is.na(HCO3) && HCO3 <= THR_REL$HCO3_SOFT) &&
    (!is.na(Lac) && Lac >= THR_REL$LAC_SOFT)

  delta_core_severe <- (!is.na(AST) && AST >= THR_REL$CORE_AST_STR) &&
    (!is.na(HCO3) && HCO3 <= THR_REL$HCO3_STR) &&
    (!is.na(Lac) && Lac >= THR_REL$LAC_STR)

  if (delta_core_severe) tags_rel <- c(tags_rel, "Lactate-core:Severe")
  else if (delta_core) tags_rel <- c(tags_rel, "Lactate-core")

  if (!is.na(Hgb) && Hgb < THR_ABS$HGB_SOFT) tags_rel <- c(tags_rel, "Anemia")
  if (!is.na(WBC) && (WBC < 4 || WBC > 12)) tags_rel <- c(tags_rel, "Leukocyte anomaly")
  if (!is.na(Na) && (Na < THR_ABS$NA_CRISIS_LOW || Na > THR_ABS$NA_CRISIS_HIGH)) tags_rel <- c(tags_rel, "Dysnatremia")
  if (!is.na(Cl) && (Cl < THR_ABS$CL_CRISIS_LOW || Cl > THR_ABS$CL_CRISIS_HIGH)) tags_rel <- c(tags_rel, "Dyschloremia")
  if (!is.na(Glu) && Glu >= THR_ABS$GLU_HIGH) tags_rel <- c(tags_rel, "Hyperglycemia")

  if (!is.na(SpO2) && SpO2 < 90) orth_abs <- c(orth_abs, "PulseOx-Low")
  if (!is.na(PF) && PF < THR_ABS$PF_REFRACTORY) orth_abs <- c(orth_abs, "Refractory-Hypoxemia")
  if (!is.na(Cr) && Cr >= THR_ABS$CR_FAILURE) orth_abs <- c(orth_abs, "Renal-Failure:Severe")
  if ((!is.na(Tb) && Tb >= THR_ABS$TB_FAILURE) ||
      (!is.na(AST) && AST >= THR_ABS$AST_FAILURE) ||
      (!is.na(ALT) && ALT >= THR_ABS$ALT_FAILURE)) orth_abs <- c(orth_abs, "Hepatic-Failure")
  if ((!is.na(Plt) && Plt <= THR_ABS$PLT_FAILURE) || (!is.na(INR) && INR >= THR_ABS$INR_FAILURE)) orth_abs <- c(orth_abs, "Coag-Failure")
  if (!is.na(Hgb) && Hgb <= THR_ABS$HGB_FAILURE) orth_abs <- c(orth_abs, "Anemia-Severe")
  if (!is.na(pH) && (pH <= THR_ABS$PH_CRISIS_LOW || pH >= THR_ABS$PH_CRISIS_HIGH)) orth_abs <- c(orth_abs, "AcidBase-Crisis")

  electrolyte_crisis_n <- 0L
  if (!is.na(Na) && (Na <= THR_ABS$NA_CRISIS_LOW || Na >= THR_ABS$NA_CRISIS_HIGH)) electrolyte_crisis_n <- electrolyte_crisis_n + 1L
  if (!is.na(K)  && (K  <= THR_ABS$K_CRISIS_LOW  || K  >= THR_ABS$K_CRISIS_HIGH))  electrolyte_crisis_n <- electrolyte_crisis_n + 1L
  if (!is.na(Cl) && (Cl <= THR_ABS$CL_CRISIS_LOW || Cl >= THR_ABS$CL_CRISIS_HIGH)) electrolyte_crisis_n <- electrolyte_crisis_n + 1L
  if (electrolyte_crisis_n >= 1L) orth_abs <- c(orth_abs, "Electrolyte-Crisis")
  if (!is.na(Glu) && Glu >= THR_ABS$GLU_SEV) orth_abs <- c(orth_abs, "Glucose-Crisis")

  flags <- list(
    core_soft = any(tags_rel == "Lactate-core"),
    core_str  = any(tags_rel == "Lactate-core:Severe"),

    shock_soft = any(tags_rel %in% c("Shock-like:Mild", "Hypotension")),
    shock_str  = ("Shock-like:Severe" %in% tags_rel),

    lac_soft = ("Hyperlactatemia" %in% tags_rel),
    lac_str  = ("Hyperlactatemia:Severe" %in% tags_rel),

    hep_soft = ("Hepatic dysfunction" %in% tags_rel),
    hep_str  = ("Hepatic dysfunction:Severe" %in% tags_rel),

    coag_soft = any(tags_rel %in% c("Coagulation dysfunction", "Coagulation dysfunction:Factor")),
    coag_str  = any(tags_rel %in% c("Coagulation dysfunction:Severe", "Coagulation dysfunction:Factor:Severe")),

    azot_soft = ("Azotemia" %in% tags_rel),
    azot_str  = ("Azotemia:Severe" %in% tags_rel),
    renal_soft = ("Renal dysfunction" %in% tags_rel),
    renal_str  = ("Renal dysfunction:Severe" %in% tags_rel)
  )

  list(
    tags_rel  = unique(tags_rel),
    rules_rel = unique(rules_rel),
    orth_abs  = unique(orth_abs),
    flags_rel = flags
  )
}

labels_from_tagone <- function(tag_out, r = NULL, THR_A_SPLIT = NULL) {
  tags <- tag_out$tags_rel
  f <- tag_out$flags_rel

  dsub <- delta_subtype_from_flags(f)
  bsub <- beta_subtype_from_flags(f)
  gsub <- if (is_gamma_rel(tags)) "G respiratory" else NA_character_

  assignment_source <- NA_character_
  a_burden_score <- NA_integer_
  a_burden_hits <- ""

  if (!is.na(dsub)) {
    subphenotype <- dsub
    assignment_source <- "D_family"
  } else if (!is.na(bsub)) {
    subphenotype <- bsub
    assignment_source <- "B_family"
  } else if (!is.na(gsub)) {
    subphenotype <- gsub
    assignment_source <- "G_family"
  } else {
    subphenotype <- "A provisional"
    assignment_source <- "A_default"

    MAP  <- if (is.null(r)) NA_real_ else get_num(r, "MAP")
    PF   <- if (is.null(r)) NA_real_ else get_num(r, "PF")
    HCO3 <- if (is.null(r)) NA_real_ else get_num(r, "HCO3")
    Hgb  <- if (is.null(r)) NA_real_ else get_num(r, "Hemoglobin")
    FiO2 <- if (is.null(r)) NA_real_ else get_num(r, "FiO2")

    hits <- character(0)
    if (!is.null(THR_A_SPLIT) && length(THR_A_SPLIT)) {
      if (!is.na(MAP)  && is.finite(THR_A_SPLIT$MAP_A2)  && MAP  <= THR_A_SPLIT$MAP_A2)  hits <- c(hits, "MAP_low")
      if (!is.na(PF)   && is.finite(THR_A_SPLIT$PF_A2)   && PF   <= THR_A_SPLIT$PF_A2)   hits <- c(hits, "PF_low")
      if (!is.na(HCO3) && is.finite(THR_A_SPLIT$HCO3_A2) && HCO3 <= THR_A_SPLIT$HCO3_A2) hits <- c(hits, "HCO3_low")
      if (!is.na(Hgb)  && is.finite(THR_A_SPLIT$HGB_A2)  && Hgb  <= THR_A_SPLIT$HGB_A2)  hits <- c(hits, "HGB_low")
      if (!is.na(FiO2) && is.finite(THR_A_SPLIT$FIO2_A2) && FiO2 >= THR_A_SPLIT$FIO2_A2) hits <- c(hits, "FIO2_high")
    }

    a_burden_score <- as.integer(length(hits))
    a_burden_hits <- paste(hits, collapse = ";")

    if (length(hits) >= THR_A_SPLIT$SCORE_CUT) {
      subphenotype <- "A2 mild burden"
      assignment_source <- "A2_rule"
    } else {
      subphenotype <- "A1 low burden"
      assignment_source <- "A1_rule"
    }
  }

  assigned_subtype <- if (startsWith(subphenotype, "D")) "D"
  else if (startsWith(subphenotype, "B")) "B"
  else if (startsWith(subphenotype, "G")) "G"
  else "A"

  ABGD <- character(0)
  if (!is.na(dsub) || startsWith(subphenotype, "D")) ABGD <- c(ABGD, "D")
  if (!is.na(gsub) || startsWith(subphenotype, "G")) ABGD <- c(ABGD, "G")
  if (!is.na(bsub) || startsWith(subphenotype, "B")) ABGD <- c(ABGD, "B")
  if (length(ABGD) == 0) ABGD <- "A"

  dom <- domains_from_flags(f, startsWith(subphenotype, "G") || !is.na(gsub))

  list(
    tags = tags,
    orth = tag_out$orth_abs,
    ABGD = unique(ABGD),
    assigned_subtype = assigned_subtype,
    subphenotype = subphenotype,
    primary_domain = dom$primary,
    secondary_domain = dom$secondary,
    assignment_source = assignment_source,
    has_respiratory_axis = !is.na(gsub),
    d3_respiratory_overlap = !is.na(dsub) && identical(dsub, "D3 shock/coagulation-overlap") && !is.na(gsub),
    a_burden_score = a_burden_score,
    a_burden_hits = a_burden_hits
  )
}

annotate_table <- function(df, THR_REL, THR_ABS, THR_A_SPLIT) {
  rows <- lapply(seq_len(nrow(df)), function(i) {
    tag_out <- tag_one(df[i, , drop = FALSE], THR_REL, THR_ABS)
    lab <- labels_from_tagone(tag_out, r = df[i, , drop = FALSE], THR_A_SPLIT = THR_A_SPLIT)

    tibble(
      state_id = df$state_id[i],
      tags_main = paste(lab$tags, collapse = "; "),
      orthogonal_types = paste(lab$orth, collapse = "; "),
      ABGD = paste(lab$ABGD, collapse = ","),
      assigned_subtype = lab$assigned_subtype,
      subphenotype = lab$subphenotype,
      primary_domain = lab$primary_domain,
      secondary_domain = lab$secondary_domain,
      assignment_source = lab$assignment_source,
      has_respiratory_axis = as.integer(isTRUE(lab$has_respiratory_axis)),
      d3_respiratory_overlap = as.integer(isTRUE(lab$d3_respiratory_overlap)),
      a_burden_score = ifelse(is.null(lab$a_burden_score), NA_integer_, as.integer(lab$a_burden_score)),
      a_burden_hits = ifelse(is.null(lab$a_burden_hits), "", lab$a_burden_hits),
      tags_top = paste(head(lab$tags, 4), collapse = ", ")
    )
  })

  bind_rows(rows)
}

read_thresholds <- function(path) {
  x <- readr::read_csv(path, show_col_types = FALSE)
  if (!all(c("threshold_name", "threshold") %in% names(x))) stop("annotation_thresholds.csv must contain threshold_name and threshold.", call. = FALSE)
  x$threshold_name <- as.character(x$threshold_name)
  x$threshold <- suppressWarnings(as.numeric(x$threshold))
  if (anyDuplicated(x$threshold_name) || any(!is.finite(x$threshold))) stop("Invalid annotation threshold table.", call. = FALSE)
  x
}
threshold_value_list <- function(tbl, names_required) {
  missing <- setdiff(names_required, tbl$threshold_name)
  if (length(missing)) stop("Missing annotation thresholds: ", paste(missing, collapse = ", "), call. = FALSE)
  z <- tbl$threshold[match(names_required, tbl$threshold_name)]
  names(z) <- names_required
  as.list(z)
}

bundle <- mimic_state_profile_bundle
rep_raw <- as.data.frame(bundle$representative_mean, check.names = FALSE)
if (!"state_id" %in% names(rep_raw)) rep_raw$state_id <- seq_len(nrow(rep_raw))
mean_profiles <- lapply(names(bundle$mean_by_stratum), function(key) {
  p <- strsplit(key, "_", fixed = TRUE)[[1L]]
  x <- as.data.frame(bundle$mean_by_stratum[[key]], check.names = FALSE)
  x$state_id <- seq_len(nrow(x)); x$age <- as.numeric(p[[1L]]); x$sex <- as.numeric(p[[2L]]); x
})
names(mean_profiles) <- names(bundle$mean_by_stratum)

thr <- read_thresholds(Sys.getenv("MIMIC_ANNOT_THRESHOLD_FORMAL", unset = resource_path("annotation_thresholds.csv")))
main_names <- REL_THRESHOLD_NAMES
THR_REL <- threshold_value_list(thr, main_names)
THR_REL$TB_MIX <- mean(c(THR_REL$TB_SOFT, THR_REL$TB_STR))
THR_REL$AST_MIX <- mean(c(THR_REL$AST_SOFT, THR_REL$AST_STR))
THR_REL$ALT_MIX <- mean(c(THR_REL$ALT_SOFT, THR_REL$ALT_STR))
THR_ABS <- threshold_value_list(thr, ABS_THRESHOLD_NAMES)

provisional <- annotate_table(rep_raw, THR_REL, THR_ABS, empty_a_split_thresholds())
a_state_ids <- provisional$state_id[provisional$assigned_subtype == "A"]
a_split_names <- c("MAP_A2", "PF_A2", "HCO3_A2", "HGB_A2", "FIO2_A2")
THR_A_SPLIT <- threshold_value_list(thr, a_split_names)
score <- thr$threshold[match("SCORE_CUT", thr$threshold_name)]
THR_A_SPLIT$SCORE_CUT <- if (length(score) && is.finite(score)) score else 2
THR_A_SPLIT$N_A_REF <- length(a_state_ids)

master_ann <- annotate_table(rep_raw, THR_REL, THR_ABS, THR_A_SPLIT)
state_master <- rep_raw %>% left_join(master_ann, by = "state_id") %>% mutate(
  annotation = sprintf("%s; subphenotype=%s; Domains: %s%s", assigned_subtype, subphenotype, primary_domain, ifelse(secondary_domain != "", paste0(" > ", secondary_domain), ""))
)
readr::write_csv(state_master, "state_annotation.csv")
readr::write_csv(thresholds_to_table(THR_REL, THR_ABS, THR_A_SPLIT), "annotation_thresholds_used.csv")

master_cols <- c("state_id", "tags_main", "orthogonal_types", "ABGD", "assigned_subtype", "subphenotype", "primary_domain", "secondary_domain", "assignment_source", "has_respiratory_axis", "d3_respiratory_overlap", "a_burden_score", "a_burden_hits", "tags_top", "annotation")
master_profiles <- rep_raw %>% rename_with(~ paste0("master_", .x), -state_id)
state_reference_profiles <- bind_rows(lapply(mean_profiles, function(raw) {
  ref_ann <- annotate_table(raw, THR_REL, THR_ABS, THR_A_SPLIT) %>% rename(reference_tags = tags_main, reference_orthogonal = orthogonal_types, reference_ABGD_check = ABGD, reference_assigned_subtype_check = assigned_subtype, reference_subphenotype_check = subphenotype, reference_primary_domain_check = primary_domain, reference_secondary_domain_check = secondary_domain, reference_assignment_source_check = assignment_source, reference_has_respiratory_axis_check = has_respiratory_axis, reference_d3_respiratory_overlap_check = d3_respiratory_overlap, reference_a_burden_score_check = a_burden_score, reference_a_burden_hits_check = a_burden_hits, reference_tags_top = tags_top)
  raw_ref <- raw %>% rename_with(~ paste0("ref_", .x), -c(state_id, age, sex))
  raw_ref %>% left_join(ref_ann, by = "state_id") %>% left_join(state_master %>% select(all_of(master_cols)), by = "state_id") %>% left_join(master_profiles, by = "state_id")
}))
readr::write_csv(state_reference_profiles, "state_reference_profiles.csv")

counts <- state_master %>% count(subphenotype, name = "n_states") %>% arrange(subphenotype)
readr::write_csv(counts, "state_annotation_counts.csv")
message("State annotation completed: ", nrow(state_master), " states.")
