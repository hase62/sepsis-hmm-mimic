
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
subphenotype_levels_order <- subphenotype_palette_levels <- c(
  "A1 low burden",
  "A2 mild burden",
  "B1 renal mild",
  "B2 renal severe",
  "G respiratory",
  "D1 lactate/shock",
  "D2 hepatic",
  "D3 shock/coagulation-overlap"
)

subphenotype_colors <- c(
  "A1 low burden"       = "#F6C65B",
  "A2 mild burden"        = "#E69F00",
  "B1 renal mild"              = "#56B4E9",
  "B2 renal severe"            = "#009E73",
  "G respiratory"             = "#CC79A7",
  "D1 lactate/shock"           = "#F0E442",
  "D2 hepatic"                = "#0072B2",
  "D3 shock/coagulation-overlap"             = "#D55E00"
)

subphenotype_linetypes <- c(
  "solid",
  "dashed",
  "dotdash",
  "dotted",
  "longdash",
  "twodash",
  "solid",
  "dashed",
  "dotted"
)
names(subphenotype_linetypes) <- subphenotype_palette_levels

subphenotype_shapes <- c(
  16,
  17,
  15,
  3,
  7,
  8,
  18,
  1,
  2
)
names(subphenotype_shapes) <- subphenotype_palette_levels

subphenotype_colors["Other/NA"] <- "#bdbdbd"
