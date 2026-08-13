
source(file.path(Sys.getenv("SEPSIS_HMM_CODE_DIR", unset = file.path(getwd(), "R")), "00_setup.R"))
source(code_path("00_packages.R"))
source(code_path("14_load_final_model.R"))
source(code_path("15_reconstruct_state_profiles.R"))
source(code_path("16_annotate_states.R"))
