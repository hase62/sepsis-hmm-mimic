# Longitudinal HMM analysis of peri-admission ventilated sepsis in MIMIC-IV

This repository contains the MIMIC-IV analysis code for a longitudinal hidden Markov model study of serial physiological subphenotypes in sepsis with peri-admission invasive mechanical ventilation. The primary cohort is restricted to qualifying sepsis ICU stays with invasive ventilation overlapping the interval from 6 hours before to 6 hours after ICU admission.

The code reproduces cohort construction, preprocessing, HMM state-number selection and fitting, deterministic clinical summarization of the fitted state space, and downstream longitudinal analyses. The repository is intended to reproduce the MIMIC-IV component of the study rather than to provide a general-purpose sepsis phenotyping model for all ICU patients with sepsis.

Publication-layout scripts for the final manuscript figures and supplementary files are not included. Individual analysis scripts may generate analysis-specific figures. Code that depends on the private Japanese hospital cohort is not included.

## Data requirements

Obtain credentialed access to MIMIC-IV v3.1 from PhysioNet. Keep the original `hosp/` and `icu/` directories intact. The only data supplied by the user are the original MIMIC-IV files; study-specific mappings, variable definitions, normalization choices, and annotation thresholds are stored in `resources/`.

Set `MIMIC_ROOT` to the directory containing `hosp/` and `icu/`.

No patient-level MIMIC-IV data are distributed with this repository.

## R requirements

Required packages are listed in `R/00_packages.R`. The principal computational dependencies are `depmixS4`, `mice`, `bestNormalize`, `data.table`, `tidyverse` components, and `duckdb`/`DBI` for the patient-day SOFA workflow. Optional packages only add selected output formats or plots.

## Repository layout

- `R/01`-`05`: cohort extraction and analytical matrix construction
- `R/06`-`10`: split-safe preprocessing and HMM state-number selection
- `R/11`-`14`: full-cohort HMM estimation
- `R/15`-`19`: state-profile reconstruction and deterministic clinical annotation
- `R/20`-`34`: state, outcome, filtering, SOFA, organ-support, and transition analyses
- `resources/`: frozen study definitions and mappings
- `check_inputs.R`: validates the MIMIC-IV directory and repository resources
- `run_pipeline.R`: convenience entry point for sequential stages
- `print_hpc_commands.R`: prints commands for parallel model-selection and multi-start fitting

## Setup

Run commands from the repository root.

```bash
export MIMIC_ROOT=/path/to/mimiciv/3.1
export SEPSIS_HMM_CODE_DIR="$PWD/R"
export SEPSIS_HMM_RESOURCE_DIR="$PWD/resources"
export SEPSIS_HMM_OUTPUT_DIR=/path/to/output
mkdir -p "$SEPSIS_HMM_OUTPUT_DIR"
```

On Windows, set the same environment variables in R or PowerShell using absolute paths.

Validate inputs before starting:

```bash
Rscript check_inputs.R "$MIMIC_ROOT"
```

## 1. Cohort and analytical matrix

```bash
Rscript run_pipeline.R prepare-data "$MIMIC_ROOT" "$SEPSIS_HMM_OUTPUT_DIR"
```

This stage creates the peri-admission invasively ventilated sepsis cohort, laboratory-anchored analytical rows, and the finalized HMM input matrix.

## 2. HMM state-number selection

The analysis evaluates 5-65 states over 20 patient-level development/held-out splits. This requires 1,220 HMM fits and is intended for parallel or batch execution.

Create one preprocessing cache for each split:

```bash
for i in $(seq 1 20); do
  Rscript "$SEPSIS_HMM_CODE_DIR/08_prepare_model_selection_split.R" "$i"
done
```

Fit each candidate state number and split:

```bash
Rscript "$SEPSIS_HMM_CODE_DIR/09_fit_model_selection_candidate.R" <K> <iteration>
```

After all jobs finish:

```bash
export MIMIC_K_MIN=5
export MIMIC_K_MAX=65
Rscript "$SEPSIS_HMM_CODE_DIR/10_select_state_dimension.R" 20
```

The primary selection rule is the one-standard-error rule applied to mean held-out three-step reconstruction error across the 20 patient-level splits.

## 3. Final HMM

Prepare the full-cohort preprocessing cache:

```bash
Rscript "$SEPSIS_HMM_CODE_DIR/11_prepare_full_cohort.R"
```

Fit 20 independent starts for the selected 60-state model, preferably in parallel:

```bash
Rscript "$SEPSIS_HMM_CODE_DIR/12_fit_full_cohort_start.R" 60 <start_id>
```

Select the valid fit with the highest full-cohort log likelihood:

```bash
Rscript "$SEPSIS_HMM_CODE_DIR/13_select_final_model.R" 60 20
```

This stage writes `model_config.R` in the analysis output directory.

## 4. State annotation

```bash
Rscript "$SEPSIS_HMM_CODE_DIR/19_run_annotation.R"
```

The eight reported subphenotypes are deterministic clinical summaries of the fitted 60-state representation using the frozen thresholds in `resources/annotation_thresholds.csv`.

## 5. Downstream analyses

```bash
Rscript "$SEPSIS_HMM_CODE_DIR/34_run_downstream_analyses.R"
```

The downstream scripts generate analysis tables and analytical figures. They do not assemble the publication-specific final figure or supplementary-file package.

## Batch command list

```bash
Rscript print_hpc_commands.R "$MIMIC_ROOT" "$SEPSIS_HMM_OUTPUT_DIR" > hmm_jobs.sh
```

Adapt the printed commands to the local scheduler or batch system.

## Reproducibility scope

This repository covers the MIMIC-IV analyses. Japanese-cohort analyses, publication-layout code, and internal development utilities are not included.
