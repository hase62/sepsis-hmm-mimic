# Frozen study resources

The public repository bundles the fixed study resources used by the analysis. End users do not create these files.

Required files:

- `antibiotics_list.txt`: antibiotic terms used to identify suspected infection.
- `mimic_item_map.csv`: MIMIC item identifiers mapped to analytical variable names; requires `itemid` and `matched_test` columns.
- `observation_variables.csv`: ordered 38-variable HMM observation set.
- `normalization_families.csv`: frozen normalization family for each HMM observation variable.
- `annotation_thresholds.csv`: frozen deterministic state-annotation thresholds; requires `threshold_name` and `threshold` columns.
