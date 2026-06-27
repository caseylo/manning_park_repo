# Data Analysis Scripts

This README gives a brief overview of the scripts used for the Manning Park vegetation resurvey analysis. The scripts are grouped by workflow stage: cleaning, exploration, and final analysis.

## Suggested reading order

For a general understanding of the analysis workflow, start with:

1. `filtered_data_for_analysis.R`
2. `categorical_glmm_model_selection.R`
3. `exploring_best_categotical_model.R`
4. `permanova_relocation_error_analysis.R`
5. `relocation_error_sr_comparisons.R`
6. bayesian_analysis_script folder


## `01_cleaning/`

### `master_species_list.R`

Builds the master species list used across the project. This includes species codes, Latin names, and functional groups.

### `historical_data_cleaning.R`

Cleans the historical vegetation data and the 2025 resurvey data.

### `cleaning_data_for_analysis.R`

Creates the main dataset for the continuous-time analysis.

### `filtered_data_for_analysis.R`

Creates the main historical vs. present dataset. This is the dataset used for the main range-shift GLMM, where time is grouped as `historical` or `present`.

### `cleaning_data_for_relocation_analysis.R`

Prepares the relocation-error dataset. It combines the 2019 plots, accurately relocated 2025 plots, and 2025 error plots.

## `02_data_exploration/`

### `categorical_glmm_model_selection.R`

Compares candidate GLMMs for the main historical vs. present analysis. This script is used to choose the final model structure.

### `continuous_glmm_model_selection.R`

Compares candidate GLMMs using survey year (time) as a continuous variable. This is used as a supporting analysis alongside the categorical model.

### `comparing_plot_and_project_sr.R`

Compares species richness across plots and survey projects. This gives a basic summary of richness patterns before modelling began to ensure BEC data was not systemically different than 2025 data in terms of richness.

### `historical_plot_relocation_distances.R`

Calculates distances between historical plot coordinates and the 2025 resurvey locations.

### `relocation_error_plot_distances_2019.R`

Calculates distances among the 2019, 2025, and error plot locations used in the relocation-error analysis to identify how far, on average, we were able to relocate plots.

## `03_analysis/`

### `exploring_best_categotical_model.R`

Runs the main historical vs. present range-shift model. This script also includes model diagnostics, functional-group analyses, relocation-confidence sensitivity checks, and final figures.

### `exploring_best_continuous_model.R`

Runs the continuous-time version of the range-shift analysis and creates the related figures.

### `permanova_relocation_error_analysis.R`

Runs the community-composition analysis for the relocation-error experiment. This includes PERMANOVA, pairwise tests, dispersion checks, and NMDS figures.

### `relocation_error_sr_comparisons.R`

Compares species richness and diversity among 2019, 2025, and error plots.

### `full_model_on_relocation_data.R`

Fits the full categorical range shift model to the relocation-error dataset. This checks whether the error plots produce similar patterns to the accurately relocated plots.
