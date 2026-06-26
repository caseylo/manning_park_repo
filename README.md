# Assessing Upslope Range Shifts with Relocation Error in E.C. Manning Provincial Park- A Vegetation Resurvey Analysis

This repository contains the code, processed data, figures, model outputs, and summary tables for a vegetation resurvey analysis in E.C. Manning Provincial Park, British Columbia. The project examines whether plant species occurrence patterns have shifted along elevation gradients over time and whether those patterns are robust to plot relocation error.

## Project overview

Montane plant communities are expected to respond to climate change, but observed range shifts can vary among species and functional groups. This project compares historical vegetation surveys with a 2025 resurvey to assess changes in species occurrence across elevation. It also includes a relocation-error analysis using 2019 plots that were resurveyed in 2025 both at the best-estimated original location and at intentionally mis-relocated comparison plots.

The main goals are to evaluate:

* evidence for upslope shifts in plant species occurrence,
* differences among herbs, shrubs, and trees,
* whether relocation error affects estimates of vegetation change.

## Scripts

The scripts are organized into three main folders:

* `scripts/01_cleaning/`: cleans and formats historical, 2025, and relocation-error datasets.
* `scripts/02_data_exploration/`: examines species richness and candidate model structures.
* `scripts/03_analysis/`: runs the AIC-selected final range-shift models, Bayesian range shift analysis, relocation-error analyses, diagnostics, and figure creation.

## Outputs

Final project outputs are stored in `outputs/`:

* `outputs/figures/`: elevation-response plots, functional-group plots, NMDS figures, diagnostics, and maps.
* `outputs/models/`: saved model objects from the main range-shift and relocation-error analyses.
* `outputs/tables/`: summary tables for richness and diversity comparisons.

## Main analyses

The main range-shift analysis models species presence as a function of time, elevation, and their interaction, while accounting for several random effects. A functional-group model compares responses among herbs, shrubs, and trees.

The relocation-error analysis compares 2019 plots, accurately relocated 2025 plots (2025correct), and intentionally relocated 2025 plots (2025error) to test whether relocation error changes community composition or diversity enough to affect interpretation of temporal change.

## Notes

This repository is intended as a project archive for supervisors, collaborators, and readers who want to understand the organization and main analyses. It is not intended as a full public data-release or step-by-step reproduction guide.

## Creator & Maintainer

Casey Lo //
University of British Columbia //
B.Sc Honours Ecology and Environmental Science
