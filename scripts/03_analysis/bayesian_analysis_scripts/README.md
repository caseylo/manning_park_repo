# Bayesian Analysis Scripts

This folder contains the Bayesian modelling scripts for the Manning Park vegetation resurvey analysis. These models were run as a follow-up to the original GLMM analysis to test how much the results change when relocation error is included through informative priors.

The general model structure is:

```r
present ~ Time * poly(Elevation_sc, 2, raw = TRUE) +
  Lat_sc +
  (poly(Elevation_sc, 2, raw = TRUE) | Species) +
  (1 | PlotNumber)
```
## General workflow

For each analysis, the workflow was:

1. Run a relocation error model using the 2025correct vs 2025error data.
2. Extract priors from the relocation error model.
3. Run the Bayesian model with those priors on the full dataset of historical vs. present data.
4. Run the same full Bayesian model without priors in order to compare.

The relocation error model is used to estimate how much plot relocation uncertainty could affect the main model terms. Those estimates are then used as priors in the full Bayesian model.

## Model sets

### 1. Original data for all Bayesian models

Script: `01_model_with_filtered_species_30_30.R`

This script uses the same 30-species dataset as the original GLMM analysis for both the relocation error model and the full Bayesian model.

* Relocation error model: 30 species, 540 observations
* Full Bayesian model with priors: 30 species, 2,040 observations
* Full Bayesian model without priors: 30 species, 2,040 observations

This is the closest Bayesian version of the original GLMM analysis.

### 2. Original full model data, broader relocation error model

Script: `02_bayesian_model_with_most_species_108_30.R`

This script uses all species available in the relocation error dataset to estimate the priors, but keeps the full Bayesian model limited to the original 30-species dataset.

* Relocation error model: 108 species, 1,944 observations
* Full Bayesian model with priors: 30 species, 2,040 observations
* Full Bayesian model without priors: 30 species, 2,040 observations

This was run because the priors from the 30-species relocation error model were not very strong. Using more species in the relocation error model gives more information for estimating relocation-error effects.

### 3. All historical species for full Bayesian models

Script: `03_bayesian_model_on_all_species_108_260.R`

This script uses the 108-species relocation error model to set priors, then runs the full Bayesian model using all species (n=260) in the broader historical dataset.

* Relocation error model: 108 species, 1,944 observations
* Full Bayesian model with priors: 260 species, 17,680 observations
* Full Bayesian model without priors: 260 species, 17,680 observations

This was run to check whether the Bayesian results are similar when the full model includes all available species instead of only the 30 species used in the original GLMM.
