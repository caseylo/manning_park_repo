#### Bayesian relocation error-informed model will 30 species from relocation error data (n=30)
#### Full model data is includes ALL species from 2025 data (n=260)

## Load libraries

library(tidyverse)
library(brms)
library(posterior)
library(ggeffects)
library(ggplot2)


set.seed(1234)

## Read in data

model_data <- read_csv("data/processed/model_data_all_species.csv")

## Read in models to avoid re-running if needed

bayesmod_error <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_error_108_spp.rds")

## Prepare the full model data

model_data <- model_data %>%
  mutate(
    Species = factor(Species),
    ProjectID = factor(ProjectID),
    PlotNumber = factor(PlotNumber),
    SiteSurveyor = factor(SiteSurveyor),
    Elevation_sc = as.numeric(scale(Elevation)),
    RelocationConfidence = factor(RelocationConfidence),
    Time_sc = as.numeric(scale(YearsSinceStart)),
    Lat_sc = as.numeric(scale(Latitude)),
    Lon_sc = as.numeric(scale(Longitude)),
    Time = factor(Time, levels = c("historical", "present")),
    Functional_group = factor(Functional_group)
  )

####

## Extract priors

error_priors <- fixef(bayesmod_error)[c(
  "TreatmentError",
  "TreatmentError:polyElevation_sc2rawEQTRUE1",
  "TreatmentError:polyElevation_sc2rawEQTRUE2"
), ]

error_priors

## Saving each mean and SD as  objects

trt_mu  <- error_priors["TreatmentError", "Estimate"]
trt_sd  <- error_priors["TreatmentError", "Est.Error"]

elev_mu <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE1", "Estimate"]
elev_sd <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE1", "Est.Error"]

elev2_mu <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE2", "Estimate"]
elev2_sd <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE2", "Est.Error"]

## Get coefficient names and set priors

get_prior(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = bernoulli(link = "logit")
)

priors_error <- c(
  set_prior(paste0("normal(", trt_mu, ", ", trt_sd, ")"),
            class = "b", coef = "Timepresent"),
  
  set_prior(paste0("normal(", elev_mu, ", ", elev_sd, ")"),
            class = "b", coef = "Timepresent:polyElevation_sc2rawEQTRUE1"),
  
  set_prior(paste0("normal(", elev2_mu, ", ", elev2_sd, ")"),
            class = "b", coef = "Timepresent:polyElevation_sc2rawEQTRUE2")
)

priors_error

####

## Run full model with priors

bayesmod_full <- brm(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = bernoulli(link = "logit"),
  prior = priors_error,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 1234
)

summary(bayesmod_full)
pp_check(bayesmod_full)
prior_summary(bayesmod_full)

####

## Bayesian with no priors

bayesmod_full_no_prior <- brm(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 1234
)

summary(bayesmod_full_no_prior)
pp_check(bayesmod_full_no_prior)
prior_summary(bayesmod_full_no_prior)

####

## Save models

saveRDS(bayesmod_full, "outputs/models/bayesian_analysis_models/bayesmod_full_108_260_spp.rds")
saveRDS(bayesmod_full_no_prior, "bayesian_analysis_scripts/bayesmod_full_no_prior_108_260_spp.rds")

####