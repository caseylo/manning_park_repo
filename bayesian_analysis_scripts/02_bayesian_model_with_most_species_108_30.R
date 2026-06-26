#### Bayesian relocation error-informed model will ALL species from relocation error data (n=108)
#### Full model data is the original species subset (n=30)

## Load libraries

library(tidyverse)
library(brms)
library(posterior)
library(ggeffects)
library(ggplot2)


set.seed(1234)

## Read in data

error_data <- read_csv("data/processed/relocation_data_long.csv")

model_data <- read_csv("data/processed/model_data_filt.csv")

## Read in models to avoid re-running if needed

#bayesmod_error <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_error_108_spp.rds")
#bayesmod_full <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_full_108_30_spp.rds")
#bayesmod_full_no_prior <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_full_no_prior_30_spp.rds")

## Check data structure
str(error_data)

error_data <- error_data  %>%
  select(-1)

####

## Making a long dataframe for the full model

## Keep only species used in full model

species_keep <- error_data %>%
  distinct(Species) %>%
  pull(Species)

length(unique(error_data$Species))
#[1] 108

## Keep only 2025 and Error plots, and all species in error data

error_data <- error_data %>%
  filter(Treatment %in% c("2025", "Error")) %>%
  filter(Species %in% species_keep) %>%
  mutate(present = 1)

## Make complete plot x species dataset

error_grid <- error_data %>%
  distinct(
    PlotNumber, Site, ProjectID, Treatment, Year,
    SiteSurveyor, RelocationConfidence,
    Longitude, Latitude, Elevation
  ) %>%
  crossing(Species = species_keep)

## Add present/absent column

error_model_data <- error_grid %>%
  left_join(
    error_data %>%
      select(PlotNumber, Site, Treatment, Species, present),
    by = c("PlotNumber", "Site", "Treatment", "Species")
  ) %>%
  mutate(
    present = replace_na(present, 0)
  )

table(error_model_data$Treatment, error_model_data$present)
length(unique(error_model_data$Species))

## Add functional groups to error data in case needed later

species_fg <- model_data %>%
  distinct(Species, Functional_group)

error_model_data <- error_model_data %>%
  left_join(species_fg, by = "Species")

## Save as csv

write_csv(error_model_data,"data/processed/error_model_data_all_species.csv")

####

## Re-run model to set priors on new error data

error_model_data <- read_csv("data/processed/error_model_data_all_species.csv")

## Check data structure
str(error_model_data)
str(model_data)

## Prepare error model data

error_model_data <- error_model_data %>%
  mutate(
    Species = factor(Species),
    PlotNumber = factor(PlotNumber),
    Treatment = factor(Treatment, levels = c("2025", "Error")),
    Elevation_sc = as.numeric(scale(Elevation)),
    Lat_sc = as.numeric(scale(Latitude))
  )

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

## Running Bayesian model

bayesmod_error <- brm(
  present ~ Treatment * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = error_model_data,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 1234)

summary(bayesmod_error)
pp_check(bayesmod_error)

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

bayesmod_full_no_prior <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_full_no_prior_30_spp.rds")

summary(bayesmod_full_no_prior)
pp_check(bayesmod_full_no_prior)
prior_summary(bayesmod_full_no_prior)

####

## Save models

saveRDS(bayesmod_error, "outputs/models/bayesian_analysis_models/bayesmod_error_all_108_spp.rds")
saveRDS(bayesmod_full, "outputs/models/bayesian_analysis_models/bayesmod_full_108_30_spp.rds")

####