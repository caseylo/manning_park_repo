#### Cleaning error data for Bayesian analysis

## Load libraries

library(tidyverse)

## Read in data

error_data <- read_csv("manning_park_data/relocation_data_long.csv")

model_data <- read_csv("manning_park_data/model_data_filt.csv")

## Check data structure
str(error_data)

error_data <- error_data  %>%
  select(-1)

####

## Making a long dataframe for the full model

## Keep only species used in full model

species_keep <- model_data %>%
  distinct(Species) %>%
  pull(Species)

## Keep only 2025 and Error plots, and only the 30 species from the full model

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

write_csv(error_model_data,"bayesian_analysis_scripts/error_model_data.csv")




