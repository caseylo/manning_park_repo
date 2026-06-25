## Load packages
library(tidyverse)

data <- read.csv("manning_park_data/model_data.csv")
names <- read.csv("manning_park_data/cleaning_data/species_code_and_latin_names_2025.csv")

names <- names %>%
  rename(Species = BC_species_code)

####

## Calculating the SD of pres/abs for Error vs. 2025

# Load data
reloc <- read.csv("manning_park_data/relocation_data_long.csv")

# Keep only 2025 and Error plots
reloc_2025_error <- reloc %>%
  filter(Treatment %in% c("2025", "Error"))

# Make observed presence/absence data
reloc_pa_obs <- reloc_2025_error %>%
  group_by(Site, Treatment, Species) %>%
  summarise(
    presence = ifelse(any(Cover > 0), 1, 0),
    .groups = "drop")

# Make sure absences are included
# This creates every Site x Treatment x Species combination
reloc_pa <- reloc_2025_error %>%
  distinct(Site, Treatment) %>%
  crossing(Species = unique(reloc_2025_error$Species))%>%
  left_join(reloc_pa_obs, by = c("Site", "Treatment", "Species")) %>%
  mutate(
    presence = replace_na(presence, 0)
  )

reloc_pa <- reloc_pa %>%
  left_join(names)

# Pair 2025 and Error observations within each Site x Species
reloc_site_species_diff <- reloc_pa %>%
  select(Site, Treatment, Functional_group, Species, presence) %>%
  pivot_wider(
    names_from = Treatment,
    values_from = presence
  ) %>%
  mutate(
    presence_difference = Error - `2025`,
    abs_presence_difference = abs(presence_difference)
  )

# Calculate SD of relocation error across Site x Species comparisons
SDerror_site_species <-  reloc_site_species_diff %>%
  group_by(Functional_group) %>%
  summarise(sd = sd(presence_difference,
  na.rm = TRUE))

SDerror_site_species
#[1] 0.2992431``

