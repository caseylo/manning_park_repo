#### Creating dataframe with binned data

## load libraries

library(tidyverse)
library(dplyr)
library(lubridate)

## Read in select data

BEC_clean <- read_csv("data/cleaned/BEC_select.csv")
species_names <- read_csv("data/cleaned/master_species_list.csv")
data_2025_clean <- read_csv("data/cleaned/data_2025_select.csv")

####

## Getting rid of space

data_2025_clean <- data_2025_clean %>%
  mutate(RelocationConfidence =
           recode(RelocationConfidence, "High/ Med-High" = "High/Med-High"))

####

## Different filters to create data for sensitivity analysis based on reocation confidence

## High, High/Med-High, Med-High, Med confidence only:

#data_2025_clean<- data_2025_clean %>%
#  filter(RelocationConfidence != "Low")
#length(unique(data_2025_clean$PlotNumber))
#[1] 49
#BEC_clean <- BEC_clean %>%
#  filter(PlotNumber %in% data_2025_clean$PlotNumber)
#length(unique(BEC_clean$PlotNumber))
#[1] 49

## High, High/Med-High, Med-High confidence only:

#medhigh <- c("High", "Med-High", "High/Med-High")
#data_2025_clean<- data_2025_clean %>%
#  filter(RelocationConfidence %in% medhigh)
#length(unique(data_2025_clean$PlotNumber))
#[1] 42
#BEC_clean <- BEC_clean %>%
#  filter(PlotNumber %in% data_2025_clean$PlotNumber)
#length(unique(BEC_clean$PlotNumber))
#[1] 42

## High, High/Med-High Confidence only:

#high <- c("High", "High/Med-High")
#data_2025_clean<- data_2025_clean %>%
#  filter(RelocationConfidence %in% high)
#length(unique(data_2025_clean$PlotNumber))
#[1] 25
#BEC_clean <- BEC_clean %>%
#  filter(PlotNumber %in% data_2025_clean$PlotNumber)
#length(unique(BEC_clean$PlotNumber))
#[1] 25

####

## Adding a Year only column

data_2025_clean <- data_2025_clean %>%
  mutate(Year = year(as.Date(Date2025))) %>%
  select(-Date2025, -DateOriginal) %>%
  relocate(Year, .after = ProjectID)# %>%
  #select(-1,-2)

## BEC data now

BEC_clean <- BEC_clean %>%
  mutate(Year = year(as.Date(Date))) %>%
  select(-Date) %>%
  relocate(Year, .after = ProjectID)

####

## Filtering out years >2000 in BEC data

BEC_clean <- BEC_clean %>%
  filter(Year <= 2000)

## Keeping only 2025 plots found in BEC data

data_2025_clean <- data_2025_clean %>%
  filter(PlotNumber %in% BEC_clean$PlotNumber)

####

## Now I need to filter out all the species with a combined occurrences less than 20

## BEC: one record per plot-species (Some species have been duplicated in plots)
## The plots causing inflated occurences:

BEC_clean %>%
  count(PlotNumber, Species) %>%
  filter(n > 1) %>%
  arrange(desc(n))

##2025: one record per plot-species
data_2025_plot_species <- data_2025_clean %>%
  distinct(PlotNumber, ProjectID, Elevation, Species)

##BEC: one record per plot-species
BEC_plot_species <- BEC_clean %>%
  distinct(PlotNumber, ProjectID, Elevation, Species)

##Count occurrences as number of unique plots
species_2025 <- data_2025_plot_species %>%
  count(Species, name = "occurences_2025")

species_BEC <- BEC_plot_species %>%
  count(Species, name = "occurences_BEC")

species_both <- inner_join(species_2025, species_BEC, by = "Species") %>%
  mutate(total_occurences = occurences_2025 + occurences_BEC)
head(species_both)

##Filtering for species with >=20 occurences and at lease >=7 in either BEC or 2025 survey

species_filtered <- species_both %>%
  filter(total_occurences >= 20, occurences_2025 >= 7, occurences_BEC >= 7)
species_filtered

#### Filter BEC and 2025 data for only species from this list:

BEC_clean_filt <- BEC_clean %>%
  filter(Species %in% species_filtered$Species)

data_2025_clean_filt <- data_2025_clean %>%
  filter(Species %in% species_filtered$Species)

####Making final species list:

final_species_list <- species_filtered %>%
  select(-occurences_2025,
         -occurences_BEC,
         -total_occurences) %>%
  mutate(Species = as.character(Species))

####

####Making 2025 data long

plot_grid <- data_2025_clean_filt %>%
  distinct(PlotNumber, ProjectID, Latitude, Longitude, Elevation, Year, SiteSurveyor, RelocationConfidence) %>%
  crossing(final_species_list)   # adds every Species to every plot

##Join observed data onto the grid + create presence/absence
data_2025_long <- plot_grid %>%
  left_join(data_2025_clean_filt,
    by = c(
      "PlotNumber",
      "ProjectID",
      "Latitude",
      "Longitude",
      "Elevation",
      "Year",
      "Species",
      "SiteSurveyor",
      "RelocationConfidence")) %>%
  mutate(present = if_else(is.na(Cover), 0L, 1L), Cover   = replace_na(Cover, 0)) %>%
  left_join(species_names, by = c("Species" = "BC_species_code"))

data_2025_long <- data_2025_long %>%
  select(-Functional_group.x,
         -Latin_name) %>%
  rename(Functional_group = Functional_group.y)

#write.csv(data_2025_long, file = "data/cleaned/data_2025_filt_long.csv")

####

## Making BEC data long

## Combining duplicated species in plots/ collapse BEC observations to one row per plot-species

BEC_obs_plot_species <- BEC_clean_filt %>%
  group_by(PlotNumber,
           ProjectID,
           Latitude,
           Longitude,
           Elevation,
           Year,
           Species,
           SiteSurveyor) %>%
  summarise(Cover = sum(Cover, na.rm = TRUE), .groups = "drop")

BEC_plot_grid <- BEC_clean_filt %>%
  distinct(PlotNumber, ProjectID, Latitude, Longitude, Elevation, Year, SiteSurveyor) %>%
  crossing(final_species_list)   # adds every Species to every plot

##Join observed data onto the grid + create presence/absence

BEC_data_long <- BEC_plot_grid %>%
  left_join(
    BEC_obs_plot_species,
    by = c(
      "PlotNumber",
      "ProjectID",
      "Latitude",
      "Longitude",
      "Elevation",
      "Year",
      "Species",
      "SiteSurveyor")) %>%
  mutate(present = if_else(is.na(Cover), 0L, 1L), Cover   = replace_na(Cover, 0)) %>%
  left_join(species_names, by = c("Species" = "BC_species_code"))

BEC_data_long <- BEC_data_long %>%
  select(-Latin_name)

#write.csv(BEC_data_long, file = "data/cleaned/BEC_data_filt_long.csv")

####

## Converting to years since first survey

min_year <- min(
  c(BEC_data_long$Year, data_2025_long$Year),
  na.rm = TRUE
)

BEC_data_long <- BEC_data_long %>%
  mutate(YearsSinceStart = Year - min_year) %>%
  mutate(SiteSurveyor = if_else(ProjectID == "BEC-Vancouver", "BEC-Vancouver", SiteSurveyor)) %>%
  mutate(RelocationConfidence = NA_character_)

data_2025_long <- data_2025_long %>%
  mutate(YearsSinceStart = Year - min_year) %>%
  mutate(ProjectID = if_else(Year == 2025, "Manning Park Project 2025", ProjectID)) %>%
  mutate(SiteSurveyor = if_else(Year == 2025, "CML RNW", SiteSurveyor))

## Combine BEC and 2025 data into 1 df

model_data <- bind_rows(
  BEC_data_long %>%
    select(Species, Functional_group, ProjectID, PlotNumber, SiteSurveyor, Year, YearsSinceStart, RelocationConfidence,
           Latitude, Longitude, Elevation, present, Cover),
  data_2025_long %>%
    select(Species, Functional_group, ProjectID, PlotNumber, SiteSurveyor, Year, YearsSinceStart, RelocationConfidence,
           Latitude, Longitude, Elevation, present, Cover))

## Making an object with the relocation confidence of each plot from 2025 to attach to the relocation condidence column of BEC plots

rc_2025 <- model_data %>%
  filter(Year == 2025, !is.na(RelocationConfidence)) %>%
  select(PlotNumber, RelocationConfidence) %>%
  distinct(PlotNumber, .keep_all = TRUE)

## Fill historical RelocationConfidence NAs using the 2025 value for the same PlotNumber
model_data <- model_data %>%
  left_join(rc_2025, by = "PlotNumber", suffix = c("", "_2025")) %>%
  mutate(
    RelocationConfidence = if_else(
      Year != 2025 & is.na(RelocationConfidence),
      RelocationConfidence_2025,
      RelocationConfidence)) %>%
  select(-RelocationConfidence_2025)

####

model_data <- model_data %>%
  mutate(Time = case_when(
    Year <= 2000 ~ "historical",
    Year == 2025 ~ "present"
  )) %>%
  mutate(Time = factor(Time, levels = c("historical", "present")))

## Creating bins for historical (1976-2000) and new data (2025)

## Save full model_data
#write.csv(model_data, file = "data/processed/model_data_filt.csv")

## Save sensitivity data

## High, Med-High, Med confidence only:
#write.csv(model_data, file = "data/processed/sensitivity_analysis_data/model_data_conf.csv")

## High, High/Med-High, Med-High confidence only
#write.csv(model_data, file = "data/processed/sensitivity_analysis_data/model_data_medhigh_conf.csv")

## High, High/Med-High Confidence only 
#write.csv(model_data, file = "data/processed/sensitivity_analysis_data/model_data_high_conf.csv")
