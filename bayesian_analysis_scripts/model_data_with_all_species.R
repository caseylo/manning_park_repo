#### Script for cleaning data for analysis

## Filtering data

library(tidyverse)
library(dplyr)
library(lubridate)

## Read in data

BEC_data <- read_csv("manning_park_data/cleaning_data/BEC_clean.csv")
data_2025 <- read_csv("manning_park_data/cleaning_data/data_2025_clean.csv")
species_names <- read_csv("manning_park_data/cleaning_data/master_species_list.csv")
model_data_filt <- read_csv("manning_park_data/model_data_filt.csv")

####

## Working with 2025 data first

## Filter out with "spp."

data_2025 <- data_2025 %>%
  filter(!grepl("SP\\.?$", Species)) %>%
  filter(ErrorPlot != "Y")

## Remove unnecessary columns in 2025 data

nm <- names(data_2025)

## helper to get a vector of column indices between two column names
col_range <- function(start_col, end_col, nm) {
  i1 <- match(start_col, nm)
  i2 <- match(end_col, nm)
  if (is.na(i1) || is.na(i2)) {
    stop(paste("Column name not found:", start_col, "or", end_col))
  }
  seq(min(i1, i2), max(i1, i2))
}

drop_idx <- c(
  col_range("2019Plot", "WaypointNumber", nm),
  col_range("SlopeGradient", "Photo", nm),
  col_range("SiteNotes", "CanopyCoverQ4", nm),
  col_range("MossCoverQ1", "VegNotes", nm)
)

data_2025 <- data_2025[, -unique(drop_idx), drop = FALSE]

## Adding up covers and moving them to a new column "Cover"

data_2025 <- data_2025 %>%
  mutate(CoverTree = AverageCoverA + AverageCoverB1 + AverageCoverB2)

data_2025 <- data_2025 %>%
  mutate(Cover = coalesce(AverageTotalB, AverageCoverC, CoverTree))

## Check for NA
which(is.na(data_2025$Cover)) 
sum(is.na(data_2025$Cover))

## Removing other cover columns

cols_to_drop <- c(
  "CoverAQ1","CoverAQ2","CoverAQ3","CoverAQ4","AverageCoverA",
  "CoverB1Q1","CoverB1Q2","CoverB1Q3","CoverB1Q4","AverageCoverB1",
  "CoverB2Q1","CoverB2Q2","CoverB2Q3","CoverB2Q4","AverageCoverB2",
  "TotalBQ1","TotalBQ2","TotalBQ3","TotalBQ4","AverageTotalB",
  "CoverCQ1","CoverCQ2","CoverCQ3","CoverCQ4","AverageCoverC",
  "CoverTree","AverageCanopyCover", "Latin_name")

data_2025_clean <- data_2025 %>%
  select(-any_of(cols_to_drop)) %>%
  select(-1)

#write.csv(data_2025_clean, file = "manning_park_data/data_2025_select.csv")

####

## Working to clean BEC data now

## Filter out "spp.", mosses, lichens, and grasses, removing unneeded columns

BEC_no_func <- BEC_data %>%
  filter(!grepl("SP\\.?$", Species)) %>%
  filter(!grepl("^(moss|lichen|grass)$", Functional_group)) %>%
  filter(PlotNumber %in% data_2025$PlotNumber) %>%
  # drop the first two unnamed columns by position
  select(-1,-2) %>%
  select(-any_of(
    c(
      "PlotRepresenting",
      "Location",
      "Ecosection",
      "LocationAccuracy",
      "Latin_name",
      "SlopeGradient",
      "Aspect",
      "MesoSlopePosition",
      "SurfaceShape",
      "SiteNotes",
      "SoilNotes",
      "VegSurveyor",
      "VegNotes",
      "SpeciesListComplete",
      "StrataCoverTree",
      "StrataCoverShrub",
      "StrataCoverMoss",
      "StrataCoverHerb",
      "Cover E",
      "Cover F")))

## Adding up covers and moving them to a new column "Cover"

BEC_clean <- BEC_no_func %>%
  mutate(Cover = rowSums(across(c(
    TotalA, TotalB, `Cover D`, `Cover C`)), na.rm = TRUE)) %>% 
  select(-any_of( #Removing unneeded columns after this
    c(
      "Cover A1",
      "Cover A2",
      "Cover A3",
      "Cover C",
      "Cover D",
      "TotalA",
      "TotalB",
      "Cover B1",
      "Cover B2")))

#write.csv(BEC_clean, file = "manning_park_data/BEC_select.csv")

####

## Make full species list WITHOUT filtering by occurrence minimum

final_species_list <- tibble(
  Species = sort(unique(c(
    as.character(BEC_clean$Species),
    as.character(data_2025_clean$Species)
  )))
) %>%
  filter(!is.na(Species), Species != "")

## Optional check: species missing from master species list
missing_from_master <- final_species_list %>%
  anti_join(species_names, by = c("Species" = "BC_species_code"))

missing_from_master

####

## Add year columns

data_2025_clean_all <- data_2025_clean %>%
  mutate(Year = year(as.Date(Date2025))) %>%
  select(-Date2025, -DateOriginal) %>%
  relocate(Year, .after = ProjectID)

BEC_clean_all <- BEC_clean %>%
  mutate(Year = year(as.Date(Date))) %>%
  select(-Date) %>%
  relocate(Year, .after = ProjectID)

####

## Make 2025 data long: every plot gets every species as 0/1

data_2025_obs_plot_species <- data_2025_clean_all %>%
  group_by(
    PlotNumber,
    ProjectID,
    Latitude,
    Longitude,
    Elevation,
    Year,
    Species,
    SiteSurveyor,
    RelocationConfidence
  ) %>%
  summarise(Cover = sum(Cover, na.rm = TRUE), .groups = "drop")

data_2025_plot_grid <- data_2025_clean_all %>%
  distinct(
    PlotNumber,
    ProjectID,
    Latitude,
    Longitude,
    Elevation,
    Year,
    SiteSurveyor,
    RelocationConfidence
  ) %>%
  crossing(final_species_list)

data_2025_long <- data_2025_plot_grid %>%
  left_join(
    data_2025_obs_plot_species,
    by = c(
      "PlotNumber",
      "ProjectID",
      "Latitude",
      "Longitude",
      "Elevation",
      "Year",
      "Species",
      "SiteSurveyor",
      "RelocationConfidence"
    )
  ) %>%
  mutate(
    present = if_else(!is.na(Cover) & Cover > 0, 1L, 0L),
    Cover = replace_na(Cover, 0)
  ) %>%
  left_join(species_names, by = c("Species" = "BC_species_code")) %>%
  select(-Latin_name)

####

#### Make BEC data long: every plot gets every species as 0/1

BEC_obs_plot_species <- BEC_clean_all %>%
  group_by(
    PlotNumber,
    ProjectID,
    Latitude,
    Longitude,
    Elevation,
    Year,
    Species,
    SiteSurveyor
  ) %>%
  summarise(Cover = sum(Cover, na.rm = TRUE), .groups = "drop")

BEC_plot_grid <- BEC_clean_all %>%
  distinct(
    PlotNumber,
    ProjectID,
    Latitude,
    Longitude,
    Elevation,
    Year,
    SiteSurveyor
  ) %>%
  crossing(final_species_list)

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
      "SiteSurveyor"
    )
  ) %>%
  mutate(
    present = if_else(!is.na(Cover) & Cover > 0, 1L, 0L),
    Cover = replace_na(Cover, 0)
  ) %>%
  left_join(species_names, by = c("Species" = "BC_species_code")) %>%
  select(-Latin_name)

####

## Converting to years since first survey

min_year <- min(
  c(BEC_data_long$Year, data_2025_long$Year),
  na.rm = TRUE)

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
           Latitude, Longitude, Elevation, present, Cover)
)

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

## Preparing the variables

model_data <- model_data %>%
  mutate(
    Species = factor(Species),
    ProjectID = factor(ProjectID),
    PlotNumber = factor(PlotNumber),
    SiteSurveyor = factor(SiteSurveyor),
    Elevation_sc = as.numeric(scale(Elevation)),
    RelocationConfidence = factor(RelocationConfidence),
    Time_sc = as.numeric(scale(YearsSinceStart))) %>%
  mutate(RelocationConfidence = 
           recode(RelocationConfidence,"High/ Med-High" = "High/Med-High"))

####

## filter for plots from 1976-2000 only
## If model_data_filt is already loaded:
plots_keep <- model_data_filt %>%
  distinct(PlotNumber)

model_data <- model_data %>%
  semi_join(plots_keep, by = "PlotNumber")

model_data <- model_data %>%
  mutate(Time = case_when(
    Year <= 2000 ~ "historical",
    Year == 2025 ~ "present"
  )) %>%
  mutate(Time = factor(Time, levels = c("historical", "present")))

## Save model_data

write.csv(model_data, file = "bayesian_analysis_scripts/model_data_all_species.csv")


####