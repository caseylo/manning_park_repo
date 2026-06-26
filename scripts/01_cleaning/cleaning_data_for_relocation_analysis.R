#### Relocation error estimates

## Cleaning data for relocation error analysis

## Load libraries

library(tidyverse)
library(dplyr)
library(lubridate)
library(tidyr)

## Read in data

BEC_data <- read_csv("data/cleaned/BEC_clean.csv")
data_2025 <- read_csv("data/cleaned/data_2025_clean.csv")
species_names <- read_csv("data/cleaned/master_species_list.csv")

####

## Working with 2025 data first

## Filter out with "spp."

data_2025 <- data_2025 %>%
  filter(!grepl("SP\\.?$", Species),
         `2019Plot` == "Y") %>%
  ## Adding up covers and moving them to a new column "Cover"
  mutate(CoverTree = AverageCoverA + AverageCoverB1 + AverageCoverB2) %>%
  mutate(Cover = coalesce(AverageTotalB, AverageCoverC, CoverTree))

unique(data_2025$PlotNumber)

## Check for NA
which(is.na(data_2025$Cover)) 
sum(is.na(data_2025$Cover))

## Removing other cover columns
data_2025_filt <- data_2025 %>%
select(-any_of(
  c( "SlopeGradient", "Aspect", "MesoSlopePosition", "SurfaceShape", "Photo",
     "SiteNotes", "CanopyCoverQ1", "CanopyCoverQ2", "CanopyCoverQ3", "CanopyCoverQ4",
     "MossCoverQ1", "MossCoverQ2", "MossCoverQ3", "MossCoverQ4", "AverageMossCover", 
     "GrassCoverQ1", "GrassCoverQ2", "GrassCoverQ3", "GrassCoverQ4", "AverageGrassCover", 
     "VegNotes","WaypointNumber","CoverAQ1","CoverAQ2","CoverAQ3","CoverAQ4","AverageCoverA",
     "CoverB1Q1","CoverB1Q2","CoverB1Q3","CoverB1Q4","AverageCoverB1",
     "CoverB2Q1","CoverB2Q2","CoverB2Q3","CoverB2Q4","AverageCoverB2",
     "TotalBQ1","TotalBQ2","TotalBQ3","TotalBQ4","AverageTotalB",
     "CoverCQ1","CoverCQ2","CoverCQ3","CoverCQ4","AverageCoverC",
     "CoverTree","AverageCanopyCover", "Latin_name"))) %>%
  mutate(Year = year(as.Date(Date2025))) %>%
  select(-Date2025, -DateOriginal) %>%
  relocate(Year, .after = ProjectID) %>%
  select(-1)

data_2025_clean <- data_2025_filt %>%
  mutate(
    Site = if_else(
      ErrorPlot == "N",
      as.character(PlotNumber),
      sub("E$", "", as.character(PlotNumber))),
    PlotNumber = if_else(
      ErrorPlot == "N",
      paste0(as.character(PlotNumber), "_2025"),
      paste0(sub("E$", "", as.character(PlotNumber)), "_Error"))) %>%
  relocate(Site, .after = PlotNumber) %>%
  mutate(SiteSurveyor = "CML RNW") %>%
  select(-`2019Plot`,-ErrorPlot) %>%
  mutate(ProjectID = "Manning Park Project 2025")

#write.csv(data_2025_clean, file = "data/cleaned/data_2025_relocation.csv")

####

## Working to clean BEC data now

## Filter out "spp.", mosses, lichens, and grasses, removing unneeded columns

BEC_clean<- BEC_data %>%
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
      "Cover F"))) %>%
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
      "Cover B2"))) %>%
  mutate(Year = year(as.Date(Date))) %>%
  select(-Date) %>%
  relocate(Year, .after = ProjectID) %>%
  mutate(
    Site = PlotNumber,
    PlotNumber = paste0(PlotNumber, "_2019")) %>%
  relocate(Site, .after = PlotNumber)

#write.csv(BEC_clean, file = "data/cleaned/BEC_relocation.csv")

####

## BEC: one record per plot-species (Some species have been duplicated in plots)
## Making sure there are no plots causing inflated occurences:

BEC_clean %>%
  count(PlotNumber, Species) %>%
  filter(n > 1) %>%
  arrange(desc(n))

####

## Double-check everything is in order

unique(BEC_clean$PlotNumber)
unique(BEC_clean$Site)
unique(data_2025_clean$PlotNumber)
unique(data_2025_clean$Site)

setdiff(names(BEC_clean), names(data_2025_clean))
setdiff(names(data_2025_clean), names(BEC_clean))

####

## Join together data set

data_all <- bind_rows(BEC_clean, data_2025_clean)

## Need to mutate cover so that andything =< 1 is changed to = 1
data_all %>%
  filter(Cover <= 1) %>%
  nrow()
#[1] 429

data_all <- data_all %>%
  mutate(Treatment = case_when(
    grepl("2019", PlotNumber) ~ "2019",
    grepl("2025", PlotNumber) ~ "2025",
    grepl("Error", PlotNumber) ~ "Error")) %>%
  relocate(RelocationConfidence, .after = Longitude) %>%
  relocate(Treatment, .after = ProjectID) %>%
  mutate(Cover = if_else(Cover == 0, 0.1, Cover),
         Cover = if_else(Cover < 1, 1.0, Cover)) %>% ## changing any cover values from 2025 <1 to =1, to match the 2019 data
  group_by(PlotNumber, Species) %>%
  mutate(Cover = sum(Cover)) %>%
  distinct(PlotNumber, Species, .keep_all = TRUE) %>%
  ungroup() %>%
  filter(Site != "14-4311") ## Filter out the plot that experienced a post 2019 fire

data_all %>%
  filter(Cover <= 1) %>%
  nrow()
#[1] 360

## Check for no duplicated Species

data_all %>%
  count(PlotNumber, Species) %>%
  filter(n > 1)

write.csv(data_all, file = "data/processed/relocation_data_long.csv")

## Turn species into columns and plots into rows

comm_matrix <- data_all %>%
  select(PlotNumber, Species, Cover) %>%
  pivot_wider(
    names_from = Species,
    values_from = Cover,
    values_fill = 0)

write.csv(comm_matrix, file = "data/processed/relocation_data_wide.csv")

