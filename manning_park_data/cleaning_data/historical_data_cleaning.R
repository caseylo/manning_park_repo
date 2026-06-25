####Cleaning original Manning Park data in R. Changes are also copied on a the "Tracking Data Tidying" Word Doc

library(tidyverse)
library(dplyr)
library(readr)

BEC_raw <- read_csv("manning_park_data/Manning_plots_allveg.csv")

####

####Removing unused columns from original BEC data

#Using any_of because I was running into errors using the select and -column function

BEC_select <- BEC_raw %>%
  select(
    -any_of(c(
      # Excel column A
      "X",
      
      # Site / admin metadata
      "FieldNumber", "FSRegionDistrict",
      "NtsMapSheet",
      
      # Coordinates / mapping
      "UTMZone", "UTMEasting", "UTMNorthing",
      "AirPhotoNum", "XCoord", "YCoord",
      "Zone", "SubZone",
      
      # Classification / site description
      "SiteSeries", "SiteModifier1", "SiteModifier2",
      "TransDistrib", "RealmClass", "MapUnit",
      "SnowCoverregime", "MoistureRegime", "NutrientRegime",
      "SuccessionalStatus", "StructuralStage", "StructuralStageMod",
      "StandAge",
      
      # Topography / hydrology
      "SurfaceTopographyType", "SurfaceTopographySize", "WaterSource",
      "Photo",
      
      # Disturbance / exposure
      "Exposure1", "Exposure2",
      "SiteDisturbance1", "SiteDisturbance2", "SiteDisturbance3",
      
      # Substrate
      "SubstrateDecWood",
      
      # Survey metadata
      "HydroGeoSystem", "HydroGeoSubSystem",
      "Temporary",
      "Layer",
      
      # Height columns
      "HeightA", "Height1", "Height2", "Height3", "Height4", "Height5",
      "Height5a", "Height5b", "Height5c", "HeightB", "Height6",
      
      # Cover extras without data
      "Cover5a", "Cover5b", "Cover5c", "Cover10"
    )),
    
    # Drop contiguous column blocks (only if present)
    -any_of(names(BEC_raw)[
      which(names(BEC_raw) == "SubstrateBedRock"):
        which(names(BEC_raw) == "SubstrateWater")
    ]),
    -any_of(names(BEC_raw)[
      which(names(BEC_raw) == "SoilSurveyor"):
        which(names(BEC_raw) == "pHMethodCodeOrganic")
    ]),
    -any_of(names(BEC_raw)[
      which(names(BEC_raw) == "Flag"):
        which(names(BEC_raw) == "DD_delayed")
    ]),
    -any_of(names(BEC_raw)[
      which(names(BEC_raw) == "Collected"):
        which(names(BEC_raw) == "Other2")
    ])
  ) %>% #Renaming some columns to make my life easier
  rename(
    `Cover A1` = Cover1,
    `Cover A2` = Cover2,
    `Cover A3` = Cover3,
    `Cover B1` = Cover4,
    `Cover B2` = Cover5,
    `Cover C`  = Cover6,
    `Cover D`  = Cover7,
    `Cover E`  = Cover8,
    `Cover F`  = Cover9
  )

##Checking all columns present are the wanted ones
names(BEC_select)

####

####Filtering out historical plots that were not resurveyed in 2025

#Read in 2025 resurvey data

data_2025 <- read_csv("manning_park_data/manning_park_resurvey_data.csv")

BEC_resurveyed <- BEC_select %>%
  filter(PlotNumber %in% data_2025$PlotNumber)

##Sanity check to make sure there are the same number of plots. There should be 10 more plots in the 2025 data because 2019 plots were resurveyed twice.

table(unique(data_2025$PlotNumber))
length(unique(data_2025$PlotNumber))
#[1] 62

table(unique(BEC_resurveyed$PlotNumber))
length(unique(BEC_resurveyed$PlotNumber))
#[1] 52


####

####Correcting species code typos, grouping subspecies, ect.

#Exploring species and number of occurrences from 2025 resurvey data

species_2025 <- as.data.frame(table(data_2025$Species))
colnames(species_2025) <- c("species_2025", "occurences_2025")

#Total species in 2025 resurvey
nrow(species_2025)
#[1] 226

write.csv(species_2025, file = "manning_park_data/cleaning_data/species_2025.csv")

#Making changes to the species codes

BEC_clean <- BEC_resurveyed %>%
  mutate(
    Species = recode(
      Species,
      "ABIELASD"   = "ABIELAS",
      "ACERGLA1"    = "ACERGLA",
      "ARNICA"      = "ARNICA SP.",
      "ACHIMIL2"    = "ACHIMIL",
      "ALNUALN2"    = "ALNUALN",
      "ANTEALPI"   = "ANTEALP",
      "ARENARIA"   = "ARENARIA SP.",
      "ANTENARIA"  = "ANTENARIA SP.",
      "ASTER"      = "ASTER SP.",
      "BRACHYT"    = "BRACHYTHECIUM SP.",
      "CALAMAGR"   = "CALAMAGROTIS SP.",
      "CAREX"      = "CAREX SP.",
      "CEPHALOZ"   = "CEPHALOZIA SP.",
      "CIRSIUM"    = "CIRSIUM SP.",
      "CLADONIA"   = "CLADONIA SP.",
      "CORNUNA"    = "CORNCAN",
      "DICRANUM"   = "DICRANUM SP.",
      "DRABA"      = "DRABA SP.",
      "EQUISETU"   = "EQUISETUM SP.",
      "ERIGERON"   = "ERIGERON SP.",
      "FESTUCA"    = "FESTUCA SP.",
      "HIERACIU"   = "HIERACIUM SP.",
      "LETHARIA"   = "LETHARIA SP.",
      "LICHEN"     = "LICHEN SP.",
      "LUPILEP2"   = "LUPILEP",
      "LOPHOZIA"   = "LOPHOZIA SP.",
      "MITELLA"    = "MITELLA SP.",
      "PICEENE"    = "PICEENG",
      "PICEENGD"   = "PICEENG",
      "PICEENED"   = "PICEENG",
      "PINUCON1"   = "PINUCON",
      "PINUCON2"   = "PINUCON",
      "PLAGIOMN"   = "PLAGIOMNIUM SP.",
      "POA_CUS2"   = "POA_CUS",
      "POACEAE_"    = "POACEAE SP.",
      "POLYGONU"   = "POLYGONUM SP.",
      "PSEUMEN1"   = "PSEUMEN",
      "PSEUMEN2"   = "PSEUMEN",
      "PSEUMEND"   = "PSEUMEN",
      "POLYTRIC"    = "POLYTRICHUM SP.",
      "POTENTIL"   = "POTENTILLA SP.",
      "POTEUNIF"   = "POTEUNI",
      "RACOMITR"   = "RACOMITRIUM SP.",
      "PETAFRIG"   = "PETAFRI",
      "PELTIGER"   = "PELTIGERA SP.",
      "ROSA"       = "ROSA SP.",
      "RHYTIDIO"   = "RHYTIDIOPSIS SP.",
      "SALIX"      = "SALIX SP.",
      "SENECIO"    = "SENECIO SP.",
      "LINOINV"    = "LONIINV",
      "SABURUB"    = "SAMBRAC",
      "POHLIA"     = "POHLIA SP.",
      "TIARTRI2"   = "TIARTRI",
      "VIOLA"      = "VIOLA SP."  
    )
  )


##Exploring species and number of occurrences from historical data

BEC_species <- as.data.frame(table(BEC_clean$Species))
colnames(BEC_species) <- c("species_BEC", "occurences_BEC")

#Total species in historical data
nrow(BEC_species)
#[1] 354

write.csv(BEC_species, file = "manning_park_data/cleaning_data/BEC_species.csv")

####

####Adding latin names and functional group next to the species column in each data frame (historical and 2025 data)

#For BEC data

BEC_species_code_and_latin_names <- read.csv("manning_park_data/cleaning_data/BEC_species_code_and_latin_names.csv")


BEC_clean <- BEC_clean %>%
  left_join(
    BEC_species_code_and_latin_names %>%
      select(BC_species_code, Latin_name, Functional_group),
    by = c("Species" = "BC_species_code")
  )

BEC_clean <- BEC_clean %>%
  relocate(Latin_name, Functional_group, .after = Species)

write.csv(BEC_clean, file = "manning_park_data/cleaning_data/BEC_clean.csv")

#For 2025 data

species_code_and_latin_names_2025 <- read.csv("manning_park_data/cleaning_data/species_code_and_latin_names_2025.csv")

data_2025 <- data_2025 %>%
  left_join(species_code_and_latin_names_2025, by = c("Species" = "BC_species_code")) %>%
  relocate(Latin_name, Functional_group, .after = Species)

data_2025 <- data_2025 %>%
  relocate(Latin_name, Functional_group, .after = Species)

####

####Sanity checks to make sure the functional groups and latin names transferred properly for 2025 data:

table(unique(data_2025$Functional_group))
#herb shrub shurb  tree 
#1     1     1     1 

#Need to mutate "shurb" to "shrub"

data_2025 <- data_2025 %>%
  mutate(Functional_group = recode(Functional_group, "shurb" = "shrub", ))

table(data_2025$Functional_group)
# herb shrub  tree 
#1063   531   204 

write.csv(data_2025, file = "manning_park_data/cleaning_data/data_2025_clean.csv")

####Sanity checks to make sure the functional groups and latin names transferred properly for 2025 data:

table(BEC_clean$Functional_group)
#grass   herb lichen   moss  shrub   tree 
#92    596     53    168    347    194 

table(data_2025$Latin_name)

####


