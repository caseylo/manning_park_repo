library("tidyverse")

####Subsetting Hist data to include only sites surveyed in resruvey

#Subset resurvey to include only plots from <2019

resurvey_data_all <- read.csv("manning_park_data/manning_park_resurvey_data.csv")

#Object with pre 2019 plots only
resurvey_data_all %>%
  filter(X2019Plot =="N") %>%
  select(PlotNumber, Species) -> resurveyed_2025

#Object with 2019 plots only
resurvey_data_all %>%
  filter(X2019Plot =="Y") %>%
  select(PlotNumber, Species) -> resurvey_2019_only

####

####Read in historical Manning_allveg data

historical_data <- read.csv("manning_park_data/historical_data/Manning_plots_allveg.csv")

#Subset to include only historical plots resurveyed in 2025 (Not including 2019 plots)

historical_resurveyed <- historical_data %>%
  filter(PlotNumber %in% resurveyed_2025$PlotNumber)

#Sanity check that this worked:
table(unique(resurveyed_2025$PlotNumber))

length(unique(resurveyed_2025$PlotNumber))
#[1] 42

table(unique(historical_resurveyed$PlotNumber))
length(unique(historical_resurveyed$PlotNumber))
#[1] 42

####

####Exploring species and number of occurrences from 2025 resurvey data

resurvey_2025_species_counts <- as.data.frame(table(resurveyed_2025$Species))
colnames(resurvey_2025_species_counts) <- c("species_2025", "occurences_2025")

#Total species in 2025 resurvey
nrow(resurvey_2025_species_counts)
#[1] 227

sum(as.numeric(resurvey_2025_species_counts$occurences_2025) >= 10, na.rm = TRUE)
#[1] 38

#Create df with species 10 or more
resurvey_common_species <- subset(resurvey_2025_species_counts, occurences_2025 >= 10)

####

####Exploring species and number of occurrences from historical data

historical_species_counts <- as.data.frame(table(historical_resurveyed$Species))
colnames(historical_species_counts) <- c("species_hist", "occurences_hist")

#Total species in historical data
nrow(historical_species_counts)
#[1] 306

sum(as.numeric(historical_species_counts$occurences_hist) >= 10, na.rm = TRUE)
#[1] 29

#Create df with species 10 or more
hist_common_species <- subset(historical_species_counts, occurences_hist >= 10)

####

####Creating a list of common species between historical and resurvey data

species_both <- merge(
  hist_common_species,
  resurvey_common_species,
  by.x = "species_hist",
  by.y = "species_2025",
  all = FALSE
)

colnames(species_both)[colnames(species_both) == "species_hist"] <- "species"

species_both

#Make CSV of new shared data
write.csv(resurvey_common_species, file = "manning_park_data/resurvey_common_species_10plus.csv", row.names = FALSE)

