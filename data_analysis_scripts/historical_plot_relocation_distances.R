#### Calculatin the differnce between 2019 and error plots and 2019 and 2025 plots

## Load libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(geosphere)

## Read data
data_long <- read_csv("manning_park_data/model_data.csv")
#data_long <- read_csv("manning_park_data/model_data_filt.csv")

###

data_long <- data_long %>%
  mutate(Time = case_when(
    Year <= 2019 ~ "historical",
    Year == 2025 ~ "present"
  )) %>%
  mutate(Time = factor(Time, levels = c("historical", "present")))

#data_long <- data_long %>%
#  filter(Year != 2019)

data_long <- data_long %>%
  filter(!PlotNumber %in% c(
    "35079", "35379", "00NCe89", "MP 288", "SK93722",
    "SK93735", "SK93721", "SK93734", "00NCe84", "MP 287",
    "SK93736", "SK93711", "SK93710", "00NCe92"
  ))

## Keep one set of coordinates per Site x Treatment
coords <- data_long %>%
  select(PlotNumber, Time, Latitude, Longitude) %>%
  distinct()

# Make treatments into separate columns
coords_wide <- coords %>%
  pivot_wider(
    names_from = Time,
    values_from = c(Latitude, Longitude))

# Calculate distances (in meters) at each site
distances <- coords_wide %>%
  mutate(
    distance = distHaversine(
      cbind(Longitude_historical, Latitude_historical),
      cbind(Longitude_present, Latitude_present)))

distances

## Summarize

distances %>%
  summarise(
    mean_distance = mean(distance, na.rm = TRUE))

## Make a table

