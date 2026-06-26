#### Calculating the difference between 2019 and error plots and 2019 and 2025 plots

## Load libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(geosphere)

## Read data
data_long <- read_csv("data/processed/relocation_data_long.csv")

## Keep one set of coordinates per Site x Treatment
coords <- data_long %>%
  select(Site, Treatment, Latitude, Longitude) %>%
  distinct()

# Make treatments into separate columns
coords_wide <- coords %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(Latitude, Longitude))

# Calculate distances (in meters) at each site
distances <- coords_wide %>%
  mutate(
    dist_2019_Error = distHaversine(
      cbind(Longitude_2019, Latitude_2019),
      cbind(Longitude_Error, Latitude_Error)),
    dist_2019_2025 = distHaversine(
      cbind(Longitude_2019, Latitude_2019),
      cbind(Longitude_2025, Latitude_2025)))

distances

## Summarize

distances %>%
  summarise(
    mean_2019_Error = mean(dist_2019_Error, na.rm = TRUE),
    mean_2019_2025 = mean(dist_2019_2025, na.rm = TRUE))

## Make a table

distance_table <- coords_wide %>%
  mutate(
    `2019_vs_Error_m` = distHaversine(
      cbind(Longitude_2019, Latitude_2019),
      cbind(Longitude_Error, Latitude_Error)),
    `2019_vs_2025_m` = distHaversine(
      cbind(Longitude_2019, Latitude_2019),
      cbind(Longitude_2025, Latitude_2025))) %>%
  select(
    Site,
    `2019_vs_Error_m`,
    `2019_vs_2025_m`)

distance_table

## Add an average at the end

distance_table <- distance_table %>%
  summarise(
    Site = "MEAN",
    `2019_vs_Error_m` = mean(`2019_vs_Error_m`, na.rm = TRUE),
    `2019_vs_2025_m` = mean(`2019_vs_2025_m`, na.rm = TRUE)
  ) %>%
  bind_rows(distance_table, .)

distance_table
