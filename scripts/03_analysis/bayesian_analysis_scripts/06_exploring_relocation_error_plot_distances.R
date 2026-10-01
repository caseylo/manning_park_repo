## Exploring distances and relocation confidence

## Relocation confidence

# Load libraries
library(tidyverse)
library(geosphere)

# Load data
all_data <- read_csv("data/processed/model_data_filt.csv")
error_data <- read_csv("data/processed/relocation_data_long.csv")

# Count plots with medium-high or high relocation confidence


unique(all_data$RelocationConfidence)
unique(error_data$RelocationConfidence)

all_data %>%
  filter(RelocationConfidence %in% c("Med-High", "High")) %>%
  summarise(
    n_plots = n_distinct(PlotNumber)
  )
#24/34

error_data %>%
  filter(Treatment == "2025") %>%
  filter(RelocationConfidence %in% c("Med-High", "High")) %>%
  summarise(
    n_plots = n_distinct(PlotNumber)
  )
#9/9

#Total 33/43 plots were confidently relocated

# Count total unique plots in each dataset

n_distinct(all_data$PlotNumber)
#[1] 34
n_distinct(error_data$PlotNumber)
#[1] 18

####

## Calculate distance between 2025correct and 2025 error plots on average

## 2025 vs Error

# Extract coordinates for each treatment
coords_2025 <- error_data %>%
  filter(Treatment == "2025") %>%
  distinct(Site, Latitude, Longitude) %>%
  select(Site, Latitude, Longitude) %>%
  rename(
    Latitude_2025 = Latitude,
    Longitude_2025 = Longitude
  )

coords_error <- error_data %>%
  filter(Treatment == "Error") %>%
  distinct(Site, Latitude, Longitude) %>%
  select(Site, Latitude, Longitude) %>%
  rename(
    Latitude_error = Latitude,
    Longitude_error = Longitude
  )

# Calculate distances in metres
distances <- inner_join(coords_2025, coords_error, by = "Site") %>%
  mutate(
    Distance_m = geosphere::distHaversine(
      p1 = as.matrix(select(., Longitude_2025, Latitude_2025)),
      p2 = as.matrix(select(., Longitude_error, Latitude_error))
    )
  )

# Summarise distances
distances_2025_error_sum <- distances %>%
  summarise(
    n_sites = n(),
    mean_distance_m = mean(Distance_m),
    median_distance_m = median(Distance_m),
    min_distance_m = min(Distance_m),
    max_distance_m = max(Distance_m),
    sd_distance_m = sd(Distance_m)
  )

print(distances_2025_error_sum)

## 2019 vs 2025

# Extract 2019 coordinates
coords_2019 <- error_data %>%
  filter(Treatment == "2019") %>%
  distinct(Site, Latitude, Longitude) %>%
  select(Site, Latitude, Longitude) %>%
  rename(
    Latitude_2019 = Latitude,
    Longitude_2019 = Longitude
  )

# Calculate distances in metres
distances_2019_2025 <- inner_join(
  coords_2019,
  coords_2025,
  by = "Site"
) %>%
  mutate(
    Distance_m = geosphere::distHaversine(
      p1 = as.matrix(select(., Longitude_2019, Latitude_2019)),
      p2 = as.matrix(select(., Longitude_2025, Latitude_2025))
    )
  )

# Summarise distances
distances_2019_2025_sum <- distances_2019_2025 %>%
  summarise(
    n_sites = n(),
    mean_distance_m = mean(Distance_m),
    median_distance_m = median(Distance_m),
    min_distance_m = min(Distance_m),
    max_distance_m = max(Distance_m),
    sd_distance_m = sd(Distance_m)
  )

print(distances_2019_2025_sum)

## 2019 vs Error

# Match sites and calculate distances in metres
distances_2019_error <- inner_join(
  coords_2019,
  coords_error,
  by = "Site"
) %>%
  mutate(
    Distance_m = geosphere::distHaversine(
      p1 = as.matrix(select(., Longitude_2019, Latitude_2019)),
      p2 = as.matrix(select(., Longitude_error, Latitude_error))
    )
  )


# Summarise distances
distances_2019_error_sum <- distances_2019_error %>%
  summarise(
    n_sites = n(),
    mean_distance_m = mean(Distance_m),
    median_distance_m = median(Distance_m),
    min_distance_m = min(Distance_m),
    max_distance_m = max(Distance_m),
    sd_distance_m = sd(Distance_m)
  )

print(distances_2019_error_sum)


## All comparsions
print(distances_2025_error_sum)

print(distances_2019_2025_sum)

print(distances_2019_error_sum)
