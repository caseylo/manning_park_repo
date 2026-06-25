#### Testing for differences in species richness between sites

## libraries
library(tidyverse)
library(dplyr)
library(vegan)
library(ggplot2)

## Read data
data_long <- read_csv("manning_park_data/relocation_data_long.csv")
comm_matrix <- read_csv("manning_park_data/relocation_data_wide.csv")

## Create metadata
metadata <- data_long %>%
  select(PlotNumber, Site, ProjectID, Treatment, Year,
         SiteSurveyor, RelocationConfidence,
         Latitude, Longitude, Elevation) %>%
  distinct() %>%
  mutate(
    Treatment = as.factor(Treatment),
    PlotNumber = as.factor(PlotNumber),
    Site = as.factor(Site))

## Species richness per plot
richness_df <- data_long %>%
  group_by(Site, Treatment) %>%
  summarise(richness = n_distinct(Species), .groups = "drop")

## Make wide. One row per plot and one column per year
richness_wide <- richness_df %>%
  pivot_wider(names_from = Treatment, values_from = richness)

## Pairwise test of differences in SR

wilcox.test(richness_wide$`2019`, richness_wide$`2025`, paired = TRUE)

wilcox.test(richness_wide$`2019`, richness_wide$Error, paired = TRUE)

wilcox.test(richness_wide$`2025`, richness_wide$Error, paired = TRUE)

####

## Testing Simpsons and Shannons diversity

#str(comm_matrix)

comm_matrix <- comm_matrix %>%
  as.data.frame() %>%
  select(-1) 

rownames(comm_matrix) <- comm_matrix$PlotNumber
comm_matrix$PlotNumber <- NULL

# Calculate diversity
div_df <- data.frame(
  PlotNumber = rownames(comm_matrix),
  Shannon = diversity(comm_matrix, index = "shannon"),
  Simpson = diversity(comm_matrix, index = "simpson"))

div_df <- div_df %>%
  left_join(metadata %>% 
            select(PlotNumber, Site, Treatment),
            by = "PlotNumber")

div_site <- div_df %>%
  group_by(Site, Treatment) %>%
  summarise(
    Shannon = mean(Shannon),
    Simpson = mean(Simpson),
    .groups = "drop"
  )

friedman.test(Shannon ~ Treatment | Site, data = div_site)
friedman.test(Simpson ~ Treatment | Site, data = div_site)

pairwise.wilcox.test(div_site$Shannon, div_site$Treatment,
                     paired = TRUE, p.adjust.method = "BH")

pairwise.wilcox.test(div_site$Simpson, div_site$Treatment,
                     paired = TRUE, p.adjust.method = "BH")


ggplot(div_site, aes(x = Treatment, y = Shannon, group = Site, colour = Site)) +
  geom_line(alpha = 0.9) +
  geom_point(size = 2) +
  theme_classic()

ggplot(div_site, aes(x = Treatment, y = Simpson, group = Site, colour = Site)) +
  geom_line(alpha = 0.9) +
  geom_point(size = 2) +
  theme_classic()
