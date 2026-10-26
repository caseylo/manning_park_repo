## Simple regression to assess warming trend at sites in Manning Park

library(tidyverse)
library(ggplot2)


# Read ClimateNA data
climate <- read_csv(
  "data/cleaned/ClimateNA_locations_manning_filtered_MAT_1976-2024Y.csv"
)

# Calculate mean annual temperature across sites
annual_mat <- climate %>%
  group_by(Year) %>%
  summarise(
    mean_MAT = mean(MAT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Year)

# Fit linear warming trend
model <- lm(mean_MAT ~ Year, data = annual_mat)

# View model results, including p-value and confidence interval
summary(model)
confint(model)

# Calculate warming rate in °C per decade
warming_per_decade <- coef(model)["Year"] * 10
warming_per_decade

## Simple MAT vs Year graph

## Simple MAT vs Year graph with mean line

siteMAT <- ggplot(climate, aes(
  x = Year,
  y = MAT,
  group = ID1,
  colour = factor(ID1)
)) +
  geom_line(linewidth = 0.3, alpha = 0.8, show.legend = FALSE) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "line",
    colour = "red",
    linewidth = 1
  ) +
  stat_summary(
    aes(group = 1),
    fun = mean,
    geom = "point",
    colour = "red",
    size = 2
  ) +
  labs(
    x = "Year",
    y = "Mean annual temperature (°C)"
  ) +
  theme_classic()

siteMAT

## Save

ggsave(
  "outputs/figures/bayesian_figures/MAT_vs_Year.png",
  width = 8,
  height = 5,
  dpi = 300
)
