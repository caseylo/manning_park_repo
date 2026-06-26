#### Analyzing relocation error

## Load libraries
library(tidyverse)
library(vegan)
library(lme4)
library(ggplot2)
library(ggbreak)
library(RColorBrewer)
library(dplyr)
library(pairwiseAdonis)
library(ggnewscale)
library(ggrepel)
library(emmeans)
library(ggeffects)
set.seed(123)

## Read data for PERMANOVA
data_long <- read_csv("data/processed/relocation_data_long.csv")

species_names <- read_csv("data/processed/master_species_list.csv")

## Create metadata frame
str(data_long)

data_long <- data_long  %>%
  select(-1)

## Structure long data to run model

## one record per site-species-treatment
data_2025_site_species <- data_long %>%
  distinct(Site, Treatment, Species, Functional_group)

## Species must occur in >1 site in BOTH 2019 and Error
species_keep_no_error <- data_long %>%
  filter(Treatment %in% c("2019", "2025")) %>%
  group_by(Species) %>%
  summarise(
    n_occ_2019 = sum(Treatment == "2019"),
    n_occ_2025 = sum(Treatment == "2025"),
    n_sites_total = n_distinct(Site),
    .groups = "drop"
  ) %>%
  filter(
    n_occ_2019 > 0,
    n_occ_2025 > 0,
    n_sites_total > 1
  ) %>%
  pull(Species)

species_keep_error <- data_long %>%
  filter(Treatment %in% c("2019", "Error")) %>%
  group_by(Species) %>%
  summarise(
    n_occ_2019 = sum(Treatment == "2019"),
    n_occ_Error = sum(Treatment == "Error"),
    n_sites_total = n_distinct(Site),
    .groups = "drop"
  ) %>%
  filter(
    n_occ_2019 > 0,
    n_occ_Error > 0,
    n_sites_total > 1
  ) %>%
  pull(Species)

## Use union of both species lists to build the full species x plot grid
species_all <- union(species_keep_no_error, species_keep_error)

species_df <- tibble(Species = as.character(species_all))

plot_grid <- data_long %>%
  distinct(
    PlotNumber, ProjectID, Site, Treatment, Latitude, Longitude,
    Elevation, Year, SiteSurveyor, RelocationConfidence
  ) %>%
  crossing(species_df)

length(unique(plot_grid$Site)) # Check to make sure this is working

data_long_long <- plot_grid %>%
  left_join(
    data_long,
    by = c(
      "PlotNumber",
      "ProjectID",
      "Treatment",
      "Site",
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
    present = if_else(is.na(Cover), 0L, 1L),
    Cover = replace_na(Cover, 0)
  ) %>%
  left_join(species_names, by = c("Species" = "BC_species_code")) %>%
  select(-Functional_group.x, -Latin_name) %>%
  rename(Functional_group = Functional_group.y)

min_year <- min(data_long_long$Year)

model_data <- data_long_long %>%
  mutate(
    YearsSinceStart = Year - min_year,
    Time = case_when(
      Year == 2019 ~ "2019",
      Year == 2025 ~ "2025"
    ),
    Time = factor(Time, levels = c("2019", "2025")))

## Prepare data for analysis
model_data <- model_data %>%
  mutate(
    Species = factor(Species),
    ProjectID = factor(ProjectID),
    PlotNumber = factor(PlotNumber),
    SiteSurveyor = factor(SiteSurveyor),
    Elevation_sc = as.numeric(scale(Elevation)),
    RelocationConfidence = factor(RelocationConfidence),
    Time_sc = as.numeric(scale(YearsSinceStart)),
    Lat_sc = as.numeric(scale(Latitude)),
    Lon_sc = as.numeric(scale(Longitude)),
    Time = factor(Time),
    Functional_group = factor(Functional_group)
  )

## data from 2019-2025
model_data_no_error <- model_data %>%
  filter(Treatment != "Error") %>%
  filter(Species %in% species_keep_no_error)

mod_no_error <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | Site)
  + Lat_sc,
  data = model_data_no_error,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod_no_error)

mod_no_error_fg <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | Site)
  + Lat_sc,
  data = model_data_no_error,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod_no_error_fg)

## data from 2019-Error
model_data_error <- model_data %>%
  filter(Treatment != "2025") %>%
  filter(Species %in% species_keep_error)

mod_error <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | Site)
  + Lat_sc,
  data = model_data_error,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod_error)

mod_error_fg <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | Site)
  + Lat_sc,
  data = model_data_error,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod_no_error)
summary(mod_no_error_fg)
summary(mod_error)
summary(mod_error_fg)

emtrends(
  mod_no_error_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

emmeans(mod_no_error_fg, pairwise ~ Time | Functional_group)

emtrends(
  mod_error_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

emmeans(mod_error_fg, pairwise ~ Time | Functional_group)

#### Visualizing elevation models

####

## Visualizing models: 2019-2025 correct
pred <- ggpredict(mod_no_error, terms = c("Elevation_sc [all]", "Time"))

## Adjust to show re-scaled elevation
mean_elev <- mean(model_data_no_error$Elevation)
sd_elev   <- sd(model_data_no_error$Elevation)
pred$Elevation_m <- pred$x * sd_elev + mean_elev

## Calculate optima and shift of predicted probability of presence
optima <- pred %>%
  group_by(group) %>%
  slice_max(predicted)

shift <- optima %>%
  select(group, Elevation_m) %>%
  pivot_wider(names_from = group, values_from = Elevation_m) %>%
  mutate(shift_m = `2025` - `2019`)

shift

## Pretty plot
no_error_plot <- ggplot(pred, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA) +
  scale_x_continuous(
    breaks = seq(
      floor(min(pred$Elevation_m) / 200) * 200,
      ceiling(max(pred$Elevation_m) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(pred$Elevation_m) / 100) * 100,
      ceiling(max(pred$Elevation_m) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(floor(min(pred$predicted) / 0.1) * 0.1, ceiling(max(pred$predicted) / 0.1 * 0.1), by = 0.1)) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence",
    colour = "Time",
    fill = "Time",
    title = "2019 - 2025 Correct Comparison") +
  geom_point(
    data = optima,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  theme_classic()

no_error_plot
#ggsave("outputs/figures/no_error_elevation_plot.png",plot = no_error_plot, width = 7, height = 4, dpi = 300)

## 2019-2025 correct with FG
pred_fg <- ggpredict(mod_no_error_fg, terms = c("Elevation_sc [all]", "Time", "Functional_group"))

## Adjust to show re-scaled elevation
mean_elev_fg <- mean(model_data_no_error$Elevation)
sd_elev_fg   <- sd(model_data_no_error$Elevation)
pred_fg$Elevation_m <- pred_fg$x * sd_elev_fg + mean_elev_fg

## Calculate optima of predicted probability of presence
optima_fg <- pred_fg %>%
  group_by(facet, group) %>%
  slice_max(predicted)

shift_fg <- optima_fg %>%
  select(facet, group, Elevation_m) %>%
  tidyr::pivot_wider(names_from = group, values_from = Elevation_m) %>%
  mutate(shift_m = `2025` - `2019`,
         label = paste0("+ ", round(shift_m, 1), " m"))

shift_fg

## Pretty plot
bestmod_fg_plot <- ggplot(pred_fg, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA) +
  scale_x_continuous(
    breaks = seq(
      floor(min(pred_fg$Elevation_m) / 200) * 200,
      ceiling(max(pred_fg$Elevation_m) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(pred_fg$Elevation_m) / 100) * 100,
      ceiling(max(pred_fg$Elevation_m) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(floor(min(pred_fg$predicted) / 0.1) * 0.1, ceiling(max(pred_fg$predicted) / 0.1 * 0.1), by = 0.1)) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence",
    colour = "Time",
    fill = "Time",
    title = "2019 - 2025 Correct Comparison") +
  geom_point(
    data = optima_fg,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  facet_wrap(~facet) +
  theme_classic()

bestmod_fg_plot
#ggsave("outputs/figures/no_error_elevation_fg_plot.png", plot = bestmod_fg_plot, width = 11, height = 4.5, dpi = 300)

####

## Visualizing models: 2019-2025 ERROR
pred_conf <- ggpredict(mod_error, terms = c("Elevation_sc [all]", "Time"))

## Adjust to show re-scaled elevation
mean_elev_conf <- mean(model_data_error$Elevation)
sd_elev_conf   <- sd(model_data_error$Elevation)
pred_conf$Elevation_m <- pred_conf$x * sd_elev_conf + mean_elev_conf

## Calculate optima of predicted probability of presence
optima_conf <- pred_conf %>%
  group_by(group) %>%
  slice_max(predicted)

## Pretty plot
conf_plot <- ggplot(pred_conf, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA) +
  scale_x_continuous(
    breaks = seq(
      floor(min(pred_conf$Elevation_m) / 200) * 200,
      ceiling(max(pred_conf$Elevation_m) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(pred_conf$Elevation_m) / 100) * 100,
      ceiling(max(pred_conf$Elevation_m) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(floor(min(pred_conf$predicted) / 0.1) * 0.1, ceiling(max(pred_conf$predicted) / 0.1 * 0.1), by = 0.1)) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence",
    colour = "Time",
    fill = "Time",
    title = "2019 - 2025 Error Comparison") +
  geom_point(
    data = optima_conf,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  theme_classic()

conf_plot
#ggsave("outputs/figures/error_elevation_plot.png",plot = conf_plot, width = 7, height = 4, dpi = 300)

####

## Confmod with FG
pred_conf_fg <- ggpredict(mod_error_fg, terms = c("Elevation_sc [all]", "Time", "Functional_group"))

## Adjust to show re-scaled elevation
mean_elev_conf_fg <- mean(model_data_error$Elevation)
sd_elev_conf_fg   <- sd(model_data_error$Elevation)
pred_conf_fg$Elevation_m <- pred_conf_fg$x * sd_elev_conf_fg + mean_elev_conf_fg

## Calculate optima of predicted probability of presence
optima_conf_fg <- pred_conf_fg %>%
  group_by(facet, group) %>%
  slice_max(predicted)

## Pretty plot
conf_fg_plot <- ggplot(pred_conf_fg, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA) +
  scale_x_continuous(
    breaks = seq(
      floor(min(pred_conf_fg$Elevation_m) / 200) * 200,
      ceiling(max(pred_conf_fg$Elevation_m) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(pred_conf_fg$Elevation_m) / 100) * 100,
      ceiling(max(pred_conf_fg$Elevation_m) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(floor(min(pred_conf_fg$predicted) / 0.1) * 0.1, ceiling(max(pred_conf_fg$predicted) / 0.1 * 0.1), by = 0.1)) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence",
    colour = "Time",
    fill = "Time",
    title = "2019 - 2025 Error Comparison") +
  geom_point(
    data = optima_conf_fg,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  facet_wrap(~facet) +
  theme_classic()

conf_fg_plot

#ggsave("outputs/figures/error_elevation_fg_plot.png", plot = conf_fg_plot, width = 11, height = 4.5, dpi = 300)

## Combining plots for easy visualization

library(patchwork)

combined_plot <- no_error_plot + conf_plot +
  plot_annotation(tag_levels = "A")

combined_plot

ggsave("outputs/figures/2019_combined_plot.png",
       plot = combined_plot,
       width = 11,
       height = 4.5,
       dpi = 300)

combined_plot2 <- bestmod_fg_plot / conf_fg_plot +
  plot_annotation(tag_levels = "A")

combined_plot2

ggsave("outputs/figures/2019_combined_plot_fg.png",
       plot = combined_plot2,
       width = 11,
       height = 8,
       dpi = 300)

## 

####

## re-scaling to compare full model with this model:
no_error_plot_scaled <- ggplot(pred, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA
  ) +
  scale_x_continuous(
    limits = c(680, 2100),
    breaks = seq(800, 1800, by = 200),
    minor_breaks = seq(700, 2100, by = 100),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  scale_y_continuous(
    limits = c(0.0, 0.8),
    breaks = seq(0.1, 0.8, by = 0.1)
  ) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence",
    colour = "Time",
    fill = "Time",
    title = "2019 - 2025 Correct Comparison"
  ) +
  geom_point(
    data = optima,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_classic()


no_error_plot_scaled

ggsave("outputs/figures/no_error_plot_scaled.png", plot = no_error_plot_scaled, width = 7, height = 4, dpi = 300)
