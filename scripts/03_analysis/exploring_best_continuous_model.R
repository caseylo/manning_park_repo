#### Analyzing model 8.4 from model selection

## Load libraries

library(tidyverse)
library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)
library(ggeffects)
library(spdep)
library(sf)
library(dplyr)
library(ggh4x)

## Read in data

model_data <- read_csv("data/processed/model_data.csv")

## Prepare data

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
    Functional_group = factor(Functional_group))

####

## Exploring best model (model 8)

bestmod <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(bestmod)

####

table(model_data$Functional_group)

## Functional group addition

fgmod <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(fgmod)

####

## Quantifying Upslope Range Shift using ggpredict

## mod8 linear elevation

predbest <- ggpredict(
  bestmod,
  terms = c("Elevation_sc [all]", "Time_sc [-2,0,2]"),
  bias_correction = F, 
)

## OR adjust to show re-scaled elevation

mean_elev <- mean(model_data$Elevation)
sd_elev   <- sd(model_data$Elevation)

## mod8

predbest$Elevation_m <- predbest$x * sd_elev + mean_elev

optima8 <- predbest %>%
  group_by(group) %>%
  slice_max(predicted)

pred_plot <- ggplot(predbest, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA
  ) +
  scale_x_continuous(
    breaks = seq(
      floor(min(predbest$Elevation_m) / 200) * 200,
      ceiling(max(predbest$Elevation_m) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(predbest$Elevation_m) / 100) * 100,
      ceiling(max(predbest$Elevation_m) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(floor(min(predbest$predicted) / 0.1) * 0.1, ceiling(max(predbest$predicted) / 0.1 * 0.1), by = 0.1)) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence (mod8)",
    colour = "Time",
    fill = "Time"
  ) +
  geom_vline(
    data = optima8,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    data = optima8,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_classic()

pred_plot

ggsave("outputs/figures/continuous_mod_elevation_curve.png",plot = pred_plot, width = 7, height = 4, dpi = 300)

####

## Visualizing functional group

## mod8 WITH functional groups

pred_fg <- ggpredict(
  fgmod,
  terms = c("Elevation_sc [all]", "Time_sc [-2,0,2]","Functional_group"),
  bias_correction = F, 
)

## Better visual

mean_elev.fg <- mean(model_data$Elevation)
sd_elev.fg   <- sd(model_data$Elevation)

##

pred_fg$Elevation_m <- pred_fg$x * sd_elev.fg + mean_elev.fg

optima.fg <- pred_fg %>%
  group_by(facet, group) %>%
  slice_max(predicted)

pred_fg_plot <- ggplot(pred_fg, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA
  ) +
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
    y = "Predicted probability of occurrence (conf_fg)",
    colour = "Time",
    fill = "Time"
  ) +
  geom_vline(
    data = optima.fg,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    data = optima.fg,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE
  ) +
  facet_wrap(~facet) +
  theme_classic()

pred_fg_plot

ggsave("outputs/figures/continuous_fg_elevation_curve.png", plot = pred_fg_plot, width = 11, height = 4.5, dpi = 300)

####