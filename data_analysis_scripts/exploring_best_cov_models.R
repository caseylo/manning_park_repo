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

model_data <- read_csv("manning_park_data/model_data.csv")

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

model_data <- model_data %>%
  filter(Cover > 0) %>%
  mutate(
    Cover = as.numeric(Cover),
    logCover = log1p(Cover))

####

covmod_c <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE)
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1|SiteSurveyor)
  + Lat_sc,
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa"))

summary(covmod_c)

## FG model

covmod_cfg <- lmer(
  logCover ~ Time_sc * Functional_group
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc
  + (1 | SiteSurveyor),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

summary(covmod_cfg)

#### IS SINGULAR

####

## Comparing filtered model

model_data2 <- read_csv("manning_park_data/model_data_filt.csv")

## Prepare data

model_data2 <- model_data2 %>%
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
    Functional_group = factor(Functional_group))

## Make sure Cover exists and is numeric
model_data2 <- model_data2 %>%
  filter(Cover > 0) %>%
  mutate(
    Cover = as.numeric(Cover),
    logCover = log1p(Cover))  # handles 0 values safely

####

covmod_b_fg <- lmer(
  logCover ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + Time * Functional_group
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc,
  data = model_data2,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(covmod_b_fg)

predcov2 <- ggpredict(cov_mod9.4, terms = c("Time"))

plot(predcov2) + theme_classic()
