#### Model Selection of categorical data. Comparing binned data

#### Logistic mixed model to evaluate change in probability or occurrence over elevation

## Load libraries

library(tidyverse)
library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)
library(ggeffects)

## Read in data

model_data <- read_csv("data/processed/model_data_filt.csv")

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
    Time = factor(Time))

## Check structure

str(model_data)

####

## Model Selection

## Fitting the Simplest Model

mod1_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) + (1 | Species),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod1_full)

####

## Adding linear elevation change as a random effect. Allows each species to have it's own elevation response

mod2_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species), # Do species differ in how much they shift
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod2_full)

####

## Adding (1 | PlotNumber), not nested
mod3_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species) # Do species differ in how much they shift
  + (1 | PlotNumber),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod3_full)

####

## Adding (1 | PlotNumber) and (1 | ProjectID)
mod3.1_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species) # Do species differ in how much they shift
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod3.1_full)

####

## Adding (1 | ProjectID) only
mod3.2_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = T) 
  + (Elevation_sc | Species) # Do species differ in how much they shift
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod3.2_full)

####

## Adding (poly(Elevation_sc,2) | Species) and taking out simpler elevation term
mod4_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species), #Do species differ in linear elevation and the shape of the curve
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod4_full)

#### 

## Adding (poly(Elevation_sc,2) | Species) and (1 | PlotNumber)
mod5_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species) #Do species differ in linear elevation and the shape of the curve
  + (1 | PlotNumber),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod5_full)

#### 

## Testing ProjectID and PlotNumber as random effects
mod5.1_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod5.1_full)

####

## Testing (1 | ProjectI) as a random effect only
mod5.2_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod5.2_full)

####

## Testing (1 | SiteSurveyor) and (1 | PlotNumber)
mod6_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + (1 | SiteSurveyor),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod6_full)

####

## Testing (1 | SiteSurveyor) only
mod6.1_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | SiteSurveyor),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod6.1_full)

####

## Testing Linear elevatoin and (1 | SiteSurveyor) and (1 | PlotNumber)
mod7_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + (1 | SiteSurveyor),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod7_full)

####

## Testing (1 | SiteSurveyor) only
mod7.1_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | SiteSurveyor),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa"))

summary(mod7.1_full)

####

## Adding  Lat lon as fixed effects

mod8_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod8_full)
ranef(mod8_full)

####

## Adding ProjectID and PlotNumber

mod8.1_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod8.1_full)

####

## Adding ProjectID only (no PlotNumber)

mod8.2_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod8.2_full)

####

## No PlotNumber no ProjectID

mod8.3_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod8.3_full)

####

## Adding  Lat only

mod8.4_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod8.4_full)

####

## Adding long only

mod8.5_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod8.5_full)

####

## Adding  Lat lon as fixed effects

mod9_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod9_full)

####

## Adding ProjectID and PlotNumber

mod9.1_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod9.1_full)

####

## Adding ProjectID only (no PlotNumber)

mod9.2_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod9.2_full)

####

## No PlotNumber no ProjectID

mod9.3_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod9.3_full)

####

## Adding  Lat only

mod9.4_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod9.4_full)

####

## Adding long only

mod9.5_full <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod9.5_full)

####


## Comparing models

AIC(mod1_full, mod2_full, mod3_full, mod3.1_full, mod3.2_full,
    mod4_full, mod5_full, mod5.1_full, mod5.2_full, mod6_full, 
    mod6.1_full, mod7_full, mod7.1_full, mod8_full, mod8.1_full,
    mod8.2_full, mod8.3_full, mod8.4_full, mod8.5_full,
    mod9_full, mod9.1_full, mod9.2_full, mod9.3_full, mod9.4_full, mod9.5_full)

## Summaries

summary(mod1_full)
summary(mod2_full)
summary(mod3_full)
summary(mod3.1_full)
summary(mod3.2_full)
summary(mod4_full)
summary(mod5_full)
summary(mod5.1_full)
summary(mod5.2_full)
summary(mod6_full)
summary(mod6.1_full)
summary(mod7_full)
summary(mod7.1_full)
summary(mod8_full)
summary(mod8.1_full)
summary(mod8.2_full)
summary(mod8.3_full)
summary(mod8.4_full)
summary(mod8.5_full)
summary(mod9_full)
summary(mod9.1_full)
summary(mod9.2_full)
summary(mod9.3_full)
summary(mod9.4_full)
summary(mod9.5_full)

####

## Visualizing raw data

ggplot(model_data, aes(x = Elevation, y = present)) + 
  geom_point() + 
  facet_wrap(~Species) + 
  geom_smooth(method = "lm", formula = y~ poly(x,2, raw = TRUE), se = FALSE)
