#### LMM Cover change

## Load libraries

library(lme4)
library(tidyverse)
library(dplyr)
library(lmerTest)
library(ggeffects)
library(ggplot2)

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
    Lon_sc = as.numeric(scale(Longitude)))

## Make sure Cover exists and is numeric
model_data <- model_data %>%
  filter(Cover > 0) %>%
  mutate(
    Cover = as.numeric(Cover),
    logCover = log1p(Cover))  # handles 0 values safely

####

## Model 1 (random intercept for Species)

cov_mod1_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (1 | Species),
  data = model_data,
  REML = FALSE
)

####

## Model 2 (Species-specific linear elevation response)

cov_mod2_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species),
  data = model_data,
  REML = FALSE
)

####

## Model 3 (add PlotNumber random intercept)

cov_mod3_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species) 
  + (1 | PlotNumber),
  data = model_data,
  REML = FALSE
)

####

## Adding ProjectID and PlotNumber

cov_mod3.1_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species) 
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  REML = FALSE
)

####

## Adding ProjectID only

cov_mod3.2_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species) 
  + (1 | ProjectID),
  data = model_data,
  REML = FALSE
)

####

## Adding quadratic elevation random effect

cov_mod4_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species),
  data = model_data,
  REML = FALSE
)

#### 

## Model 5 Adding (1 | PlotNumber)

cov_mod5_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber),
  data = model_data,
  REML = FALSE
)

####

## Adding (1 | ProjectID) and (1 | PlotNumber)

cov_mod5.1_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Adding (1 | ProjectID) only

cov_mod5.2_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | ProjectID),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#### 

## Replacing ProjectID with (1|SiteSurveyor) in the linear elevation model

cov_mod6_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + (1 | SiteSurveyor),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#### 

## Adding (1 | SiteSurveyor) only

cov_mod6.1_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | SiteSurveyor),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#### 

## Replacing ProjectID with (1|SiteSurveyor) in the quadratic elevation model

cov_mod7_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + (1 | SiteSurveyor),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#### 

## Adding (1 | ProjectID) only

cov_mod7.1_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | SiteSurveyor),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

#### 

## Replacing ProjectID with lat_sc lon_sc in linear elevation model

cov_mod8_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc,
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Adding ProjectID to lat_sc lon_sc

cov_mod8.1_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | ProjectID),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Adding PlotNumber to lat_sc lon_sc

cov_mod8.2_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | PlotNumber),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Adding PlotNumber and ProjectID to lat_sc lon_sc

cov_mod8.3_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Lat only

cov_mod8.4_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc,
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Long only

cov_mod8.5_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lon_sc,
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Adding lat long to poly elevation model

cov_mod9_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc,
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding ProjectID to lat_sc lon_sc

cov_mod9.1_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc
  + (1 | ProjectID),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Adding PlotNumber to lat_sc lon_sc

cov_mod9.2_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | PlotNumber),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Adding PlotNumber and ProjectID to lat_sc lon_sc

cov_mod9.3_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## lat only

cov_mod9.4_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc,
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Long only

cov_mod9.5_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lon_sc,
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Replacing ProjectID with SiteSurveyor

cov_mod10_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc
  + (1 | SiteSurveyor),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Lat only

cov_mod10.1_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc
  + (1 | SiteSurveyor),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####

## Long only

cov_mod10.2_full <- lmer(
  logCover ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lon_sc
  + (1 | SiteSurveyor),
  data = model_data,
  REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)

####
## AIC tables for comparison

AIC(cov_mod1_full, cov_mod2_full, cov_mod3_full,cov_mod3.1_full, cov_mod3.2_full,
    cov_mod4_full, cov_mod5_full, cov_mod5.1_full, cov_mod5.2_full, cov_mod6_full,
    cov_mod6.1_full, cov_mod7_full, cov_mod7.1_full, cov_mod8_full,
    cov_mod8.1_full, cov_mod8.2_full, cov_mod8.3_full, cov_mod8.4_full, cov_mod8.5_full,
    cov_mod9_full, cov_mod9.1_full, cov_mod9.2_full, cov_mod9.3_full,
    cov_mod9.4_full, cov_mod9.5_full, cov_mod10_full, cov_mod10.1_full, cov_mod10.2_full)

## Summaries

summary(cov_mod1_full)
summary(cov_mod2_full)
summary(cov_mod3_full)
summary(cov_mod3.1_full)
summary(cov_mod3.2_full)
summary(cov_mod4_full)
summary(cov_mod5_full)
summary(cov_mod5.1_full)
summary(cov_mod5.2_full)
summary(cov_mod6_full)
summary(cov_mod6.1_full)
summary(cov_mod7_full)
summary(cov_mod7.1_full)
summary(cov_mod8_full)
summary(cov_mod8.1_full)
summary(cov_mod8.2_full)
summary(cov_mod8.3_full)
summary(cov_mod8.4_full)
summary(cov_mod8.5_full)
summary(cov_mod9_full)
summary(cov_mod9.1_full)
summary(cov_mod9.2_full)
summary(cov_mod9.3_full)
summary(cov_mod9.4_full)
summary(cov_mod9.5_full)
summary(cov_mod10_full)
summary(cov_mod10.1_full)
summary(cov_mod10.2_full)

## Variance and Correlation

VarCorr(cov_mod1_full)
VarCorr(cov_mod2_full)
VarCorr(cov_mod3_full)
VarCorr(cov_mod3.1_full)
VarCorr(cov_mod3.2_full)
VarCorr(cov_mod4_full)
VarCorr(cov_mod5_full)
VarCorr(cov_mod5.1_full)
VarCorr(cov_mod5.2_full)
VarCorr(cov_mod6_full)
VarCorr(cov_mod6.1_full)
VarCorr(cov_mod7_full)
VarCorr(cov_mod7.1_full)
VarCorr(cov_mod8_full)
VarCorr(cov_mod8.1_full)
VarCorr(cov_mod8.2_full)
VarCorr(cov_mod8.3_full)
VarCorr(cov_mod9_full)


####

## Visualizing raw cover data

ggplot(model_data, aes(x = Elevation, y = Cover, colour = Year)) + 
  geom_point() + 
  facet_wrap(~Species, scales = "free") + 
  geom_smooth(method = "lm", formula = y~ poly(x,2, raw = TRUE), se = FALSE)
