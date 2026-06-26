#### RUNNING THE FIRST ANALYSIS!!!!!

#### Logistic mixed model to evaluate change in probability or occurrence over elevation

## Load libraries

library(tidyverse)
library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)
library(ggeffects)

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

## Check structure

str(model_data)

####

#### Model Selection

## Fitting the Simplest Model
mod1_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) + (1 | Species),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Adding linear elevation change as a random effect. Allows each species to have it's own elevation response
mod2_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species), # Do species differ in how much they shift
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Adding (1 | PlotNumber), not nested
mod3_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species) # Do species differ in how much they shift
  + (1 | PlotNumber),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Adding (1 | PlotNumber) and (1 | ProjectID)
mod3.1_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species) # Do species differ in how much they shift
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Adding (1 | ProjectID) only
mod3.2_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = T) 
  + (Elevation_sc | Species) # Do species differ in how much they shift
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Adding (poly(Elevation_sc,2) | Species) and taking out simpler elevation term
mod4_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species), #Do species differ in linear elevation and the shape of the curve
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

#### 

## Adding (poly(Elevation_sc,2) | Species) and (1 | PlotNumber)
mod5_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species) #Do species differ in linear elevation and the shape of the curve
  + (1 | PlotNumber),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

#### 

## Testing ProjectID and PlotNumber as random effects
mod5.1_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Testing (1 | ProjectI) as a random effect only
mod5.2_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Testing (1 | SiteSurveyor) and (1 | PlotNumber)
mod6_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + (1 | SiteSurveyor),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Testing (1 | SiteSurveyor) only
mod6.1_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | SiteSurveyor),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Testing Linear elevatoin and (1 | SiteSurveyor) and (1 | PlotNumber)
mod7_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + (1 | SiteSurveyor),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Testing (1 | SiteSurveyor) only
mod7.1_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | SiteSurveyor),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)

####

## Adding  Lat lon as fixed effects

mod8_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding ProjectID and PlotNumber

mod8.1_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

#### DID NOT CONVERGE

####

## Adding ProjectID only (no PlotNumber)

mod8.2_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## No PlotNumber no ProjectID

mod8.3_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding  Lat only

mod8.4_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding long only

mod8.5_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding  Lat lon as fixed effects to linear model

mod9_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding ProjectID and PlotNumber

mod9.1_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | PlotNumber)
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding ProjectID only (no PlotNumber)

mod9.2_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc
  + (1 | ProjectID),
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## No PlotNumber no ProjectID

mod9.3_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + Lat_sc + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding  Lat only

mod9.4_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Adding long only

mod9.5_full <- glmer(
  present ~ Time_sc * poly(Elevation_sc, 2, raw = TRUE) 
  + (Elevation_sc | Species)
  + (1 | PlotNumber)
  + Lon_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

####

## Comparing model AIC values

AIC(mod1_full, mod2_full, mod3_full, mod3.1_full, mod3.2_full,
    mod4_full, mod5_full, mod5.1_full, mod5.2_full, mod6_full, 
    mod6.1_full, mod7_full, mod7.1_full, mod8_full, mod8.1_full,
    mod8.2_full, mod8.3_full, mod8.4_full, mod8.5_full,
    mod9_full, mod9.1_full, mod9.2_full, mod9.3_full, mod9.4_full, mod9.5_full)

## Model summaries

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

## Moran's I test for spatial autocorrelation in GLMM residuals

#Extract Pearson residuals
model_data$resid_pearson <- residuals(mod8_full, type = "pearson")

## Collapse to ONE row per plot (because plots are the spatial units)
## This avoids identical points (two timepoints per plot at same coords)

plot_resids <- model_data %>%
  group_by(PlotNumber) %>%
  summarise(
    Lat_sc = first(Lat_sc),
    Lon_sc = first(Lon_sc),
    resid_pearson = mean(resid_pearson, na.rm = TRUE),
    .groups = "drop"
  )

coords <- cbind(plot_resids$Lon_sc, plot_resids$Lat_sc)

# Start with k=4; if the graph is disconnected, increase k
k <- 10
nb <- knn2nb(knearneigh(coords, k = k))

## Spatial weights
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

## Moran's I test
moran_res <- moran.test(plot_resids$resid_pearson, listw, zero.policy = TRUE)
moran_res

####

## Visualizing raw data

model_data %>%
  group_by(ProjectID) %>%
  summarise(mean_elev = mean(Elevation),
            sd_elev = sd(Elevation),
            n = n())
ggplot(model_data, aes(x = Elevation, y = present)) + 
  geom_point() + 
  facet_wrap(~Species) + 
  geom_smooth(method = "lm", formula = y~ poly(x,2, raw = TRUE), se = FALSE)

####

## Testing correlation between Survey Year and Elevation

## Visualize Year*Elevation correlation

ggplot(model_data, aes(x = Elevation, fill = ProjectID)) +
  geom_density(alpha = 0.3) +
  ylab("Density of projects occurences across elevation") +
  xlab("Elevation (m)") +
  theme_classic()

ggplot(model_data_hist, aes(x = Year, y = Elevation, colour = ProjectID)) +
  geom_point(
    size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  ylab("Elevation (m)") +
  theme_classic()

model_data_hist <- model_data %>%
  filter(Year < 2025) ## removing 2025 surveys

cor.test(model_data_hist$Year, model_data_hist$Elevation)

####

## Making graphs and results info

####

## Quantifying Upslope Range Shift using ggpredict

## mod8

pred8 <- ggpredict(
  mod8_full,
  terms = c("Elevation_sc [all]", "Time_sc [-2,0,2]"),
  bias_correction = TRUE, 
)

plot(pred8) + theme_classic()

## OR adjust to show re-scaled elevation

mean_elev <- mean(model_data$Elevation)
sd_elev   <- sd(model_data$Elevation)

## mod8

pred8$Elevation_m <- pred8$x * sd_elev + mean_elev

optima8 <- pred8 %>%
  group_by(group) %>%
  slice_max(predicted)

ggplot(pred8, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA
  ) +
  scale_x_continuous(breaks = seq(floor(min(pred8$Elevation_m) / 100) * 100, ceiling(max(pred8$Elevation_m) / 100) * 100, by = 100)) +
  scale_y_continuous(breaks = seq(floor(min(pred8$predicted) / 0.1) * 0.1, ceiling(max(pred8$predicted) / 0.1 * 0.1), by = 0.1)) +
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

####

## Extracting Species-Specific Responses:

ranef(mod3_full)
ranef_slopes <- ranef(mod3_full)$Species[, "Elevation_sc"]

hist(ranef_slopes)
abline(v = 0, col = "red")

