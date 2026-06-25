#### Bayesian relocation error-informed model will ALL species (n=136)

## Load libraries

library(tidyverse)
library(brms)
library(posterior)
library(ggeffects)
library(ggplot2)


set.seed(1234)

## Read in data

model_data <- read_csv("bayesian_analysis_scripts/model_data_all_species.csv")

## Read in models to avoid re-running if needed

bayesmod_error <- readRDS("bayesian_analysis_scripts/bayesmod_error_all_species.rds")

## Prepare the full model data

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
    Time = factor(Time, levels = c("historical", "present")),
    Functional_group = factor(Functional_group)
  )

####

## Extract priors

error_priors <- fixef(bayesmod_error)[c(
  "TreatmentError",
  "TreatmentError:polyElevation_sc2rawEQTRUE1",
  "TreatmentError:polyElevation_sc2rawEQTRUE2"
), ]

error_priors

## Saving each mean and SD as  objects

trt_mu  <- error_priors["TreatmentError", "Estimate"]
trt_sd  <- error_priors["TreatmentError", "Est.Error"]

elev_mu <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE1", "Estimate"]
elev_sd <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE1", "Est.Error"]

elev2_mu <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE2", "Estimate"]
elev2_sd <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE2", "Est.Error"]

## Get coefficient names and set priors

get_prior(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = bernoulli(link = "logit")
)

priors_error <- c(
  set_prior(paste0("normal(", trt_mu, ", ", trt_sd, ")"),
            class = "b", coef = "Timepresent"),
  
  set_prior(paste0("normal(", elev_mu, ", ", elev_sd, ")"),
            class = "b", coef = "Timepresent:polyElevation_sc2rawEQTRUE1"),
  
  set_prior(paste0("normal(", elev2_mu, ", ", elev2_sd, ")"),
            class = "b", coef = "Timepresent:polyElevation_sc2rawEQTRUE2")
)

priors_error

####

## Run full model with priors

bayesmod_full <- brm(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = bernoulli(link = "logit"),
  prior = priors_error,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 1234
)

summary(bayesmod_full)
pp_check(bayesmod_full)
prior_summary(bayesmod_full)

####

## Bayesian with no priors

bayesmod_full_no_prior <- brm(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 1234
)

summary(bayesmod_full_no_prior)
pp_check(bayesmod_full_no_prior)
prior_summary(bayesmod_full_no_prior)

####

## Save models
saveRDS(bayesmod_error, "bayesian_analysis_scripts/bayesmod_error_all_species.rds")
saveRDS(bayesmod_full, "bayesian_analysis_scripts/bayesmod_full_all_species.rds")

####

## Visualize results with predict

pred <- predict_response(
  bayesmod_full,
  terms = c("Elevation_sc [all]", "Time"),
  type = "fixed",
  interval = "confidence"
)

head(pred)

## Adjust to show re-scaled elevation
mean_elev <- mean(model_data$Elevation)
sd_elev   <- sd(model_data$Elevation)
pred$Elevation_m <- pred$x * sd_elev + mean_elev

## Calculate optima and shift of predicted probability of presence
optima <- pred%>%
  group_by(group) %>%
  slice_max(predicted)

shift <- optima %>%
  select(group, Elevation_m) %>%
  pivot_wider(names_from = group, values_from = Elevation_m) %>%
  mutate(shift_m = present - historical)

shift
#46.2 m

## Pretty plot

bayes_plot <- ggplot(pred, aes(x = Elevation_m, y = predicted, colour = group)) +
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
    fill = "Time") +
  geom_vline(
    data = optima,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
    data = optima,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  annotate("text", x = Inf, y = Inf, label = "+46.2 m", hjust = 1.5, vjust = 4, size = 5) +
  theme_classic()

bayes_plot
