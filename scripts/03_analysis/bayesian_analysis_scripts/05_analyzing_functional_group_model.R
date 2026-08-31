#### Bayesian functional group interaction model

## Load libraries

library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(ggeffects)
library(ggplot2)

set.seed(1234)

## Read in data

error_model_data <- read_csv("data/processed/error_model_data.csv")
model_data <- read_csv("data/processed/model_data_filt.csv")

## Check data structure

str(error_model_data)
str(model_data)

## Prepare error model data

error_model_data <- error_model_data %>%
  mutate(
    Species = factor(Species),
    PlotNumber = factor(PlotNumber),
    Treatment = factor(Treatment, levels = c("2025", "Error")),
    Functional_group = factor(Functional_group),
    Elevation_sc = as.numeric(scale(Elevation)),
    Lat_sc = as.numeric(scale(Latitude))
  )

## Prepare full model data

model_data <- model_data %>%
  mutate(
    Species = factor(Species),
    ProjectID = factor(ProjectID),
    PlotNumber = factor(PlotNumber),
    SiteSurveyor = factor(SiteSurveyor),
    RelocationConfidence = factor(RelocationConfidence),
    Time = factor(Time, levels = c("historical", "present")),
    Functional_group = factor(Functional_group),
    Elevation_sc = as.numeric(scale(Elevation)),
    Time_sc = as.numeric(scale(YearsSinceStart)),
    Lat_sc = as.numeric(scale(Latitude)),
    Lon_sc = as.numeric(scale(Longitude))
  )

####

## Run relocation error model with functional group interaction

bayesmod_error_fg <- brm(
  present ~ Treatment * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = error_model_data,
  family = bernoulli(link = "logit"),
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 1234
)

summary(bayesmod_error_fg)
pp_check(bayesmod_error_fg)

## Extract error-model coefficients

error_priors_fg <- fixef(bayesmod_error_fg)
error_priors_fg

## Get full model coefficient names

get_prior(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = bernoulli(link = "logit")
)

## Identify functional group levels

levels(model_data$Functional_group)

## Identify treatment and treatment-by-elevation coefficients

grep(
  "^TreatmentError$|^TreatmentError:polyElevation",
  rownames(error_priors_fg),
  value = TRUE
)

## Identify functional group interaction coefficients

grep(
  "TreatmentError.*Functional_group",
  rownames(error_priors_fg),
  value = TRUE
)

## Extract treatment and treatment-by-elevation priors

trt_mu <- error_priors_fg["TreatmentError", "Estimate"]
trt_sd <- error_priors_fg["TreatmentError", "Est.Error"]

elev_mu <- error_priors_fg["TreatmentError:polyElevation_sc2rawEQTRUE1", "Estimate"]
elev_sd <- error_priors_fg["TreatmentError:polyElevation_sc2rawEQTRUE1", "Est.Error"]

elev2_mu <- error_priors_fg["TreatmentError:polyElevation_sc2rawEQTRUE2", "Estimate"]
elev2_sd <- error_priors_fg["TreatmentError:polyElevation_sc2rawEQTRUE2", "Est.Error"]

## Set priors for treatment and treatment-by-elevation effects

priors_error_fg <- c(
  set_prior(
    paste0("normal(", trt_mu, ", ", trt_sd, ")"),
    class = "b",
    coef = "Timepresent"
  ),
  set_prior(
    paste0("normal(", elev_mu, ", ", elev_sd, ")"),
    class = "b",
    coef = "Timepresent:polyElevation_sc2rawEQTRUE1"
  ),
  set_prior(
    paste0("normal(", elev2_mu, ", ", elev2_sd, ")"),
    class = "b",
    coef = "Timepresent:polyElevation_sc2rawEQTRUE2"
  )
)

## Extract functional group interaction priors

fg_error_priors <- error_priors_fg[
  grep(
    "TreatmentError.*Functional_group",
    rownames(error_priors_fg)
  ),
]

fg_error_priors

## Add functional group interaction priors

for (parameter in rownames(fg_error_priors)) {
  
  prior_mean <- fg_error_priors[parameter, "Estimate"]
  prior_sd <- fg_error_priors[parameter, "Est.Error"]
  
  full_parameter <- gsub(
    "TreatmentError",
    "Timepresent",
    parameter
  )
  
  priors_error_fg <- c(
    priors_error_fg,
    set_prior(
      paste0("normal(", prior_mean, ", ", prior_sd, ")"),
      class = "b",
      coef = full_parameter))}

## View priors

priors_error_fg

####

## Run full model with functional group interaction and priors

bayesmod_full_fg <- brm(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = bernoulli(link = "logit"),
  prior = priors_error_fg,
  chains = 4,
  iter = 4000,
  warmup = 1000,
  cores = 4,
  seed = 1234
)

summary(bayesmod_full_fg)
pp_check(bayesmod_full_fg)
prior_summary(bayesmod_full_fg)

## Run full model without informative priors

bayesmod_full_fg_no_prior <- brm(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
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

summary(bayesmod_full_fg_no_prior)
pp_check(bayesmod_full_fg_no_prior)
prior_summary(bayesmod_full_fg_no_prior)

####

## Save models

saveRDS(bayesmod_error_fg, "outputs/models/bayesian_analysis_models/bayesmod_error_fg_30_spp.rds")
saveRDS(bayesmod_full_fg, "outputs/models/bayesian_analysis_models/bayesmod_full_fg_30_spp.rds")
saveRDS(bayesmod_full_fg_no_prior, "outputs/models/bayesian_analysis_models/bayesmod_full_fg_no_prior_30_spp.rds")

####

## Model diagnostics

summary(bayesmod_full_fg)
summary(bayesmod_error_fg)
summary(bayesmod_full_fg_no_prior)

## Plot models

plot(bayesmod_full_fg)
plot(bayesmod_error_fg)
plot(bayesmod_full_fg_no_prior)

## Posterior predictive checks

pp_check(bayesmod_full_fg, ndraws = 1000)
pp_check(bayesmod_error_fg, ndraws = 1000)
pp_check(bayesmod_full_fg_no_prior, ndraws = 1000)

## Posterior predictive check of the mean

pp_check(bayesmod_full_fg, type = "stat", stat = "mean", ndraws = 1000)
pp_check(bayesmod_error_fg, type = "stat", stat = "mean", ndraws = 1000)
pp_check(bayesmod_full_fg_no_prior, type = "stat", stat = "mean", ndraws = 1000)

## Bayes R2

bayes_R2(bayesmod_full_fg)
bayes_R2(bayesmod_error_fg)
bayes_R2(bayesmod_full_fg_no_prior)

## Check divergences

nuts_full_fg <- nuts_params(bayesmod_full_fg)
nuts_error_fg <- nuts_params(bayesmod_error_fg)
nuts_full_fg_no_prior <- nuts_params(bayesmod_full_fg_no_prior)

sum(nuts_full_fg$Value[nuts_full_fg$Parameter == "divergent__"])
sum(nuts_error_fg$Value[nuts_error_fg$Parameter == "divergent__"])
sum(nuts_full_fg_no_prior$Value[nuts_full_fg_no_prior$Parameter == "divergent__"])

####

## Plotting priors and posterior distributions

## Extract error-model coefficients
error_priors_fg <- fixef(bayesmod_error_fg)

## Identify treatment-related coefficients
fg_prior_parameters <- rownames(error_priors_fg)[
  grepl("^TreatmentError", rownames(error_priors_fg))
]

## Match error-model parameters to full-model parameters
posterior_parameters_fg <- str_replace(
  fg_prior_parameters,
  "^TreatmentError",
  "Timepresent"
)

## Extract posterior draws
posterior_draws_fg <- as_draws_df(bayesmod_full_fg) %>%
  select(all_of(paste0("b_", posterior_parameters_fg)))

## Rename posterior columns to match prior parameters
names(posterior_draws_fg) <- posterior_parameters_fg

## Generate prior draws
prior_draws_fg <- map_dfc(
  fg_prior_parameters,
  ~ rnorm(
    12000,
    mean = error_priors_fg[.x, "Estimate"],
    sd = error_priors_fg[.x, "Est.Error"]))

names(prior_draws_fg) <- posterior_parameters_fg

## Convert to long format
prior_draws_fg <- prior_draws_fg %>%
  pivot_longer(
    everything(),
    names_to = "Parameter",
    values_to = "Estimate") %>%
  mutate(Distribution = "Prior")

posterior_draws_fg <- posterior_draws_fg %>%
  pivot_longer(
    everything(),
    names_to = "Parameter",
    values_to = "Estimate") %>%
  mutate(Distribution = "Posterior")

## Combine
prior_posterior_fg <- bind_rows(
  prior_draws_fg,
  posterior_draws_fg)

## Rename parameters for plotting
prior_posterior_fg <- prior_posterior_fg %>%
  mutate(
    Parameter = case_when(
      Parameter == "Timepresent" ~ "Timepresent",
      Parameter == "Timepresent:polyElevation_sc2rawEQTRUE1" ~
        "Timepresent × Elevation",
      Parameter == "Timepresent:polyElevation_sc2rawEQTRUE2" ~
        "Timepresent × Elevation²",
      Parameter == "Timepresent:Functional_groupshrub" ~
        "Timepresent × Shrub",
      Parameter == "Timepresent:Functional_grouptree" ~
        "Timepresent × Tree",
      Parameter == "Timepresent:polyElevation_sc2rawEQTRUE1:Functional_groupshrub" ~
        "Timepresent × Elevation × Shrub",
      Parameter == "Timepresent:polyElevation_sc2rawEQTRUE2:Functional_groupshrub" ~
        "Timepresent × Elevation² × Shrub",
      Parameter == "Timepresent:polyElevation_sc2rawEQTRUE1:Functional_grouptree" ~
        "Timepresent × Elevation × Tree",
      Parameter == "Timepresent:polyElevation_sc2rawEQTRUE2:Functional_grouptree" ~
        "Timepresent × Elevation² × Tree",
      TRUE ~ Parameter
    )
  )

## Plot
ggplot(
  prior_posterior_fg,
  aes(x = Estimate, fill = Distribution, color = Distribution)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  facet_wrap(~ Parameter, scales = "free", ncol = 3) +
  xlab("Coefficient estimate") +
  ylab("Density") +
  theme_classic()

####

## Visualizing historical vs present by functional group

## Get elevation scaling parameters

elevation_mean <- mean(model_data$Elevation)
elevation_sd <- sd(model_data$Elevation)

## Create prediction data

newdata_fg <- expand_grid(
  Elevation = seq(
    min(model_data$Elevation),
    max(model_data$Elevation),
    length.out = 1000),
  Time = factor(
    c("historical", "present"),
    levels = c("historical", "present")),
  Functional_group = levels(model_data$Functional_group)) %>%
  mutate(
    Elevation_sc = (Elevation - elevation_mean) / elevation_sd,
    Lat_sc = 0)

## Generate posterior predictions

posterior_predictions_fg <- posterior_epred(
  bayesmod_full_fg,
  newdata = newdata_fg,
  re_formula = NA)

## Calculate posterior medians and 95% credible intervals

prediction_summary_fg <- apply(
  posterior_predictions_fg,
  2,
  quantile,
  probs = c(0.025, 0.5, 0.975)) %>%
  t() %>%
  as.data.frame() %>%
  rename(
    lower = `2.5%`,
    estimate = `50%`,
    upper = `97.5%`)

## Add predictions to new data

prediction_data_fg <- bind_cols(
  newdata_fg,
  prediction_summary_fg)

## Plot historical vs present curves by functional group

ggplot(
  prediction_data_fg,
  aes(
    x = Elevation,
    y = estimate,
    linetype = Time,
    colour = Time,
    fill = Time)) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper),
    alpha = 0.2,
    colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ Functional_group) +
  xlab("Elevation (m)") +
  ylab("Predicted probability of occurrence") +
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
  scale_colour_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4")) +
  scale_fill_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4")) +
  scale_linetype_manual(
    values = c(
      "historical" = "solid",
      "present" = "solid")) +
  theme_classic()

####

## Calculating range shifts by functional group

## Find the elevation of maximum predicted probability for each posterior draw

range_shift_draws_fg <- map_dfr(
  1:nrow(posterior_predictions_fg),
  function(i) {
    
    predictions <- posterior_predictions_fg[i, ]
    
    map_dfr(
      levels(model_data$Functional_group),
      function(fg) {
        
        historical_predictions <- predictions[
          newdata_fg$Time == "historical" &
            newdata_fg$Functional_group == fg
        ]
        
        present_predictions <- predictions[
          newdata_fg$Time == "present" &
            newdata_fg$Functional_group == fg
        ]
        
        historical_elevations <- newdata_fg$Elevation[
          newdata_fg$Time == "historical" &
            newdata_fg$Functional_group == fg
        ]
        
        present_elevations <- newdata_fg$Elevation[
          newdata_fg$Time == "present" &
            newdata_fg$Functional_group == fg
        ]
        
        historical_optimum <- historical_elevations[
          which.max(historical_predictions)
        ]
        
        present_optimum <- present_elevations[
          which.max(present_predictions)
        ]
        
        tibble(
          Functional_group = fg,
          historical_optimum = historical_optimum,
          present_optimum = present_optimum,
          shift_m = present_optimum - historical_optimum
        )
      }
    )
  }
)

## Summarize Bayesian range shifts by functional group

range_shift_summary_fg <- range_shift_draws_fg %>%
  group_by(Functional_group) %>%
  summarise(
    median_historical_optimum = median(historical_optimum),
    lower_95_historical_optimum = quantile(historical_optimum, 0.025),
    upper_95_historical_optimum = quantile(historical_optimum, 0.975),
    median_present_optimum = median(present_optimum),
    lower_95_present_optimum = quantile(present_optimum, 0.025),
    upper_95_present_optimum = quantile(present_optimum, 0.975),
    median_shift = median(shift_m),
    lower_95_shift = quantile(shift_m, 0.025),
    upper_95_shift = quantile(shift_m, 0.975),
    probability_positive = mean(shift_m > 0)
  )

range_shift_summary_fg

## Visualize posterior distribution of range shifts

ggplot(
  range_shift_draws_fg,
  aes(x = shift_m)
) +
  geom_density() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  facet_wrap(~ Functional_group) +
  xlab("Range shift (m)") +
  ylab("Posterior density") +
  theme_classic()