#### Exploring the Bayesian model

## Load libraries
library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(ggeffects)
library(ggplot2)

## Read in models
bayesmod_full <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_full_30_30_spp.rds")
bayesmod_error <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_error_30_spp.rds")
bayesmod_full_no_prior <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_full_no_prior_30_spp.rds")

## Model summaries
summary(bayesmod_full)
summary(bayesmod_error)
summary(bayesmod_full_no_prior)

## Plot models

plot(bayesmod_full)
plot(bayesmod_error)
plot(bayesmod_full_no_prior)

## pp checks
pp_check(bayesmod_full, ndraws = 1000)
pp_check(bayesmod_error, ndraws = 1000)
pp_check(bayesmod_full_no_prior, ndraws = 1000)

## Posterior predictive check of the mean (ratio of present vs absent)
pp_check(bayesmod_full, type = "stat", stat = "mean", ndraws = 1000)
pp_check(bayesmod_error, type = "stat", stat = "mean", ndraws = 1000)
pp_check(bayesmod_full_no_prior, type = "stat", stat = "mean", ndraws = 1000)

## Bayes R^2
bayes_R2(bayesmod_full)
bayes_R2(bayesmod_error)
bayes_R2(bayesmod_full_no_prior)

## Check divergences

nuts_full <- nuts_params(bayesmod_full)
nuts_error <- nuts_params(bayesmod_error)
nuts_full_no_prior <- nuts_params(bayesmod_full_no_prior)

sum(nuts_full$Value[nuts_full$Parameter == "divergent__"])
#[1] 0
sum(nuts_error$Value[nuts_error$Parameter == "divergent__"])
#[1] 0
sum(nuts_full_no_prior$Value[nuts_full_no_prior$Parameter == "divergent__"])
#[1] 0

####

## Plotting priors and posterior distributions

## Extract priors from relocation error model

error_priors <- fixef(bayesmod_error)[c(
  "TreatmentError",
  "TreatmentError:polyElevation_sc2rawEQTRUE1",
  "TreatmentError:polyElevation_sc2rawEQTRUE2"),]

error_priors

## Save each mean and SD as objects

trt_mu <- error_priors["TreatmentError", "Estimate"]
trt_sd <- error_priors["TreatmentError", "Est.Error"]

elev_mu <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE1", "Estimate"]
elev_sd <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE1", "Est.Error"]

elev2_mu <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE2", "Estimate"]
elev2_sd <- error_priors["TreatmentError:polyElevation_sc2rawEQTRUE2", "Est.Error"]

## Generate draws from the prior distributions

prior_draws <- tibble(
  Timepresent = rnorm(12000, mean = trt_mu, sd = trt_sd),
  `Timepresent:polyElevation_sc2rawEQTRUE1` = rnorm(12000, mean = elev_mu, sd = elev_sd),
  `Timepresent:polyElevation_sc2rawEQTRUE2` = rnorm(12000, mean = elev2_mu, sd = elev2_sd))

## Extract posterior draws from full model

variables(bayesmod_full)

posterior_draws <- as_draws_df(bayesmod_full) %>%
  select(
    b_Timepresent,
    `b_Timepresent:polyElevation_sc2rawEQTRUE1`,
    `b_Timepresent:polyElevation_sc2rawEQTRUE2`  )

## Convert prior draws to long format

prior_draws <- prior_draws %>%
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Estimate") %>%
  mutate(Distribution = "Prior")

## Convert posterior draws to long format

posterior_draws <- posterior_draws %>%
  pivot_longer(
    cols = everything(),
    names_to = "Parameter",
    values_to = "Estimate") %>%
  mutate(Distribution = "Posterior")

## Combine prior and posterior draws

prior_posterior <- bind_rows(prior_draws, posterior_draws)

## Rename parameters

prior_posterior <- prior_posterior %>%
  mutate(
    Parameter = case_when(
      Parameter %in% c(
        "Timepresent", "b_Timepresent") ~ "Timepresent",
      Parameter %in% c(
        "Timepresent:polyElevation_sc2rawEQTRUE1",
        "b_Timepresent:polyElevation_sc2rawEQTRUE1") ~ "Timepresent × Elevation",
      Parameter %in% c(
        "Timepresent:polyElevation_sc2rawEQTRUE2",
        "b_Timepresent:polyElevation_sc2rawEQTRUE2") ~ "Timepresent × Elevation²", TRUE ~ Parameter))

## Plot prior and posterior distributions

prior_posterior_dist <- 
  ggplot(prior_posterior,
  aes(x = Estimate, fill = Distribution, color = Distribution)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ Parameter) +
  xlab("Coefficient estimate") +
  ylab("Density") +
  theme_classic()

prior_posterior_dist

####

## Visualizing historical vs. present

model_data <- read_csv("data/processed/model_data_filt.csv")

## Get elevation scaling parameters

elevation_mean <- mean(model_data$Elevation)
elevation_sd <- sd(model_data$Elevation)

## Create prediction data

newdata <- expand_grid(
  Elevation = seq(
    min(model_data$Elevation),
    max(model_data$Elevation),
    length.out = 1000),
  Time = factor(
    c("historical", "present"),
    levels = c("historical", "present"))) %>%
  mutate(
    Elevation_sc = (Elevation - elevation_mean) / elevation_sd,
    Lat_sc = 0)

## Generate posterior predictions

posterior_predictions <- posterior_epred(
  bayesmod_full,
  newdata = newdata,
  re_formula = NA)

## Calculate posterior medians and 95% credible intervals

prediction_data <- bind_cols(
  newdata,
  t(apply(
    posterior_predictions,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975))) %>%
    as.data.frame() %>%
    rename(
      lower = `2.5%`,
      estimate = `50%`,
      upper = `97.5%`))

## Plot historical vs. present curves

bayes_mod <- ggplot(
  prediction_data,
  aes(
    x = Elevation,
    y = estimate,
    colour = Time)) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper,
      fill = Time),
    alpha = 0.2,
    colour = NA) +
  geom_line(
    linewidth = 1) +
  xlab("Elevation (m)") +
  ylab("Predicted probability of occurrence") +
  scale_colour_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4")) +
  scale_fill_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4")) +
  theme_classic()

bayes_mod

ggsave("outputs/figures/bayesian_figures/bayesmod_elevation_curve.png",plot = bestmod_plot, width = 7, height = 4, dpi = 300)


####

## Generate alternative predictions

freq_curve <- ggpredict(
  bayesmod_full,
  terms = c("Elevation_sc [all]", "Time")) %>%
  mutate(
    Elevation_m = x * elevation_sd + elevation_mean)

## Plot Bayesian and alternative prediction curves

ggplot() +
  geom_ribbon(
    data = prediction_data,
    aes(
      x = Elevation,
      ymin = lower,
      ymax = upper,
      fill = Time),
    alpha = 0.15) +
  geom_line(
    data = prediction_data,
    aes(
      x = Elevation,
      y = estimate,
      colour = Time),
    linewidth = 1) +
  geom_line(
    data = freq_curve,
    aes(
      x = Elevation_m,
      y = predicted,
      colour = group),
    linewidth = 1) +
  xlab("Elevation (m)") +
  ylab("Predicted probability of occurrence") +
  labs(
    colour = "Time",
    fill = "Time") +
  scale_colour_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4")) +
  scale_fill_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4")) +
  theme_classic()

####

## Calculating the range shift

## Find elevation of maximum predicted probability for each posterior draw

range_shift_draws <- map_dfr(
  1:nrow(posterior_predictions),
  function(i) {
    
    predictions <- posterior_predictions[i, ]
    
    historical_predictions <- predictions[
      newdata$Time == "historical"]
    
    present_predictions <- predictions[
      newdata$Time == "present"]
    
    historical_optimum <- newdata$Elevation[
      newdata$Time == "historical"
    ][which.max(historical_predictions)]
    
    present_optimum <- newdata$Elevation[
      newdata$Time == "present"
    ][which.max(present_predictions)]
    
    tibble(
      historical_optimum = historical_optimum,
      present_optimum = present_optimum,
      shift_m = present_optimum - historical_optimum)
  })

## Summarize Bayesian range shift

range_shift_summary <- range_shift_draws %>%
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
    probability_positive = mean(shift_m > 0)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Statistic",
    values_to = "Value")

range_shift_summary

## Visualize posterior distribution of range shift

ggplot(
  range_shift_draws,
  aes(x = shift_m)) +
  geom_density() +
  geom_vline(
    xintercept = 0,
    linetype = "dashed") +
  xlab("Range shift (m)") +
  ylab("Posterior density") +
  theme_classic()

## Calculate the posterior max

historical_density <- density(range_shift_draws$historical_optimum)
present_density <- density(range_shift_draws$present_optimum)

historical_max <- historical_density$x[
  which.max(historical_density$y)]

present_max <- present_density$x[
  which.max(present_density$y)]

## Plot posterior distributions with optima

bayes_optima <- ggplot(
  range_shift_draws) +
  geom_density(
    aes(
      x = historical_optimum,
      colour = "Historical"),
    linewidth = 1) +
  geom_density(
    aes(
      x = present_optimum,
      colour = "Present"),
    linewidth = 1) +
  annotate(
    "segment",
    x = historical_max,
    xend = historical_max,
    y = 0,
    yend = max(historical_density$y),
    colour = "#00BFC4",
    linetype = "dashed",
    linewidth = 0.8) +
  annotate(
    "point",
    x = historical_max,
    y = max(historical_density$y),
    colour = "#00BFC4",
    size = 3) +
  annotate(
    "segment",
    x = present_max,
    xend = present_max,
    y = 0,
    yend = max(present_density$y),
    colour = "#F8766D",
    linetype = "dashed",
    linewidth = 0.8) +
  annotate(
    "point",
    x = present_max,
    y = max(present_density$y),
    colour = "#F8766D",
    size = 3) +
  xlab("Elevation of predicted optimum (m)") +
  ylab("Posterior density") +
  scale_colour_manual(
    values = c(
      "Historical" = "#00BFC4",
      "Present" = "#F8766D")) +
  scale_x_continuous(
    breaks = seq(
      floor(min(range_shift_draws$historical_optimum) / 200) * 200,
      ceiling(max(range_shift_draws$present_optimum) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(range_shift_draws$historical_optimum) / 100) * 100,
      ceiling(max(range_shift_draws$present_optimum) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE))+
  scale_y_continuous(
    breaks = seq(0, 0.008, by = 0.002),
    minor_breaks = seq(0, 0.008, by = 0.001),
    expand = expansion(mult = c(0, 0.08)),
    guide = guide_axis(minor.ticks = TRUE)) +
  theme_classic()

bayes_optima

ggsave("ouputs/figures/bayesian_figures/bayesian_optima.png", plot = bayes_optima, width = 8, height = 6, dpi = 300)

####

## Comparing with the model with no priors

## Generate posterior predictions

posterior_predictions_np <- posterior_epred(
  bayesmod_full_no_prior,
  newdata = newdata,
  re_formula = NA)

## Calculate posterior medians and 95% credible intervals

prediction_data_np <- bind_cols(
  newdata,
  t(apply(
    posterior_predictions_np,
    2,
    quantile,
    probs = c(0.025, 0.5, 0.975))) %>%
    as.data.frame() %>%
    rename(
      lower = `2.5%`,
      estimate = `50%`,
      upper = `97.5%`))

## Plot model with no priors

ggplot(
  prediction_data_np,
  aes(
    x = Elevation,
    y = estimate,
    colour = Time)) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper,
      fill = Time),
    alpha = 0.2,
    colour = NA) +
  geom_line(
    linewidth = 1) +
  xlab("Elevation (m)") +
  ylab("Predicted probability of occurrence") +
  scale_colour_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4")) +
  scale_fill_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4")) +
  theme_classic()


