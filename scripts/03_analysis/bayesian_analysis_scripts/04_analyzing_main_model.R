#### Exploring the Bayesian model

## Load libraries
library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(ggeffects)
library(ggplot2)
library(patchwork)

set.seed(1234)

## Read in models
bayesmod_full <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_full_30_30_spp.rds")
bayesmod_error <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_error_30_spp.rds")
bayesmod_full_no_prior <- readRDS("outputs/models/bayesian_analysis_models/bayesmod_full_no_prior_30_spp.rds")

## Model summaries
summary(bayesmod_full)
summary(bayesmod_error)
summary(bayesmod_full_no_prior)

## Plot models
#bayesmod_full_plot <- plot(bayesmod_full) #+ labs(title = "Bayesian full model")
#bayesmod_error_plot <- plot(bayesmod_error) + labs(title = "Bayesian relocation error model")
#bayesmod_full_no_prior_plot <- plot(bayesmod_full_no_prior) + labs(title = "Bayesian full model without informative priors")

## Posterior predictive checks
pp_full <- pp_check(bayesmod_full, ndraws = 100) + 
  labs(x = "Proportion of species present", y = "Density")
pp_error <- pp_check(bayesmod_error, ndraws = 100) +
  labs(title = "Posterior predictive check: relocation error model", x = "Proportion of species present", y = "Density")
pp_full_no_prior <- pp_check(bayesmod_full_no_prior, ndraws = 100) +
  labs(title = "Posterior predictive check: full model without informative priors", x = "Proportion of species present", y = "Density")

pp_full
pp_error
pp_full_no_prior

## Save posterior predictive checks
ggsave("outputs/figures/bayesian_figures/pp_check_full.png", pp_full, width = 8, height = 6, dpi = 300)
ggsave("outputs/figures/bayesian_figures/pp_check_error.png", pp_error, width = 8, height = 6, dpi = 300)
ggsave("outputs/figures/bayesian_figures/pp_check_full_no_prior.png", pp_full_no_prior, width = 8, height = 6, dpi = 300)

## Posterior predictive checks of the mean
pp_mean_full <- pp_check(bayesmod_full, type = "stat", stat = "mean", ndraws = 1000) + 
  labs(x = "Proportion of species present")
pp_mean_error <- pp_check(bayesmod_error, type = "stat", stat = "mean",ndraws = 1000) + 
  labs(title = "Posterior predictive check of the mean: relocation error model", x = "Proportion of species present")
pp_mean_full_no_prior <- pp_check(bayesmod_full_no_prior, type = "stat", stat = "mean", ndraws = 1000) + 
  labs(title = "Posterior predictive check of the mean: full model without informative priors", x = "Proportion of species present")

pp_mean_full
pp_mean_error
pp_mean_full_no_prior

## Save mean posterior predictive checks
ggsave("outputs/figures/bayesian_figures/pp_mean_full.png", pp_mean_full, width = 8, height = 6, dpi = 300)
ggsave("outputs/figures/bayesian_figures/pp_mean_error.png", pp_mean_error, width = 8, height = 6, dpi = 300)
ggsave("outputs/figures/bayesian_figures/pp_mean_full_no_prior.png", pp_mean_full_no_prior, width = 8, height = 6, dpi = 300)

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
    `b_Timepresent:polyElevation_sc2rawEQTRUE2`)

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
  labs(tag = "A") +
  theme_classic()

prior_posterior_dist

#ggsave("outputs/figures/bayesian_figures/prior_posterior_dist.png", prior_posterior_dist, width = 8, height = 6, dpi = 300)

## Make a zoomed in version of the prior/ posterior graph and combine as panels

panel_b <- ggplot(
  prior_posterior,
  aes(
    x = Estimate,
    fill = Distribution,
    colour = Distribution)) +
  geom_density(alpha = 0.3) +
  facet_wrap(~ Parameter) +
  coord_cartesian(xlim = c(-1, 1)) +
  xlab("Coefficient estimate") +
  ylab("Density") +
  labs(tag = "B") +
  theme_classic()

panel_b

## Combine panels

prior_posterior_combined <- prior_posterior_dist / panel_b
prior_posterior_combined

ggsave("outputs/figures/bayesian_figures/prior_posterior_dist_zoom.png", prior_posterior_combined, width = 7, height = 9, dpi = 300)


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

bayesmod_plot <- ggplot(
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
  geom_vline(
    xintercept = historical_median,
    colour = "#00BFC4",
    linetype = "dashed",
    linewidth = 0.8) +
  annotate(
    "point",
    x = historical_median,
    y = approx(
      prediction_data$Elevation[
        prediction_data$Time == "historical"],
      prediction_data$estimate[
        prediction_data$Time == "historical"],
      xout = historical_median)$y,
    colour = "#00BFC4",
    size = 3) +
  geom_vline(
    xintercept = present_median,
    colour = "#F8766D",
    linetype = "dashed",
    linewidth = 0.8) +
  annotate(
    "point",
    x = present_median,
    y = approx(
      prediction_data$Elevation[
        prediction_data$Time == "present"],
      prediction_data$estimate[
        prediction_data$Time == "present"],
      xout = present_median)$y,
    colour = "#F8766D",
    size = 3) +
  xlab("Elevation (m)") +
  ylab("Predicted probability of occurrence") +
  scale_colour_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4"),
    labels = c(
      "present" = "Present",
      "historical" = "Historical")) +
  scale_fill_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4"),
    guide = "none") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = "+55 m",
    hjust = 1.1,
    vjust = 1.5,
    size = 4) +
  scale_x_continuous(
    breaks = seq(600, 1800, by = 100),
    minor_breaks = seq(600, 1800, by = 50),
    guide = guide_axis(minor.ticks = TRUE)) +
  theme_classic()

bayesmod_plot

ggsave("outputs/figures/bayesian_figures/bayesmod_elevation_curve.png",plot = bayesmod_plot, width = 7, height = 4, dpi = 300)


####

## Generate alternative predictions

freq_curve <- ggpredict(
  bayesmod_full,
  terms = c("Elevation_sc [all]", "Time")) %>%
  mutate(
    Elevation_m = x * elevation_sd + elevation_mean)

## Plot Bayesian and alternative prediction curves

gg_predict_curve <- ggplot() +
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

gg_predict_curve

ggsave("outputs/figures/bayesian_figures/gg_predict_bayesmod_elevation_curve.png",plot = gg_predict_curve, width = 7, height = 4, dpi = 300)


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

density_shift <- density(range_shift_draws$shift_m)

density_data <- tibble(
  x = density_shift$x,
  y = density_shift$y
)

ci_density <- density_data %>%
  filter(
    x >= 6.97,
    x <= 127)

median_density <- approx(
  density_data$x,
  density_data$y,
  xout = 54.6)$y

bayes_range_shift <- ggplot(
  density_data,
  aes(x = x, y = y)) +
  geom_ribbon(
    data = ci_density,
    aes(
      ymin = 0,
      ymax = y),
    alpha = 0.2) +
  geom_line() +
  annotate(
    "segment",
    x = 55,
    xend = 55,
    y = 0,
    yend = median_density,
    linetype = "dashed") +
  annotate(
    "point",
    x = 55,
    y = median_density,
    size = 3) +
  annotate(
    "text",
    x = 55,
    y = median_density,
    label = "55 m",
    vjust = -1,
    size = 4) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed") +
  xlab("Range shift (m)") +
  ylab("Posterior density") +
  #labs(title = "Posterior distribution of elevation range shift") +
  scale_x_continuous(
    limits = c(-50, 250),
    breaks = seq(
      floor(min(range_shift_draws$shift_m) / 100) * 100,
      ceiling(max(range_shift_draws$shift_m) / 100) * 100,
      by = 100),
    minor_breaks = seq(
      floor(min(range_shift_draws$shift_m) / 50) * 50,
      ceiling(max(range_shift_draws$shift_m) / 50) * 50,
      by = 50),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    guide = guide_axis(minor.ticks = TRUE)) +
  theme_classic()

bayes_range_shift

#ggsave("outputs/figures/bayesian_figures/bayes_range_shift.png", plot = bayes_range_shift, width = 7, height = 4, dpi = 300)

## Calculate posterior medians

historical_median <- median(range_shift_draws$historical_optimum)
present_median <- median(range_shift_draws$present_optimum)

historical_density <- density(range_shift_draws$historical_optimum)
present_density <- density(range_shift_draws$present_optimum)

## Plot posterior distributions with optima

bayes_optima <- ggplot(
  range_shift_draws) +
  geom_density(
    aes(
      x = historical_optimum,
      colour = "Historical"),
    linewidth = 1,
    key_glyph = "path") +
  geom_density(
    aes(
      x = present_optimum,
      colour = "Present"),
    linewidth = 1,
    key_glyph = "path") +
  annotate(
    "segment",
    x = historical_median,
    xend = historical_median,
    y = 0,
    yend = approx(
      historical_density$x,
      historical_density$y,
      xout = historical_median)$y,
    colour = "#00BFC4",
    linetype = "dashed",
    linewidth = 0.8) +
  annotate(
    "point",
    x = historical_median,
    y = approx(
      historical_density$x,
      historical_density$y,
      xout = historical_median)$y,
    colour = "#00BFC4",
    size = 3) +
  annotate(
    "segment",
    x = present_median,
    xend = present_median,
    y = 0,
    yend = approx(
      present_density$x,
      present_density$y,
      xout = present_median)$y,
    colour = "#F8766D",
    linetype = "dashed",
    linewidth = 0.8) +
  annotate(
    "point",
    x = present_median,
    y = approx(
      present_density$x,
      present_density$y,
      xout = present_median)$y,
    colour = "#F8766D",
    size = 3) +
  xlab("Elevation of predicted optimum (m)") +
  ylab("Posterior density") +
  scale_colour_manual(
    name = "Time",
    values = c(
      "Historical" = "#00BFC4",
      "Present" = "#F8766D")) +
  annotate(
    "text",
    x = historical_median,
    y = approx(
      historical_density$x,
      historical_density$y,
      xout = historical_median)$y,
    label = paste0(
      round(historical_median, 0),
      " m"),
    vjust = -1,
    colour = "#00BFC4",
    size = 4) +
  annotate(
    "text",
    x = present_median + 30,
    y = approx(
      present_density$x,
      present_density$y,
      xout = present_median)$y,
    label = paste0(
      round(present_median, 0),
      " m"),
    vjust = -1,
    colour = "#F8766D",
    size = 4) +
  scale_x_continuous(
    breaks = seq(
      floor(min(range_shift_draws$historical_optimum) / 200) * 200,
      ceiling(max(range_shift_draws$present_optimum) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(range_shift_draws$historical_optimum) / 100) * 100,
      ceiling(max(range_shift_draws$present_optimum) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(
    breaks = seq(0, 0.008, by = 0.002),
    minor_breaks = seq(0, 0.008, by = 0.001),
    expand = expansion(mult = c(0, 0.08)),
    guide = guide_axis(minor.ticks = TRUE)) +
  theme_classic()

bayes_optima

ggsave("outputs/figures/bayesian_figures/bayesian_optima.png", plot = bayes_optima, width = 8, height = 6, dpi = 300)

####

## Combine range shift curve and range shift distribution

bayes_range_shift2 <- bayes_range_shift +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    guide = guide_axis(minor.ticks = TRUE))
bayes_range_shift2

bayes_combined <- bayesmod_plot / bayes_range_shift2 +
  plot_annotation(tag_levels = "A")

bayes_combined

ggsave("outputs/figures/bayesian_figures/range_shift_combined.png", plot = bayes_combined, width = 7, height = 8, dpi = 300)


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


####

## Calculate range shift for model with no priors

range_shift_draws_np <- map_dfr(
  1:nrow(posterior_predictions_np),
  function(i) {
    
    predictions <- posterior_predictions_np[i, ]
    
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

range_shift_summary_np <- range_shift_draws_np %>%
  summarise(
    median_historical_optimum = median(historical_optimum),
    lower_95_historical_optimum = quantile(
      historical_optimum, 0.025),
    upper_95_historical_optimum = quantile(
      historical_optimum, 0.975),
    median_present_optimum = median(present_optimum),
    lower_95_present_optimum = quantile(
      present_optimum, 0.025),
    upper_95_present_optimum = quantile(
      present_optimum, 0.975),
    median_shift = median(shift_m),
    lower_95_shift = quantile(
      shift_m, 0.025),
    upper_95_shift = quantile(
      shift_m, 0.975),
    probability_positive = mean(shift_m > 0)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Statistic",
    values_to = "Value")

range_shift_summary_np


####

## Calculate posterior medians

historical_median_np <- median(
  range_shift_draws_np$historical_optimum)

present_median_np <- median(
  range_shift_draws_np$present_optimum)

median_shift_np <- median(
  range_shift_draws_np$shift_m)

lower_95_shift_np <- quantile(
  range_shift_draws_np$shift_m, 0.025)

upper_95_shift_np <- quantile(
  range_shift_draws_np$shift_m, 0.975)


####

## Plot historical vs. present curves for model with no priors

bayesmod_plot_np <- ggplot(
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
  geom_vline(
    xintercept = historical_median_np,
    colour = "#00BFC4",
    linetype = "dashed",
    linewidth = 0.8) +
  annotate(
    "point",
    x = historical_median_np,
    y = approx(
      prediction_data_np$Elevation[
        prediction_data_np$Time == "historical"],
      prediction_data_np$estimate[
        prediction_data_np$Time == "historical"],
      xout = historical_median_np)$y,
    colour = "#00BFC4",
    size = 3) +
  geom_vline(
    xintercept = present_median_np,
    colour = "#F8766D",
    linetype = "dashed",
    linewidth = 0.8) +
  annotate(
    "point",
    x = present_median_np,
    y = approx(
      prediction_data_np$Elevation[
        prediction_data_np$Time == "present"],
      prediction_data_np$estimate[
        prediction_data_np$Time == "present"],
      xout = present_median_np)$y,
    colour = "#F8766D",
    size = 3) +
  xlab("Elevation (m)") +
  ylab("Predicted probability of occurrence") +
  scale_colour_manual(
    name = "Time",
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4"),
    labels = c(
      "present" = "Present",
      "historical" = "Historical")) +
  scale_fill_manual(
    values = c(
      "present" = "#F8766D",
      "historical" = "#00BFC4"),
    guide = "none") +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = paste0(
      "+", round(median_shift_np, 0), " m"),
    hjust = 1.1,
    vjust = 1.5,
    size = 4) +
  scale_x_continuous(
    breaks = seq(600, 1800, by = 100),
    minor_breaks = seq(600, 1800, by = 50),
    guide = guide_axis(minor.ticks = TRUE)) +
  theme_classic()

bayesmod_plot_np

ggsave(
  "outputs/figures/bayesian_figures/bayesmod_elevation_curve_no_prior.png",
  plot = bayesmod_plot_np,
  width = 7,
  height = 4,
  dpi = 300)


####

## Visualize posterior distribution of range shift

density_shift_np <- density(
  range_shift_draws_np$shift_m)

density_data_np <- tibble(
  x = density_shift_np$x,
  y = density_shift_np$y)

ci_density_np <- density_data_np %>%
  filter(
    x >= lower_95_shift_np,
    x <= upper_95_shift_np)

median_density_np <- approx(
  density_data_np$x,
  density_data_np$y,
  xout = median_shift_np)$y

bayes_range_shift_np <- ggplot(
  density_data_np,
  aes(x = x, y = y)) +
  geom_ribbon(
    data = ci_density_np,
    aes(
      ymin = 0,
      ymax = y),
    alpha = 0.2) +
  geom_line() +
  annotate(
    "segment",
    x = median_shift_np,
    xend = median_shift_np,
    y = 0,
    yend = median_density_np,
    linetype = "dashed") +
  annotate(
    "point",
    x = median_shift_np,
    y = median_density_np,
    size = 3) +
  annotate(
    "text",
    x = median_shift_np,
    y = median_density_np,
    label = paste0(
      round(median_shift_np, 0), " m"),
    vjust = -1,
    size = 4) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed") +
  xlab("Range shift (m)") +
  ylab("Posterior density") +
  scale_x_continuous(
    limits = c(-50, 250),
    breaks = seq(-50, 250, by = 50),
    minor_breaks = seq(-50, 250, by = 25),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    guide = guide_axis(minor.ticks = TRUE)) +
  theme_classic()

bayes_range_shift_np


####

## Combine elevation curve and range shift distribution

bayes_combined_np <- bayesmod_plot_np / bayes_range_shift_np +
  plot_annotation(tag_levels = "A")

bayes_combined_np

ggsave(
  "outputs/figures/bayesian_figures/range_shift_combined_no_prior.png",
  plot = bayes_combined_np,
  width = 7,
  height = 8,
  dpi = 300)
