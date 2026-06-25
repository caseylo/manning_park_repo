#### Analyzing model 8.4 from model selection

## Load libraries

library(tidyverse)
library(dplyr)
library(lme4)
library(tidyr)
library(ggplot2)
library(ggeffects)
library(DHARMa)
library(spdep)
library(sf)
library(dplyr)
library(ggh4x)
library(emmeans)
library(lmerTest)

## Read in data

model_data <- read_csv("manning_park_data/model_data_filt.csv")

## Sensitivity analysis data

model_data_conf <- read_csv("manning_park_data/sensitivity_analysis_data/model_data_conf.csv")
model_data_medhigh_conf <- read_csv("manning_park_data/sensitivity_analysis_data/model_data_medhigh_conf.csv")

# No data in this dataset
#model_data_high_conf <- read_csv("manning_park_data/sensitivity_analysis_data/model_data_high_conf.csv")

####


## Best model (all data)

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
    Functional_group = factor(Functional_group))

## Check structure

str(model_data)

## Check number of plots

length(unique(model_data$PlotNumber))
table(unique(model_data$PlotNumber))
#[1] 34

####

## Running the best model on full data

bestmod <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(bestmod)

## Evaluate model fit 
sim_res <- simulateResiduals(bestmod)

png("figures/model_outliers.png", width = 15, height = 6, units = "in", res = 300)
plot(sim_res)
dev.off()

testDispersion(sim_res)
testUniformity(sim_res)
testOutliers(sim_res)

plotResiduals(sim_res, model_data$Elevation_sc)
plotResiduals(sim_res, model_data$Lat_sc)


overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  ratio <- Pearson.chisq / rdf
  pval <- pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE)
  c(chisq = Pearson.chisq, ratio = ratio, p = pval)
}

overdisp_fun(bestmod)

## Moran's I test for spatial autocorrelation in GLMM residuals

#Extract Pearson residuals
model_data$resid_pearson <- residuals(bestmod, type = "pearson")

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
k <- 8
nb <- knn2nb(knearneigh(coords, k = k))

## Spatial weights
listw <- nb2listw(nb, style = "W", zero.policy = TRUE)

## Moran's I test
moran_res <- moran.test(plot_resids$resid_pearson, listw, zero.policy = TRUE)
moran_res

####

## Examine with Functional group interaction

bestmod_fg <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(bestmod_fg)
summary(bestmod_fg)

## Pairwise significance of linear elevation

emtrends(
  bestmod_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

emmeans(bestmod_fg, pairwise ~ Time | Functional_group)


####

## SENSITIVITY ANALYSIS: High, High/Med-High, Med-High, and Med confidence:

model_data_conf <- model_data_conf %>%
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

confmod <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data_conf,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(confmod)

## With Functional_group

confmod_fg <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc, 2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data_conf,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(confmod_fg)
emmeans(confmod_fg, pairwise ~ Time | Functional_group)
emtrends(
  confmod_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

####

## SENSITIVITY ANALYSIS: High, High/Med-High, Med-High confidence: 

model_data_medhigh_conf <- model_data_medhigh_conf %>%
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

medhighmod <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data_medhigh_conf,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(medhighmod)

## With Functional group 

medhighmod_fg <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | PlotNumber)
  + Lat_sc,
  data = model_data_medhigh_conf,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(medhighmod_fg)
emmeans(medhighmod_fg, pairwise ~ Time | Functional_group)
emtrends(
  medhighmod_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

####

## All summaries together

summary(bestmod)
summary(confmod)
summary(medhighmod)

summary(bestmod_fg)
summary(confmod_fg)
summary(medhighmod_fg)

emtrends(
  bestmod_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

emmeans(bestmod_fg, pairwise ~ Time | Functional_group)

emtrends(
  confmod_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

emmeans(confmod_fg, pairwise ~ Time | Functional_group)

emtrends(
  medhighmod_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

emmeans(medhighmod_fg, pairwise ~ Time | Functional_group)

####

## Visualizing models: Best model
pred <- ggpredict(bestmod, terms = c("Elevation_sc [all]", "Time"))

## Adjust to show re-scaled elevation
mean_elev <- mean(model_data$Elevation)
sd_elev   <- sd(model_data$Elevation)
pred$Elevation_m <- pred$x * sd_elev + mean_elev

## Calculate optima and shift of predicted probability of presence
optima <- pred %>%
  group_by(group) %>%
  slice_max(predicted)

shift <- optima %>%
  select(group, Elevation_m) %>%
  pivot_wider(names_from = group, values_from = Elevation_m) %>%
  mutate(shift_m = present - historical)

shift

## 1. Difference at each peak (same elevation comparison)
peak_probs <- optima %>%
  select(group, Elevation_m, predicted) %>%
  pivot_wider(names_from = group, values_from = predicted)
peak_probs

peak_probs$present[2] - peak_probs$historical[1]
#[1] 0.1110458

## 2. Difference at 1700 m
target_elev <- 1700

prob_1700 <- pred %>%
  group_by(group) %>%
  slice_min(abs(Elevation_m - target_elev), n = 1) %>%  # closest value
  ungroup() %>%
  select(group, Elevation_m, predicted) %>%
  pivot_wider(names_from = group, values_from = predicted) %>%
  mutate(
    location = "1700 m",
    diff_present_minus_historical = present - historical
  )

prob_1700
#diff =  0.232

## Pretty plot
bestmod_plot <- ggplot(pred, aes(x = Elevation_m, y = predicted, colour = group)) +
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
  annotate("text", x = Inf, y = Inf, label = "+63.2 m", hjust = 1.5, vjust = 4, size = 5) +
  theme_classic()

bestmod_plot

ggsave("figures/bestmod_elevation_curve.png",plot = bestmod_plot, width = 7, height = 4, dpi = 300)

####

## Bestmod with FG
pred_fg <- ggpredict(bestmod_fg, terms = c("Elevation_sc [all]", "Time", "Functional_group"))

## Adjust to show re-scaled elevation
mean_elev_fg <- mean(model_data$Elevation)
sd_elev_fg   <- sd(model_data$Elevation)
pred_fg$Elevation_m <- pred_fg$x * sd_elev_fg + mean_elev_fg

## Calculate optima of predicted probability of presence
optima_fg <- pred_fg %>%
  group_by(facet, group) %>%
  slice_max(predicted)

shift_fg <- optima_fg %>%
  select(facet, group, Elevation_m) %>%
  tidyr::pivot_wider(names_from = group, values_from = Elevation_m) %>%
  mutate(shift_m = present - historical,
         label = paste0("+ ", round(shift_m, 1), " m"))

shift_fg

shift_fg <- shift_fg %>%
  mutate(
    x = Inf,
    y = Inf)

## Difference in probability of occurrence at peak probabilities of occurence
optima_fg %>%
  group_by(facet) %>%
  summarise(diff_peak = predicted[group == "present"] - predicted[group == "historical"])

#facet diff_peak
#<fct>     <dbl>
#1 herb     0.111 
#2 shrub    0.146 
#3 tree     0.0392

## Differece at 1700m
pred_fg %>%
  group_by(facet, group) %>%
  slice_min(abs(Elevation_m - 1700), n = 1, with_ties = FALSE) %>%
  summarise(predicted = predicted, .groups = "drop") %>%
  group_by(facet) %>%
  summarise(diff_1700 = predicted[group == "present"] - predicted[group == "historical"])

#A tibble: 3 × 2
#facet diff_1700
#<fct>     <dbl>
# 1 herb     0.411 
#2 shrub    0.189 
#3 tree    -0.0124

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
    fill = "Time") +
  geom_vline(
    data = optima_fg,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
    data = optima_fg,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  facet_wrap(~facet) +
  geom_text(
    data = shift_fg,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.1,
    vjust = 1.5,
    size = 4.5
  ) +
  theme_classic()

bestmod_fg_plot

ggsave("figures/bestmod_fg_elevation_curve.png", plot = bestmod_fg_plot, width = 11, height = 4.5, dpi = 300)

####

## Visualizing models: confmod
pred_conf <- ggpredict(confmod, terms = c("Elevation_sc [all]", "Time"))

## Adjust to show re-scaled elevation
mean_elev_conf <- mean(model_data_conf$Elevation)
sd_elev_conf   <- sd(model_data_conf$Elevation)
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
    y = "Predicted probability of occurrence (confmod)",
    colour = "Time",
    fill = "Time") +
  geom_vline(
    data = optima_conf,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
    data = optima_conf,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  theme_classic()

conf_plot

ggsave("figures/confmod_elevation_curve.png",plot = conf_plot, width = 7, height = 4, dpi = 300)


####

## Confmod with FG
pred_conf_fg <- ggpredict(bestmod_fg, terms = c("Elevation_sc [all]", "Time", "Functional_group"))

## Adjust to show re-scaled elevation
mean_elev_conf_fg <- mean(model_data_conf$Elevation)
sd_elev_conf_fg   <- sd(model_data_conf$Elevation)
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
    y = "Predicted probability of occurrence (confmod)",
    colour = "Time",
    fill = "Time") +
  geom_vline(
    data = optima_conf_fg,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
    data = optima_conf_fg,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  facet_wrap(~facet) +
  theme_classic()

conf_fg_plot

ggsave("figures/confmod_fg_elevation_curve.png", plot = conf_fg_plot, width = 11, height = 4.5, dpi = 300)


####

## Visualizing models: medhighconf
pred_medhigh <- ggpredict(medhighmod, terms = c("Elevation_sc [all]", "Time"))

## Adjust to show re-scaled elevation
mean_elev_medhigh <- mean(model_data_medhigh_conf$Elevation)
sd_elev_medhigh  <- sd(model_data_medhigh_conf$Elevation)
pred_medhigh$Elevation_m <- pred_medhigh$x * sd_elev_medhigh+ mean_elev_medhigh

## Calculate optima of predicted probability of presence
optima_medhigh <- pred_medhigh %>%
  group_by(group) %>%
  slice_max(predicted)

## Pretty plot
medhigh_plot <- ggplot(pred_medhigh, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA) +
  scale_x_continuous(
    breaks = seq(
      floor(min(pred_medhigh$Elevation_m) / 200) * 200,
      ceiling(max(pred_medhigh$Elevation_m) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(pred_medhigh$Elevation_m) / 100) * 100,
      ceiling(max(pred_medhigh$Elevation_m) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(floor(min(pred_medhigh$predicted) / 0.1) * 0.1, ceiling(max(pred_medhigh$predicted) / 0.1 * 0.1), by = 0.1)) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence (medhighmod)",
    colour = "Time",
    fill = "Time") +
  geom_vline(
    data = optima_medhigh,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
    data = optima_medhigh,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  theme_classic()

medhigh_plot

ggsave("figures/medhighmod_elevation_curve.png",plot = medhigh_plot, width = 7, height = 4, dpi = 300)

####

## medhighmod with FG
pred_medhigh_fg <- ggpredict(bestmod_fg, terms = c("Elevation_sc [all]", "Time", "Functional_group"))

## Adjust to show re-scaled elevation
mean_elev_medhigh_fg <- mean(model_data_medhigh_conf$Elevation)
sd_elev_medhigh_fg   <- sd(model_data_medhigh_conf$Elevation)
pred_medhigh_fg$Elevation_m <- pred_medhigh_fg$x * sd_elev_medhigh_fg + mean_elev_medhigh_fg

## Calculate optima of predicted probability of presence
optima_medhigh_fg <- pred_medhigh_fg %>%
  group_by(facet, group) %>%
  slice_max(predicted)

## Pretty plot
medhigh_fg_plot <- ggplot(pred_medhigh_fg, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA) +
  scale_x_continuous(
    breaks = seq(
      floor(min(pred_medhigh_fg$Elevation_m) / 200) * 200,
      ceiling(max(pred_medhigh_fg$Elevation_m) / 200) * 200,
      by = 200),
    minor_breaks = seq(
      floor(min(pred_medhigh_fg$Elevation_m) / 100) * 100,
      ceiling(max(pred_medhigh_fg$Elevation_m) / 100) * 100,
      by = 100),
    guide = guide_axis(minor.ticks = TRUE)) +
  scale_y_continuous(breaks = seq(floor(min(pred_medhigh_fg$predicted) / 0.1) * 0.1, ceiling(max(pred_medhigh_fg$predicted) / 0.1 * 0.1), by = 0.1)) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence (medhighmod)",
    colour = "Time",
    fill = "Time") +
  geom_vline(
    data = optima_medhigh_fg,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
    data = optima_medhigh_fg,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE) +
  facet_wrap(~facet) +
  theme_classic()

medhigh_fg_plot

ggsave("figures/medhighmod_fg_elevation_curve.png", plot = medhigh_fg_plot, width = 11, height = 4.5, dpi = 300)



## Playing around
library(plotly)
library(dplyr)

plot_ly(
  model_data_conf,
  x = ~Lon_sc,
  y = ~Lat_sc,
  z = ~Elevation,
  color = ~Time,
  colors = c("red","blue"),
  type = "scatter3d",
  mode = "markers"
) %>%
  layout(
    scene = list(
      xaxis = list(title = "Longitude"),
      yaxis = list(title = "Latitude"),
      zaxis = list(title = "Elevation (m)")
    )
  )

pred3d <- ggpredict(
  confmod_fg,
  terms = c("Elevation_sc [all]", "Lat_sc [all]", "Time")
)

####

## Rescaling for visualization

## Pretty plot

bestmod_plot_scaled <- ggplot(pred, aes(x = Elevation_m, y = predicted, colour = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(
    aes(ymin = conf.low, ymax = conf.high, fill = group),
    alpha = 0.2,
    colour = NA
  ) +
  scale_x_continuous(
    breaks = seq(800, 1800, by = 200),
    minor_breaks = seq(700, 2100, by = 100),
    guide = guide_axis(minor.ticks = TRUE)
  ) +
  scale_y_continuous(
    breaks = seq(0.1, 0.8, by = 0.1)
  ) +
  coord_cartesian(
    xlim = c(680, 2100),
    ylim = c(0.0, 0.8)
  ) +
  labs(
    x = "Elevation (m)",
    y = "Predicted probability of occurrence",
    colour = "Time",
    fill = "Time",
    title = "Full model (1976-2000) - 2025"
  ) +
  geom_vline(
    data = optima,
    aes(xintercept = Elevation_m, colour = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(
    data = optima,
    aes(x = Elevation_m, y = predicted, colour = group),
    size = 2.5,
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = 2050,   # place explicitly instead of Inf for consistency
    y = 0.75,
    label = "+63.2 m",
    hjust = 1,
    size = 5
  ) +
  theme_classic()

bestmod_plot_scaled

ggsave("figures/bestmod_elevation_curve_scaled.png",plot = bestmod_plot_scaled, width = 7, height = 4, dpi = 300)

## Join 2019 to full model
library(patchwork)
library(png)
library(grid)

## Load PNG
img <- readPNG("figures/bestmod_elevation_curve_scaled.png")
img2 <- readPNG("figures/no_error_plot_scaled.png")
img_grob <- rasterGrob(img, interpolate = TRUE)
img_grob2 <- rasterGrob(img2, interpolate = TRUE)

## Combine
combined_plot <- wrap_elements(img_grob2) / wrap_elements(img_grob)

combined_plot

ggsave("figures/combined_scaled_plot.png",plot = combined_plot, width = 7, height = 8, dpi = 300)
