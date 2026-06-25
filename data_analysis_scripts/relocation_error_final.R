#### Analyzing relocation error

## Load libraries

library(tidyverse)
library(vegan)
library(lme4)
library(ggplot2)
library(ggeffects)
library(ggbreak)
library(RColorBrewer)
library(dplyr)
library(pairwiseAdonis)
library(ggnewscale)
library(ggrepel)
library(emmeans)
set.seed(123)

## Read data for PERMANOVA

comm_matrix <- read_csv("manning_park_data/relocation_data_wide.csv")
data_long <- read_csv("manning_park_data/relocation_data_long.csv")

species_names <- read_csv("manning_park_data/cleaning_data/master_species_list.csv")

## Create metadata frame
str(data_long)

data_long <- data_long  %>%
  select(-1)

metadata <- data_long %>%
  select(PlotNumber, Site, ProjectID, Treatment, Year,
         SiteSurveyor, RelocationConfidence,
         Latitude, Longitude, Elevation) %>%
  distinct() %>%
  mutate(
    Treatment = as.factor(Treatment),
    PlotNumber = as.factor(PlotNumber),
    Site = as.factor(Site))

## Structure matrix data
comm_matrix <- comm_matrix %>%
  as.data.frame() %>%
  select(-1) 

rownames(comm_matrix) <- comm_matrix$PlotNumber
comm_matrix$PlotNumber <- NULL

####

## Structure long data to run model

## 2025: one record per plot-species
data_2025_site_species <- data_long %>%
  distinct(Site, Treatment, Species, Functional_group)

## Count occurrences as number of unique plots
species_filter <- data_2025_site_species %>%
  filter(Treatment %in% c("2019", "2025", "Error")) %>%
  group_by(Species) %>%
  summarise(
    n_sites_2019_2025 = n_distinct(Site[Treatment %in% c("2019", "2025")]),
    n_sites_2019_Error = n_distinct(Site[Treatment %in% c("2019", "Error")]),
    .groups = "drop"
  ) %>%
  filter(n_sites_2019_2025 > 1 | n_sites_2019_Error > 1)

species_keep <- species_filter$Species

data_filtered <- data_2025_site_species %>%
  filter(Species %in% species_keep)

species_2025_count <- data_filtered %>%
  count(Species, name = "occurences_2025") %>%
  mutate(Species = as.character(Species))

#species_2025_count <- data_2025_plot_species %>%
  #count(Species, name = "occurences_2025") %>%
  #mutate(Species = as.character(Species)) %>%
  #filter(occurences_2025 >= 4)

species_2025 <- species_2025_count %>%
  select(-occurences_2025)

plot_grid <- data_long %>%
  distinct(PlotNumber, ProjectID, Site, Treatment, Latitude, Longitude, Elevation, Year, SiteSurveyor,
           RelocationConfidence) %>%
  crossing(species_2025)

length(unique(plot_grid$Site)) #Check to make sure this is working

data_long_long <- plot_grid %>%
  left_join(
    data_long,
    by = c(
      "PlotNumber",
      "ProjectID",
      "Treatment",
      "Site",
      "Latitude",
      "Longitude",
      "Elevation",
      "Year",
      "Species",
      "SiteSurveyor",
      "RelocationConfidence")) %>%
  mutate(present = if_else(is.na(Cover), 0L, 1L), Cover   = replace_na(Cover, 0)) %>%
  left_join(species_names, by = c("Species" = "BC_species_code")) %>%
  select(-Functional_group.x,
         -Latin_name) %>%
  rename(Functional_group = Functional_group.y)

min_year <- min(data_long_long$Year)

model_data <- data_long_long %>%
  mutate(YearsSinceStart = Year - min_year) %>%
  mutate(Time = case_when(
    Year == 2019 ~ "2019",
    Year == 2025 ~ "2025")) %>%
  mutate(Time = factor(Time, levels = c("2019", "2025")))

## Making an occurrences excel to look at if there was surveyor bias in identifying shrubs
species_occurrence_by_time <- data_2025_site_species %>%
  count(Species, Functional_group, Treatment, name = "occurrences") %>%
  tidyr::pivot_wider(
    names_from = Treatment,
    values_from = occurrences,
    values_fill = 0)

#write.csv(species_occurrence_by_time, "manning_park_data/relocation_error_occurrences.csv", row.names = FALSE) 

## Prepare data for analysis

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

## data from 2019-2025
model_data_no_error <- model_data %>%
  filter(Treatment != "Error")

mod_no_error <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | Site)
  + Lat_sc,
  data = model_data_no_error,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod_no_error)

mod_no_error_fg <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | Site)
  + Lat_sc,
  data = model_data_no_error,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod_no_error_fg)

## data from 2019-Error
model_data_error <- model_data %>%
  filter(Treatment != "2025")

mod_error <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) 
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | Site)
  + Lat_sc,
  data = model_data_error,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod_error)

mod_error_fg <- glmer(
  present ~ Time * poly(Elevation_sc, 2, raw = TRUE) * Functional_group
  + (poly(Elevation_sc,2, raw = TRUE) | Species)
  + (1 | Site)
  + Lat_sc,
  data = model_data_error,
  family = binomial,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))

summary(mod_no_error)
summary(mod_no_error_fg)
summary(mod_error)
summary(mod_error_fg)

emtrends(
  mod_no_error_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

emmeans(mod_no_error_fg, pairwise ~ Time | Functional_group)

emtrends(
  mod_error_fg,
  pairwise ~ Time | Functional_group,
  var = "Elevation_sc",
  at = list(Elevation_sc = 0))

emmeans(mod_error_fg, pairwise ~ Time | Functional_group)


## Start analysis!!

## PERMANOVA for Bray-Curtis Dissimilarity. Asking if the centroids of groups differ in multivariate space.
## Are communities in different treatments(2019, 2025, error) compositionally different from each other?

## Full model

perm_mod_bray <- adonis2(
  comm_matrix ~ Treatment,
  data = metadata,
  method = "bray",
  strata = metadata$Site,
  permutations = 10000)

perm_mod_bray

## Pairwise analysis of PERMANOVA

perm_mod_bray_pairwise <- pairwise.adonis2(
  comm_matrix ~ Treatment,
  data = as.data.frame(metadata),
  method = "bray",
  strata = 'Site',
  nperm = 10000)

perm_mod_bray_pairwise

####

## Make distance matrix for dispersion test
dist_mat_bray <- vegdist(comm_matrix, method = "bray")
summary(dist_mat_bray)

## Test homogeneity of multivariate dispersion
disp_bray <- betadisper(dist_mat_bray, metadata$Treatment)
anova(disp_bray)
permutest(disp_bray, permutations = 10000)

####

## Exploring herbs, shrubs, trees only

## herbs first
herb_species <- data_long %>%
  filter(Functional_group == "herb") %>%
  distinct(Species) %>%
  pull(Species)

comm_matrix_herb <- comm_matrix[, colnames(comm_matrix) %in% herb_species]

perm_mod_bray_pairwise_herb <- pairwise.adonis2(
  comm_matrix_herb ~ Treatment,
  data = as.data.frame(metadata),
  method = "bray",
  strata = 'Site',
  nperm = 10000
)
perm_mod_bray_pairwise_herb

## shrubs

shrub_species <- data_long %>%
  filter(Functional_group == "shrub") %>%
  distinct(Species) %>%
  pull(Species)

comm_matrix_shrub <- comm_matrix[, colnames(comm_matrix) %in% shrub_species]

perm_mod_bray_pairwise_shrub <- pairwise.adonis2(
  comm_matrix_shrub ~ Treatment,
  data = as.data.frame(metadata),
  method = "bray",
  strata = 'Site',
  nperm = 10000
)
perm_mod_bray_pairwise_shrub

## trees

tree_shrub_species <- data_long %>%
  filter(Functional_group %in% c("tree", "shrub")) %>%
  distinct(Species) %>%
  pull(Species)

comm_matrix_tree_shrub <- comm_matrix[, colnames(comm_matrix) %in% tree_shrub_species]

perm_mod_bray_pairwise_tree_shrub <- pairwise.adonis2(
  comm_matrix_tree_shrub ~ Treatment,
  data = as.data.frame(metadata),
  method = "bray",
  strata = 'Site',
  nperm = 10000
)
perm_mod_bray_pairwise_tree_shrub

####

## PERMANOVA for Jaccard Dissimilarity. Asking if the centroids of groups differ in multivariate space.
## Are communities in different treatments(2019, 2025, error) compositionally different from each other?

perm_mod_jaccard <- adonis2(
  comm_matrix ~ Treatment,
  data = metadata,
  method = "jaccard",
  binary = TRUE,
  strata = metadata$Site,
  permutations = 10000)

perm_mod_jaccard

perm_mod_jaccard_pairwise <- pairwise.adonis2(
  comm_matrix ~ Treatment,
  data = as.data.frame(metadata),
  method = "jaccard",
  binary = TRUE,
  strata = 'Site',
  nperm = 10000)

perm_mod_jaccard_pairwise

## Make distance matrix for dispersion test
dist_mat_jaccard <- vegdist(comm_matrix, method = "jaccard", binary = TRUE)

## Test homogeneity of multivariate dispersion
disp_jaccard <- betadisper(dist_mat_jaccard, metadata$Treatment)
anova(disp_jaccard)
permutest(disp_jaccard, permutations = 10000)

####

## Exploring herbs, shrubs, trees only

## herbs first

perm_mod_jaccard_pairwise_herb <- pairwise.adonis2(
  comm_matrix_herb ~ Treatment,
  data = as.data.frame(metadata),
  method = "jaccard",
  strata = 'Site',
  nperm = 10000
)
perm_mod_jaccard_pairwise_herb

## shrubs

perm_mod_jaccard_pairwise_shrub <- pairwise.adonis2(
  comm_matrix_shrub ~ Treatment,
  data = as.data.frame(metadata),
  method = "jaccard",
  strata = 'Site',
  nperm = 10000
)
perm_mod_jaccard_pairwise_shrub

## trees

perm_mod_jaccard_pairwise_tree_shrub <- pairwise.adonis2(
  comm_matrix_tree_shrub ~ Treatment,
  data = as.data.frame(metadata),
  method = "jaccard",
  strata = 'Site',
  nperm = 10000
)
perm_mod_jaccard_pairwise_tree_shrub

####

####

## Plot an NMDS for Bray-Curtis Dissimilarity

nmds_bray <- metaMDS(comm_matrix, distance = "bray", k = 2)

## Extract NMDS scores
scores_df_bray <- as.data.frame(scores(nmds_bray, display = "sites"))
scores_df_bray$PlotNumber <- rownames(scores_df_bray)

scores_df_bray <- scores_df_bray %>%
  left_join(metadata, by = "PlotNumber") %>%
  mutate(Treatment = factor(Treatment, levels = c("2025", "2019", "Error"))) %>%
  arrange(Site, Treatment)

## Calculate treatment centroids
centroids_bray <- scores_df_bray %>%
  group_by(Treatment) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2))

####

nmds_bray_plot <- ggplot(scores_df_bray, aes(x = NMDS1, y = NMDS2)) +
  geom_path(
    aes(group = Site),
    color = "grey60",
    alpha = 1,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(aes(color = Treatment, shape = Treatment), size = 3.2) +
  geom_text(
    data = scores_df_bray %>% filter(Treatment == "2019"),
    aes(label = Site),
    size = 3,
    vjust = -1.5
  ) +
  scale_color_manual(
    name = "Plot type",
    breaks = c("2019", "2025", "Error"),
    values = c(
      "2019" = "black",
      "2025" = "#1f78b4",
      "Error" = "#e31a1c"
    )
  ) +
  scale_shape_discrete(name = "Plot type", breaks = c("2019", "2025", "Error")) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black")
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Bray-Curtis",
    hjust = -0.1,
    vjust = 1.2,
    size = 5
  )
nmds_bray_plot

ggsave("figures/nmds_bray.png", plot = nmds_bray_plot, width = 7, height = 5, dpi = 300)

## Plot NMDS with treatment ellipses and centroids

nmds_bray_centeroid <- ggplot(scores_df_bray, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Treatment, shape = Treatment), size = 3) +
  stat_ellipse(
    aes(color = Treatment, group = Treatment),
    linewidth = 1,
    linetype = "solid",
    level = 0.95) +
  geom_point(
    data = centroids_bray,
    aes(x = NMDS1, y = NMDS2, color = Treatment),
    size = 4,
    shape = 4,
    stroke = 1,
    inherit.aes = FALSE) +
  scale_color_manual(
    name = "Plot type",
    breaks = c("2019", "2025", "Error"),
    values = c(
      "2019" = "black",
      "2025" = "#1f78b4",
      "Error" = "#e31a1c" )) +
  scale_shape_discrete(
    name = "Plot type",
    breaks = c("2019", "2025", "Error")) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black"))

nmds_bray_centeroid

ggsave("figures/nmds_bray_centeroid.png", plot = nmds_bray_centeroid, width = 7, height = 5, dpi = 300)

####

## Plot an NMDS for Jaccard Dissimilarity

nmds_jaccard <- metaMDS(comm_matrix, distance = "jaccard", binary = TRUE, k = 2)

## Extract NMDS scores
scores_df_jaccard <- as.data.frame(scores(nmds_jaccard, display = "sites"))
scores_df_jaccard$PlotNumber <- rownames(scores_df_jaccard)

scores_df_jaccard <- scores_df_jaccard %>%
  left_join(metadata, by = "PlotNumber") %>%
  mutate(Treatment = factor(Treatment, levels = c("2025", "2019", "Error"))) %>%
  arrange(Site, Treatment)

## Calculate treatment centroids
centroids_jaccard <- scores_df_jaccard %>%
  group_by(Treatment) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2))

####

#### NMDS with lines between each site treatment
nmds_jaccard_plot <- ggplot(scores_df_jaccard, aes(x = NMDS1, y = NMDS2)) +
  geom_path(
    aes(group = Site),
    color = "grey60",
    alpha = 1,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_point(aes(color = Treatment, shape = Treatment), size = 3.2) +
  geom_text(
    data = scores_df_jaccard %>% filter(Treatment == "2019"),
    aes(label = Site),
    size = 3,
    vjust = -1.5
  ) +
  scale_color_manual(
    name = "Plot type",
    breaks = c("2019", "2025", "Error"),
    values = c(
      "2019" = "black",
      "2025" = "#1f78b4",
      "Error" = "#e31a1c"
    )
  ) +
  scale_shape_discrete(name = "Plot type", breaks = c("2019", "2025", "Error")) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black")
  ) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Jaccard",
    hjust = -0.1,
    vjust = 1.2,
    size = 5
  )
nmds_jaccard_plot

ggsave("figures/nmds_jaccard.png", plot = nmds_jaccard_plot, width = 7, height = 5, dpi = 300)

## Plot NMDS with treatment ellipses and centroids

nmds_jaccard_centeroid <- ggplot(scores_df_jaccard, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Treatment, shape = Treatment), size = 3) +
  stat_ellipse(
    aes(color = Treatment, group = Treatment),
    linewidth = 1,
    linetype = "solid",
    level = 0.95) +
  geom_point(
    data = centroids_jaccard,
    aes(x = NMDS1, y = NMDS2, color = Treatment),
    size = 4,
    shape = 4,
    stroke = 1,
    inherit.aes = FALSE) +
  scale_color_manual(
    name = "Plot type",
    breaks = c("2019", "2025", "Error"),
    values = c(
      "2019" = "black",
      "2025" = "#1f78b4",
      "Error" = "#e31a1c" )) +
  scale_shape_discrete(
    name = "Plot type",
    breaks = c("2019", "2025", "Error")) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black"))

nmds_jaccard_centeroid

ggsave("figures/nmds_jaccard_centeroid.png", plot = nmds_jaccard_centeroid, width = 7, height = 5, dpi = 300)

####

#### Playing around

## Extract species scores from the NMDS
species_scores_bray <- as.data.frame(scores(nmds_bray, display = "species"))
species_scores_bray$BC_species_code <- rownames(species_scores_bray)

## Join functional group information
species_scores_bray <- species_scores_bray %>%
  left_join(species_names, by = "BC_species_code", "Functional_group")

## Make functional group centeroid
fg_centroids_bray <- species_scores_bray %>%
  group_by(Functional_group) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop")

## Pretty plot

nmds_bray_plot_max <- ggplot(scores_df_bray, aes(x = NMDS1, y = NMDS2)) +
  ## Species labels (no legend)
  geom_text_repel(
    data = species_scores_bray,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = BC_species_code,
      color = Functional_group),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey70",
    alpha = 0.8,
    show.legend = FALSE) + 
  ## Dummy points ONLY for legend
  geom_point(
    data = species_scores_bray,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Functional_group),
    shape = 16,
    size = 3,
    alpha = 0,            # invisible on plot
    show.legend = TRUE) +
  scale_color_manual(
    name = "Functional group",
    values = c(
      "herb" = "forestgreen",
      "shrub" = "darkorange",
      "tree" = "darkblue")) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 4, shape = 16))) +
  ggnewscale::new_scale_color() +
  ## ---- Plot structure SECOND (on top) ----
geom_path(
  aes(group = Site),
  color = "grey60",
  alpha = 1,
  linewidth = 0.8,
  show.legend = FALSE) +
  geom_point(
    aes(color = Treatment, shape = Treatment),
    size = 3.2) +
  geom_text(
    data = scores_df_bray %>% filter(Treatment == "2019"),
    aes(label = Site),
    size = 3,
    vjust = -1.5) +
  scale_color_manual(
    name = "Plot type",
    breaks = c("2019", "2025", "Error"),
    values = c(
      "2019" = "black",
      "2025" = "#1f78b4",
      "Error" = "#e31a1c" )) +
  scale_shape_discrete(
    name = "Plot type",
    breaks = c("2019", "2025", "Error")) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black")) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Bray-Curtis",
    hjust = -0.1,
    vjust = 1.2,
    size = 5
  )

nmds_bray_plot_max

ggsave("figures/nmds_bray_with_species.png", plot = nmds_bray_plot_max, width = 16, height = 8, dpi = 300)

nmds_bray_plot_functional <- ggplot(scores_df_bray, aes(x = NMDS1, y = NMDS2)) +
   ## Functional group ellipses (no legend, no labels)
  stat_ellipse(
    data = species_scores_bray,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Functional_group,
      group = Functional_group
    ),
    linewidth = 1.1,
    linetype = "solid",
    level = 0.95,
    alpha = 0.8,
    show.legend = FALSE) +
  ## Dummy points ONLY for legend
  geom_point(
    data = species_scores_bray,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Functional_group),
    shape = 16,
    size = 3,
    alpha = 0,
    show.legend = TRUE) +
  scale_color_manual(
    name = "Functional group",
    values = c(
      "herb" = "forestgreen",
      "shrub" = "darkorange",
      "tree" = "darkblue")) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 4, shape = 16)))+
  
  ggnewscale::new_scale_color() +
  ##Plot structure on top
  geom_path(
    aes(group = Site),
    color = "grey60",
    alpha = 1,
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
      aes(color = Treatment, shape = Treatment),
      size = 3.2) +
  geom_text(
    data = scores_df_bray %>% filter(Treatment == "2019"),
    aes(label = Site),
    size = 3,
    vjust = -1.5) +
  scale_color_manual(
    name = "Plot type",
    breaks = c("2019", "2025", "Error"),
    values = c(
      "2019" = "black",
      "2025" = "#1f78b4",
      "Error" = "#e31a1c")) +
  scale_shape_discrete(
    name = "Plot type",
    breaks = c("2019", "2025", "Error")) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black")) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Bray-Curtis",
    hjust = -0.1,
    vjust = 1.2,
    size = 5
  )

nmds_bray_plot_functional

ggsave("figures/nmds_bray_plot_functional.png", plot = nmds_bray_plot_functional, width = 11, height = 7, dpi = 300)


####

species_scores_jaccard <- as.data.frame(scores(nmds_jaccard, display = "species")) %>%
  rownames_to_column("BC_species_code") %>%
  left_join(species_names, by = "BC_species_code")

## Make functional group centeroid
fg_centroids_jaccard <- species_scores_jaccard %>%
  group_by(Functional_group) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop")

nmds_jaccard_plot_max <- ggplot(scores_df_jaccard, aes(x = NMDS1, y = NMDS2)) +
## Species labels (no legend)
  geom_text_repel(
    data = species_scores_jaccard,
    aes(
      x = NMDS1,
      y = NMDS2,
      label = BC_species_code,
      color = Functional_group),
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.color = "grey70",
    alpha = 0.8,
    show.legend = FALSE) + 
    ## Dummy points ONLY for legend
  geom_point(
    data = species_scores_jaccard,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Functional_group),
    shape = 16,
    size = 3,
    alpha = 0,            # invisible on plot
    show.legend = TRUE) +
  scale_color_manual(
    name = "Functional group",
    values = c(
      "herb" = "forestgreen",
      "shrub" = "darkorange",
      "tree" = "darkblue")) +
    guides(
      color = guide_legend(
        override.aes = list(alpha = 1, size = 4, shape = 16))) +
  ggnewscale::new_scale_color() +
    ## ---- Plot structure SECOND (on top) ----
  geom_path(
    aes(group = Site),
    color = "grey60",
    alpha = 1,
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
    aes(color = Treatment, shape = Treatment),
      size = 3.2) +
  geom_text(
    data = scores_df_jaccard %>% filter(Treatment == "2019"),
    aes(label = Site),
    size = 3,
    vjust = -1.5) +
  scale_color_manual(
    name = "Plot type",
    breaks = c("2019", "2025", "Error"),
    values = c(
      "2019" = "black",
      "2025" = "#1f78b4",
      "Error" = "#e31a1c" )) +
  scale_shape_discrete(
    name = "Plot type",
    breaks = c("2019", "2025", "Error")) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black"))

nmds_jaccard_plot_max

ggsave("figures/nmds_jaccard_with_species.png", plot = nmds_jaccard_plot_max, width = 16, height = 8, dpi = 300)

## With functional group centeroids and no species

nmds_jaccard_plot_functional <- ggplot(scores_df_jaccard, aes(x = NMDS1, y = NMDS2)) +
  ## Functional group ellipses (no legend, no labels)
  stat_ellipse(
    data = species_scores_jaccard,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Functional_group,
      group = Functional_group
    ),
    linewidth = 1.1,
    linetype = "solid",
    level = 0.95,
    alpha = 0.8,
    show.legend = FALSE) +
  ## Dummy points ONLY for legend
  geom_point(
    data = species_scores_jaccard,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Functional_group),
    shape = 16,
    size = 3,
    alpha = 0,
    show.legend = TRUE) +
  scale_color_manual(
    name = "Functional group",
    values = c(
      "herb" = "forestgreen",
      "shrub" = "darkorange",
      "tree" = "darkblue")) +
  guides(
    color = guide_legend(
      override.aes = list(alpha = 1, size = 4, shape = 16)))+
  
  ggnewscale::new_scale_color() +
  ##Plot structure on top
  geom_path(
    aes(group = Site),
    color = "grey60",
    alpha = 1,
    linewidth = 0.8,
    show.legend = FALSE) +
  geom_point(
    aes(color = Treatment, shape = Treatment),
    size = 3.2) +
  geom_text(
    data = scores_df_jaccard %>% filter(Treatment == "2019"),
    aes(label = Site),
    size = 3,
    vjust = -1.5) +
  scale_color_manual(
    name = "Plot type",
    breaks = c("2019", "2025", "Error"),
    values = c(
      "2019" = "black",
      "2025" = "#1f78b4",
      "Error" = "#e31a1c")) +
  scale_shape_discrete(
    name = "Plot type",
    breaks = c("2019", "2025", "Error")) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "white", colour = "black")) +
  annotate(
    "text",
    x = -Inf,
    y = Inf,
    label = "Jaccard",
    hjust = -0.1,
    vjust = 1.2,
    size = 5
  )

nmds_jaccard_plot_functional

ggsave("figures/nmds_jaccard_plot_functional.png", plot = nmds_jaccard_plot_functional, width = 11, height = 7, dpi = 300)


#### Visualizing elevation models

####

## Visualizing models: 2019-2025 correct
pred <- ggpredict(mod_no_error, terms = c("Elevation_sc [all]", "Time"))

## Adjust to show re-scaled elevation
mean_elev <- mean(model_data_no_error$Elevation)
sd_elev   <- sd(model_data_no_error$Elevation)
pred$Elevation_m <- pred$x * sd_elev + mean_elev

## Calculate optima and shift of predicted probability of presence
optima <- pred %>%
  group_by(group) %>%
  slice_max(predicted)

shift <- optima %>%
  select(group, Elevation_m) %>%
  pivot_wider(names_from = group, values_from = Elevation_m) %>%
  mutate(shift_m = `2025` - `2019`)

shift

## Pretty plot
no_error_plot <- ggplot(pred, aes(x = Elevation_m, y = predicted, colour = group)) +
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
    fill = "Time",
    title = "2019 - 2025 Correct Comparison") +
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
  theme_classic()

no_error_plot
ggsave("figures/no_error_elevation_plot.png",plot = no_error_plot, width = 7, height = 4, dpi = 300)

## 2019-2025 correct with FG
pred_fg <- ggpredict(mod_no_error_fg, terms = c("Elevation_sc [all]", "Time", "Functional_group"))

## Adjust to show re-scaled elevation
mean_elev_fg <- mean(model_data_no_error$Elevation)
sd_elev_fg   <- sd(model_data_no_error$Elevation)
pred_fg$Elevation_m <- pred_fg$x * sd_elev_fg + mean_elev_fg

## Calculate optima of predicted probability of presence
optima_fg <- pred_fg %>%
  group_by(facet, group) %>%
  slice_max(predicted)

shift_fg <- optima_fg %>%
  select(facet, group, Elevation_m) %>%
  tidyr::pivot_wider(names_from = group, values_from = Elevation_m) %>%
  mutate(shift_m = `2025` - `2019`,
         label = paste0("+ ", round(shift_m, 1), " m"))

shift_fg

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
    fill = "Time",
    title = "2019 - 2025 Correct Comparison") +
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
  theme_classic()

bestmod_fg_plot
ggsave("figures/no_error_elevation_fg_plot.png", plot = bestmod_fg_plot, width = 11, height = 4.5, dpi = 300)

####

## Visualizing models: 2019-2025 ERROR
pred_conf <- ggpredict(mod_error, terms = c("Elevation_sc [all]", "Time"))

## Adjust to show re-scaled elevation
mean_elev_conf <- mean(model_data_error$Elevation)
sd_elev_conf   <- sd(model_data_error$Elevation)
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
    fill = "Time",
    title = "2019 - 2025 Error Comparison") +
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
ggsave("figures/error_elevation_plot.png",plot = conf_plot, width = 7, height = 4, dpi = 300)

####

## Confmod with FG
pred_conf_fg <- ggpredict(mod_error_fg, terms = c("Elevation_sc [all]", "Time", "Functional_group"))

## Adjust to show re-scaled elevation
mean_elev_conf_fg <- mean(model_data_error$Elevation)
sd_elev_conf_fg   <- sd(model_data_error$Elevation)
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
    fill = "Time",
    title = "2019 - 2025 Error Comparison") +
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
ggsave("figures/error_elevation_fg_plot.png", plot = conf_fg_plot, width = 11, height = 4.5, dpi = 300)
