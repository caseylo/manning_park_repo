#### Analyzing relocation error with PERMANOVA

## Load libraries

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggbreak)
library(RColorBrewer)
library(dplyr)
library(pairwiseAdonis)
library(ggnewscale)
library(ggrepel)

set.seed(123)

## Read data for PERMANOVA

comm_matrix <- read_csv("data/processed/relocation_data_wide.csv")
data_long <- read_csv("data/processed/relocation_data_long.csv")

species_names <- read_csv("data/cleaned/master_species_list.csv")

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

## Saving main permanovas

#### Save relocation error PERMANOVA models

saveRDS(perm_mod_bray, file = "outputs/models/relocation_error_permanova/perm_mod_bray.rds")
saveRDS(perm_mod_bray_pairwise, file = "outputs/models/relocation_error_permanova/perm_mod_bray_pairwise.rds")

saveRDS(perm_mod_jaccard, file = "outputs/models/relocation_error_permanova/perm_mod_jaccard.rds")
saveRDS(perm_mod_jaccard_pairwise, file = "outputs/models/relocation_error_permanova/perm_mod_jaccard_pairwise.rds")

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

ggsave("outputs/figures/nmds_bray.png", plot = nmds_bray_plot, width = 7, height = 5, dpi = 300)

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

ggsave("outputs/figures/nmds_bray_centeroid.png", plot = nmds_bray_centeroid, width = 7, height = 5, dpi = 300)

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

ggsave("outputs/figures/nmds_jaccard.png", plot = nmds_jaccard_plot, width = 7, height = 5, dpi = 300)

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

ggsave("outputs/figures/nmds_jaccard_centeroid.png", plot = nmds_jaccard_centeroid, width = 7, height = 5, dpi = 300)

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

ggsave("outputs/figures/nmds_bray_with_species.png", plot = nmds_bray_plot_max, width = 16, height = 8, dpi = 300)

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

ggsave("outputs/figures/nmds_bray_plot_functional.png", plot = nmds_bray_plot_functional, width = 11, height = 7, dpi = 300)

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

ggsave("outputs/figures/nmds_jaccard_with_species.png", plot = nmds_jaccard_plot_max, width = 16, height = 8, dpi = 300)

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

ggsave("outputs/figures/nmds_jaccard_plot_functional.png", plot = nmds_jaccard_plot_functional, width = 11, height = 7, dpi = 300)
