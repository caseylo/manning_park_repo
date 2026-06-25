#### Analyzing relocation error

## Load libraries

library(tidyverse)
library(vegan)
library(emmeans)
library(ggplot2)
library(ggbreak)
library(RColorBrewer)
library(dplyr)
library(pairwiseAdonis)
set.seed(123)

## Read data

comm_matrix <- read_csv("manning_park_data/relocation_data_wide.csv")
data_long <- read_csv("manning_park_data/relocation_data_long.csv")

## Create metadata frame
str(data_long)

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

## Start analysis!!

## Compute Bray_curtis distance matrix
dist_mat_bray <- vegdist(comm_matrix, method = "bray")

## Convert distance matrix to long format
dist_df_bray <- as.data.frame(as.matrix(dist_mat_bray)) %>%
  rownames_to_column("Plot1") %>%
  pivot_longer(-Plot1, names_to = "Plot2", values_to = "Distance") %>%
  filter(Plot1 != Plot2)

## Attach metadata for both plots
dist_df_bray <- dist_df_bray %>%
  left_join(metadata, by = c("Plot1" = "PlotNumber")) %>%
  rename(Site1 = Site, Treatment1 = Treatment) %>%
  left_join(metadata, by = c("Plot2" = "PlotNumber")) %>%
  rename(Site2 = Site, Treatment2 = Treatment)

## Define comparison type and keep within-site unique pairs
dist_df_bray <- dist_df_bray %>%
  mutate(
    Comparison = case_when(
      (Treatment1 == "2019" & Treatment2 == "2025") |
        (Treatment1 == "2025" & Treatment2 == "2019") ~ "2019_vs_2025",
      
      (Treatment1 == "2019" & Treatment2 == "Error") |
        (Treatment1 == "Error" & Treatment2 == "2019") ~ "2019_vs_Error",
      
      (Treatment1 == "2025" & Treatment2 == "Error") |
        (Treatment1 == "Error" & Treatment2 == "2025") ~ "2025_vs_Error")) %>%
  filter(
    Site1 == Site2,
    Plot1 < Plot2)

## Format to wide shape for 2019 comparisons
dist_wide_relocation_bray <- dist_df_bray %>%
  filter(Comparison %in% c("2019_vs_2025", "2019_vs_Error")) %>%
  select(Site1, Comparison, Distance) %>%
  pivot_wider(names_from = Comparison, values_from = Distance) %>%
  mutate(delta = `2019_vs_Error` - `2019_vs_2025`)

## Check normality of paired differences
hist(dist_wide_relocation_bray$delta)
qqnorm(dist_wide_relocation_bray$delta)
qqline(dist_wide_relocation_bray$delta)
shapiro.test(dist_wide_relocation_bray$delta)

## Paired t-test

t.test(dist_wide_relocation_bray$`2019_vs_Error`,
       dist_wide_relocation_bray$`2019_vs_2025`,
       paired = TRUE)

wilcox.test(dist_wide_relocation_bray$`2019_vs_Error`,
            dist_wide_relocation_bray$`2019_vs_2025`,
            paired = TRUE)

## Visualize

ggplot(dist_wide_relocation_bray,
       aes(x = `2019_vs_2025`, y = `2019_vs_Error`, colour = Site1)) +
  geom_point(size = 3) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed") +
  coord_cartesian(xlim = c(0.3, NA), ylim = c(0.3, NA)) +
  expand_limits(x = 0, y = 0) +
  labs(x = "2019 vs. 2025 plot dissimilarity", y = "2019 vs. error plot dissimilarity") +
  theme_classic()

####

## Compute Jaccard distance matrix

dist_mat_jaccard <- vegdist(comm_matrix, method = "jaccard", binary = TRUE)

## Convert distance matrix to long format
dist_df_jaccard <- as.data.frame(as.matrix(dist_mat_jaccard)) %>%
  rownames_to_column("Plot1") %>%
  pivot_longer(-Plot1, names_to = "Plot2", values_to = "Distance") %>%
  filter(Plot1 != Plot2)

## Attach metadata for both plots
dist_df_jaccard <- dist_df_jaccard %>%
  left_join(metadata, by = c("Plot1" = "PlotNumber")) %>%
  rename(Site1 = Site, Treatment1 = Treatment) %>%
  left_join(metadata, by = c("Plot2" = "PlotNumber")) %>%
  rename(Site2 = Site, Treatment2 = Treatment)

## Define comparison type and keep within-site unique pairs
dist_df_jaccard <- dist_df_jaccard %>%
  mutate(
    Comparison = case_when(
      (Treatment1 == "2019" & Treatment2 == "2025") |
        (Treatment1 == "2025" & Treatment2 == "2019") ~ "2019_vs_2025",
      
      (Treatment1 == "2019" & Treatment2 == "Error") |
        (Treatment1 == "Error" & Treatment2 == "2019") ~ "2019_vs_Error",
      
      (Treatment1 == "2025" & Treatment2 == "Error") |
        (Treatment1 == "Error" & Treatment2 == "2025") ~ "2025_vs_Error")) %>%
  filter(
    Site1 == Site2,
    Plot1 < Plot2)

## Format to wide shape for 2019 comparisons
dist_wide_relocation_jaccard <- dist_df_jaccard %>%
  filter(Comparison %in% c("2019_vs_2025", "2019_vs_Error")) %>%
  select(Site1, Comparison, Distance) %>%
  pivot_wider(names_from = Comparison, values_from = Distance) %>%
  mutate(delta = `2019_vs_Error` - `2019_vs_2025`)

## Check normality of paired differences
hist(dist_wide_relocation_jaccard$delta)
qqnorm(dist_wide_relocation_jaccard$delta)
qqline(dist_wide_relocation_jaccard$delta)
shapiro.test(dist_wide_relocation_jaccard$delta)

## Paired t-test

t.test(dist_wide_relocation_jaccard$`2019_vs_Error`,
       dist_wide_relocation_jaccard$`2019_vs_2025`,
       paired = TRUE)

wilcox.test(dist_wide_relocation_jaccard$`2019_vs_Error`,
            dist_wide_relocation_jaccard$`2019_vs_2025`,
            paired = TRUE)

## Visualize

ggplot(dist_wide_relocation_jaccard,
       aes(x = `2019_vs_2025`, y = `2019_vs_Error`, colour = Site1)) +
  geom_point(size = 3) +
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed") +
  coord_cartesian(xlim = c(0.3, NA), ylim = c(0.3, NA)) +
  expand_limits(x = 0, y = 0) +
  labs(x = "2019 vs. 2025 plot dissimilarity", y = "2019 vs. error plot dissimilarity") +
  theme_classic()

####

## PERMANOVA for Bray-Curtis Dissimilarity. Asking if the centroids of groups differ in multivariate space.
## Are communities in different treatments(2019, 2025, error) compositionally different from each other?

perm_mod_bray <- adonis2(
  comm_matrix ~ Treatment,
  data = metadata,
  method = "bray",
  strata = metadata$Site,
  permutations = 10000)

perm_mod_bray

perm_mod_bray_func <- adonis2(
  comm_matrix ~ Treatment,
  data = metadata,
  method = "bray",
  strata = metadata$Site,
  permutations = 10000)

perm_mod_bray_func

## Make distance matrix for dispersion test
dist_mat_bray <- vegdist(comm_matrix, method = "bray")
summary(dist_mat_bray)

## Test homogeneity of multivariate dispersion
disp_bray <- betadisper(dist_mat_bray, metadata$Treatment)
anova(disp_bray)
permutest(disp_bray, permutations = 10000)

##Subset for only 2025 and error
keep <- metadata$Treatment %in% c("2025", "Error")
meta_2025_error <- metadata[keep, ]
comm_2025_error <- comm_matrix[meta_2025_error$PlotNumber, ]

adonis2(comm_2025_error ~ Treatment,
        data = meta_2025_error,
        method = "bray",
        strata = meta_2025_error$Site,
        permutations = 10000)

## 2019 and 2025 only
meta_2019_2025 <- metadata[metadata$Treatment %in% c("2019", "2025"), ]
comm_2019_2025 <- comm_matrix[meta_2019_2025$PlotNumber, ]

adonis2(comm_2019_2025 ~ Treatment,
        data = meta_2019_2025,
        method = "bray",
        strata = meta_2019_2025$Site,
        permutations = 10000)

## 2019 and error
meta_2019_error <- metadata[metadata$Treatment %in% c("2019", "Error"), ]
comm_2019_error <- comm_matrix[meta_2019_error$PlotNumber, ]

adonis2(comm_2019_error ~ Treatment,
        data = meta_2019_error,
        method = "bray",
        strata = meta_2019_error$Site,
        permutations = 10000)

####


## Pairwise analysis of PERMANOVA

#NOT AVAILABLE???

pairwise.adonis2(
  comm_matrix ~ Treatment,
  data = as.data.frame(metadata),
  method = "bray",
  strata = 'Site',
  nperm = 10000)

####

## Plot an NMDS for Bray-Curtis Dissimilarity

nmds_bray2 <- metaMDS(comm_matrix, distance = "bray", k = 2)

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

plot(nmds_bray2)
points(nmds_bray2, display = "species")
text(nmds_bray2, display = "species")

####

nmds_bray <- ggplot(scores_df_bray, aes(x = NMDS1, y = NMDS2)) +
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
  )
nmds_bray

ggsave("figures/nmds_bray.png", plot = nmds_bray, width = 7, height = 5, dpi = 300)

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

pairwise.adonis2(
  comm_matrix ~ Treatment,
  data = as.data.frame(metadata),
  method = "jaccard",
  binary = TRUE,
  strata = 'Site',
  nperm = 10000)

## Make distance matrix for dispersion test
dist_mat_jaccard <- vegdist(comm_matrix, method = "jaccard", binary = TRUE)

## Test homogeneity of multivariate dispersion
disp_jaccard <- betadisper(dist_mat_jaccard, metadata$Treatment)
anova(disp_jaccard)
permutest(disp_jaccard, permutations = 10000)

## 2025 and error only
adonis2(comm_2025_error ~ Treatment,
        data = meta_2025_error,
        method = "jaccard",
        binary = TRUE,
        strata = meta_2025_error$Site,
        permutations = 10000)

## 2019 and 2025 only
adonis2(comm_2019_2025 ~ Treatment,
        data = meta_2019_2025,
        method = "jaccard",
        binary = TRUE,
        strata = meta_2019_2025$Site,
        permutations = 10000)

## 2019 and error only
adonis2(comm_2019_error ~ Treatment,
        data = meta_2019_error,
        method = "jaccard",
        binary = TRUE,
        strata = meta_2019_error$Site,
        permutations = 10000)

####

## Pairwise analysis of PERMANOVA

#NOT AVAILABLE???

#pairwise.adonis2(
#  comm_matrix ~ Treatment,
#  data = metadata,
#  method = "jaccard",
#  binary = TRUE,
#  strata = metadata$Site)

####

## Plot an NMDS for Bray-Curtis Dissimilarity

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
nmds_jaccard <- ggplot(scores_df_jaccard, aes(x = NMDS1, y = NMDS2)) +
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
  )
nmds_jaccard

ggsave("figures/nmds_jaccard.png", plot = nmds_jaccard, width = 7, height = 5, dpi = 300)

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