#### PERMANOVA analysis of all plots

#### Comparing 2025 resurvey plots with correctly relocated 2025 plots

## Load libraries

library(tidyverse)
library(vegan)
library(ggplot2)
library(ggrepel)

set.seed(1234)


## Read data

error_model_data <- read_csv("data/processed/error_model_data.csv")

model_data <- read_csv("data/processed/model_data_filt.csv")


#### Prepare data

## 2025 resurvey plots from model_data

data_2025 <- model_data %>%
  filter(Time == "present") %>%
  mutate(
    Treatment = "Historical plots",
    Dataset = "2025Resurvey"
  )


## Correctly relocated 2025 plots from error_model_data

data_2025_correct <- error_model_data %>%
  filter(Treatment == "2025") %>%
  mutate(
    Treatment = "Relocation error plots",
    Dataset = "2025Correct"
  )


## Identify species shared between datasets

shared_species <- intersect(
  unique(data_2025$Species),
  unique(data_2025_correct$Species)
)


## Combine datasets

combined_data <- bind_rows(
  data_2025,
  data_2025_correct
) %>%
  filter(Species %in% shared_species) %>%
  mutate(
    PlotID = paste(Dataset, PlotNumber, sep = "_"),
    Treatment = factor(
      Treatment,
      levels = c("Historical plots", "Relocation error plots")
    )
  )


## Check number of plots in each group

combined_data %>%
  distinct(PlotID, Treatment) %>%
  count(Treatment)


#### Create metadata

metadata <- combined_data %>%
  select(PlotID, PlotNumber, Treatment) %>%
  distinct() %>%
  mutate(
    Treatment = as.factor(Treatment),
    PlotNumber = as.factor(PlotNumber)
  )


#### Create community matrix

comm_matrix <- combined_data %>%
  group_by(PlotID, Species) %>%
  summarise(
    present = max(present),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Species,
    values_from = present,
    values_fill = 0
  ) %>%
  as.data.frame()

rownames(comm_matrix) <- comm_matrix$PlotID
comm_matrix$PlotID <- NULL


## Ensure metadata matches community matrix order

metadata <- metadata %>%
  slice(match(rownames(comm_matrix), PlotID))


## Check community matrix and metadata

dim(comm_matrix)

table(metadata$Treatment)

head(metadata)


## Confirm that the rows match

stopifnot(identical(rownames(comm_matrix), metadata$PlotID))


#### PERMANOVA for Bray-Curtis dissimilarity

## Test whether community composition differs between groups

perm_mod_bray <- adonis2(
  comm_matrix ~ Treatment,
  data = metadata,
  method = "bray",
  permutations = 10000
)

perm_mod_bray


#### Test homogeneity of multivariate dispersion for Bray-Curtis

## Calculate Bray-Curtis distance matrix

dist_mat_bray <- vegdist(
  comm_matrix,
  method = "bray"
)

summary(dist_mat_bray)


## Calculate distance to group centroid

disp_bray <- betadisper(
  dist_mat_bray,
  metadata$Treatment
)


## Test differences in dispersion

anova(disp_bray)

permutest(
  disp_bray,
  permutations = 10000
)


#### PERMANOVA for Jaccard dissimilarity

## Test whether community composition differs between groups

perm_mod_jaccard <- adonis2(
  comm_matrix ~ Treatment,
  data = metadata,
  method = "jaccard",
  binary = TRUE,
  permutations = 10000
)

perm_mod_jaccard


#### Test homogeneity of multivariate dispersion for Jaccard

## Calculate Jaccard distance matrix

dist_mat_jaccard <- vegdist(
  comm_matrix,
  method = "jaccard",
  binary = TRUE
)

summary(dist_mat_jaccard)


## Calculate distance to group centroid

disp_jaccard <- betadisper(
  dist_mat_jaccard,
  metadata$Treatment
)


## Test differences in dispersion

anova(disp_jaccard)

permutest(
  disp_jaccard,
  permutations = 10000
)


#### NMDS for Bray-Curtis dissimilarity

## Run NMDS

nmds_bray <- metaMDS(
  comm_matrix,
  distance = "bray",
  k = 2,
  autotransform = FALSE
)


## Check NMDS stress

nmds_bray$stress


## Extract NMDS scores

scores_df_bray <- as.data.frame(
  scores(nmds_bray, display = "sites")
)

scores_df_bray$PlotID <- rownames(scores_df_bray)

scores_df_bray <- scores_df_bray %>%
  left_join(metadata, by = "PlotID") %>%
  arrange(Treatment, PlotNumber) %>%
  mutate(
    PlotLabel = str_remove(as.character(PlotNumber), "_2025$")
  )


## Calculate group centroids

centroids_bray <- scores_df_bray %>%
  group_by(Treatment) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )


#### Plot Bray-Curtis NMDS

nmds_bray_plot <- ggplot(
  scores_df_bray,
  aes(x = NMDS1, y = NMDS2)
) +
  stat_ellipse(
    aes(color = Treatment, group = Treatment),
    linewidth = 1,
    linetype = "solid",
    level = 0.95
  ) +
  geom_point(
    aes(color = Treatment),
    size = 3.2
  ) +
  geom_text_repel(
    aes(
      label = PlotLabel,
      color = Treatment
    ),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_point(
    data = centroids_bray,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Treatment
    ),
    size = 4,
    shape = 4,
    stroke = 1.2,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c(
      "Historical plots",
      "Relocation error plots"
    ),
    values = c(
      "Historical plots" = "#1f78b4", 
      "Relocation error plots" = "#e31a1c"
    )
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    legend.title = element_blank()
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


#### NMDS for Jaccard dissimilarity

## Run NMDS

nmds_jaccard <- metaMDS(
  comm_matrix,
  distance = "jaccard",
  binary = TRUE,
  k = 2,
  autotransform = FALSE
)


## Check NMDS stress

nmds_jaccard$stress


## Extract NMDS scores

scores_df_jaccard <- as.data.frame(
  scores(nmds_jaccard, display = "sites")
)

scores_df_jaccard$PlotID <- rownames(scores_df_jaccard)

scores_df_jaccard <- scores_df_jaccard %>%
  left_join(metadata, by = "PlotID") %>%
  arrange(Treatment, PlotNumber) %>%
  mutate(
    PlotLabel = str_remove(as.character(PlotNumber), "_2025$")
  )


## Calculate group centroids

centroids_jaccard <- scores_df_jaccard %>%
  group_by(Treatment) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  )


#### Plot Jaccard NMDS

nmds_jaccard_plot <- ggplot(
  scores_df_jaccard,
  aes(x = NMDS1, y = NMDS2)
) +
 # stat_ellipse(
 #   aes(color = Treatment, group = Treatment),
 #   linewidth = 1,
 #   linetype = "solid",
 #   level = 0.95
 # ) +
  geom_point(
    aes(color = Treatment),
    size = 3.2
  ) +
  geom_text_repel(
    aes(
      label = PlotLabel,
      color = Treatment
    ),
    size = 3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  geom_point(
    data = centroids_jaccard,
    aes(
      x = NMDS1,
      y = NMDS2,
      color = Treatment
    ),
    size = 4,
    shape = 4,
    stroke = 1.2,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name = NULL,
    breaks = c(
      "Historical plots",
      "Relocation error plots"
    ),
    values = c(
      "Historical plots" = "#1f78b4",
      "Relocation error plots" = "#e31a1c"
    )
  ) +
  theme_classic() +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.background = element_rect(
      fill = "white",
      colour = "black"
    ),
    legend.title = element_blank()
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

ggsave(
  "outputs/figures/bayesian_figures/nmds_jaccard_2025_resurvey_correct.png",
  plot = nmds_jaccard_plot,
  width = 9,
  height = 7,
  dpi = 300
)



#### Save PERMANOVA models

saveRDS(
  perm_mod_bray,
  file = "outputs/models/perm_mod_bray_2025_resurvey_correct.rds"
)

saveRDS(
  perm_mod_jaccard,
  file = "outputs/models/perm_mod_jaccard_2025_resurvey_correct.rds"
)

saveRDS(
  disp_bray,
  file = "outputs/models/disp_bray_2025_resurvey_correct.rds"
)

saveRDS(
  disp_jaccard,
  file = "outputs/models/disp_jaccard_2025_resurvey_correct.rds"
)