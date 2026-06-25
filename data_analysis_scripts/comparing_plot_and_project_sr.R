#### Comparing  species richness for each project WITH moss, lichens, grasses

## Load libaries

library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)

####

## Filtering out all genus level species ID

BEC_clean <- read_csv("manning_park_data/cleaning_data/BEC_clean.csv")

data_2025_clean <- read_csv("manning_park_data/cleaning_data/data_2025_clean.csv")

BEC_filt <- BEC_clean %>%
  filter(!grepl("SP\\.?$", Species))

data_2025_filt <- data_2025_clean %>%
  filter(!grepl("SP\\.?$", Species))

## Checking to make sure it worked

BEC_species <- as.data.frame(table(BEC_filt$Species))
colnames(BEC_species) <- c("species_BEC", "occurences_BEC")
BEC_species

species_2025 <- as.data.frame(table(data_2025_filt$Species))
colnames(species_2025) <- c("species_2025", "occurences_2025")
species_2025

####

## Comparing average species richness per plot per project in BEC data:

richness_BEC <- BEC_filt %>%
  group_by(ProjectID, PlotNumber) %>%
  summarise(
    species_richness = n_distinct(Species),
    .groups = "drop"
  ) %>%
  mutate(dataset = "BEC")

## In 2025 data:

richness_2025 <- data_2025_filt %>%
  group_by(ProjectID, PlotNumber) %>%
  summarise(species_richness = n_distinct(Species), .groups = "drop") %>%
  mutate(dataset = "2025") %>%
  filter(PlotNumber %in% richness_BEC$PlotNumber) ##Filtering out the 2019 error plots that were surveyed in 2025

## Combining both datasets

richness_all <- bind_rows(richness_BEC, richness_2025)

####

## Comparing average richness by ProjectID

richness_summary <- richness_all %>%
  group_by(dataset, ProjectID) %>%
  summarise(
    mean_richness = mean(species_richness),
    sd_richness   = sd(species_richness),
    n_plots       = n(),
    .groups = "drop")

richness_summary

####

## Statistical comparison of average species/ project using Mann Whitney U-test

richness_tests <- richness_all %>%
  group_by(ProjectID) %>%
  do(tidy(wilcox.test(species_richness ~ dataset, data = .)))

## So I can read the p-values better lol

richness_tests <- richness_tests %>%
  mutate(
    p.value = ifelse(
      p.value < 0.001,
      "<0.001",
      formatC(p.value, format = "f", digits = 3)))

richness_tests

## Violin plot to confirm distribution of SR of plots is not normal

ggplot(richness_all, aes(x = dataset, y = species_richness, fill = dataset)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 1) +
  facet_wrap(~ ProjectID, scales = "free_y") +
  labs(x = "Dataset", y = "Species richness per plot with functional groups") +
  theme(legend.position = "none")

####

##Visualize the difference in SR between sites

ggplot(richness_all, aes(x = dataset, y = species_richness, fill = dataset)) +
  geom_boxplot() +
  facet_wrap( ~ ProjectID) +
  labs(x = "Dataset", y = "Species richness per plot with functional groups")

####

## Comparing  species richness for each project WITHOUT moss, lichens, grasses

## Filtering out all genus level species ID and species with functional groups moss, lichen, or grass

BEC_filt_no_func <- BEC_clean %>%
  filter(!grepl("SP\\.?$", Species)) %>%
  filter(!grepl("^(moss|lichen|grass)$", Functional_group))

## Checking to make sure it worked

BEC_species_no_func <- as.data.frame(table(BEC_filt_no_func$Species))
colnames(BEC_species_no_func) <- c("species_BEC", "occurences_BEC")
BEC_species_no_func

#### 

##Comparing average species richness per plot per project in BEC data:

richness_BEC_no_func <- BEC_filt_no_func %>%
  group_by(ProjectID, PlotNumber) %>%
  summarise(
    species_richness = n_distinct(Species),
    .groups = "drop"
  ) %>%
  mutate(dataset = "BEC")

##Combining both datasets

richness_all_no_func <- bind_rows(richness_BEC_no_func, richness_2025)

####

## Comparing average richness by ProjectID

richness_summary_no_func <- richness_all_no_func %>%
  group_by(dataset, ProjectID) %>%
  summarise(
    mean_richness = mean(species_richness),
    sd_richness   = sd(species_richness),
    n_plots       = n(),
    .groups = "drop")

richness_summary_no_func

####

## Statistical comparison of average species/ project using Mann Whitney U-test

richness_tests_no_func <- richness_all_no_func %>%
  group_by(ProjectID) %>%
  do(tidy(wilcox.test(species_richness ~ dataset, data = .)))

## So I can read the p-values better lol

richness_tests_no_func <- richness_tests_no_func %>%
  mutate(
    p.value = ifelse(
      p.value < 0.001,
      "<0.001",
      formatC(p.value, format = "f", digits = 3)))

## Violin plot to confirm distribution of SR of plots is not normal

ggplot(richness_all_no_func, aes(x = dataset, y = species_richness, fill = dataset)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_jitter(width = 0.12, alpha = 0.6, size = 1) +
  facet_wrap(~ ProjectID, scales = "free_y") +
  labs(x = "Dataset", y = "Species richness per plot") +
  theme(legend.position = "none") +
  theme_classic()

####

## Visualize the difference in SR between sites

ggplot(richness_all_no_func, aes(x = dataset, y = species_richness, fill = dataset)) +
  geom_boxplot() +
  facet_wrap( ~ ProjectID) +
  labs(x = "Dataset", y = "Species richness per plot") +
  theme_classic()

####

## Saving results to csv files so I can paste tables into the results words doc

write.csv(richness_tests, "results/richness_tests_table.csv", row.names = FALSE)

write.csv(richness_tests_no_func, "results/richness_tests_table_no_func.csv", row.names = FALSE)

## Making the summary tables wider for with/without functinoal groups

## With functional groups
richness_summary_wide <- richness_summary %>%
  pivot_wider(
    id_cols = ProjectID,
    names_from = dataset,
    values_from = c(mean_richness, sd_richness, n_plots),
    names_glue = "{.value}_{dataset}")

## Without functional groups
richness_summary_wide_no_func <- richness_summary_no_func %>%
  pivot_wider(
    id_cols = ProjectID,
    names_from = dataset,
    values_from = c(mean_richness, sd_richness, n_plots),
    names_glue = "{.value}_{dataset}")

#### 

## Saving results to csv files so I can paste tables into the results words doc

write.csv(richness_summary_wide, "results/richness_summary_table.csv", row.names = FALSE)

write.csv(richness_summary_wide_no_func, "results/richness_summary_table_no_func.csv", row.names = FALSE)

####

## Running an ANOVA/ Kruskal Wallis on data (without functional groups) to test for significant differences between projects and within projects:

richness_all_no_func <- richness_all_no_func %>%
  mutate(
    ProjectID = as.factor(ProjectID),
    dataset   = as.factor(dataset))

richness_all_no_func <- richness_all_no_func %>%
  mutate(group = interaction(ProjectID, dataset, sep = "_"))

global_kruskal <- tidy(kruskal.test(species_richness ~ group, data = richness_all_no_func))
global_kruskal

summary(global_kruskal)
