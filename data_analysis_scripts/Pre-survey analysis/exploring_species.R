##Exploring species and number of occurrences

data <- read.csv("manning_park_data/Manning_plots_allveg.csv")

species_counts <- as.data.frame(table(data$Species))
colnames(species_counts) <- c("Species", "Occurrences")

sum(as.numeric(species_counts$Occurrences) >= 10, na.rm = TRUE)
#[1] 90


##Create df with species 10 or more
common_species <- subset(species_counts, Occurrences >= 10)

write.csv(common_species, "common_species_10plus.csv", row.names = FALSE)


##Making a histogram of species occurrences
library(ggplot2)

ggplot(common_species) +
  geom_histogram(aes(x = Occurrences), stat = "count")

