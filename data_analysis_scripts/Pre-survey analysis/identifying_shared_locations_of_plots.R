#load data
all_data <- read.csv("manning_park_data/Manning_plots_allveg.csv")
all_data

##Subset and clean data

#Subset unique plots
plots_d <- all_data[!duplicated(all_data$PlotNumber),c("PlotNumber", "ProjectID", "Latitude", "Longitude", "Date", "Elevation")]
plots_d

list(plots$PlotNumber)

library(dyplr)
library(tidyverse)

duplicates <- plots_d %>%
  group_by(Latitude, Longitude) %>%
  filter(n() > 1) %>%
  arrange(Latitude, Longitude)

write.csv(duplicates, "Manning_plot_duplicates.csv", row.names = FALSE)
print(duplicates)


##Identifying plots with the same overall data

e85 <- all_data %>%
  filter(PlotNumber %in% c("essf85", "00NCe85"))
print(e85)

subset_data <- plots_d %>%
  filter(PlotNumber %in% c("essf85", "00NCe85")) %>%
  arrange(PlotNumber)  # Ensure consistent row order

# Remove the PlotNumber column
compare_data <- subset_data %>% select(-PlotNumber)

# Check if all other columns are equal
identical_rows <- all(compare_data[1, ] == compare_data[2, ])
print(identical_rows)

#Show if columns share data

if (!identical_rows) {
  differences <- names(compare_data)[compare_data[1, ] != compare_data[2, ]]
  print(differences)
}