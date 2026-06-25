# Basic exploratory Analysis of Manning Data
library(ggplot2)

#load data
all_data <- read.csv("manning_park_data/Manning_plots_allveg.csv")
all_data

##Subset and clean data

#Subset unique plots
plots <- all_data[!duplicated(all_data$PlotNumber),c("PlotNumber", "LocationAccuracy", "Latitude", "Longitude", "Date", "Elevation")]
plots

#Convert date to year column
install.packages("lubridate")
library(lubridate)

plots$Date <- as.Date(plots$Date)
plots$Year <- format(plots$Date, "%Y")

plots

##AVERAGE location accuracy: 

# Filter for years before 2015
plots_before_2015 <- subset(plots, Year < 2015)

# Calculate mean of LocationAccuracy, ignoring NA values
mean_location_accuracy <- mean(plots_before_2015$LocationAccuracy, na.rm = TRUE)

# View the result
mean_location_accuracy

plots_unique <- plots_before_2015[!duplicated(plots_before_2015[, c("Latitude", "Longitude")]), ]

# Calculate mean LocationAccuracy
mean_location_accuracy2 <- mean(plots_unique$LocationAccuracy, na.rm = TRUE)

# View result
mean_location_accuracy2

filtered_plots <- subset(plots, Year < 2015 & LocationAccuracy != 1000)

# Remove duplicates based on Latitude and Longitude
unique_plots2 <- filtered_plots[!duplicated(filtered_plots[, c("Latitude", "Longitude")]), ]

# Calculate mean LocationAccuracy
mean_location_accuracy3 <- mean(unique_plots2$LocationAccuracy, na.rm = TRUE)

# View result
mean_location_accuracy3

unique_plots <- filtered_plots[!duplicated(filtered_plots[, c("Latitude", "Longitude")]), ]

# Calculate median LocationAccuracy (including 1000s)
median_location_accuracy <- median(unique_plots$LocationAccuracy, na.rm = TRUE)

# View result
median_location_accuracy


#Subset 2019 plots

# Create a subset of unique plots with selected columns
plots <- all_data[!duplicated(all_data$PlotNumber), 
                  c("PlotNumber", "Latitude", "Longitude", "Date", "Elevation", "SpeciesListComplete")]

# Convert Date column to Date format if needed
plots$Date <- as.Date(plots$Date)

# Subset to include only plots from the year 2019
plots_2019 <- subset(plots, format(Date, "%Y") == "2019")
plots_2019
# Save to CSV
write.csv(plots_2019, "plots_2019.csv", row.names = FALSE)

##Exploring plots and years

range(as.numeric(plots$Year), na.rm = TRUE)
#Range is from 1976-2019

#Histogram of plots and years
ggplot(plots) +
  geom_histogram(aes(x = Year), stat = "count")

#number of plots 2015 and before
sum(as.numeric(plots$Year) <= 2015, na.rm = TRUE)
#[1] 63

#number of plots 2007 and before
sum(as.numeric(plots$Year) <= 2007, na.rm = TRUE)
#[1] 57

#number of plots 2000 and before
sum(as.numeric(plots$Year) <= 2000, na.rm = TRUE)
#[1] 50

#number of plots 1993 and before
sum(as.numeric(plots$Year) <= 1993, na.rm = TRUE)
#[1] 29

#number of plots from each year
sum(plots$Year == 2019, na.rm = T)
#[1] 52
sum(plots$Year == 2015, na.rm = T)
#[1] 6
sum(plots$Year == 2007, na.rm = T)
#[1] 7
sum(plots$Year == 2000, na.rm = T)
#[1] 21
sum(plots$Year == 1993, na.rm = T)
#[1] 7
sum(plots$Year == 1979, na.rm = T)
#[1] 7
sum(plots$Year == 1976, na.rm = T)
#[1] 15

##Exploring plots and elevation

range(as.numeric(plots$Elevation), na.rm = TRUE)
#Elevation ranges between 640-2336m


##Correlation between Year and Elevation

plots$Year <- as.numeric(plots$Year)
plots$Elevation <- as.numeric(plots$Elevation)

cor(plots$Year, plots$Elevation, use = "complete.obs")
#[1] 0.4445787

cor.test(plots$Year, plots$Elevation)
#p-value = 5.744e-07, significant positive correlation b/w elev. and year

plot(plots$Year, plots$Elevation,
     main = "Elevation vs Year",
     xlab = "Year", ylab = "Elevation",
     pch = 19, col = "blue")
abline(lm(Elevation ~ Year, data = plots), col = "red")


##Saving a table with PlotNumber, Year, and Elevation only

plots_subset <- plots[, c("PlotNumber", "Year", "Elevation")]
write.csv(plots_subset, "plots_year_elevation.csv", row.names = FALSE)




##Summary of brief analysis

#1. Years range from 1976-2019
#2. Elevation ranges from 640-2336m
#3. There are 63 plots from and before 2015, 57 from and before 2007, and 50 from and before 2000
#4. There is a significant weak positive correlation between year and elevation
#5. I saved a csv file with only plot, year, and elevation for easy reference
