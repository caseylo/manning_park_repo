##Subsetting Data to look at individual years on a map


#load data
all_data <- read.csv("manning_park_data/Manning_plots_allveg.csv")
all_data

##Subset and clean data

#Subset unique plots
plots <- all_data[!duplicated(all_data$PlotNumber),c("PlotNumber", "Latitude", "Longitude", "Date", "Elevation")]
plots

#Convert date to year column
install.packages("lubridate")
library(lubridate)

plots$Date <- as.Date(plots$Date)
plots$Year <- format(plots$Date, "%Y")

##Check which plots have the same Lat/Long

same_lat_long <- plots[duplicated(plots[c("Latitude", "Longitude")]) | 
                         duplicated(plots[c("Latitude", "Longitude")], fromLast = TRUE), ]

same_lat_long

write.csv(same_lat_long, "8_duplicated_location_plots.csv", row.names = FALSE)
#There are 8 plots total that share the same Lat/long


##Adding a Year and Elevation label to  the PlotNumber

plots$PlotLabel <- paste(plots$PlotNumber, plots$Year, plots$Elevation_m, sep = " - ")

plots$Elevation_m <- paste(plots$Elevation, "m")

write.csv(plots, "PlotNumber_Year_Elevation_Label.csv", row.names = FALSE)




##For Loop to subset each year:
# Get unique years from the Year column
years <- unique(plots$Year)

# Loop through each year
for (yr in years) {
  # Create a filtered data frame for the current year
  temp_df <- plots[format(plots$Year) == yr, ]
  
  # Assign the data frame to a variable like plots_1976, plots_1977, etc.
  assign(paste0("plots_", yr), temp_df)
  
  # Print elevation range for the current year
  elev_range <- range(as.numeric(temp_df$Elevation), na.rm = TRUE)
  print(paste("Year:", yr, "- Elevation range:", elev_range[1], "to", elev_range[2]))
  
  # Save to CSV
  write.csv(temp_df, paste0("plots_", yr, ".csv"), row.names = FALSE)
}

#[1] "Year: 2019 - Elevation range: 1046 to 2336"
#[1] "Year: 2018 - Elevation range: 1055 to 1055"
#[1] "Year: 2015 - Elevation range: 906 to 1151"
#[1] "Year: 2007 - Elevation range: 2100 to 2225"
#[1] "Year: 2000 - Elevation range: 685 to 1840"
#[1] "Year: 1993 - Elevation range: 686 to 869"
#[1] "Year: 1979 - Elevation range: 640 to 1210"
#[1] "Year: 1976 - Elevation range: 1250 to 1921"

