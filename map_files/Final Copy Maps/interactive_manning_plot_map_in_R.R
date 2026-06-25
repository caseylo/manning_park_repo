##Making an interactive plot map of Manning Park of Historical Plots

library(leaflet)
library(lubridate)
library(readr)
library(dplyr)

#Load location, plot number, year, and elevation data

data <- read_csv("manning_park_data/ClimateNA_locations_manning.csv")
model_data <- read_csv("manning_park_data/model_data.csv")
BEC_clean <- read_csv("manning_park_data/cleaning_data/BEC_clean.csv")

##Create leaflet map

data_filtered <- data %>%
  filter(ID1 %in% model_data$PlotNumber)

plots_needed <- c("14-4311", "14-4316", "14-4317")

bec_subset <- BEC_clean %>%
  filter(PlotNumber %in% plots_needed) %>%
  select(PlotNumber, Latitude, Longitude, Date, Elevation) %>%
  mutate(Year = year(as.Date(Date))) %>%
  mutate(ID1 = PlotNumber,
         ID2 = Year,
         lat = Latitude,
         long = Longitude,
         el = Elevation)

bec_subset <- bec_subset %>%
  select(any_of(names(data_filtered)))

data_filtered <- data_filtered %>%
  bind_rows(bec_subset) %>%
  distinct(ID1, .keep_all = TRUE) %>%
  mutate(ID2 = if_else(ID1 == "00NCe91" & is.na(ID2), 2000, ID2))


# Create popup/label text
data_filtered <- data_filtered %>%
  mutate(
    label_text = paste0(
      "Plot: ", ID1, "<br>",
      "Year: ", ID2, "<br>",
      "Elevation: ", el, " m<br>",
      "Latitude:  ", lat, "<br>",
      "Longitude:  ", long, "<br>"))


# Create a color palette for the 'Year' column
year_palette <- colorFactor(palette = "Set1", domain = data_filtered$ID2)

#Create the map with colored points (Use this, no topo map but basic map)

leaflet(data_filtered) %>%
  addTiles() %>%  # grey topo style
  addCircleMarkers(
    ~long, ~lat,
    radius = 3,
    color = ~year_palette(ID2),
    stroke = TRUE,
    fillOpacity = 0.8,
    label = ~as.character(ID1),
    popup = ~lapply(label_text, htmltools::HTML)
  ) %>%
  addLegend(
    "bottomright",
    pal = year_palette,
    values = ~ID2,
    title = "Year",
    opacity = 1
  )

##Light coloured topomap, not very detailed but shows general topo):
leaflet(data_filtered) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addProviderTiles("OpenTopoMap", options = providerTileOptions(opacity = 0.25)) %>%  # subtle topo texture
  addCircleMarkers(
    ~long, ~lat,
    radius = 3,
    color = ~year_palette(ID2),
    stroke = TRUE,
    fillOpacity = 0.8,
    label = ~as.character(ID1),
    popup = ~lapply(label_text, htmltools::HTML)
  ) %>%
  addLegend(
    "bottomright",
    pal = year_palette,
    values = ~ID2,
    title = "Year",
    opacity = 1
  )


## Circle shapes with park boundary
library(osmdata)
library(sf)

manning_boundary <- opq(
  "E.C. Manning Provincial Park",
  timeout = 120
) %>%
  add_osm_feature(key = "boundary", value = "protected_area") %>%
  osmdata_sf()
park_poly <- manning_boundary$osm_multipolygons

park_poly <- st_transform(park_poly, 4326)

leaflet(data_filtered) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addProviderTiles("OpenTopoMap", options = providerTileOptions(opacity = 0.25)) %>%
  addPolygons(
    data = park_poly,
    color = "black",
    weight = 2,
    fill = FALSE
  ) %>%
  addCircleMarkers(
    ~long, ~lat,
    radius = 4.7,
    color = ~year_palette(ID2),
    stroke = TRUE,
    fillOpacity = 0.8,
    label = ~as.character(ID1),
    popup = ~lapply(label_text, htmltools::HTML)
  ) %>%
  addLegend(
    "bottomright",
    pal = year_palette,
    values = ~ID2,
    title = "Year",
    opacity = 1
  )

####

## X site shapes
yrs <- sort(unique(data_filtered$ID2))

park_poly <- manning_boundary$osm_multipolygons %>%
  st_transform(4326)

m <- leaflet(data_filtered) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addProviderTiles("OpenTopoMap", options = providerTileOptions(opacity = 0.25)) %>%
  addPolygons(
    data = park_poly,
    color = "darkgreen",
    weight = 2,
    fill = FALSE,
    opacity = 0.8
  ) %>%
  addCircleMarkers(
    ~long, ~lat,
    radius = 8,
    opacity = 0,
    fillOpacity = 0,
    popup = ~lapply(label_text, htmltools::HTML)
  )

for (y in yrs) {
  dat_y <- data_filtered[data_filtered$ID2 == y, ]
  
  m <- m %>%
    addLabelOnlyMarkers(
      data = dat_y,
      ~long, ~lat,
      label = ~htmltools::HTML(
        paste0("<span style='color:", year_palette(y), ";'>✕</span>")
      ),
      labelOptions = labelOptions(
        noHide = TRUE,
        textOnly = FALSE,
        direction = "center",
        style = list(
          "background" = "transparent",
          "border" = "none",
          "box-shadow" = "none",
          "font-size" = "25px",
          "font-weight" = "bold",
          "padding" = "0px",
          "margin" = "0px"
        )
      )
    )
}

m %>%
  addLegend(
    "bottomright",
    pal = year_palette,
    values = ~ID2,
    title = "Survey Year",
    opacity = 1
  ) %>%
  addScaleBar(position = "topright")
