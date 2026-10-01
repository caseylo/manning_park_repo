## Making a better map

library(leaflet)
library(tidyverse)
library(sf)
library(osmdata)
library(ggplot2)
library(rnaturalearth)
library(base64enc)
library(htmltools)
library(mapview)
library(webshot)

# Read data
data_all <- read_csv("data/processed/model_data_filt.csv")
error_data <- read_csv("data/processed/error_model_data.csv")


# Prepare plot locations

present_plots <- data_all %>%
  filter(Time == "historical") %>%
  select(PlotNumber, Latitude, Longitude) %>%
  distinct() %>%
  drop_na(Latitude, Longitude)

error_plots <- error_data %>%
  select(PlotNumber, Treatment, Latitude, Longitude) %>%
  distinct() %>%
  drop_na(Latitude, Longitude) %>%
  mutate(Location = if_else(Treatment == "Error", "Error", "Correct"))

correct_plots <- error_plots %>%
  filter(Location == "Correct")

mislocated_plots <- error_plots %>%
  filter(Location == "Error")


# Manning Park boundary

manning_boundary <- opq(
  "E.C. Manning Provincial Park",
  timeout = 120
) %>%
  add_osm_feature(key = "boundary", value = "protected_area") %>%
  osmdata_sf()

park_poly <- manning_boundary$osm_multipolygons %>%
  st_transform(4326)


# British Columbia inset with ONE five-point star

north_america <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

bc <- ne_states(
  country = "Canada",
  returnclass = "sf"
) %>%
  filter(name == "British Columbia")

# Combine park polygons and create one point
manning_point <- st_sfc(
  st_point(c(-120.8, 49.2)),
  crs = 4326
)

bc_inset <- ggplot() +
  geom_sf(
    data = north_america,
    fill = "grey90",
    colour = "grey65",
    linewidth = 0.25
  ) +
  geom_sf(
    data = bc,
    fill = "#DCEEF2",
    colour = "grey30",
    linewidth = 0.4
  ) +
  geom_sf_text(
    data = manning_point,
    label = "\u2605",
    size = 6,
    colour = "black"
  ) +
  coord_sf(
    xlim = c(-139, -113),
    ylim = c(48, 60),
    expand = FALSE
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    plot.margin = margin(3, 3, 3, 3)
  )

inset_file <- tempfile(fileext = ".png")

ggsave(
  inset_file,
  plot = bc_inset,
  width = 3,
  height = 2.5,
  dpi = 300,
  bg = "white"
)

inset_image <- base64encode(
  readBin(inset_file, "raw", n = file.info(inset_file)$size)
)

inset_html <- paste0(
  '<div style="background:white;padding:7px;',
  'border:1px solid #888;border-radius:4px;">',
  '<img src="data:image/png;base64,',
  inset_image,
  '" style="width:500px;display:block;">',
  '</div>'
)


# Simple legend

legend_html <- paste0(
  '<div style="background:white;padding:14px 18px;',
  'border:2px solid #888;border-radius:6px;',
  'font-family:Arial,sans-serif;font-size:29px;',
  'line-height:1.6;min-width:260px;">',
  
  '<b>Resurvey Plots</b><br>',
  '<span style="color:#0072B2;font-size:38px;',
  '-webkit-text-stroke:2px white;paint-order:stroke fill;">▲</span> ',
  'Resurvey Plot<br>',
  
  '<b>Relocation Error Plots</b><br>',
  '<span style="color:black;font-size:38px;">○</span> Correct<br>',
  '<span style="color:black;font-size:38px;">×</span> Error',
  
  '</div>'
)


# Compass

compass_html <- paste0(
  '<div style="background:white;padding:16px 20px;',
  'border:2px solid #888;border-radius:5px;',
  'text-align:center;font-family:Arial,sans-serif;',
  'width:75px;">',
  '<div style="font-size:32px;font-weight:bold;">N</div>',
  '<div style="font-size:64px;line-height:1.1;">▲</div>',
  '</div>'
)
# Create map

manning_map <- leaflet(options = leafletOptions(
  zoomSnap = 0,
  zoomDelta = 0.0000001
)) %>%

  # Basemap
  addProviderTiles("Esri.WorldGrayCanvas") %>%
  addProviderTiles(
    "OpenTopoMap",
    options = providerTileOptions(opacity = 0.25)
  ) %>%
  
  # Manning Park boundary
  addPolygons(
    data = park_poly,
    color = "black",
    weight = 6,
    fill = FALSE
  ) %>%
  
  # Resurvey plots: blue triangles with white outline
  addLabelOnlyMarkers(
    data = present_plots,
    lng = ~Longitude,
    lat = ~Latitude,
    label = ~htmltools::HTML(
      paste0(
        "<span style='color:#0072B2;",
        "-webkit-text-stroke:1.5px white;",
        "paint-order:stroke fill;",
        "font-size:40px;font-weight:bold;'>▲</span>"
      )
    ),
    labelOptions = labelOptions(
      noHide = TRUE,
      textOnly = FALSE,
      direction = "center",
      style = list(
        "background" = "transparent",
        "border" = "none",
        "box-shadow" = "none",
        "padding" = "0px",
        "margin" = "0px"
      )
    )
  ) %>%
  
  # Correct relocation plots: open circles
  addCircleMarkers(
    data = correct_plots,
    lng = ~Longitude,
    lat = ~Latitude,
    radius = 16,
    color = "black",
    fillColor = "white",
    fillOpacity = 1,
    weight = 1.5
  ) %>%
  
  # Relocation error plots: black Xs
  addLabelOnlyMarkers(
    data = mislocated_plots,
    lng = ~Longitude,
    lat = ~Latitude,
    label = ~htmltools::HTML(
      "<span style='color:black;font-size:40px;font-weight:bold;'>×</span>"
    ),
    labelOptions = labelOptions(
      noHide = TRUE,
      textOnly = FALSE,
      direction = "center",
      style = list(
        "background" = "transparent",
        "border" = "none",
        "box-shadow" = "none",
        "padding" = "0px",
        "margin" = "0px"
      )
    )
  ) %>%
  
  # Legend, BC inset, compass and scale bar
  addControl(
    html = HTML(legend_html),
    position = "bottomright"
  ) %>%
  addControl(
    html = HTML(inset_html),
    position = "topright"
  ) %>%
  addControl(
    html = HTML(compass_html),
    position = "topleft"
  ) %>%
  addScaleBar(
    position = "bottomleft",
    options = scaleBarOptions(
      maxWidth = 300,
      metric = TRUE,
      imperial = FALSE
   )
  ) %>%
  
  setView(
    lng = -120.86,
    lat = 49.120,
    zoom = 12.4999999)

manning_map

manning_map <- htmlwidgets::prependContent(
  manning_map,
  htmltools::tags$style(
    htmltools::HTML("
      .leaflet-control-scale-line {
        font-size: 24px !important;
        font-weight: bold !important;
        line-height: 1.5 !important;
        padding: 8px 12px !important;
        border: 4px solid black !important;
        border-top: none !important;
        box-sizing: content-box !important;
        background: white !important;
        color: black !important;
      }
    ")
  )
)

mapview::mapshot(
  manning_map,
  file = "outputs/figures/bayesian_figures/Manning_Park_Map.png",
  remove_controls = FALSE,
  delay = 5,
  vwidth = 1800,
  vheight = 1200,
  zoom = 1.3
)


