#### Making a better plot map of Manning Park sites

## Read in data

data_all <- read_csv("data/processed/model_data_filt.csv")
error_data <- read_csv("data/processed/error_model_data.csv")

# Load packages

library(tidyverse)
library(sf)
library(ggplot2)
library(ggspatial)
library(osmdata)
library(cowplot)


# Prepare historical survey years

original_years <- data_all %>%
  filter(Time == "historical") %>%
  select(
    PlotNumber,
    OriginalSurveyYear = Year
  ) %>%
  distinct()


# Prepare present-day plots

present_plots <- data_all %>%
  filter(Time == "present") %>%
  select(
    PlotNumber,
    Latitude,
    Longitude
  ) %>%
  distinct() %>%
  left_join(
    original_years,
    by = "PlotNumber"
  ) %>%
  mutate(
    Location = "Correct"
  ) %>%
  filter(
    !is.na(Latitude),
    !is.na(Longitude),
    !is.na(OriginalSurveyYear)
  )


# Prepare relocation error plots

error_plots <- error_data %>%
  select(
    PlotNumber,
    Treatment,
    Latitude,
    Longitude
  ) %>%
  distinct() %>%
  filter(
    !is.na(Latitude),
    !is.na(Longitude)
  ) %>%
  mutate(
    OriginalSurveyYear = 2025,
    Location = if_else(
      Treatment == "Error",
      "Error",
      "Correct"
    )
  )


# Convert plots to spatial data

present_sf <- st_as_sf(
  present_plots,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)

error_sf <- st_as_sf(
  error_plots,
  coords = c("Longitude", "Latitude"),
  crs = 4326
)


# Get Manning Park boundary from OpenStreetMap

manning_boundary <- opq(
  "E.C. Manning Provincial Park",
  timeout = 120
) %>%
  add_osm_feature(
    key = "boundary",
    value = "protected_area"
  ) %>%
  osmdata_sf()

park_poly <- manning_boundary$osm_multipolygons

park_poly <- st_transform(
  park_poly,
  4326
)


# Check the park boundary

park_poly


# Set map extent

all_coordinates <- bind_rows(
  present_plots %>%
    select(
      Longitude,
      Latitude
    ),
  error_plots %>%
    select(
      Longitude,
      Latitude
    )
)

xlim <- range(
  all_coordinates$Longitude,
  na.rm = TRUE
)

ylim <- range(
  all_coordinates$Latitude,
  na.rm = TRUE
)

xpad <- diff(xlim) * 0.10
ypad <- diff(ylim) * 0.10

xlim <- xlim + c(-xpad, xpad)
ylim <- ylim + c(-ypad, ypad)


# Create main map

main_map <- ggplot() +
  
  annotation_map_tile(
    type = "https://tile.openstreetmap.org/{z}/{x}/{y}.png",
    zoom = 10,
    cachedir = "map_tiles",
    progress = "none"
  ) +
  
  geom_sf(
    data = park_poly,
    fill = NA,
    colour = "black",
    linewidth = 1
  ) +
  
  geom_sf(
    data = present_sf,
    aes(
      colour = factor(OriginalSurveyYear)
    ),
    shape = 16,
    size = 3.5,
    alpha = 0.85
  ) +
  
  geom_sf(
    data = error_sf %>%
      filter(Location == "Correct"),
    shape = 21,
    size = 5,
    fill = "white",
    colour = "black",
    stroke = 1.2
  ) +
  
  geom_sf(
    data = error_sf %>%
      filter(Location == "Error"),
    shape = 4,
    size = 5,
    colour = "black",
    stroke = 1.3
  ) +
  
  scale_colour_viridis_d(
    option = "turbo",
    name = "Original survey year"
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.8
  ) +
  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    style = north_arrow_fancy_orienteering(
      text_size = 10
    )
  ) +
  
  coord_sf(
    xlim = xlim,
    ylim = ylim,
    expand = FALSE
  ) +
  
  labs(
    x = NULL,
    y = NULL
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(
      5, 5, 5, 5
    )
  )


# Create relocation error legend

relocation_legend <- ggplot() +
  
  geom_point(
    aes(
      x = 1,
      y = 2
    ),
    shape = 21,
    size = 5,
    fill = "white",
    colour = "black",
    stroke = 1.2
  ) +
  
  geom_point(
    aes(
      x = 1,
      y = 1
    ),
    shape = 4,
    size = 5,
    colour = "black",
    stroke = 1.3
  ) +
  
  annotate(
    "text",
    x = 1.3,
    y = 2,
    label = "Correct relocation",
    hjust = 0,
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = 1.3,
    y = 1,
    label = "Relocation error",
    hjust = 0,
    size = 3.5
  ) +
  
  annotate(
    "text",
    x = 1,
    y = 2.7,
    label = "Relocation Error Plots",
    fontface = "bold",
    hjust = 0,
    size = 3.8
  ) +
  
  xlim(
    0.8,
    4.5
  ) +
  
  ylim(
    0.5,
    3
  ) +
  
  theme_void()


# Get British Columbia

bc <- rnaturalearth::ne_states(
  country = "Canada",
  returnclass = "sf"
) %>%
  filter(
    name_en == "British Columbia"
  )


# Create British Columbia inset

bc_inset <- ggplot() +
  
  geom_sf(
    data = bc,
    fill = "grey90",
    colour = "grey30",
    linewidth = 0.4
  ) +
  
  geom_sf(
    data = park_poly,
    fill = "black",
    colour = "black",
    alpha = 0.8
  ) +
  
  coord_sf(
    xlim = c(-139, -114),
    ylim = c(48, 61),
    expand = FALSE
  ) +
  
  theme_void()


# Combine map and inset

final_map <- ggdraw(
  main_map
) +
  
  draw_plot(
    bc_inset,
    x = 0.70,
    y = 0.68,
    width = 0.25,
    height = 0.25
  ) +
  
  draw_plot(
    relocation_legend,
    x = 0.68,
    y = 0.08,
    width = 0.27,
    height = 0.18
  )


# Display map

final_map


# Save high-resolution map

ggsave(
  "Manning_Park_map.png",
  final_map,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)