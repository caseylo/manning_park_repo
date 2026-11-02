############################################################
# FIGURE 1
# Manning Park Study Area
############################################################

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggplot2)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(osmdata)
library(patchwork)
library(bcdata)

############################################################
# READ DATA
############################################################

data_all <- read_csv(
  "data/processed/model_data_filt.csv"
)

error_data <- read_csv(
  "data/processed/error_model_data.csv"
)

############################################################
# HISTORICAL PLOTS
############################################################

historical_sf <- data_all %>%
  filter(Time=="historical") %>%
  select(
    PlotNumber,
    Latitude,
    Longitude
  ) %>%
  distinct() %>%
  drop_na() %>%
  st_as_sf(
    coords = c("Longitude","Latitude"),
    crs = 4326
  )

############################################################
# RELOCATION PLOTS
############################################################

relocation_sf <- error_data %>%
  filter(Treatment=="Error") %>%
  select(
    PlotNumber,
    Latitude,
    Longitude
  ) %>%
  distinct() %>%
  drop_na() %>%
  st_as_sf(
    coords = c("Longitude","Latitude"),
    crs = 4326
  )

############################################################
# PARK BOUNDARY
############################################################

park_poly <- bcdc_get_data(
  "WHSE_TANTALIS.TA_PARK_ECORES_PA_SVW"
) %>%
  filter(
    str_detect(
      toupper(PROTECTED_LANDS_NAME),
      "MANNING"
    )
  ) %>%
  st_transform(4326)

############################################################
# BUFFER
############################################################

park_buffer <- park_poly %>%
  st_transform(3005) %>%
  st_buffer(5000) %>%
  st_transform(4326)

map_bbox <- st_bbox(park_buffer)

map_extent <- ext(
  map_bbox["xmin"],
  map_bbox["xmax"],
  map_bbox["ymin"],
  map_bbox["ymax"]
)

############################################################
# DEM
############################################################

dem <- rast(
  "data/spatial/manning_dem.tif"
)

dem <- project(
  dem,
  "EPSG:4326"
)

dem_crop <- crop(
  dem,
  ext(vect(park_buffer))
)

############################################################
# HILLSHADE
############################################################

slope <- terrain(
  dem_crop,
  "slope",
  unit="radians"
)

aspect <- terrain(
  dem_crop,
  "aspect",
  unit="radians"
)

hillshade <- shade(
  slope,
  aspect,
  angle = 45,
  direction = 315
)

############################################################
# CONTOURS
############################################################

contours_minor <- as.contour(
  dem_crop,
  levels = seq(
    500,
    2500,
    by = 250
  )
)

contours_major <- as.contour(
  dem_crop,
  levels = seq(
    500,
    2500,
    by = 250
  )
)

# contours_major is a SpatVector

contours_major <- crop(
  contours_major,
  vect(st_as_sfc(map_bbox))
)

# label locations
contour_labels <- centroids(
  contours_major,
  inside = TRUE
)

############################################################
# OSM WATER FEATURES
############################################################
library(osmdata)

#set_overpass_url("https://overpass.kumi.systems/api/interpreter")
set_overpass_url("https://lz4.overpass-api.de/api/interpreter")
#set_overpass_url("https://overpass.openstreetmap.ru/api/interpreter")


# bbox_poly <- st_as_sfc(
#   st_bbox(park_buffer)
# )
water_query <- opq(bbox = map_bbox)

lakes <- water_query %>%
  add_osm_feature(
    key = "natural",
    value = "water"
  ) %>%
  osmdata_sf()

lakes_sf <- lakes$osm_polygons
lakes_sf <- st_crop(
  lakes_sf,
  map_bbox
)

############################################################
# OSM RIVERS
############################################################

# rivers <- water_query %>%
#   add_osm_feature(
#     key = "waterway"
#   ) %>%
#   osmdata_sf()
# 
# rivers_sf <- rivers$osm_lines
# 
# rivers_sf <- st_crop(
#   rivers_sf,
#   map_bbox
# )

############################################################
# OSM ROADS
############################################################

roads <- opq(bbox = map_bbox) %>%
  add_osm_feature(
    key = "highway"
  ) %>%
  osmdata_sf()

roads_sf <- roads$osm_lines

############################################################
# MAJOR ROADS ONLY
############################################################

major_roads <- roads_sf %>%
  filter(
    highway %in%
      c(
        "motorway",
        "trunk",
        "primary",
        "secondary"
      )
  )

major_roads <- st_crop(
  major_roads,
  map_bbox
)
############################################################
# LOCATOR MAP
############################################################

north_america <- ne_countries(
  scale="medium",
  returnclass="sf"
)

bc <- ne_states(
  country="Canada",
  returnclass="sf"
) %>%
  filter(
    name=="British Columbia"
  )
park_point <- st_centroid(
  st_union(park_poly)
)
bc <- ne_states(
  country = "Canada",
  returnclass = "sf"
) %>%
  filter(name == "British Columbia")

bc_line <- st_boundary(bc)

bc_line <- st_intersection(
  bc_line,
  st_as_sfc(map_bbox)
)
##################

locator_map <- ggplot() +
  
  geom_sf(
    data=north_america,
    fill="grey95",
    colour="grey75",
    linewidth=0.15
  ) +
  
  geom_sf(
    data=bc,
    fill="#CFE8F3",
    colour="grey35",
    linewidth=0.3
  ) +
  
  # geom_sf(
  #   data = park_poly,
  #   fill = "#FEE0D2",
  #   colour = "red",
  #   linewidth = 0.3
  # ) +
  # geom_sf(
  #   data = park_point,
  #   shape = 21,
  #   fill = "red",
  #   colour = "black",
  #   size = 4
  # )+
  geom_sf_text(
    data = park_point,
    label = "\u2605",
    colour = "darkred",
    fontface = "bold",
    size = 5
  )+
  
  coord_sf(
    xlim=c(-140,-113),
    ylim=c(48,60),
    expand=FALSE
  ) +
  
  theme_void() +
  
  theme(
    panel.background=
      element_rect(
        fill="white",
        colour="black"
      )
  )

############################################################
# MAIN MAP
############################################################

study_map <- ggplot() +
  
  geom_spatraster(
    data = hillshade
  ) +
  
  scale_fill_gradient(
    low = "grey55",
    high = "white",
    guide = "none"
  ) +
  
  geom_sf(
    data = lakes_sf,
    fill = "#BFDFF5",
    colour = NA,
    alpha = 0.8
  ) +
  
  geom_spatvector(
    data = contours_minor,
    colour = "grey75",
    linewidth = 0.12,
    alpha = 0.6
  ) +
  
  geom_spatvector(
    data = contours_major,
    colour = "grey55",
    linewidth = 0.20,
    alpha = 0.65
  ) +
  
#  geom_sf_text(
#    data = contour_labels,
#    aes(label = paste0(level, " m")),
#    size = 2.8,
#    colour = "grey20",
#    fontface = "bold"
#  ) +
  
  geom_sf(
    data=park_poly,
    fill=NA,
    colour="black",
    linewidth=0.9
  ) +
  geom_hline(
    yintercept = 49,
    colour = "black",
    linewidth = 0.9,
    linetype = 1
  )+

  annotate(
    "text",
    x = -120.95,
    y = 48.99,
    label = "Canada–USA Border",
    size = 4,
    colour = "black",
    hjust = 0
  )+
  
  geom_sf(
    data=major_roads,
    colour="#C0392B",
    linewidth=0.7
  ) +  
  geom_sf(
    data = historical_sf,
    shape = 21,
    fill = "#0072B2",
    colour = "white",
    stroke = 0.7,
    size = 2.8
  ) +
  
  geom_sf(
    data=relocation_sf,
    shape=17,
    colour="#D55E00",
    size=2.8
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.15,
    pad_x = unit(0.1, "cm"),
    pad_y = unit(0.1, "cm"),
    line_width = 0.8,
    text_cex = 0.7
  )+

annotation_north_arrow(
  location = "tl",
  pad_x = unit(0.15, "cm"),
  pad_y = unit(0.15, "cm"),
  height = unit(1, "cm"),
  width = unit(1, "cm"),
  style = north_arrow_orienteering()
) +
  
  coord_sf(
    xlim = c(map_bbox["xmin"], map_bbox["xmax"]),
    ylim = c(map_bbox["ymin"], map_bbox["ymax"]),
    expand = FALSE
  )+
  
  labs(
    x=NULL,
    y=NULL
  ) +
  
  theme_bw() +
  
  theme(
    panel.grid=element_blank(),
    legend.position="none"
  )

############################################################
# MANUAL LEGEND
############################################################

# Legend box dimensions — bottom left, above scale bar
legend_xmin <- map_bbox["xmin"] + 0.009 * (map_bbox["xmax"] - map_bbox["xmin"])
legend_xmax <- map_bbox["xmin"] + 0.30 * (map_bbox["xmax"] - map_bbox["xmin"])

legend_ymin <- map_bbox["ymin"] + 0.05 * (map_bbox["ymax"] - map_bbox["ymin"])
legend_ymax <- map_bbox["ymin"] + 0.19 * (map_bbox["ymax"] - map_bbox["ymin"])

# Positions within the legend box
legend_x_point <- legend_xmin + 0.10 * (legend_xmax - legend_xmin)
legend_x_text  <- legend_xmin + 0.25 * (legend_xmax - legend_xmin)

legend_y_hist  <- legend_ymin + 0.70 * (legend_ymax - legend_ymin)
legend_y_reloc <- legend_ymin + 0.30 * (legend_ymax - legend_ymin)


study_map <- study_map +
  
  # Legend box
  annotate(
    "rect",
    xmin = legend_xmin,
    xmax = legend_xmax,
    ymin = legend_ymin,
    ymax = legend_ymax,
    fill = "white",
    colour = "black"
  ) +
  
  # Historical plots
  annotate(
    "point",
    x = legend_x_point,
    y = legend_y_hist,
    shape = 21,
    fill = "#0072B2",
    colour = "white",
    stroke = 0.8,
    size = 3.5
  ) +
  
  # Relocation plots
  annotate(
    "point",
    x = legend_x_point,
    y = legend_y_reloc,
    shape = 17,
    colour = "#D55E00",
    size = 3.5
  ) +
  
  # Historical label
  annotate(
    "text",
    x = legend_x_text,
    y = legend_y_hist,
    label = "Historical plots",
    hjust = 0,
    size = 3
  ) +
  
  # Relocation label
  annotate(
    "text",
    x = legend_x_text,
    y = legend_y_reloc,
    label = "Relocation plots",
    hjust = 0,
    size = 3
  )

############################################################
# ADD BC INSET
############################################################

locator_map <- locator_map +
  theme(
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.6
    )
  )

figure1 <- study_map +
  inset_element(
    locator_map,
    left   = 0.70,
    bottom = 0.68,
    right  = 0.98,
    top    = 0.98
  )

############################################################
# VIEW
############################################################

figure1

############################################################
# EXPORT
############################################################

ggsave(
  "outputs/figures/Figure1_ManningPark.tiff",
  figure1,
  width = 210,
  height = 120,
  units = "mm",
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "outputs/figures/Figure1_ManningPark.png",
  figure1,
  width = 210,
  height = 120,
  units = "mm"
)

