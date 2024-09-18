#remotes::install_github("EMODnet/EMODnetWCS")
#remotes::install_github("EMODnet/EMODnetWFS")

#library(EMODnetWFS)
library(EMODnetWCS)
library(ows4R)
library(sf)
#library(ggsn)
library(ggplot2)
library(terra)
library(tidyterra)
library(grid)
library(svglite)

# Set stations
stations <- data.frame(
  label = c('Station 51', 'Station 130'),
  latitude = c(51.53166115, 51.21634923),
  longitude = c(3.18280367, 2.85125683)
  ) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Define area bounding box
bbox <- st_as_sfc("MULTIPOLYGON (((1.8 51, 1.8 52, 3.8 52, 3.8 51, 1.8 51)))")
st_crs(bbox) <- 4326

# Get coastline
wfs <- WFSClient$new("https://geo.vliz.be/geoserver/MarineRegions/wfs", "2.0.0")

countries <- wfs$getFeatures("MarineRegions:worldcountries_esri_2014",
       cql_filter = "territory IN ('Belgium', 'France', 'Netherlands')") %>%
       st_cast(to = "GEOMETRYCOLLECTION")

# Turn GEOMETRYCOLLECTION into multipolygons
geoms <- lapply(countries$the_geom, `[`)
mp <- lapply(geoms, function(x) sf::st_multipolygon(x = x))
sfc_mp <- sf::st_sfc(mp)
countries$mp <- sfc_mp
countries <- sf::st_set_geometry(countries, sfc_mp)
countries$mp <- NULL
rm(geoms);rm(mp);rm(sfc_mp)

# Crop countries to bounding box
st_crs(countries) <- 4326
countries <- countries %>% st_intersection(bbox)

# Get belgian part of the north sea
be_eez <- wfs$getFeatures("MarineRegions:eez_boundaries",
  cql_filter = "line_id IN (3730, 3729, 128, 127, 126, 125, 124)") %>%
  st_cast(to = "MULTILINESTRING") %>%
  st_cast(to = "LINESTRING")

st_crs(be_eez) <- 4326

# Get bathymetry
wcs <- emdn_init_wcs_client(service = "bathymetry")

cov <- emdn_get_coverage(wcs,
                         coverage_id = "emodnet__mean",
                         bbox = c(xmin = 1.8,
                                  ymin = 51,
                                  xmax = 3.8,
                                  ymax = 52)
)

# Mask bathymetry
cov_masked <- terra::mask(cov, countries, inverse = TRUE)

# Remove bathymetry values above 0
cov_masked[cov_masked > 0] <- NA

# Get the coordinates from the 'stations' sf object
station_coords <- st_coordinates(stations)

# Add a column for x and y coordinates to the stations data
stations$X <- station_coords[, 1]
stations$Y <- station_coords[, 2]

# Function to create a custom north arrow
north_arrow_custom <- function(x, y, size = 0.1) {
  annotation_custom(
    grob = textGrob("N", x = x, y = y + 0.1, gp = gpar(fontsize = 10)),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
  ) +
    annotation_custom(
      grob = linesGrob(arrow = arrow(type = "closed", length = unit(size, "inches")), gp = gpar(col = "black")),
      xmin = x, xmax = x, ymin = y, ymax = y + 0.1
    )
}

# Custom scale bar function
scale_bar_custom <- function(x, y, length_km, label_dist = 1, height = 0.02) {
  geom_segment(aes(x = x, xend = x + length_km / 100, y = y, yend = y), color = "black", size = 0.5) +
    geom_segment(aes(x = x + length_km / 100, xend = x + 2 * length_km / 100, y = y, yend = y), color = "black", size = 0.5) +
    annotate("text", x = x + length_km / 200, y = y - height, label = paste0(label_dist, " km"), size = 3, hjust = 0.5)
}

# Plot with adjustments
map <- ggplot() +
  # Add Bathymetry
  geom_spatraster(data = cov_masked) +
  scale_fill_hypso_c(
    palette = "wiki-2.0_bathy",
    name = "Depth (m)",
    breaks = c(0, -20, -40, -60, -80),
    labels = c("0", "20", "40", "60", "80")
  ) +

  # Add countries and boundaries background
  geom_sf(mapping = aes(), data = countries) +

  # Add stations and label them using the extracted coordinates
  geom_sf(mapping = aes(), data = stations, size = 1) +
  geom_text(aes(x = X, y = Y, label = label),
            data = stations,
            size = 3, nudge_x = 0, nudge_y = 0.05, fontface = "bold") +

  # Add a frame around the whole plot
  geom_rect(aes(xmin = 1.8, xmax = 3.8, ymin = 51, ymax = 52),
            fill = NA, color = "black", size = 0.2) +

  # Add custom north arrow
  north_arrow_custom(x = 3.5, y = 51.8, size = 0.2) +

  # Add custom scale bar (10 km)
  scale_bar_custom(x = 1.9, y = 51.1, length_km = 10, label_dist = 10, height = 0.01) +

  # Set theme
  theme_void() +
  theme(text = element_text(size = 8),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Sampling locations in the BPNS")

# Save the updated map
ggsave("figures/environmental/BPNS_map.svg", map,
       width = 10, height = 8, unit = "cm", dpi = 600)
