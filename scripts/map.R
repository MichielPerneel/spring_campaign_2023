#remotes::install_github("EMODnet/EMODnetWCS")

library(EMODnetWCS)
library(ows4R)
library(sf)
library(ggplot2)
library(terra)
library(tidyterra)
library(svglite)
library(ggspatial)

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

# Plot with adjustments
plot.new()
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
  geom_sf(mapping = aes(), data = stations, size = 0.5) +
  geom_text(aes(x = X, y = Y, label = label),
            data = stations,
            size = 2.5, nudge_x = 0, nudge_y = 0.05) +

  # Add a frame around the whole plot
  geom_rect(aes(xmin = 1.8, xmax = 3.8, ymin = 51, ymax = 52),
            fill = NA, color = "black", size = 0.2) +

  # Add north arrow and scale bar
  annotation_north_arrow(location = "tl", which_north = "true",
                         height = unit(0.5, "cm"), width = unit(0.5, "cm"),
                         style = north_arrow_minimal()) +
  annotation_scale(location = "br", bar_cols = c("black", "white"),
                   line_width = 0.5, height = unit(0.15, "cm"),
                   unit_category = "metric") +

  # Set blank theme
  theme_minimal() +
  theme(text = element_text(size = 8)) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.caption = element_text(hjust = 0.5))

# Save the updated map
ggsave("figures/environmental/BPNS_map.svg", map,
       width = 10, height = 8, unit = "cm", dpi = 600)
