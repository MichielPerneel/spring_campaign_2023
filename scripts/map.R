# Load required libraries
library(dplyr)
library(tidyr)
library(sf)
library(terra)
library(downloader)
library(directlabels)
library(ggplot2)
library(rgl)
library(ncdf4)
library(mapdata)
library(geojsonio)
library(lattice)
library(reshape2)
library(XML)
library(raster)

# Load BPNS EEZ shapefile (download from https://www.marineregions.org/gazetteer.php?p=details&id=26567)
BPNS_EEZ <- st_read("./data/raw/eez_iho/eez_iho.shp")
# Convert to data frame for plotting
BPNS_EEZ <- st_as_sf(BPNS_EEZ)

# Get the spatial extent for the BPNS
xmin <- extent(BPNS_EEZ)@xmin
ymin <- extent(BPNS_EEZ)@ymin
xmax <- extent(BPNS_EEZ)@xmax
ymax <- extent(BPNS_EEZ)@ymax

# Define a function to read in raster data from the EMODnet bathymetry WCS
getbathymetry <- function(name = "emodnet:mean",
                          resolution = "0.2km",
                          xmin = 2, xmax = 3.5,
                          ymin = 51, ymax = 52) {
  bbox <- paste(xmin, ymin, xmax, ymax, sep = ",")

  con <- paste("https://ows.emodnet-bathymetry.eu/wcs?service=wcs&version=1.0.0&request=getcoverage&coverage=",
               name, "&crs=EPSG:4326&BBOX=", bbox,
               "&format=image/tiff&interpolation=nearest&resx=0.00208333&resy=0.00208333",
               sep = "")

  print(con)

  nomfich <- paste(name, "img.tiff", sep = "_")
  nomfich <- tempfile(nomfich)
  download(con, nomfich, quiet = TRUE, mode = "wb")

  # Use terra to read the downloaded raster file
  img <- rast(nomfich)
  img[img == 0] <- NA
  img[img < 0] <- 0
  names(img) <- paste(name)
  return(img)
}

# Get the bathymetry data for the BPNS using terra
bathy_img <- getbathymetry(name = "emodnet:mean",
                           resolution = "0.2km",
                           xmin = 2, xmax = 3.5,
                           ymin = 51, ymax = 52)

# Convert the raster to a data frame for ggplot2
bathy_df <- as.data.frame(bathy_img, xy = TRUE)

# Define sampling stations (replace with your actual coordinates)
stations <- data.frame(
  Name = c('Station 51', 'Station 130'),
  Latitude = c(51.53166115, 51.21634923),
  Longitude = c(3.18280367, 2.85125683)
)

# Create bathymetry map with sampling points
map <- ggplot() +
  geom_raster(data = bathy_df, aes(x = x, y = y, fill=bathy_df$`emodnet:mean`), alpha = 0.75) +
  scale_fill_gradient(low = "white", high = "darkblue", name = "Depth (m)") +
  geom_sf(data = BPNS_EEZ_fort, colour = "black", fill = NA) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  ggtitle("Bathymetry of the Belgian Part of the North Sea with Sampling Points") +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

# Add sampling stations as points on the map
map <- map +
  geom_point(data = stations, aes(x = Longitude, y = Latitude),
             color = "red", size = 3) +
  geom_text(data = stations, aes(x = Longitude, y = Latitude, label = Name),
            vjust = -1, hjust = 0.5, color = "black")

# Display the map
print(map)

# Optional: Add isobaths for more detailed visualization of depth contours
map_iso <- map +
  geom_contour(data = bathy_df, aes(x = x, y = y, z = bathy_df$`emodnet:mean`),
               colour = "gray30", na.rm = TRUE)

# Show the map with isobaths
print(map_iso)
