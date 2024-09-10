# Install packages not found online
# Quick installer function
install_if_needed <- function(pkg, tar_file = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!is.null(tar_file)) {
      install.packages(tar_file, repos = NULL, type = "source")
    } else {
      install.packages(pkg)
    }
  }
}

# Paths to the .tar.gz files
tar_files <- list(
  "FME" = "scripts/R_packages/FME_1.3.6.3.tar.gz",
  "deSolve" = "scripts/R_packages/deSolve_1.40.tar.gz",
  "dtWad" = "scripts/R_packages/dtWad_0.0.1.tar.gz",
  "dtBioG" = "scripts/R_packages/dtBioG_0.0.1.tar.gz",
  "FRRF" = "scripts/R_packages/FRRF_1.0.tar.gz"
)

# Install packages from .tar.gz files if not installed
for (pkg in names(tar_files)) {
  install_if_needed(pkg, tar_file = tar_files[[pkg]])
}

# Load the libraries
library(tidyverse)
library(FME)
library(deSolve)
library(ggsci)
library(plot3D)
library(sf)
library(ncdf4)
library(dtWad)
library(dtBioG)
library(FRRF)
library(ggplot2)

# From the raw LabSTAF files, we need to create a dataframe to work with
## First, create a function that reads in the FRFF text data and stores additional info from the file header
read_LABstaf <- function(file, dir){
  # Split the filename by ".txt" to remove the extension
  specs <- strsplit(file, ".txt")[[1]]

  # Further split the result by underscores
  specs <- strsplit(specs, "_")[[1]]

  # Extract the station from the filename
  station <- specs[2]  # Assuming the second element is the station

  # Call the readFRRF function
  FRRF_data <- readFRRF(file=file, dir=dir, txt = "delim")

  # Add the station to the FRRF data
  FRRF_data$station <- station

  # Ensure the date column is in POSIXct format
  FRRF_data$date <- as.POSIXct(FRRF_data$date, format = "%b %d, %Y %H:%M")

  return(FRRF_data)
}

## Second, let's specify the directories on which to run the data
DIRS <- c("data/raw/LabSTAF/text_files/51", "data/raw/LabSTAF/text_files/130")

## Third, we'll run the data collection function on all files in the directories
LabSTAF <- NULL
for (dir in DIRS){
  file_list <- list.files(dir, pattern=".txt")
  for (file in file_list)
    LabSTAF <- rbind(LabSTAF,read_LABstaf(file, dir))
}
# Ensure column names are unique
colnames(LabSTAF) <- make.unique(colnames(LabSTAF))

# Restandardize, thereby converting JVPII estimates to mmol e/m3/hour
station_51 <- standardizeFRRF(subset(LabSTAF, subset = station==51), convJVPII = 3.6)
station_130 <- standardizeFRRF(subset(LabSTAF, subset = station==130), convJVPII = 3.6)

attributes(station_130)$aLHII_0
head(attributes(station_130)$processing)
# Create a plot for station 130
# Convert datetime to a factor for discrete coloring
station_130$datetime_factor <- as.factor(format(station_130$date, "%Y-%m-%d %H:%M"))

# Create the plot using a discrete color scheme
ggplot(station_130, aes(x = E, y = a_LHII, color = datetime_factor)) +
  geom_point(size = 3) + # Adjust the size of the points
  scale_color_manual(values = ggplot2::scale_colour_hue(c = 40, l = 70)(length(unique(station_130$datetime_factor)))) +
  labs(title = "a_LHII vs E Colored by Datetime",
       x = "E",
       y = "a_LHII",
       color = "Datetime") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

## WORKS UNTIL HERE
## TO DO:
# - standardize
# - Fit the PI regression to the labstaf data
# - Calculate PP from the fit by combining it with chl a data

# Function to create and save publication-ready plots
create_jvpii_plots <- function(data, output_dir) {
  # Ensure the output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Loop over each unique station
  for (station in unique(data$station)) {
    # Filter the data for the current station
    station_data <- data %>% filter(station == !!station)

    # Create the plot
    p <- ggplot(station_data, aes(x = date, y = JVPII)) +
      geom_line(size = 1, color = "blue") +
      geom_point(size = 2, color = "red") +
      labs(title = paste("JVPII over Time for Station", station),
           x = "Datetime",
           y = "JVPII",
           caption = "Data Source: LabSTAF") +
      scale_x_datetime(date_labels = "%b %d, %H:%M", date_breaks = "1 hour") +
      theme_minimal(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.title.x = element_text(size = 12),
            axis.title.y = element_text(size = 12),
            axis.text = element_text(size = 10),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    # Save the plot as a PNG file
    ggsave(filename = paste0(output_dir, "/JVPII_Station_", station, ".png"),
           plot = p,
           width = 10, height = 6, dpi = 300)

    # Print a message indicating the plot has been saved
    print(paste("Plot saved for station:", station))
  }
}

# Create an output directory to save plots
output_directory <- "figures/FRRF"

# Call the function to create and save plots
create_jvpii_plots(FRRF, output_directory)
