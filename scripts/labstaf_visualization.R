# Load necessary libraries
library(ggplot2)
library(dplyr)

# Function to create plots for JVPIIm over time for each station
plot_jvpiim_over_time <- function(data, output_dir = "figures") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Get unique stations
  stations <- unique(data$Station)

  # Loop over each station and create a plot
  for (station in stations) {
    # Filter data for the current station
    station_data <- data %>% filter(Station == station)

    # Create the plot with a smoother
    p <- ggplot(station_data, aes(x = Datetime, y = JVPIIm)) +
      geom_line(color = "blue", size = 1) +
      geom_point(color = "red", size = 2) +
      geom_smooth(method = "loess", color = "darkgreen", size = 1, se = FALSE) +
      theme_minimal(base_size = 15) +
      labs(
        title = paste("JVPIIm Over Time - Station", station),
        x = "Time",
        y = "JVPIIm (µmol photons / m³ / s)"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_datetime(date_labels = "%b %d, %H:%M", date_breaks = "1 hour")

    # Save the plot as a PNG file
    plot_filename <- paste0(output_dir, "/JVPIIm_Station_", station, ".png")
    ggsave(plot_filename, plot = p, width = 10, height = 6)
  }
}

# Function to create plots for Alpha over time for each station
plot_alpha_over_time <- function(data, output_dir = "figures") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Get unique stations
  stations <- unique(data$Station)

  # Loop over each station and create a plot
  for (station in stations) {
    # Filter data for the current station
    station_data <- data %>% filter(Station == station)

    # Create the plot with a smoother
    p <- ggplot(station_data, aes(x = Datetime, y = Alpha)) +
      geom_line(color = "blue", size = 1) +
      geom_point(color = "red", size = 2) +
      geom_smooth(method = "loess", color = "darkgreen", size = 1, se = FALSE) +
      theme_minimal(base_size = 15) +
      labs(
        title = paste("Alpha Over Time - Station", station),
        x = "Time",
        y = "Alpha"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_datetime(date_labels = "%b %d, %H:%M", date_breaks = "1 hour")

    # Save the plot as a PNG file
    plot_filename <- paste0(output_dir, "/Alpha_Station_", station, ".png")
    ggsave(plot_filename, plot = p, width = 10, height = 6)
  }
}

# Function to create plots for Ek over time for each station
plot_Ek_over_time <- function(data, output_dir = "figures") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Get unique stations
  stations <- unique(data$Station)

  # Loop over each station and create a plot
  for (station in stations) {
    # Filter data for the current station
    station_data <- data %>% filter(Station == station)

    # Create the plot with a smoother
    p <- ggplot(station_data, aes(x = Datetime, y = Ek)) +
      geom_line(color = "blue", size = 1) +
      geom_point(color = "red", size = 2) +
      geom_smooth(method = "loess", color = "darkgreen", size = 1, se = FALSE) +
      theme_minimal(base_size = 15) +
      labs(
        title = paste("Ek Over Time - Station", station),
        x = "Time",
        y = "Ek"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_datetime(date_labels = "%b %d, %H:%M", date_breaks = "1 hour")

    # Save the plot as a PNG file
    plot_filename <- paste0(output_dir, "/Ek_Station_", station, ".png")
    ggsave(plot_filename, plot = p, width = 10, height = 6)
  }
}

# Function to create plots for rPm over time for each station
plot_rPm_over_time <- function(data, output_dir = "figures") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Get unique stations
  stations <- unique(data$Station)

  # Loop over each station and create a plot
  for (station in stations) {
    # Filter data for the current station
    station_data <- data %>% filter(Station == station)

    # Create the plot with a smoother
    p <- ggplot(station_data, aes(x = Datetime, y = rPm)) +
      geom_line(color = "blue", size = 1) +
      geom_point(color = "red", size = 2) +
      geom_smooth(method = "loess", color = "darkgreen", size = 1, se = FALSE) +
      theme_minimal(base_size = 15) +
      labs(
        title = paste("rPm Over Time - Station", station),
        x = "Time",
        y = "rPm"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_datetime(date_labels = "%b %d, %H:%M", date_breaks = "1 hour")

    # Save the plot as a PNG file
    plot_filename <- paste0(output_dir, "/rPm_Station_", station, ".png")
    ggsave(plot_filename, plot = p, width = 10, height = 6)
  }
}

# Function to create plots for GOPIIm over time for each station
plot_GOPIIm_over_time <- function(data, output_dir = "figures") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Get unique stations
  stations <- unique(data$Station)

  # Loop over each station and create a plot
  for (station in stations) {
    # Filter data for the current station
    station_data <- data %>% filter(Station == station)

    # Create the plot with a smoother
    p <- ggplot(station_data, aes(x = Datetime, y = GOPIIm)) +
      geom_line(color = "blue", size = 1) +
      geom_point(color = "red", size = 2) +
      geom_smooth(method = "loess", color = "darkgreen", size = 1, se = FALSE) +
      theme_minimal(base_size = 15) +
      labs(
        title = paste("GOPIIm Over Time - Station", station),
        x = "Time",
        y = "GOPIIm"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_datetime(date_labels = "%b %d, %H:%M", date_breaks = "1 hour")

    # Save the plot as a PNG file
    plot_filename <- paste0(output_dir, "/GOPIIm_Station_", station, ".png")
    ggsave(plot_filename, plot = p, width = 10, height = 6)
  }
}

# Function to create plots for PP over time for each station
plot_PP_over_time <- function(data, output_dir = "figures") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Get unique stations
  stations <- unique(data$Station)

  # Loop over each station and create a plot
  for (station in stations) {
    # Filter data for the current station
    station_data <- data %>% filter(Station == station)

    # Create the plot with a smoother
    p <- ggplot(station_data, aes(x = Datetime, y = PP)) +
      geom_line(color = "blue", size = 1) +
      geom_point(color = "red", size = 2) +
      geom_smooth(method = "loess", color = "darkgreen", size = 1, se = FALSE) +
      theme_minimal(base_size = 15) +
      labs(
        title = paste("PP Over Time - Station", station),
        x = "Time",
        y = "PP mg C / m / h"
      ) +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_datetime(date_labels = "%b %d, %H:%M", date_breaks = "1 hour")

    # Save the plot as a PNG file
    plot_filename <- paste0(output_dir, "/PP_Station_", station, ".png")
    ggsave(plot_filename, plot = p, width = 10, height = 6)
  }
}

final_data <-read.csv("data/raw/LabSTAF/labstaf_combined_data.csv")
final_data$Datetime <- as.POSIXct(final_data$Datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Call the functions to generate plots
plot_jvpiim_over_time(final_data, output_dir = "figures/FRRF/")
plot_alpha_over_time(final_data,  output_dir = "figures/FRRF/")
plot_Ek_over_time(final_data,  output_dir = "figures/FRRF/")
plot_rPm_over_time(final_data,  output_dir = "figures/FRRF/")
plot_GOPIIm_over_time(final_data,  output_dir = "figures/FRRF/")
plot_PP_over_time(final_data, output_dir = "figures/FRRF/")
