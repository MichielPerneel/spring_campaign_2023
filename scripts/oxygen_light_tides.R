#------ Load Required Libraries ------#
library(dplyr)
library(ggplot2)
library(minpack.lm)
library(scales)  # For standardization
library(ggpubr)  # For combining plots
library(lubridate)
library(cowplot)  # For improved multi-panel plots

#------ Load and preprocess data ------#

# Load datasets
env <- read.csv("data/samples_env.csv", header = TRUE)
o2_data <- read.csv("data/analysis/O2_smoother_station_130.csv", header = TRUE)

# Convert Date columns to POSIXct
env$Date <- as.POSIXct(env$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
o2_data$Date <- as.POSIXct(o2_data$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Retain only station 130 data
station <- 130
env <- env %>%
  filter(StationPrefix == station)

# Interpolate oxygen to match sample times
oxygen_interpolated <- approx(x = o2_data$Date, y = o2_data$y, xout = env$Date, method = "linear", rule = 2)
env$O2 <- oxygen_interpolated$y

# Define colors for light phases
light_colors <- c(
  "Night" = "black",
  "Astronomical twilight" = "darkblue",
  "Nautical twilight" = "blue",
  "Civil twilight" = "lightblue",
  "Day" = "lightyellow"
)

# Date formatting as in other figures
generate_breaks <- function(start_date, end_date) {
  seq(
    from = floor_date(start_date, unit = "day"), # Round down to start of the day
    to = ceiling_date(end_date, unit = "day"),  # Round up to end of the day
    by = "4 hours"
  ) %>%
    .[format(., "%H:%M") %in% c("00:00", "08:00", "12:00", "16:00", "20:00")] # Keep only desired times
}

# Define a custom function for date labels
custom_date_format <- function(x) {
  ifelse(format(x, "%H:%M") == "00:00", format(x, "%d-%m"), format(x, "%H:%M"))
}

date_breaks <- generate_breaks(min(env$Date), max(env$Date))

#------ Visualization: O₂ Concentration, Tidal Dynamics, and Light Intensity ------#
# Plot O₂ Concentration
o2_plot <- ggplot(env, aes(x = Date)) +
  geom_line(aes(y = O2), size = 1) +
  labs(y = "O2 Concentration (µmol/L)",
       x = "Date") +
  ylim(250, 310) +
  theme_minimal() +
  scale_x_datetime(
    breaks = date_breaks,
    labels = custom_date_format
  )

# Plot Tidal Height
tide_plot <- ggplot(env, aes(x = Date)) +
  geom_rect(data = env,
            aes(xmin = Date, xmax = lead(Date, default = last(Date)), ymin = -Inf, ymax = Inf, fill = day_moment), alpha = 0.3) +
  scale_fill_manual(values = light_colors, guide = guide_legend(title = "Light Phase")) +
  geom_line(aes(y = sea_surface_height_above_sea_level), size = 1) +
  labs(y = "Tidal Height (m)",
       x = "Date") +
  ylim(-3, 3) +
  theme_minimal() +
  scale_x_datetime(
    breaks = date_breaks,
    labels = custom_date_format
  )

# Combine Plots
final_plot <- plot_grid(o2_plot, tide_plot, ncol = 1, align = "v", rel_heights = c(1, 1))
print(final_plot)

# Save plot as svg
ggsave("figures/environmental/o2_light_tide_130.svg", final_plot, width = 13, height = 10, dpi = 800, units = "cm")

# Repeat for station 51

#------ Load and preprocess data ------#
# Load datasets
env <- read.csv("data/samples_env.csv", header = TRUE)
o2_data <- read.csv("data/analysis/O2_smoother_station_51.csv", header = TRUE)

# Convert Date columns to POSIXct
env$Date <- as.POSIXct(env$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
o2_data$Date <- as.POSIXct(o2_data$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Retain only station 130 data
station <- 51
env <- env %>%
  filter(StationPrefix == station)

# Interpolate oxygen to match sample times
oxygen_interpolated <- approx(x = o2_data$Date, y = o2_data$y, xout = env$Date, method = "linear", rule = 2)
env$O2 <- oxygen_interpolated$y

date_breaks <- generate_breaks(min(env$Date), max(env$Date))

#------ Visualization: O₂ Concentration, Tidal Dynamics, and Light Intensity ------#
# Plot O₂ Concentration
o2_plot <- ggplot(env, aes(x = Date)) +
  geom_line(aes(y = O2), size = 1) +
  labs(y = "O2 Concentration (µmol/L)",
       x = "Date") +
  ylim(250, 310) +
  theme_minimal() +
  scale_x_datetime(
    breaks = date_breaks,
    labels = custom_date_format
  )

# Plot Tidal Height
tide_plot <- ggplot(env, aes(x = Date)) +
  geom_rect(data = env,
            aes(xmin = Date, xmax = lead(Date, default = last(Date)), ymin = -Inf, ymax = Inf, fill = day_moment), alpha = 0.3) +
  scale_fill_manual(values = light_colors, guide = guide_legend(title = "Light Phase")) +
  geom_line(aes(y = sea_surface_height_above_sea_level), size = 1) +
  labs(y = "Tidal Height (m)",
       x = "Date") +
  ylim(-3, 3) +
  theme_minimal() +
  scale_x_datetime(
    breaks = date_breaks,
    labels = custom_date_format
  )

# Combine Plots
final_plot <- plot_grid(o2_plot, tide_plot, ncol = 1, align = "v", rel_heights = c(1, 1))
print(final_plot)

# Save plot as svg
ggsave("figures/environmental/o2_light_tide_51.svg", final_plot, width = 13, height = 10, dpi = 800, units = "cm")
