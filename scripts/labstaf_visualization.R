library(ggplot2)
library(dplyr)
library(ggpubr)

# Function to create multipanel plots for a given parameter over time for both stations
plot_parameter <- function(data_51, data_130, parameter, y_label, output_dir = "figures/STAF/") {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  data_51 <- data_51 %>%
    filter(!is.na(!!sym(parameter)), !is.na(Date), is.finite(!!sym(parameter)))
  data_130 <- data_130 %>%
    filter(!is.na(!!sym(parameter)), !is.na(Date), is.finite(!!sym(parameter)))

  # Print min and max values with timestamps
  if (nrow(data_51) > 0) {
    min_51 <- data_51 %>% filter(!!sym(parameter) == min(!!sym(parameter), na.rm = TRUE))
    max_51 <- data_51 %>% filter(!!sym(parameter) == max(!!sym(parameter), na.rm = TRUE))
    cat(sprintf("\n[Station 51] %s:\n", parameter))
    cat(sprintf("  Min: %.3f at %s\n", min_51[[parameter]][1], format(min_51$Date[1], "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("  Max: %.3f at %s\n", max_51[[parameter]][1], format(max_51$Date[1], "%Y-%m-%d %H:%M:%S")))
  }

  if (nrow(data_130) > 0) {
    min_130 <- data_130 %>% filter(!!sym(parameter) == min(!!sym(parameter), na.rm = TRUE))
    max_130 <- data_130 %>% filter(!!sym(parameter) == max(!!sym(parameter), na.rm = TRUE))
    cat(sprintf("\n[Station 130] %s:\n", parameter))
    cat(sprintf("  Min: %.3f at %s\n", min_130[[parameter]][1], format(min_130$Date[1], "%Y-%m-%d %H:%M:%S")))
    cat(sprintf("  Max: %.3f at %s\n", max_130[[parameter]][1], format(max_130$Date[1], "%Y-%m-%d %H:%M:%S")))
  }

  # Colour scheme for light phases
  light_colors <- c(
    "Night" = "#d9d9d9",
    "Astronomical twilight" = "#ffb347",
    "Nautical twilight" = "#ffc870",
    "Civil twilight" = "#ffe0a3",
    "Day" = "#ffffb3"
  )

  first_break_51 <- min(data_51$Date, na.rm = TRUE)
  # Change first break to 08:00 that day
  first_break_51 <- as.POSIXct(paste0(format(first_break_51, "%Y-%m-%d"), " 08:00:00"), tz = "UTC")
  first_break_130 <- min(data_130$Date, na.rm = TRUE)
  # Change first break to 08:00 that day
  first_break_130 <- as.POSIXct(paste0(format(first_break_130, "%Y-%m-%d"), " 08:00:00"), tz = "UTC")

  y_limits <- range(c(data_51[[parameter]], data_130[[parameter]]), na.rm = TRUE)

  # Left panel (station 51)
  p1 <- ggplot(data_51, aes(x = Date, y = !!sym(parameter))) +
    geom_rect(data = data_51,
              aes(xmin = Date, xmax = lead(Date, default = last(Date)), ymin = -Inf, ymax = Inf, fill = day_moment), alpha = 0.6) +
    scale_fill_manual(values = light_colors, guide = guide_legend(title = "Light phase")) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, color = "black") +
    theme_minimal(base_size=10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    ) +
    scale_x_datetime(
      date_labels = "%b %d\n%H:%M",
      breaks = seq(first_break_51, max(data_51$Date, na.rm = TRUE), by = "4 hours"),
      minor_breaks = seq(first_break_51, max(data_51$Date, na.rm = TRUE), by = "2 hours")
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(y = y_label) +
    ggtitle("Station 51")

  # Right panel (station 130)
  p2 <- ggplot(data_130, aes(x = Date, y = !!sym(parameter))) +
    geom_rect(data = data_130,
              aes(xmin = Date, xmax = lead(Date, default = last(Date)), ymin = -Inf, ymax = Inf, fill = day_moment), alpha = 0.6) +
    scale_fill_manual(values = light_colors, guide = guide_legend(title = "Light phase")) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "loess", se = TRUE, color = "black") +
    theme_minimal(base_size=10) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(size = 10),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8)
    ) +
    scale_x_datetime(
      date_labels = "%b %d\n%H:%M",
      breaks = seq(first_break_130, max(data_130$Date, na.rm = TRUE), by = "4 hours"),
      minor_breaks = seq(first_break_130, max(data_130$Date, na.rm = TRUE), by = "2 hours")
    ) +
    coord_cartesian(ylim = y_limits) +
    ylab(NULL) +
    ggtitle("Station 130")

  # Join both panels in a single plot
  final_plot <- ggarrange(p1, p2, ncol = 2, align = "hv", common.legend = TRUE, legend = "none")

  # Save the plot as a PNG file
  plot_filename <- paste0(output_dir, parameter, "_Multipanel.png")
  ggsave(plot_filename, plot = final_plot, width = 12, height = 5.5, dpi = 1200, units = "cm")

  # Save the plot as a SVG file
  plot_filename <- paste0(output_dir, parameter, "_Multipanel.svg")
  ggsave(plot_filename, plot = final_plot, width = 12, height = 5.5, dpi = 1200, units = "cm")

  # Save the plot as a PDF file
  plot_filename <- paste0(output_dir, parameter, "_Multipanel.pdf")
  ggsave(plot_filename, plot = final_plot, width = 12, height = 5.5, dpi = 1200, units = "cm")
}

# Read the csv parsed data
station_130 <- read.csv("data/raw/LabSTAF/labstaf_csv_data_130.csv")
station_51 <- read.csv("data/raw/LabSTAF/labstaf_csv_data_51.csv")
metadata <- read.csv("data/samples_env.csv")

# Fix column names
colnames(station_51)[1:22] <- c("Station", "Alpha", "Beta", "Ek", "EkBeta", "rPm",
                                "JVPIIm", "GOPIIm", "Fo", "Fm", "Fv", "Fv_Fm", "Fv_Fmc",
                                "F_prime", "Fm_prime", "Fq_prime", "Fq_Fm_prime",
                                "Fq_Fmc_prime", "NPQ", "NSV", "SigmaPII_rhofit", "PP")
colnames(station_130)[1:22] <- c("Station", "Alpha", "Beta", "Ek", "EkBeta", "rPm",
                                 "JVPIIm", "GOPIIm", "Fo", "Fm", "Fv", "Fv_Fm", "Fv_Fmc",
                                 "F_prime", "Fm_prime", "Fq_prime", "Fq_Fm_prime",
                                 "Fq_Fmc_prime", "NPQ", "NSV", "SigmaPII_rhofit", "PP")

# Merge metadata to bring in Date and day_moment
station_51 <- station_51 %>%
  left_join(metadata, by = "Station") %>%
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  # Remove rows with very big negative values
  filter(JVPIIm > -1000, GOPIIm > -1000, PP > -1000) %>%
  # Sort according to Date
  arrange(Date) %>%
  mutate(
    GOPIIm = ifelse(day_moment == "Night", 0, GOPIIm),
    JVPIIm = ifelse(day_moment == "Night", 0, JVPIIm),
    PP = ifelse(day_moment == "Night", 0, PP)
  )

station_130 <- station_130 %>%
  left_join(metadata, by = "Station") %>%
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
  # Remove rows with very big negative values
  filter(JVPIIm > -1000, GOPIIm > -1000, PP > -1000) %>%
  # Sort according to Date
  arrange(Date) %>%
  mutate(
    GOPIIm = ifelse(day_moment == "Night", 0, GOPIIm),
    JVPIIm = ifelse(day_moment == "Night", 0, JVPIIm),
    PP = ifelse(day_moment == "Night", 0, PP)
  )

# Define biological parameters and their corresponding labels
parameters <- list(
  list("JVPIIm", "JVPIIm (µmol photons / m³ / s)"),
  list("Alpha", "Alpha"),
  list("Beta", "Beta"),
  list("Ek", "Ek (µmol photons / m² / s)"),
  list("EkBeta", "EkBeta (µmol photons / m² / s)"),
  list("rPm", "rPm"),
  list("GOPIIm", "GOPIIm (mmol O2 / m³ / h)"),
  list("PP", "PP (mg C/m²/h)"),
  list("NPQ", "Non-Photochemical Quenching (NPQ)"),
  list("NSV", "Non-Photochemical Quenching (NSV)"),
  list("Fo", "Fo (Minimum Fluorescence)"),
  list("Fm", "Fm (Maximum Fluorescence)"),
  list("Fv", "Fv (Variable Fluorescence)"),
  list("Fv_Fm", "Fv/Fm (PSII Efficiency)"),
  list("F_prime", "F' (Instantaneous Fluorescence)"),
  list("Fm_prime", "Fm' (Max Fluorescence under Light)"),
  list("Fq_prime", "Fq' (Fluorescence Quantum Yield)"),
  list("Fq_Fm_prime", "Fq'/Fm' (Effective PSII Quantum Yield)"),
  list("Fq_Fmc_prime", "Fq'/Fmc' (Alternative Yield Normalization)"),
  list("SigmaPII_rhofit", "SigmaPII (mmol photons / m² / s)")
)

# Generate multipanel plots for each parameter with station-specific shading
for (param in parameters) {
  plot_parameter(station_51, station_130, param[[1]], param[[2]])
}
