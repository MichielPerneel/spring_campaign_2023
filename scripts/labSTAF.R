library(ggplot2)
library(dplyr)
library(ggpubr)
library(stringr)

# Function to create multipanel plots for a given parameter over time for both stations
plot_parameter_long <- function(df, parameter, y_label, output_dir = "figures/STAF_updated/") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Subset data
  data_sub <- df %>%
    filter(Parameter == parameter, is.finite(Value), !is.na(datetime), !is.na(Value))

  # Print min/max for each station
  for (st in unique(data_sub$station)) {
    st_data <- data_sub %>% filter(station == st)
    if (nrow(st_data) > 0) {
      min_val <- st_data %>% filter(Value == min(Value, na.rm = TRUE))
      max_val <- st_data %>% filter(Value == max(Value, na.rm = TRUE))
      cat(sprintf("\n[Station %s] %s:\n", st, parameter))
      cat(sprintf("  Min: %.3f at %s\n", min_val$Value[1], format(min_val$datetime[1], "%Y-%m-%d %H:%M:%S")))
      cat(sprintf("  Max: %.3f at %s\n", max_val$Value[1], format(max_val$datetime[1], "%Y-%m-%d %H:%M:%S")))
    }
  }

  # Colour scheme for light phases
  light_colors <- c(
    "Night" = "#d9d9d9",
    "Astronomical twilight" = "#ffb347",
    "Nautical twilight" = "#ffc870",
    "Civil twilight" = "#ffe0a3",
    "Day" = "#ffffb3"
  )

  # Create plot per station
  p <- ggplot(data_sub, aes(x = diel_index, y = Value)) +
    geom_rect(
      data = data_sub %>% filter(station == 130, !is.na(lead(diel_index))),
      aes(xmin = diel_index,
          xmax = lead(diel_index),
          ymin = -Inf, ymax = Inf,
          fill = day_moment),
      alpha = 0.6,
      show.legend = FALSE
    ) +
    ggtitle(y_label) +
    scale_fill_manual(values = light_colors) +
    geom_point(aes(color = factor(station)), alpha = 0.7, size = 1) +
    geom_smooth(aes(color = factor(station)), method = "loess", se = TRUE) +
    scale_color_manual(values = c("51" = "#DE9120", "130" = "#0173B2"), name = "Station") +
    labs(
      y = NULL, x = NULL
    ) +
    scale_x_continuous(
      breaks = seq(min(data_sub$diel_index, na.rm = TRUE), max(data_sub$diel_index, na.rm = TRUE), by = 4),
      labels = function(x) {
        times <- (x - 1) %% 24 + 8  # Start at 08:00
        sprintf("%02d:00", times %% 24)
      }) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 10),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Save plot
  ggsave(paste0(output_dir, parameter, "_Multipanel.png"), p, width = 4.5, height = 5, dpi = 1200, units = "cm")
  ggsave(paste0(output_dir, parameter, "_Multipanel.svg"), p, width = 4.5, height = 5, dpi = 1200, units = "cm")
  ggsave(paste0(output_dir, parameter, "_Multipanel.pdf"), p, width = 4.5, height = 5, dpi = 1200, units = "cm")
}

# Read the csv parsed data
station_130 <- read_csv2("data/raw/LabSTAF/all_130_NF_final.csv")
station_51 <- read_csv2("data/raw/LabSTAF/all_051_NF_final.csv")

station_all <- bind_rows(station_130, station_51)

station_all <- station_all %>%
  group_by(station) %>%
  mutate(nPSII = Ka * ((cFo.1 * 1e-6)/(SigmaPII * 1e-18)) * cPEC,
          PP = GOPIIm * 12,
          cYII.1 = 1 - (cYNPQ.1 + cYNO.1),
          cYII.2 = 1 - (cYNPQ.2 + cYNO.2),
          cYII.3 = 1 - (cYNPQ.3 + cYNO.3),
          cYII.4 = 1 - (cYNPQ.4 + cYNO.4),
          cYII.5 = 1 - (cYNPQ.5 + cYNO.5),
          cYII.6 = 1 - (cYNPQ.6 + cYNO.6),
          cYII.7 = 1 - (cYNPQ.7 + cYNO.7),
          cYII.8 = 1 - (cYNPQ.8 + cYNO.8),
          cYII.9 = 1 - (cYNPQ.9 + cYNO.9),
          cYII.10 = 1 - (cYNPQ.10 + cYNO.10),
          cYII.11 = 1 - (cYNPQ.11 + cYNO.11),
          cYII.12 = 1 - (cYNPQ.12 + cYNO.12),
          SUM.2 = cYII.2 + cYNPQ.2 + cYNO.2,
         datetime = lubridate::ymd_hms(datetime) - lubridate::hours(1),
         sample_id = as.numeric(str_extract(Note, "[0-9]+$")),
         diel_index = case_when(
           station == 51  ~ sample_id + 3,
           station == 130 ~ sample_id,
           TRUE ~ NA_real_
         )) %>%
  ungroup()

photophys_long <- station_all %>%
  # Select columns to keep
  select(station, datetime, hours, Note, diel_index, Alpha, AlphaPII, Rho, SigmaPII, sig.Ek,
         sig.PAR_500, 'Fv/Fm', cFqFm.Ek, cFqFm.PAR_500, cNPQ.Ek, cNPQ.PAR_500, GOPIIm,
         Ek, cNSV.1, cNSV.Ek, cNSV.PAR_500, cFqFm.1, cFqFm.2, JVPIIm, JVPII.Ek,
         JVPII.PAR_500, nPSII, PP, cYII.1) %>%
  pivot_longer(cols = c(Alpha, AlphaPII, Rho, SigmaPII, sig.Ek, sig.PAR_500,
                        'Fv/Fm', cFqFm.Ek, cFqFm.PAR_500, cNPQ.Ek, cNPQ.PAR_500,
                        Ek,cNSV.1, cNSV.Ek, cNSV.PAR_500, cFqFm.1, cFqFm.2, JVPIIm,
                        JVPII.Ek, JVPII.PAR_500, nPSII, PP, cYII.1, GOPIIm),
               names_to = "Parameter",
               values_to = "Value") %>%
  # Fix values in the Note column, 130_01 should be 130_1 for example, 051_01 should be 51_1
  mutate(
    Note = str_replace(Note, "^0?(51|130)_0*([1-9][0-9]*)$", "\\1_\\2"),
    Parameter = str_replace(Parameter, "Fv/Fm", "FvFm") # This otherwise makes directories when plotting
  )

metadata <- read.csv("data/samples_env.csv")

# Merge metadata to bring in Date and day_moment
photophys_long <- photophys_long %>%
  left_join(metadata, by = c("Note" = "Station"))

# Define biological parameters and their corresponding labels
parameters <- list(
  list("JVPIIm", "JVPIIm"),
  list("JVPII.Ek", "JVPII(Ek)"),
  list("Ek", "Ek"),
  list("PP", "PPm"),
  list("SigmaPII", "SigmaPII"),
  list("FvFm", "Fv/Fm"),
  list("cFqFm.Ek", "Fq'/Fmc'(Ek)"),
  list("cNSV.PAR_500", "NPQ_NSV(500)"),
  list("GOPIIm", "GOPIIm"),
  list("nPSII", "[PSII]"),
  list("sig.Ek", "sigmaPSII'(Ek)"),
  list("cNPQ.PAR_500", "NPQ(500)")
)

# Generate multipanel plots for each parameter with station-specific shading
for (param in parameters) {
  plot_parameter_long(photophys_long, param[[1]], param[[2]])
}
