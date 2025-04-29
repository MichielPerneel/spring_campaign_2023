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
env <- read.csv("data/samples_env_carbon.csv", header = TRUE)
o2_data <- read.csv("data/analysis/O2_smoother_station_130.csv", header = TRUE)
phaeo_KO <- read.csv("data/analysis/130_phaeo_KO_TPL_standardized_expression.csv",
                     header = TRUE, check.names = FALSE, row.names = 1)
photosynthesis_KOs <- read.table("data/annotation/functional_eggnog/photosynthesis_KOs.txt")
sun_data <- read.csv("data/analysis/sun_data.csv", header = TRUE, row.names = 1)

# Only retain photosynthesis-related KOs
phaeo_KO <- phaeo_KO[rownames(phaeo_KO) %in% photosynthesis_KOs$V1,]

# Convert Date columns to POSIXct
env$Date <- as.POSIXct(env$Date, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
o2_data$Date <- as.POSIXct(o2_data$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
sun_data$Date <- as.POSIXct(sun_data$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Retain only station 130 data
station <- 130
env <- env %>%
  filter(StationPrefix == station)
# Rename Station column to sample for merging
names(env)[names(env) == "Station"] <- "sample"

# Interpolate oxygen to match sample times
oxygen_interpolated <- approx(x = o2_data$Date, y = o2_data$y, xout = env$Date, method = "linear", rule = 2)
env$O2 <- oxygen_interpolated$y

# Merge sun altitude data
env <- env %>%
  left_join(sun_data %>% select(Date, SunAltitude), by = "Date")

#------ Define Tidal Cycle Based on Observations ------#

# Define observed high tides at station 130 (Ostend)
tide_times <- as.POSIXct(c("2023-04-20 01:45:00", "2023-04-20 14:12:00", "2023-04-21 02:26:00"), tz = "UTC")

# Find the most recent high tide for each sample
env$tide_ref <- sapply(env$Date, function(dt) max(tide_times[tide_times <= dt]))

# Compute hours since last high tide
env$hours_since_high_tide <- as.numeric(difftime(env$Date, env$tide_ref, units = "hours"))

# Normalize to tidal cycle (12.42h) and compute tidal phase
env$tide_phase <- cos((2 * pi / 12.42) * env$hours_since_high_tide)

# Scale tidal effect to match observed O2 amplitude
observed_O2_range <- max(env$O2, na.rm = TRUE) - min(env$O2, na.rm = TRUE)
env$O2_tidal <- env$tide_phase * (observed_O2_range * 0.4) + median(env$O2, na.rm = TRUE)  # Assume tide contributes ~40% of variation

#------ Define Diel Cycle Based on Sunlight Data ------#

# Normalize SunAltitude (0 at night, peaks at midday)
env$diel_phase <- env$SunAltitude / max(env$SunAltitude, na.rm = TRUE)

# Scale diel effect to match observed O2 amplitude
env$O2_diel <- env$diel_phase * (observed_O2_range * 0.6) + median(env$O2, na.rm = TRUE) # Assume diel contributes ~60% of variation

#------ Process phaeocystis photosynthesis gene expression ------#

# Process Phaeocystis KO expression (convert to long format)
phaeo_KO_long <- phaeo_KO %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(KO, TPL, -sample) %>%
  left_join(env %>% select(sample, O2_diel, O2_tidal, Salinity), by = "sample")

# Compute correlation of each KO with diel O2 model
correlations <- phaeo_KO_long %>%
  group_by(KO) %>%
  summarise(cor_with_diel_O2 = cor(TPL, O2_diel, use = "complete.obs")) %>%
  arrange(desc(cor_with_diel_O2))

# Select top 10 KOs with strongest positive correlation
top_10_KOs <- correlations$KO[1:10]

# Add category column for coloring (top 10 vs others)
phaeo_KO_long <- phaeo_KO_long %>%
  mutate(category = ifelse(KO %in% top_10_KOs, KO, "Other"))

# Standardize expression values
scale_z <- function(x) {(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)}

phaeo_KO_long <- phaeo_KO_long %>%
  group_by(KO) %>%
  mutate(TPL_z = scale(TPL)) %>%
  ungroup() %>%
  mutate(
    O2_diel_z = scale_z(O2_diel),
    O2_tidal_z = scale_z(O2_tidal),
    Salinity_z = scale_z(Salinity)
  )

# Merge env and phaeo_KO_long into new dataframe
env_2 <- env %>% right_join(phaeo_KO_long, by = "sample")

# Compute sum of expression per sample
phaeo_KO_sum <- phaeo_KO_long %>%
  group_by(sample) %>%
  summarise(photo_TPL = sum(TPL)) %>%
  ungroup() %>%
  left_join(env, by = "sample")

#------ Standardize (Z-Score Normalize) for Middle & Bottom Panels ------#
env <- env %>%
  mutate(
    O2_diel_z = scale_z(O2_diel),
    O2_tidal_z = scale_z(O2_tidal),
    photo_TPL_z = scale_z(photo_TPL),
    photo_TPL_lag2 = lag(photo_TPL, n = 2),  # 2-hour lag for gene expression
    photo_TPL_lag2_z = scale_z(photo_TPL_lag2),
    Salinity_z = scale_z(Salinity)
  )

#------ Visualizations ------#

# Observed O2 and Hypothetical Contributions with Dual Y-Axis
p0 <- ggplot(env, aes(x = Date)) +
  geom_line(aes(y = O2, color = "Observed O2"), size = 1) +
  geom_line(aes(y = O2_diel, color = "Diel Contribution"), linetype = "dashed", size = 1) +
  geom_line(aes(y = O2_tidal, color = "Tidal Contribution"), linetype = "dotted", size = 1) +
  labs(title = "Observed O2 and Hypothetical Contributions",
       x = "Time") +
  scale_color_manual(values = c("Observed O2" = "black", "Diel Contribution" = "blue", "Tidal Contribution" = "red")) +
  theme_minimal()

# Diel O2 Model vs Phaeocystis Expression
p1 <- ggplot(env_2, aes(x = Date)) +
  geom_line(aes(y = O2_diel_z, color = "Diel O2 Model"), linetype = "dashed", size = 1) +
  geom_line(data = subset(env_2, category == "Other"), aes(y = TPL_z, group = KO), color = "grey", alpha = 0.3) +
  geom_line(data = subset(env_2, category != "Other"), aes(y = TPL_z, color = KO, group = KO), size = 1) +
  scale_color_manual(values = c(setNames(rainbow(10), top_10_KOs), "Diel O2 Model" = "blue", "Other" = "grey")) +
  geom_line(aes(y = O2_diel_z, color = "Diel O2 Model"), linetype = "dashed", size = 1) +
  labs(title = "Phaeocystis Photosynthesis KO Expression vs. Diel O2 Contribution",
       y = "Z-Score (Standardized)", x = "Time") +
  scale_color_manual(values = c("Diel O2 Model" = "blue", "Phaeocystis Gene Expression" = "green")) +
  theme_minimal()

# Tidal O2 Model vs Salinity
p2 <- ggplot(env, aes(x = Date)) +
  geom_line(aes(y = O2_tidal_z, color = "Tidal O2 Model"), linetype = "dashed", size = 1) +
  geom_line(aes(y = Salinity_z, color = "Salinity"), size = 1) +
  labs(title = "Salinity vs. Tidal O2 Contribution",
       y = "Z-Score (Standardized)", x = "Time") +
  scale_color_manual(values = c("Tidal O2 Model" = "red", "Salinity" = "orange")) +
  theme_minimal()

# Combine all plots into one final figure
final_plot <- plot_grid(p0, p1, p2, ncol = 1, rel_heights = c(1.2, 1, 1))
print(final_plot)
