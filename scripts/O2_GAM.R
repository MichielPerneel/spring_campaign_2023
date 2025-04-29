library(vegan)
library(dplyr)
library(funrar)
library(tidyr)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(reshape2)
library(mgcv)
library(lme4)
library(car)
library(gridExtra)
library(tibble)

#------ Load and preprocess data ------#
sun_data <- read.csv("data/analysis/sun_data.csv", header = TRUE, row.names = 1)
o2_data <- read.csv("data/analysis/O2_smoother_station_130.csv", header = TRUE)
carbon <- read.csv("data/samples_env_carbon.csv", header = TRUE, row.names = 1)
phaeo_KO <- read.csv("data/analysis/130_phaeo_KO_TPL_standardized_expression.csv", header = TRUE, check.names = FALSE, row.names = 1)
photosynthesis_KOs <- read.table("data/annotation/functional_eggnog/photosynthesis_KOs.txt")
# Only retain photosynthesis KOs
phaeo_KO <- phaeo_KO[rownames(phaeo_KO) %in% photosynthesis_KOs$V1,]
zooscan <- read.csv("data/raw/zooscan_data.csv", header = TRUE, sep = ";", row.names = 1)

# Ensure correct column names
o2_data <- rename(o2_data, Date = Date.Time)

# Convert Phaeocystis KO dataframe to long format
phaeo_KO_long <- phaeo_KO %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  gather(KO, TPL, -sample)

# Process carbon and environmental data
carbon <- carbon %>%
  filter(StationPrefix == "130") %>%
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")) %>%
  select(-StationPrefix) %>%
  rownames_to_column("sample")

# Convert Date columns to POSIXct
carbon$Date <- as.POSIXct(carbon$Date, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
o2_data$Date <- as.POSIXct(o2_data$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
sun_data$Date <- as.POSIXct(sun_data$Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Process oxygen interpolation
oxygen_interpolated <- approx(x = o2_data$Date, y = o2_data$y, xout = carbon$Date, method = "linear", rule = 2)
oxygen_interpolated <- data.frame(Date = carbon$Date, O2 = oxygen_interpolated$y)

sun_data <- sun_data %>%
  select(Station, Date, SunAzimuth, SunAltitude) %>%
  mutate(Time_since_start = as.numeric(difftime(Date, min(Date), units = "secs"))) %>%
  select(-Date) %>%
  rename(sample = Station)

# Process zooscan data
zooscan <- zooscan %>%
  filter(grepl("130", Station)) %>%
  # Station should be rownames instead of the dates now
  tibble::rownames_to_column("Date") %>%
  select(-Date, -Month, -Year, -Day, -Hour, -Minute) %>%
  # Rename Station to sample
  rename(sample = Station) %>%
  na.omit()

## Heatmap of Phaeocystis photosynthesis KOs expression over time
phaeo_KO_plot <- phaeo_KO %>%
  t() %>%
  melt() %>%
  # Rename columns
  rename(c("sample" = "Var1", "KO" = "Var2", "TPL" = "value")) %>%
  merge(carbon, by = "sample")

## Plot the heatmap
ggplot(phaeo_KO_plot, aes(x = Date, y = KO, fill = TPL)) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Phaeocystis photosynthesis KOs expression over time")

## Plot the sums of KO expression over time
## Group by Date and sum the TPL values
phaeo_KO_sum <- phaeo_KO_long %>%
  group_by(sample) %>%
  summarise(TPL = sum(TPL)) %>%
  merge(carbon, by = "sample")

# Plot the sum of KO expression over time
ggplot(phaeo_KO_sum, aes(x = Date, y = TPL)) +
  geom_line() +
  ggtitle("Sum of Phaeocystis photosynthesis KOs expression over time")

# Disentangle the tidal influence on O2 and the contribution of Phaeocystis PP (estimated by gene expression of photosynthesis KOs)
## Test: Simulate a test dataset of two overlayed sine functions
n <- 100
x <- seq(0, 4 * pi, length.out = n)
set.seed(1)
noise <- rnorm(n, mean = 0, sd = 0.3)
y1 <- 2*sin(x) + noise
y2 <- sin(x / 2) + noise
y_avg <- (y1 + y2) / 2
plot(x, y_avg, col = "green", lwd = 2, type="l",
     xlab = "X", ylab = "Y", main = "Sine Waves with Variation and Their Average")

#Define a sine fitting function
fit_function <- function(x, A1, A2, phi1, phi2) {
  (A1 * sin(x + phi1) + A2 * sin(x / 2 + phi2)) / 2
}

#Fit the model
fit <- nls(y_avg ~ fit_function(x, A1, A2, phi1, phi2),
           start = list(A1 = 1, A2 = 1, phi1 = 0, phi2 = 0))

#Reconstruct the original functions
A1_est <- coef(fit)["A1"]
A2_est <- coef(fit)["A2"]
phi1_est <- coef(fit)["phi1"]
phi2_est <- coef(fit)["phi2"]
y1_est <- A1_est * sin(x + phi1_est)
y2_est <- A2_est * sin(x / 1.89 + phi2_est)

#Plot the reconstructed functions
plot(x, y_avg, type = "l", col = "black", lwd = 2)
lines(x, y1_est, col = "blue", lwd = 2, lty = 2)
lines(x, y2_est, col = "red", lwd = 2, lty = 2)

## Now with biological data
