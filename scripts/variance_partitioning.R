library(vegan)
library(plyr)
library(dplyr)
library(funrar)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(reshape2)
library(mgcv)
library(lme4)
library(car)

#----- GAM of O2 vs P. globosa photosynthesis and N. scintillans biomass ------#
# Load environmental data
env <- read.csv("data/samples_env_carbon.csv", header = TRUE, row.names = 1)
sun_data <- read.csv("data/analysis/sun_data.csv", header = TRUE, row.names = 1) # This was generated in the script "phaeocystis_metabolism.R"

# Collect the tidal samples
tidal <- env %>%
    filter(StationPrefix == "130") %>%
    select(sea_surface_height_above_sea_level, surface_baroclinic_sea_water_velocity, Depth) %>%
    mutate(across(everything(), as.numeric))

# Collect the diel samples
diel <- sun_data %>%
    select(Station, Date, SunAzimuth, SunAltitude) %>%
    mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")) %>%
    mutate(Time_since_start = as.numeric(difftime(Date, min(Date), units = "secs"))) %>%
    select(-Date) %>%
    tibble::column_to_rownames("Station") %>%
    mutate(across(everything(), as.numeric))

# Collect the environmental samples
env_subset <- env %>%
  select(StationPrefix, Temperature, Salinity, NH4, NOX,PO4, Si) %>%
  filter(StationPrefix == "130") %>%
  select(-StationPrefix) %>%
  mutate(across(everything(), as.numeric))

# Last column of env_subset contains NAs. Remove sample 130_25 from all datasets
env_subset <- env_subset[-which(rownames(env_subset) == "130_25"),]
tidal <- tidal[-which(rownames(tidal) == "130_25"),]
diel <- diel[-which(rownames(diel) == "130_25"),]

oxygen_smoother <- read.csv("data/analysis/O2_smoother_station_130.csv", header = TRUE)
# Preprocess the oxygen smoother data
oxygen_smoother <- oxygen_smoother %>%
  rename(Date = Date.Time, O2 = y, O2_min = ymin, O2_max = ymax, O2_se = se) %>%
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

# Preprocess the carbon response variable
carbon <- env %>%
  select(Date, StationPrefix, pCO2, DIC) %>%
  filter(StationPrefix == "130") %>%
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")) %>%
  select(-StationPrefix) %>%
  .[-which(rownames(.) == "130_25"),]

# Interpolate oxygen values for each sample time in carbon
oxygen_interpolated <- approx(x = oxygen_smoother$Date,
                              y = oxygen_smoother$O2,
                              xout = carbon$Date,
                              method = "linear",
                              rule = 2)
## Transform the interpolated values into a dataframe
oxygen_interpolated <- data.frame(Date = carbon$Date, O2 = oxygen_interpolated$y)

##-------------------------- Phaeocystis Metabolism --------------------------##
phaeo_KO <- read.csv("data/analysis/130_phaeo_KO_TPL_standardized_expression.csv", header = TRUE, row.names = 1, check.names = FALSE)
## Read list of photosynthesis KOs
photosynthesis_KOs <- read.table("data/annotation/functional_eggnog/photosynthesis_KOs.txt")

## Remove all zero columns and rows
phaeo_KO <- phaeo_KO[, colSums(phaeo_KO) != 0]
phaeo_KO <- phaeo_KO[rowSums(phaeo_KO) != 0,]

## Only retain the photosynthesis related KOs
phaeo_KO <- phaeo_KO[rownames(phaeo_KO) %in% photosynthesis_KOs$V1,]

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

# Plot the sum of all photosynthesis KOs expression over time
ggplot(phaeo_KO_plot, aes(x = Date, y = TPL)) +
  geom_point() +
  geom_smooth(method = "loess") +
  ggtitle("Sum of Phaeocystis photosynthesis KOs expression over time")

##--------------------- Zooplankton Community Composition --------------------##

# Load ZooScan data as a zooplankton community response variable
zooscan <- read.csv("data/raw/zooscan_data.csv", header = TRUE, sep = ";", row.names = 1)
zooscan <- zooscan %>%
    filter(grepl("130", Station)) %>%
    # Station should be rownames instead of the dates now
    tibble::rownames_to_column("Date") %>%
    select(-Date, -Month, -Year, -Day, -Hour, -Minute) %>%
    tibble::column_to_rownames("Station") %>%
    na.omit() %>%
    .[-which(rownames(.) == "130_25"),]

#-#------------------------------- GAM model ---------------------------------##
# Calculate first PC of Phaeocystis photosynthesis KO expression
phaeo_KO_pca <- prcomp(t(phaeo_KO), scale = FALSE)
phaeo_KO_pca <- phaeo_KO_pca$x[,1]

# Combine PC value and N. scintillans abundances into a single dataframe
O2_test <- data.frame(phaeo_KO_pca)
## Inner join with N. scintillans abundances from zooscan data and environmental data
O2_test <- O2_test %>%
    merge(zooscan, by = "row.names") %>%
    rename(sample = "Row.names") %>%
    merge(carbon, by = "sample") %>%
    merge(oxygen_interpolated, by = "Date") %>%
    merge(tidal, by = "sample") %>%
    merge(diel, by = "sample")

# Plot PC axes against date and time
## PCA of Phaeocystis KO expression
ggplot(O2_test, aes(x = Date, y = phaeo_KO_pca)) +
  geom_point() +
  geom_smooth(method = "loess") +
  ggtitle("PCA1 of Phaeocystis KO expression over time")

## PCA of zooscan data
ggplot(O2_test, aes(x = Date, y = Noctiluca)) +
  geom_point() +
  geom_smooth(method = "loess") +
  ggtitle("PCA1 of zooscan data over time")

vif_model <- lm(O2 ~ phaeo_KO_pca + Noctiluca, data = O2_test)
vif_values <- vif(vif_model)
print(vif_values)

# Fit a GAM model
gam_O2 <- gam(O2 ~ s(phaeo_KO_pca, bs = "cs") +
                s(Noctiluca, bs = "cs"),
                #s(Time_since_start, bs = "cc"),
              data = O2_test, method = "REML")

summary(gam_O2)
plot(gam_O2, pages = 1)
