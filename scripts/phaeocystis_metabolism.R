library(RCM)
library(corrplot)
library(gclus)
library(argparse)
library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(phyloseq)
library(ggplot2)
library(oce)

genus_name <- 'Phaeocystis'

# Define file paths
genus_data_file <- "data/analysis/phaeocystis_metabolic_functions_counts.csv"
env_data_file <- "data/samples_env.csv"
pp_data_file <- "data/raw/LabSTAF/labstaf_combined_data.csv"
zooplankton_data_file <- 'data/analysis/zooplankton_counts.csv'

# Read the data
genus_data <- read.csv(genus_data_file, row.names = 1, check.names=FALSE)
env_data <- read.csv(env_data_file, row.names = 1)
# Delete columns Oxygen and Fluorescence, as they contain a lot of NaN values
env_data <- env_data[, -c(12, 13)]
# Only retain samples in station 130
env_data <- env_data[env_data$StationPrefix == 130, ]
# Read primary production data (PP)
pp_data <- read.csv(pp_data_file)
# Combine Station and Sample columns in PP_data into one column
pp_data$Station <- paste(pp_data$Station, pp_data$Sample, sep = "_")
# Read zooplankton data
zooplankton_data <- read.csv(zooplankton_data_file)

# Merge PP data with environmental data based on the Station column
env_data <- env_data %>% rownames_to_column("Station")
env_data <- env_data %>%
  inner_join(pp_data %>% select(Station, PP), by = "Station") %>%
  inner_join(zooplankton_data, by = "Station")
# Restore row names from Station column after merging
env_data <- env_data %>% column_to_rownames("Station")

# Create directory to save plots
plot_dir <- file.path("figures/metatranscriptomics/phaeocystis_rcm")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Print the dimensions of the environmental data
cat("Environmental data dimensions:", dim(env_data), "\n")
# Print the dimensions of the genus data
cat("Genus data dimensions:", dim(genus_data), "\n")

# Remove rows with NaN values in the environmental data
#env_data <- env_data[complete.cases(env_data), ]

# Preprocess the metabolic data
zero_cols <- colSums(genus_data == 0) == nrow(genus_data)
cat("Columns with all zeros in the genus data:", sum(zero_cols), "\n")
genus_data <- genus_data[, !zero_cols]

# Print the dimensions of the environmental data after removing NaN values
cat("Environmental data dimensions after removing NaN values:", dim(env_data), "\n")
cat("Genus data dimensions after preprocessings:", dim(genus_data), "\n")

# Remove unoverlapping samples
genus_data <- genus_data[rownames(genus_data) %in% rownames(env_data), ]
env_data <- env_data[rownames(env_data) %in% rownames(genus_data), ]
cat("Samples retained in both dataframes:", nrow(genus_data), "\n")
cat("Samples retained:", rownames(genus_data), "\n")

# Extract time from the date column
env_data$Date <- as.POSIXct(env_data$Date, format = "%Y-%m-%d %H:%M:%S")
env_data$Time <- format(env_data$Date, "%H:%M:%S")
# Convert the Time column to minutes since midnight
env_data$Time <- as.numeric(difftime(as.POSIXct(env_data$Time, format = "%H:%M:%S"),
                                     as.POSIXct("00:00:00", format = "%H:%M:%S"),
                                     units = "mins"))
# Encode the day_moment column as a factor and dummy code it
## Group Astronomical, Civil, and Nautical twilight as twilight
env_data$day_moment <- ifelse(env_data$day_moment %in% c("Astronomical twilight", "Civil twilight", "Nautical twilight"), "Twilight", env_data$day_moment)
env_data$day_moment <- factor(env_data$day_moment)
# Extract azimuth and altitude from the sunAngle function
sun_positions <- sunAngle(env_data$Date, env_data$Latitude, env_data$Longitude)
# Add azimuth and altitude as new columns to your environmental data
env_data$SunAzimuth <- sun_positions$azimuth
env_data$SunAltitude <- sun_positions$altitude
# Create a new binary variable for day (1) and night (0) based on sun altitude
env_data$Sun <- ifelse(env_data$SunAltitude > 0, 1, 0)
# Set altitude of the sun to 0 when it is below the horizon
env_data$SunAltitude <- ifelse(env_data$SunAltitude < 0, 0, env_data$SunAltitude)

## Write the sun data to a new file
sun_export <- env_data %>%
  # Row names to column
  rownames_to_column("Station") %>%
  select(Station, Date, day_moment, SunAzimuth, SunAltitude, Sun)
write.csv(sun_export, "data/analysis/sun_data.csv")

# Remove uninformative columns
env_data <- env_data %>%
  select(-surface_baroclinic_sea_water_velocity, -sea_surface_height_above_sea_level,
         -StationPrefix, -StationSuffix, -Latitude, -Longitude,
         -Date, -day_length, -Conductivity, -NOX)

# Convert to phyloseq object
OTU <- otu_table(as.matrix(genus_data), taxa_are_rows = FALSE)
ENV <- sample_data(env_data)
physeq <- phyloseq(OTU, ENV)

# Unconstrained analysis, conditioned on sequencing provider
GenusRCM <- RCM(physeq, k = 2, round = TRUE)

# Constrained analysis, conditioned on sequencing provider
## Define covariates and optionally confounders
covariates <- c("Salinity", "Temperature", "Depth", "TEP", "PP",  "Time",
                "SunAzimuth", "SunAltitude", "NH4", "NO3", "NO2", "PO4", "Si", "Total_Zooplankton_Count")

## Day_moment does not really work
covariates <- c("Salinity", "Temperature", "Depth", "TEP", "PP",
                "SunAltitude", "NH4", "NO3", "NO2", "PO4", "Si", "Total_Zooplankton_Count")

GenusRCMconstr <- RCM(physeq,
                      k = 2, round = TRUE, responseFun = "linear",
                      covariates = covariates)

GenusRCMconstrNonParam <- RCM(physeq,
                              k = 2, round = TRUE, responseFun = "nonparametric",
                              covariates = covariates)

# Save plots
plot_path <- function(name) file.path(plot_dir, paste0(name, ".png"))
plot_path_svg <- function(name) file.path(plot_dir, paste0(name, ".svg"))

png(plot_path("GenusRCM_species_samples"))
print(plot(GenusRCM, plotType = c("species", "samples")))
dev.off()

png(plot_path("GenusRCMconstr_species_samples"))
print(plot(GenusRCMconstr, plotType = c("species", "samples")))
dev.off()

png(plot_path("GenusRCMconstr_variables"))
print(plot(GenusRCMconstr, plotType = "variables"))
dev.off()

svg(plot_path_svg("GenusRCMconstr_variables"))
print(plot(GenusRCMconstr, plotType = "variables"))
dev.off()

png(plot_path("GenusRCMconstr_variables_samples"))
print(plot(GenusRCMconstr, plotType = c("variables", "samples"), inflVar = "psi"))
dev.off()

for (var in c("Temperature", "Salinity", "Depth", "TEP", "PP", "NO3", "NO2", "PO4", "Si", "Total_Zooplankton_Count", "Time", "day_moment", "SunAzimuth", "SunAltitude", "Sun", "Shannon", "Deviance")) {
  png(plot_path(paste0("GenusRCM_species_samples_", var)))
  print(plot(GenusRCM, plotType = c("species", "samples"), samColour = var))
  dev.off()
}

for (var in c("Temperature", "Salinity", "Depth", "TEP", "PP", "NO3", "NO2", "PO4", "Si", "Total_Zooplankton_Count",  "Time", "day_moment", "SunAzimuth", "SunAltitude", "Sun", "Shannon", "Deviance")) {
  png(plot_path(paste0("GenusRCMconstr_species_samples_", var)))
  print(plot(GenusRCMconstr, plotType = c("species", "samples"), samColour = var))
  dev.off()
}

png(plot_path("GenusRCMconstrNonParam_species_samples"))
print(plot(GenusRCMconstrNonParam, plotType = c("species", "samples")))
dev.off()

png(plot_path("GenusRCMconstrNonParam_species_variables"))
print(plot(GenusRCMconstrNonParam, plotType = "variables"))
dev.off()

png(plot_path("GenusRCMconstr_triplot"))
print(plot(GenusRCMconstr))
dev.off()

png(plot_path("GenusRCMconstrNonParam_response_function"))
print(plotRespFun(GenusRCMconstrNonParam, taxa = NULL, subdivisions = 50L, Palette = "Set1", angle = 90, yLocSam = -20, axisTitleSize = 16, axisLabSize = 11, legendTitleSize = 18, legendLabSize = 12, samShape = "Diagnosis", labSize = 5))
dev.off()

for (type in c("response", "runs")) {
  png(plot_path(paste0("residualPlot_", type)))
  print(residualPlot(GenusRCMconstr, whichTaxa = type, numTaxa = 9))
  dev.off()

  png(plot_path(paste0("residualPlot_", type, "_Pearson")))
  print(residualPlot(GenusRCMconstr, whichTaxa = type, resid = "Pearson", numTaxa = 9))
  dev.off()

  png(plot_path(paste0("residualPlot_", type, "_Deviance")))
  print(residualPlot(GenusRCMconstr, whichTaxa = type, resid = "Deviance", numTaxa = 9))
  dev.off()
}

# Extract most important modules contributing to the separation of the samples
genusCoords = extractCoord(GenusRCM)
moduleSignals = rowSums(genusCoords$species[, c("end1", "end2")]^2)
sortedModules = taxa_names(GenusRCM$physeq)[order(moduleSignals, decreasing = TRUE)]
print(sortedModules[1:50])

# Plot environmental data
plot_dir <- file.path("figures/functional_analysis/genus_counts_rcm_2/environmental")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

## Create marginal distribution plots for the environmental data
env_data_long <- env_data %>%
  select(c("Temperature", "Salinity", "Depth", "TEP", "NO3", "NO2", "PO4", "Si", "Total_Zooplankton_Count", "Time", "day_moment", "SunAzimuth", "SunAltitude", "Sun", "PP")) %>%
  gather(key = "Variable", value = "Value")

# Histogram plot
png(file.path(plot_dir, "env_data_marginal_histogram.png"))
print(ggplot(env_data_long, aes(x = Value)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_minimal() +
  labs(title = "Marginal Distribution of Environmental Variables"))
dev.off()

# Density plot
png(file.path(plot_dir, "env_data_marginal_density.png"))
print(ggplot(env_data_long, aes(x = Value)) +
  geom_density(fill = "green", alpha = 0.5) +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_minimal() +
  labs(title = "Marginal Distribution (Density) of Environmental Variables"))
dev.off()
