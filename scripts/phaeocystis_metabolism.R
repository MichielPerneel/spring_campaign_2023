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

genus_name <- 'Phaeocystis'

# Define file paths
genus_data_file <- "data/analysis/phaeocystis_bin_tpm.csv"
env_data_file <- "data/samples_env.csv"

# Read the data
genus_data <- read.csv(genus_data_file, row.names = 1)
env_data <- read.csv(env_data_file, row.names = 1, sep = ";")

# Create directory to save plots
plot_dir <- file.path("figures/metatranscriptomics/phaeocystis_rcm")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

############################        TO DO        ###############################
# - Create a KEGG KO ID x samples matrix containing count sums for transcripts annotated to Phaeocystis
# - Correctly read env and metabolic data
# - proceed from here
################################################################################

# Print the dimensions of the environmental data
cat("Environmental data dimensions:", dim(env_data), "\n")
# Print the dimensions of the genus data
cat("Genus data dimensions:", dim(genus_data), "\n")

# Remove rows with NaN values in the environmental data
  env_data <- env_data[complete.cases(env_data), ]

  # Preprocess the metabolic data
  zero_cols <- colSums(genus_data == 0) == nrow(genus_data)
  cat("Columns with all zeros in the genus data:", sum(zero_cols), "\n")
  genus_data <- genus_data[, !zero_cols]

  y_genus <- genus_data

  # Print the dimensions of the environmental data after removing NaN values
  cat("Environmental data dimensions after removing NaN values:", dim(env_data), "\n")
  cat("Genus data dimensions after preprocessings:", dim(y_genus), "\n")

  # Remove unoverlapping samples
  y_genus <- y_genus[rownames(y_genus) %in% rownames(env_data), ]
  env_data <- env_data[rownames(env_data) %in% rownames(y_genus), ]
  cat("Samples retained in both dataframes:", nrow(y_genus), "\n")
  cat("Samples retained:", rownames(y_genus), "\n")

  # Convert to phyloseq object
  OTU <- otu_table(as.matrix(y_genus), taxa_are_rows = FALSE)
  ENV <- sample_data(env_data[, -c(1)])
  physeq <- phyloseq(OTU, ENV)

  # Unconstrained analysis, conditioned on sequencing provider
  GenusRCM <- RCM(physeq, k = 2, round = TRUE, confounders = c("sequencing_provider"))

  # Constrained analysis, conditioned on sequencing provider
  GenusRCMconstr <- RCM(physeq, k = 2, round = TRUE, covariates = c("salinity", "nitrate", "nitrite", "phosphate", "silicate"), responseFun = "linear", confounders = c("sequencing_provider", "day_length", "temperature"))
  GenusRCMconstrNonParam <- RCM(physeq, k = 2, round = TRUE, covariates = c("salinity", "nitrate", "nitrite", "phosphate", "silicate"), responseFun = "nonparametric", confounders = c("sequencing_provider", "day_length", "temperature"))

  # Save plots
  plot_path <- function(name) file.path(plot_dir, paste0(name, ".png"))

  png(plot_path("GenusRCM_species_samples"))
  print(plot(GenusRCM, plotType = c("species", "samples"), samShape = 'station'))
  dev.off()

  png(plot_path("GenusRCMconstr_species_samples"))
  print(plot(GenusRCMconstr, plotType = c("species", "samples"), samShape = 'station'))
  dev.off()

  png(plot_path("GenusRCMconstr_variables"))
  print(plot(GenusRCMconstr, plotType = "variables"))
  dev.off()

  png(plot_path("GenusRCMconstr_variables_samples"))
  print(plot(GenusRCMconstr, plotType = c("variables", "samples"), inflVar = "psi", samShape = 'station'))
  dev.off()

  for (var in c("day_length", "temperature", "salinity", "nitrate", "nitrite", "silicate", "phosphate", "Shannon", "Deviance")) {
    png(plot_path(paste0("GenusRCM_species_samples_", var)))
    print(plot(GenusRCM, plotType = c("species", "samples"), samColour = var, samShape = 'station'))
    dev.off()
  }

  for (var in c("day_length", "temperature", "salinity", "nitrate", "nitrite", "silicate", "phosphate", "Shannon", "Deviance")) {
    png(plot_path(paste0("GenusRCMconstr_species_samples_", var)))
    print(plot(GenusRCMconstr, plotType = c("species", "samples"), samColour = var, samShape = 'station'))
    dev.off()
  }

  png(plot_path("GenusRCMconstrNonParam_species_samples"))
  print(plot(GenusRCMconstrNonParam, plotType = c("species", "samples"), samShape = 'station'))
  dev.off()

  png(plot_path("GenusRCMconstrNonParam_species_variables"))
  print(plot(GenusRCMconstrNonParam, plotType = "variables"))
  dev.off()

  png(plot_path("GenusRCMconstr_triplot"))
  print(plot(GenusRCMconstr, samShape = 'station'))
  dev.off()

  png(plot_path("GenusRCMconstrNonParam_response_function"))
  print(plotRespFun(GenusRCMconstrNonParam, taxa = NULL, subdivisions = 50L, Palette = "Set1", angle = 90, yLocSam = -20, axisTitleSize = 16, axisLabSize = 11, legendTitleSize = 18, legendLabSize = 12, samShape = "Diagnosis", labSize = 5))
  dev.off()

  # Below code does not work yet, I have colinear vectors creating an error in the solver function
  #for (var in c("nitrate", "silicate", "phosphate")) {
  #  png(plot_path(paste0("GenusRCMconstr_samples_inflVar_", var)))
  #  print(plot(GenusRCMconstr, plotType = c("variables", "samples"), inflVar = var))
  #  dev.off()
  #}

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

  # Close the log file
  sink()

  cat("Finished processing genus:", genus_name, "\n")
}

# Plot environmental data
plot_dir <- file.path("figures/functional_analysis/genus_counts_rcm_2/environmental")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

## Create marginal distribution plots for the environmental data
env_data_long <- env_data %>%
  select(c("day_length", "salinity", "temperature", "nitrate", "nitrite", "phosphate", "silicate")) %>%
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

# Ensure that rownames of env_data are dates or can be converted to dates
# Assuming that rownames are dates, otherwise convert them accordingly
env_data$date <- as.Date(env_data$date)

# Check if station column exists in the env_data
if ("station" %in% colnames(env_data)) {
  stations <- unique(env_data$station)

  for (station in stations) {
    # Filter data for the current station
    env_data_station <- env_data %>% filter(station == !!station)

    # Create a time series plot for the current station
    png(file.path(plot_dir, sprintf("day_length_temperature_lag_station_%s.png", station)))
    print(ggplot(env_data_station, aes(x = Date)) +
      geom_line(aes(y = day_length, color = "Day Length"), size = 1) +
      geom_line(aes(y = temperature, color = "Temperature"), size = 1, linetype = "dashed") +
      scale_color_manual(values = c("Day Length" = "blue", "Temperature" = "red")) +
      theme_minimal() +
      labs(title = paste("Day Length and Temperature Over Time - Station", station),
           x = "Date",
           y = "Value",
           color = "Variable") +
      theme(legend.position = "bottom"))
    dev.off()
  }
} else {
  cat("Station column not found in environmental data.")
}
