# WGCNA pipeline adapted from Natalie Cohen
# (https://github.com/cnatalie/METZYME/blob/master/WGCNA.R)

library(WGCNA)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(flashClust)
library(DESeq2)
library(clusterProfiler)
library(svglite)
library(reshape2)
library(tidyr)
library(forcats)
library(lubridate)

options(stringsAsFactors = FALSE)
allowWGCNAThreads()
enableWGCNAThreads()

# Define the genus we'll be working with
genus <- "Phaeocystis"

# Read in transcript counts
samples_genes_matrix <- read.csv(
  paste0(
    "data/annotation/taxonomy_eukprot/130/genus_bins/",
    genus, "_transcript_expression_sum.csv"),
    check.names = FALSE
)

# Make directories for the results
dir.create(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus), recursive = TRUE)
dir.create(paste0("data/analysis/WGCNA_130/transcripts/", genus), recursive = TRUE)
dir.create(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus), recursive = TRUE)
dir.create(paste0("data/analysis/WGCNA_130/transcripts/", genus), recursive = TRUE)

# Set row names and column names
rownames(samples_genes_matrix) <- samples_genes_matrix[, 1]
samples_genes_matrix <- samples_genes_matrix[, -1]

# Set values below 1 to 0 (lower detection limit)
samples_genes_matrix[samples_genes_matrix < 1] <- 0
cat("Dimensions of Gene matrix before processing:", dim(samples_genes_matrix), "\n")
# (optional: Remove transcripts that are expressed in one sample only or have sums across samples below 1
samples_genes_matrix <- samples_genes_matrix[rowSums(samples_genes_matrix > 0) >= 2 & rowSums(samples_genes_matrix) >= 10, ]
cat("Dimensions of Gene matrix after  processing (removing transcripts expressed in one sample only or sums across samples <10):", dim(samples_genes_matrix), "\n")

# Transform to datExpr matrix
datExpr <- as.matrix(t(samples_genes_matrix))
cat("Dimensions of datExpr matrix:", dim(datExpr), "\n")

# Check for outliers in the data
gsg <- goodSamplesGenes(datExpr, verbose = 3)
cat("All genes OK:", gsg$allOK, "\n")

# Filter out constant genes
constant_genes <- apply(datExpr, 2, function(x) var(x) == 0)
datExpr <- datExpr[, !constant_genes]
gsg <- goodSamplesGenes(datExpr, verbose = 3)
cat("All genes OK after filtering constant genes:", gsg$allOK, "\n")

# Pick a random subset of genes (for development purposes)
#datExpr <- datExpr[, sample(ncol(datExpr), 5000)]

# log2 transform, add pseudocount
datExpr <- log2(datExpr + 1)
cat("Dimensions of datExpr matrix after log2 transformation:", dim(datExpr), "\n")

# Read in environmental data
data_env <- read.csv("data/samples_env.csv", stringsAsFactors = FALSE, sep = ",", header = TRUE, row.names = 1)
# Read primary production data (PP)
pp_data <- read.csv("data/raw/LabSTAF/labstaf_combined_data.csv")
# Combine Station and Sample columns in PP_data into one column
pp_data$Station <- paste(pp_data$Station, pp_data$Sample, sep = "_")
# Read zooplankton data
zooplankton_data <- read.csv("data/analysis/zooplankton_counts.csv")

# Merge PP data with environmental data based on the Station column
data_env <- data_env %>% rownames_to_column("Station")
data_env <- data_env %>%
  inner_join(pp_data %>% select(Station, PP), by = "Station") %>%
  inner_join(zooplankton_data, by = "Station")
# Restore row names from Station column after merging
data_env <- data_env %>% column_to_rownames("Station")

# Remove uninteresting columns
columns_to_remove <- c(
  'StationPrefix', 'StationSuffix', 'Latitude', 'Longitude', 'NOX',
  'Conductivity', 'Depth',
  'sea_surface_height_above_sea_level', 'surface_baroclinic_sea_water_velocity'
  )

data_env <- data_env[, !names(data_env) %in% columns_to_remove]

# Only keep the rows that are present in both datasets
data_env <- data_env[rownames(datExpr), ]

# Order rows of datExpr to match the order of the environmental data
datExpr <- datExpr[rownames(data_env), ]

# Ensure the datasets align correctly
stopifnot(all(rownames(data_env) == rownames(datExpr)))

# Get the metadata from data_env
metadata <- data_env[, c("Date", "day_moment")]

# Get the environmental parameters from data_env
env_params <- data_env[, c(
  "day_length", "Temperature", "Salinity",
  "Oxygen", "Fluorescence", "NH4", "NO2", "NO3",
  "PO4", "Si", "TEP", "PP", "Total_Zooplankton_Count"
)]

#------------- 1. Dendrogram and trait heatmap showing outliers ------------#
cat("Calculating sample network adjacency matrix...\n")
A <- adjacency(t(datExpr), type = "signed")
k <- as.numeric(apply(A, 2, sum)) - 1
Z.k <- scale(k)
thresholdZ.k <- -2.5
outlierColor <- ifelse(Z.k < thresholdZ.k, "red", "black")
sampleTree <- flashClust(as.dist(1 - A), method = "average")
traitColors <- data.frame(numbers2colors(env_params, signed = FALSE))
dimnames(traitColors)[[2]] <- paste(names(env_params))
datColors <- data.frame(outlierC = outlierColor, traitColors)

cat("Plotting sample dendrogram and trait heatmap...\n")
svg(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus, "/sample_dendrogram_and_trait_heatmap.svg"))
plotDendroAndColors(sampleTree,
                    groupLabels = names(datColors),
                    colors = datColors,
                    main = "Sample dendrogram and trait heatmap")
dev.off()

# Remove outlying samples from the expression matrix and trait matrix
samples_to_remove <- Z.k < thresholdZ.k | is.na(Z.k)
datExpr <- datExpr[!samples_to_remove, ]
env_params <- env_params[!samples_to_remove, ]
metadata <- metadata[!samples_to_remove, ]

#------------- 2. Network construction and module detection  ------------#
cat("Choosing a set of soft-thresholding powers...\n")
powers <- c(seq(1, 20, by = 1), seq(20, 60, by = 5))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 2)

# Convert to a dataframe for plotting
sft_data <- as.data.frame(sft$fitIndices)
sft_data$Power <- sft_data[, 1]
sft_data$ScaledR2 <- -sign(sft_data[, 3]) * sft_data[, 2]
sft_data$MeanConnectivity <- sft_data[, 5]
sft_data$Labels <- powers

# Plot 1: Scale-free topology fit index
cat("Plotting Scale-free topology fit index...\n")
plot1 <- ggplot(sft_data, aes(x = Power, y = ScaledR2, label = Labels)) +
  geom_point() +
  geom_text(colour = "grey", nudge_y = 0.05) +
  geom_hline(yintercept = c(0.80, 0.90), colour = "red", linetype = "dashed") +
  labs(x = "Soft Threshold (power)", y = "Scale Free Topology Model Fit, signed R^2", title = "Scale independence") +
  theme_minimal()
ggsave(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus, "/scale_independence.svg"), plot = plot1, width = 10, height = 5, dpi = 600)

# Plot 2: Mean connectivity
cat("Plotting Mean connectivity...\n")
plot2 <- ggplot(sft_data, aes(x = Power, y = MeanConnectivity, label = Labels)) +
  geom_point() +
  geom_text(colour = "grey", nudge_y = 0.10) +
  labs(x = "Soft Threshold (power)", y = "Mean Connectivity", title = "Mean connectivity") +
  theme_minimal()
ggsave(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus, "/mean_connectivity.svg"), plot = plot2, width = 10, height = 5, dpi = 600)

# Choose soft thresholding power based on the plot above
softPower <- sft$powerEstimate
cat("Chosen soft threshold power:", softPower, "\n")

# Calculate the adjacency matrix
adjacency <- adjacency(datExpr, power = softPower, type = "signed")

# Translate the adjacency into topological overlap matrix, calculate the corresponding dissimilarity
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

# Each branch corresponds to a transcript id, branches grouping together densely are co-expressed KOs
geneTree <- flashClust(as.dist(dissTOM), method = "average")
minModuleSize <- 70
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 4, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)

cat("Plotting Gene dendrogram and module colors...\n")
svg(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus, "/gene_dendrogram_and_module_colors.svg"))
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# Merge modules whose expression profiles are very similar
MEList <- moduleEigengenes(datExpr, colors = dynamicColors, softPower = softPower)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- flashClust(as.dist(MEDiss), method = "average")

cat("Plotting Clustering of module eigengenes...\n")
svg(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus, "/clustering_of_module_eigengenes.svg"))
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h = 0.4, col = "red")
dev.off()

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.4, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

cat("Plotting Merged dynamic colors...\n")
svg(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus, "/merged_dynamic_colors.svg"))
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors <- mergedColors
colorOrder <- c("grey", standardColors(100))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs

#-------- 3. Relating modules to traits and finding important genes -------#
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
cat("Number of genes:", nGenes, "\n")
cat("Number of samples:", nSamples, "\n")

# Recalculate module eigengenes with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
if ("MEgrey" %in% colnames(MEs0)) {
  MEs0 <- subset(MEs0, select = -MEgrey)
}
MEs <- orderMEs(MEs0)

# Check correlation between module eigengenes and environmental parameters
moduleTraitCor <- cor(MEs, env_params, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

cat("Creating Module-trait relationships heatmap...\n")
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

svg(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus, "/module_trait_relationships.svg"), width = 9, height = 10)
par(mar = c(8, 12.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(env_params),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.6,
               zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

# Save the module eigengene expression data as csv
write.csv(MEs, file = paste0("data/analysis/WGCNA_130/transcripts/", genus, "/module_eigengenes.csv"), row.names = TRUE)

# Plot the expression of the eigengene of each module at every hour
cat("Plotting module eigengene expression per month and station...\n")
MElong <- MEs %>% as.data.frame() %>% rownames_to_column("sample")
MElong <- melt(MElong, id.vars = "sample", variable.name = "module", value.name = "ME_expression")
metadata2 <- metadata %>% rownames_to_column("sample")
MElong <- MElong %>% left_join(metadata2, by = "sample")
MElong$Date <- ymd_hms(MElong$Date)  # Parse both date and time
MElong$hour <- round_date(MElong$Date, "hour")
MElong$hour <- as.factor(format(MElong$hour, "%Y-%m-%d %H:%M:%S"))

date_range <- seq(from = ymd_hms("2023-04-20 09:00:00"), to = ymd_hms("2023-04-21 10:00:00"), by = "1 hour")
date_labels <- format(date_range, "%Y-%m-%d %H:%M:%S")
MElong <- MElong %>% complete(module, hour = date_labels) %>% arrange(module, hour)

p <- ggplot(MElong, aes(x = hour, y = module, fill = ME_expression)) +
  geom_tile() +
  scale_fill_distiller(palette = "Spectral", na.value = "white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = "bottom") +
  labs(x = "Hour", y = "Module", fill = "Module Eigengene Expression") +
  scale_x_discrete(limits = date_labels) +
  guides(fill = guide_colorbar(barwidth = 10, barheight = 1, title.position = "top", title.hjust = 0.5))

ggsave(paste0("figures/metatranscriptomics/WGCNA_130/transcripts/", genus, "/module_eigengene_expression_per_hour.svg"), plot = p, width = 10, height = 12)

#----------------------- 4. Module content ----------------------#
cat("Finding important genes in each module...\n")
Transcript_ModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))

# Rank modules' most correlated modules by decreasing order of correlation
# Construct a dataframe that contains the module name, IDs, and membership
for (module in names(MEs)) {
  cat("Processing module:", module, "\n")
  module_cor <- select(Transcript_ModuleMembership, module)
  module_cor <- module_cor[order(-module_cor[, 1]), , drop = FALSE]
  module_cor <- module_cor %>% as.data.frame() %>% rownames_to_column("transcript_id")
  cat("Number of transcripts in module", module, ":", nrow(module_cor), "\n")
  write.table(module_cor, file = paste0("data/analysis/WGCNA_130/transcripts/", genus, "/", module, "_content.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
}

cat("Saving all transcript IDs...\n")
write.table(rownames(Transcript_ModuleMembership), file = paste0("data/analysis/WGCNA_130/transcripts/", genus, "/transcript_id_list.txt"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)