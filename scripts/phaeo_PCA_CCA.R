library(ggplot2)
library(vegan)
library(dplyr)

# Set file paths
tpl_file <- "/data/gent/vo/001/gvo00125/vsc43619/spring_campaign_2023/data/quantification/130/130_tpl.csv"
annotation_file <- "/data/gent/vo/001/gvo00125/vsc43619/spring_campaign_2023/data/annotation/taxonomy_eukprot/130/genus_bins/Phaeocystis_transcriptome_bin.csv"
metadata_file <- "/data/gent/vo/001/gvo00125/vsc43619/spring_campaign_2023/data/samples_env.csv"
output_file <- "/data/gent/vo/001/gvo00125/vsc43619/spring_campaign_2023/figures/metatranscriptomics/pca_cca.pdf"

# Load data
cat("Loading data...\n")
tpl_data <- read.csv(tpl_file, row.names = 1)
annotations <- read.csv(annotation_file)
metadata <- read.csv(metadata_file)

# Filter TPL data for Phaeocystis transcripts
cat("Filtering TPL data for Phaeocystis transcripts...\n")
phaeo_transcripts <- annotations$query_id
tpl_data <- tpl_data[rownames(tpl_data) %in% phaeo_transcripts, ]

# Normalize TPL data by total Phaeocystis activity
cat("Normalizing TPL data by library size...\n")
phaeo_sums <- colSums(tpl_data)
tpl_data_normalized <- sweep(tpl_data, 2, phaeo_sums, FUN = "/") * 1e6

# Match metadata with expression data
cat("Matching metadata with normalized TPL data...\n")
metadata <- metadata[metadata$SampleID %in% colnames(tpl_data_normalized), ]
tpl_data_normalized <- tpl_data_normalized[, metadata$SampleID]

# PCA (Panel A)
cat("Performing PCA...\n")
pca_result <- prcomp(t(tpl_data_normalized), scale. = TRUE)  # Transpose: samples as rows, genes as columns
pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleID <- rownames(pca_scores)

# Merge PCA scores with metadata for plotting
pca_plot_data <- merge(pca_scores, metadata, by = "SampleID")

# Plot PCA
pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = Station)) +
    geom_point(size = 3) +
    labs(
        title = "PCA of Phaeocystis Transcript Expression",
        x = paste0("PC1 (", round(100 * pca_result$sdev[1]^2 / sum(pca_result$sdev^2), 1), "% variance)"),
        y = paste0("PC2 (", round(100 * pca_result$sdev[2]^2 / sum(pca_result$sdev^2), 1), "% variance)")
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

# CCA (Panel B)
cat("Performing CCA...\n")
# Subset environmental variables
env_vars <- c("Temperature", "Salinity", "NH4", "NO2", "NO3", "PO4", "Si", "TEP", "Oxygen")
env_data <- metadata[, env_vars, drop = FALSE]
expr_cca <- as.matrix(tpl_data_normalized)

# Run CCA
cca_result <- cca(expr_cca ~ ., data = env_data)

# Extract CCA scores
cca_sites <- as.data.frame(scores(cca_result, display = "sites"))
cca_sites$SampleID <- rownames(cca_sites)
cca_env <- as.data.frame(scores(cca_result, display = "bp"))

# Merge CCA scores with metadata for plotting
cca_plot_data <- merge(cca_sites, metadata, by = "SampleID")

# Plot CCA
cca_plot <- ggplot() +
    geom_point(data = cca_plot_data, aes(x = CCA1, y = CCA2, color = Station), size = 3) +
    geom_segment(data = cca_env, aes(x = 0, y = 0, xend = CCA1, yend = CCA2), arrow = arrow(length = unit(0.2, "cm")), color = "red") +
    geom_text(data = cca_env, aes(x = CCA1, y = CCA2, label = rownames(cca_env)), hjust = 1.1, vjust = 1.1, color = "red") +
    labs(
        title = "Panel B: CCA of P. globosa Expression and Environmental Conditions",
        x = "CCA1",
        y = "CCA2"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

# Combine plots and save
cat("Saving supplementary figure...\n")
pdf(output_file, width = 12, height = 6)
gridExtra::grid.arrange(pca_plot, cca_plot, ncol = 2)
dev.off()

cat("Supplementary figure saved to", output_file, "\n")