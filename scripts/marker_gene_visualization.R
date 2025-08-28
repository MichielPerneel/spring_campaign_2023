# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(tibble)
library(ggfortify)
library(stats)
library(purrr)
library(factoextra)

# Load the data
DIC <- read.csv('data/analysis/DIC_smoother_station_130.csv')
O2 <- read.csv('data/analysis/O2_smoother_station_130.csv')
meta <- read.csv('data/samples_env.csv')

marker_genes <- read.csv('data/analysis/significant_derivative_markers.csv')
TPL <- read.csv('data/quantification/130/130_tpl.csv', check.names = FALSE)
TPM <- read.csv('data/quantification/130/130_tpm.csv', check.names = FALSE)
clusters <- read.csv('data/analysis/phaeocystis_O2_clusters.csv')
annotation_file <- read.csv('data/annotation/taxonomy_eukprot/130/genus_bins/Phaeocystis_transcriptome_bin.csv')
phaeo_total_tpl <- read.csv('/Users/michiel/gitlab/spring_campaign_2023/data/analysis/phaeocystis_bin_tpl_130.csv')
phaeo_total_tpm <- read.csv('/Users/michiel/gitlab/spring_campaign_2023/data/analysis/phaeocystis_bin_tpm_130.csv')

# Process the DIC and O2 data
DIC$Date <- as.POSIXct(DIC$Date)
O2$Date <- as.POSIXct(O2$Date)
meta$Date <- as.POSIXct(meta$Date)
meta <- meta %>% filter(StationPrefix == 130)

# Intrepolate O2 and DIC data to match metadata
meta <- meta %>%
  mutate(DIC = approx(x = DIC$Date, y = DIC$y, xout = Date, method = "linear", rule = 2)$y,
         O2 = approx(x = O2$Date, y = O2$y, xout = Date, method = "linear", rule = 2)$y)

# Determine scaling factors for O2
o2_min <- min(O2$y, na.rm = TRUE)
o2_max <- max(O2$y, na.rm = TRUE)
dic_min <- min(DIC$y, na.rm = TRUE)
dic_max <- max(DIC$y, na.rm = TRUE)

# Function to scale O2 values to align with DIC scale
scale_o2_to_dic <- function(o2_value) {
  ((o2_value - o2_min) / (o2_max - o2_min)) * (dic_max - dic_min) + dic_min
}

# Function to reverse scale for secondary axis
reverse_scale_dic_to_o2 <- function(dic_value) {
  ((dic_value - dic_min) / (dic_max - dic_min)) * (o2_max - o2_min) + o2_min
}

# Add scaled O2 values to O2 dataframe
O2 <- O2 %>%
  mutate(y_scaled = scale_o2_to_dic(y))

# Process the TPL data
marker_TPL <- TPL %>%
  pivot_longer(-target_id, names_to = "sample", values_to = "TPL") %>%
  filter(target_id %in% marker_genes$query_id) %>%
  left_join(meta, by = c("sample" = "Station")) %>%
  mutate(TPL_standardized = TPL / sum(TPL, na.rm = TRUE),
         TPL_standardized_z = scale(TPL_standardized))

marker_TPM <- TPM %>%
  pivot_longer(-target_id, names_to = "sample", values_to = "TPM") %>%
  filter(target_id %in% marker_genes$query_id) %>%
  left_join(meta, by = c("sample" = "Station")) %>%
  mutate(TPM_standardized = TPM / sum(TPM, na.rm = TRUE),
         TPM_standardized_z = scale(TPM_standardized))

# Separate positively and negatively correlated genes
marker_genes <- marker_genes %>%
  mutate(correlation = ifelse(dO2_resid_dt_coef > 0, "Positive", "Negative"))

# Highlight top 3 markers for each correlation type
top_positive_markers <- marker_genes %>%
  filter(correlation == "Positive") %>%
  arrange(desc(dO2_resid_dt_coef)) %>%
  slice(1:3) %>%
  pull(query_id)

top_negative_markers <- marker_genes %>%
  filter(correlation == "Negative") %>%
  arrange(dO2_resid_dt_coef) %>%
  slice(1:3) %>%
  pull(query_id)

# Assign colors to top markers
positive_colors <- c("darkgreen", "#0d629e", "#0d9e66")
negative_colors <- c("darkred", "firebrick", "tomato")

# Assign unique colors for the top 3 markers
marker_TPL <- marker_TPL %>%
  mutate(color = case_when(
    target_id %in% top_positive_markers ~ positive_colors[match(target_id, top_positive_markers)],
    target_id %in% top_negative_markers ~ negative_colors[match(target_id, top_negative_markers)],
    TRUE ~ "gray"
  ))

marker_TPM <- marker_TPM %>%
  mutate(color = case_when(
    target_id %in% top_positive_markers ~ positive_colors[match(target_id, top_positive_markers)],
    target_id %in% top_negative_markers ~ negative_colors[match(target_id, top_negative_markers)],
    TRUE ~ "gray"
  ))


# Merge clusters with TPL data
cluster_TPL <- TPL %>%
  pivot_longer(-target_id, names_to = "sample", values_to = "TPL") %>%
  inner_join(clusters, by = c("target_id" = "Transcript")) %>%
  left_join(meta, by = c("sample" = "Station")) %>%
  mutate(TPL_standardized = TPL / sum(TPL, na.rm = TRUE),
         TPL_standardized_z = scale(TPL_standardized))

# Calculate mean trajectory for each cluster, mean trajectory is mean expression of all transcripts in cluster
cluster_means <- cluster_TPL %>%
  group_by(Cluster, Date) %>%
  # Calculate mean TPL expression per cluster per sample
  dplyr::summarize(mean_TPL = mean(TPL_standardized, na.rm=TRUE), .groups = "drop")

# Add standardized log of TPL expression
cluster_TPL <- cluster_TPL %>%
  mutate(TPL_standardized_log = log10(TPL_standardized + 1e-6))

cluster_means <- cluster_means %>%
  mutate(mean_TPL_log = log10(mean_TPL + 1e-6))

# Calculate correlation of each transcript with cluster trajectory
representative_transcripts <- cluster_TPL %>%
  left_join(cluster_means, by = c("Cluster", "Date")) %>%
  group_by(Cluster, target_id) %>%
  dplyr::summarize(
    correlation = cor(TPL_standardized_log, mean_TPL_log, use = "complete.obs"),
    .groups = "drop"
  )

# Identify the top representative transcripts
top_representative_transcripts <- representative_transcripts %>%
  group_by(Cluster) %>%
  arrange(desc(correlation)) %>%
  ungroup()

# Save the results to a CSV file
output_csv <- "data/analysis/identified_clusters_transcripts_correlation.csv"
write.csv(top_representative_transcripts, output_csv, row.names = FALSE)

## Now we have this dataset, we can run MWU enrichment on the selected clusters
## using submit_cluster_MWU.pbs on the HPC

# Plot within-cluster sum of squares to justify number of clusters that we used
# Compute elbow curve (pseudo-inertia) values for different k
# Transpose matrix to cluster genes (rows = genes, columns = samples)
phaeo_transcripts <- annotation_file$query_id
tpl_phaeo <- TPL[TPL$target_id %in% phaeo_transcripts, ]
# Set target_id as rownames
rownames(tpl_phaeo) <- tpl_phaeo$target_id
tpl_phaeo <- tpl_phaeo[, -1]  # Remove target_id column
# Turn total TPL dataframe into a named numeric vector
tpl_totals <- setNames(phaeo_total_tpl$TPL, phaeo_total_tpl$sample)
# Ensure column order in tpl_phaeo matches tpl_totals
tpl_phaeo <- tpl_phaeo[, names(tpl_totals)]
normalized_data <- sweep(tpl_phaeo, 2, tpl_totals, FUN = "/") * 1e6
df <- as.matrix(normalized_data)
df <- scale(df)  # Important: scale data
# Subset top N most variable transcripts
top_n <- 5000
var_genes <- apply(df, 1, var)
top_genes <- order(var_genes, decreasing = TRUE)[1:top_n]
df <- df[top_genes, ]
wss_plot <- fviz_nbclust(df, kmeans, method = "wss", k.max = 6, nstart = 25) +
  labs(title = "Elbow Method for Optimal k (WSS)")
ggsave("figures/metatranscriptomics/wss_elbow_plot.png", plot = wss_plot, width = 6, height = 4, dpi = 800)
sil_plot <- fviz_nbclust(df, kmeans, method = "silhouette", k.max = 6, nstart = 25) +
  labs(title = "Average Silhouette Width for Optimal k")
ggsave("figures/metatranscriptomics/silhouette_plot.png", plot = sil_plot, width = 6, height = 4, dpi = 800)

# Merge clusters with TPL data
cluster_TPM <- TPM %>%
  pivot_longer(-target_id, names_to = "sample", values_to = "TPM") %>%
  inner_join(clusters, by = c("target_id" = "Transcript")) %>%
  left_join(meta, by = c("sample" = "Station")) %>%
  mutate(TPM_standardized = TPM / sum(TPM, na.rm = TRUE),
         TPM_standardized_z = scale(TPM_standardized))

# Calculate mean trajectory for each cluster, mean trajectory is mean expression of all transcripts in cluster
cluster_TPM_means <- cluster_TPM %>%
  group_by(Cluster, Date) %>%
  # Calculate mean TPM expression per cluster per sample
  dplyr::summarize(mean_TPM = mean(TPM_standardized, na.rm=TRUE), .groups = "drop")

# Add standardized log of TPL expression
cluster_TPM <- cluster_TPM %>%
  mutate(TPM_standardized_log = log10(TPM_standardized + 1e-6))

cluster_TPM_means <- cluster_TPM_means %>%
  mutate(mean_TPM_log = log10(mean_TPM + 1e-6))

# Add metadata to phaeo tpl sums
phaeo_total_tpl <- phaeo_total_tpl %>%
  dplyr::rename(sample = sample, phaeo_total_TPL = TPL) %>%
  left_join(meta, by = c("sample" = "Station"))

## VISUALIZE ##
# Top Left Panel: DIC and O2 smoother plot with secondary y-axis
plot_dic_o2 <- ggplot() +
  geom_line(data = DIC, aes(x = Date, y = y, color = "DIC"), size = 1) +
  geom_line(data = O2, aes(x = Date, y = y_scaled, color = "O2"), size = 1) +
  scale_color_manual(values = c("DIC" = "blue", "O2" = "green")) +
  scale_y_continuous(
    name = "DIC (µmol/kg)",
    sec.axis = sec_axis(~ reverse_scale_dic_to_o2(.), name = "O2 (µmol/kg)")
  ) +
  labs(title = "DIC and O2 Concentration Over Time",
       x = "Time",
       color = "Parameter") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.y.right = element_text(color = "green"),
    axis.title.y.left = element_text(color = "blue")
  )

# Top Right Panels: Positively and Negatively Correlated Marker Genes
plot_marker_positive_TPM <- ggplot(marker_TPM %>% filter(target_id %in% marker_genes$query_id[marker_genes$correlation == "Positive"]),
                               aes(x = Date, y = TPM_standardized, group = target_id, color = color)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_line(size = 1) +
  scale_color_identity() +
  labs(title = "Positively Correlated Marker Genes",
       x = "Time",
       y = "Phaeocystis TPM Fraction") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_marker_negative_TPM <- ggplot(marker_TPM %>% filter(target_id %in% marker_genes$query_id[marker_genes$correlation == "Negative"]),
                               aes(x = Date, y = TPM_standardized, group = target_id, color = color)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_line(size = 1) +
  scale_color_identity() +
  labs(title = "Negatively Correlated Marker Genes",
       x = "Time",
       y = "Phaeocystis TPM Fraction") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_marker_positive <- ggplot(marker_TPL %>% filter(target_id %in% marker_genes$query_id[marker_genes$correlation == "Positive"]),
                               aes(x = Date, y = TPL_standardized, group = target_id, color = color)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_line(size = 1) +
  scale_color_identity() +
  labs(title = "Positively Correlated Marker Genes",
       x = "Time",
       y = "Phaeocystis TPL Fraction") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_marker_negative <- ggplot(marker_TPL %>% filter(target_id %in% marker_genes$query_id[marker_genes$correlation == "Negative"]),
                               aes(x = Date, y = TPL_standardized, group = target_id, color = color)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_line(size = 1) +
  scale_color_identity() +
  labs(title = "Negatively Correlated Marker Genes",
       x = "Time",
       y = "Phaeocystis TPL Fraction") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Middle Barplot: Total Phaeo TPL sums
plot_total_tpl <- ggplot(phaeo_total_tpl, aes(x = Date, y = phaeo_total_TPL)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  labs(
    x = NULL,
    y = "Total Phaeocystis TPL"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),  # Remove x-axis labels to save space
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12, hjust = 0.5)  # Center-align title
  )

# Create cluster expression dynamics plot
plot_clusters <- ggplot() +
  geom_line(data = cluster_means, aes(x = Date, y = mean_TPL_log, group = Cluster, color = as.factor(Cluster)),
            size = 1.2) +
  scale_color_brewer(palette = "Set3", name = "Cluster") +
  labs(
    title = "Cluster Expression Dynamics",
    x = "Time",
    y = "Log10(avg(Standardized TPL))"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Create cluster TPM expression dynamics plot
plot_TPM_clusters <- ggplot() +
  geom_line(data = cluster_TPM_means, aes(x = Date, y = mean_TPM_log, group = Cluster, color = as.factor(Cluster)),
            size = 1.2) +
  scale_color_brewer(palette = "Set3", name = "Cluster") +
  labs(
    title = "Cluster Expression Dynamics",
    x = "Time",
    y = "Log10(avg(Standardized TPM))"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine the barplot and the cluster TPM expression dynamics plot
combined_cluster_plot_TPM <- (plot_total_tpl / plot_TPM_clusters) +
  plot_layout(heights = c(1, 3))
# Save the plot
ggsave("figures/metatranscriptomics/cluster_expression_TPM.svg",
       plot = combined_cluster_plot_TPM, width = 7, height = 5)

# Create an empty placeholder panel for whitespace
empty_panel <- ggplot() + theme_void()

# Combine the barplot and cluster expression dynamics plot **vertically**
combined_cluster_plot <- (plot_total_tpl / plot_clusters) +
  plot_layout(heights = c(1, 3))

# Add empty panel to the right
adjusted_cluster_plot <- combined_cluster_plot | empty_panel +
  plot_layout(widths = c(1, 1))

# Combine all panels
combined_plot <- (plot_dic_o2 + plot_marker_positive + plot_marker_negative) / adjusted_cluster_plot

# Save the plot
ggsave("figures/metatranscriptomics/marker_gene_expression_multiplot.svg",
       plot = combined_plot, width = 16, height = 12)
ggsave("figures/metatranscriptomics/marker_gene_expression_multiplot.png",
       plot = combined_plot, width = 16, height = 12, dpi = 800)
