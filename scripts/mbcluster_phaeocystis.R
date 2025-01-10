# Load necessary libraries
if (!requireNamespace("MBCluster.Seq", quietly = TRUE)) {
    install.packages("MBCluster.Seq", repos = "http://cran.us.r-project.org")
}
library(MBCluster.Seq)

# Command-line arguments for file paths
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: Rscript mbcluster_nearest_dic.R <tpl_file> <annotation_file> <dic_file> <metadata_file> <output_file>")
}

tpl_file <- args[1]           # Path to the TPL data (expression counts)
annotation_file <- args[2]    # Path to the Phaeocystis annotation file
dic_file <- args[3]           # Path to the DIC smoother file
metadata_file <- args[4]      # Path to the sample metadata file
output_file <- args[5]        # Output file for clustering results

# Load expression data
cat("Loading TPL data...\n")
tpl_data <- read.csv(tpl_file, row.names = 1)

# Load Phaeocystis annotations
cat("Loading Phaeocystis annotations...\n")
annotations <- read.csv(annotation_file)

# Filter TPL data for Phaeocystis transcripts
cat("Filtering for Phaeocystis-associated transcripts...\n")
phaeo_transcripts <- annotations$query_id
tpl_data <- tpl_data[rownames(tpl_data) %in% phaeo_transcripts, ]

# Filter out lowly expressed transcripts
cat("Filtering lowly expressed transcripts...\n")
min_reads <- 5
tpl_data <- tpl_data[rowMeans(tpl_data) >= min_reads, ]
if (nrow(tpl_data) < 2 || ncol(tpl_data) < 2) {
    stop("Insufficient data after filtering: At least two transcripts and two samples are required.")
}

# Normalize by total Phaeocystis activity
cat("Normalizing by library size...\n")
phaeo_sums <- colSums(tpl_data)
normalized_data <- sweep(tpl_data, 2, phaeo_sums, FUN = "/") * 1e6

# Load DIC smoother and metadata
cat("Loading DIC smoother and metadata...\n")
dic_data <- read.csv(dic_file)
metadata <- read.csv(metadata_file)

# Filter metadata for Station 130
metadata <- metadata[metadata$StationPrefix == 130, ]

# Convert dates to POSIXct for nearest timepoint matching
metadata$Date <- as.POSIXct(metadata$Date, format = "%Y-%m-%d %H:%M:%S")
dic_data$Date.Time <- as.POSIXct(dic_data$Date.Time, format = "%Y-%m-%d %H:%M:%S")

# Match nearest DIC value to each sample
cat("Matching nearest DIC values...\n")
metadata$DIC <- sapply(metadata$Date, function(sample_time) {
    # Find the nearest DIC timestamp
    diffs <- abs(difftime(dic_data$Date.Time, sample_time, units = "secs"))
    nearest_index <- which.min(diffs)
    dic_data$y[nearest_index]  # Return the matched DIC value
})

# Define low- and high-DIC groups
median_dic <- median(metadata$DIC, na.rm = TRUE)  # Use median DIC as a threshold
metadata$DIC_group <- ifelse(metadata$DIC < median_dic, "Low", "High")

# Prepare treatment variable
cat("Preparing treatment variable...\n")
Treatment <- ifelse(metadata$DIC_group == "Low", 1, 2)  # 1 for Low-DIC, 2 for High-DIC
GeneID <- rownames(normalized_data)

# Create RNASeq.Data object
cat("Creating RNASeq.Data object...\n")
RNASeq_data <- RNASeq.Data(
    Count = as.matrix(normalized_data),
    Normalize = NULL,
    Treatment = Treatment,
    GeneID = GeneID
)

# Initialize clustering with K-means
cat("Initializing clustering with K-means...\n")
num_clusters <- 10  # Adjust the number of clusters as needed
kmeans_centers <- KmeansPlus.RNASeq(RNASeq_data, nK = num_clusters)$centers

# Perform clustering with Cluster.RNASeq
cat("Running MBCluster.Seq clustering...\n")
cluster_result <- Cluster.RNASeq(
    data = RNASeq_data,
    model = "nbinom",   # Negative binomial model
    centers = kmeans_centers,
    method = "EM"       # Expectation-Maximization algorithm
)

# Save clustering results
cat("Saving clustering results...\n")
clusters <- data.frame(Transcript = RNASeq_data$GeneID, Cluster = cluster_result$cluster)
write.csv(clusters, output_file, row.names = FALSE)

cat("Clustering completed successfully.\n")