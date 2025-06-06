# Load necessary libraries
if (!requireNamespace("MBCluster.Seq", quietly = TRUE)) {
    install.packages("MBCluster.Seq", repos = "http://cran.us.r-project.org")
}
library(MBCluster.Seq)

# Command-line arguments for file paths
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: Rscript mbcluster_phaeocystis.R <tpl_file> <annotation_file> <O2_file> <metadata_file> <output_file>")
}

tpl_file <- args[1]           # Path to the TPL data (expression counts)
annotation_file <- args[2]    # Path to the Phaeocystis annotation file
O2_file <- args[3]           # Path to the O2 smoother file
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

# Normalize by total Phaeocystis activity
cat("Normalizing by library size...\n")
phaeo_sums <- colSums(tpl_data)
normalized_data <- sweep(tpl_data, 2, phaeo_sums, FUN = "/") * 1e6

# Load O2 smoother and metadata
cat("Loading O2 smoother and metadata...\n")
O2_data <- read.csv(O2_file)
metadata <- read.csv(metadata_file)

# Filter metadata for Station 130
metadata <- metadata[metadata$StationPrefix == 130, ]

# Convert dates to POSIXct for nearest timepoint matching
metadata$Date <- as.POSIXct(metadata$Date, format = "%Y-%m-%d %H:%M:%S")
O2_data$Date <- as.POSIXct(O2_data$Date, format = "%Y-%m-%d %H:%M:%S")

# Match nearest O2 value to each sample
cat("Matching nearest O2 values...\n")
metadata$O2 <- sapply(metadata$Date, function(sample_time) {
    # Find the nearest O2 timestamp
    diffs <- abs(difftime(O2_data$Date, sample_time, units = "secs"))
    nearest_index <- which.min(diffs)
    O2_data$y[nearest_index]  # Return the matched O2 value
})

# Define low- and high-O2 groups
median_O2 <- median(metadata$O2, na.rm = TRUE)  # Use median O2 as a threshold
metadata$O2_group <- ifelse(metadata$O2 < median_O2, "Low", "High")

# Prepare treatment variable
cat("Preparing treatment variable...\n")
Treatment <- ifelse(metadata$O2_group == "Low", 1, 2)  # 1 for Low-O2, 2 for High-O2
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
num_clusters <- 2 # Adjust the number of clusters if needed
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