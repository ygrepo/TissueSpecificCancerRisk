#!/usr/bin/env Rscript

library(here)
library(dplyr)

# --- 1. Load Functions ---
# We need to load the same helper functions your
# main script uses. 'umap_clustering' is in one of these.
# Make sure 'src/util.R' or 'src/clustering.R' is sourced
# inside one of these files, or add it here.
# Based on your files, 'util.R' is the one.
source(here("src/util.R"))
source(here("src/clustering.R"))

# --- 2. Define File Paths ---
# This is the data file you are passing to --cnv_data
cnv_input_file <- "data/allele_specific_cn/B2HET16-hscn.csv" 
# This is the new cluster file we will create
output_cluster_file <- "data/B216_clone_assignments.csv"

# --- 3. Load Data ---
message(paste("Loading CNV data from:", cnv_input_file))
if (!file.exists(cnv_input_file)) {
  stop("Could not find CNV file. Make sure the path is correct.")
}
cnv_data <- read.csv(cnv_input_file)

# --- 4. Set Clustering Parameters ---
# These are the default parameters from your plotHeatmap function
message("Setting clustering parameters...")
ncells <- length(unique(cnv_data$cell_id))
pctcells <- 0.05
min_pts <- max(round(pctcells * ncells), 2)
umap_metric <- "euclidean"
seed <- NULL

# --- 5. Run Clustering ---
# This is the same function plotHeatmap would call
message("Running umap_clustering... This may take a moment.")
clustering_results <- umap_clustering(
  cnv_data,
  minPts = min_pts,
  field = "copy",
  umapmetric = umap_metric,
  seed = seed
)

# --- 6. Extract and Save Clusters ---
# The 'clustering_results' object contains the tree and the clusters.
# We just want the 'clustering' data frame.
message("Extracting clone assignments...")
clone_assignments <- clustering_results$clustering %>%
  dplyr::select(cell_id, clone_id)

# Save the assignments to our new file
write.csv(clone_assignments, output_cluster_file, row.names = FALSE, quote = FALSE)

message(paste("\nSuccess! Clone assignments saved to:", output_cluster_file))
message("The first few rows of your new file:")
print(head(clone_assignments))