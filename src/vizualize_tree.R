#!/usr/bin/env Rscript

library(ggplot2)
library(ggtree)
library(ape)    # For the read.tree function
library(here)
library(dplyr)
library(ComplexHeatmap) 

# --- 1. Load Custom Functions ---
source(here("src/heatmap_plot.R"))
source(here("src/col_palettes.R"))

# --- 2. Set File Paths & Parameters ---
treefile <- "trees/B216/tree.newick"
output_file <- "output/figures/B216_tree_only.pdf"

# Plot dimensions
plot_width <- 10

# Parameters (from plotHeatmap defaults)
clone_pal <- NULL
ladderize <- TRUE
normalize_tree <- FALSE # Set as needed
branch_length <- 1     # Only used if normalize_tree is TRUE

# --- 3. Read The Provided Tree ---
if (!file.exists(treefile)) {
  stop("Tree file not found: ", treefile)
}
message(paste("Reading tree from:", treefile))
mytree <- ape::read.tree(file = treefile)
cat("Tree has", length(mytree$tip.label), "tips\n")

# --- 4. Mimic plotHeatmap Logic ---
# Normalize tree if requested
normalize_tree <- FALSE
if (normalize_tree == TRUE) {
  # Assuming 'format_tree' is in your sourced files
  mytree <- format_tree(mytree, branch_length) 
}

# Create the dummy 'clusters' dataframe, just like plotHeatmap does.
# We use the tree's tip labels as the cell_id's.
message("Creating dummy cluster dataframe...")
clusters <- data.frame(
  cell_id = mytree$tip.label, 
  clone_id = "0"
)

# --- 5. Create the ggplot Tree Object ---
message("Creating ggtree plot...")
tree_ggplot <- make_tree_ggplot(
  tree = mytree,
  clones = as.data.frame(clusters), # Pass the dataframe to the 'clones' argument
  clone_pal = clone_pal,
  ladderize = ladderize
)
message("Creating tree...")
tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = 1)

message(paste("Saving tree to:", output_file))
pdf_width_inches <- 1.5 
pdf_height_inches <- 5 # Keep it tall to see the tips

message(paste("Saving tree to:", output_file))
pdf(
  output_file,
  width = pdf_width_inches,
  height = pdf_height_inches
)

# This is the 'save' command for ComplexHeatmap objects
draw(tree_hm)

# Close the PDF file
dev.off()

message("Tree saved successfully.")