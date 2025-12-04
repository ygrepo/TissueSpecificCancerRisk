#!/usr/bin/env Rscript

library(ape)    # For reading and cutting the tree
library(dplyr)
library(here)

# --- 1. Define File Paths ---
tree_file <- "trees/B216/tree.newick"
output_cluster_file <- "data/B216_clone_assignments.csv"

# --- 2. Load Tree ---
message(paste("Loading tree from:", tree_file))
if (!file.exists(tree_file)) {
  stop("Could not find tree file. Make sure the path is correct.")
}
mytree <- ape::read.tree(tree_file)

# --- 3. Cut The Tree into 2 Groups ---
# We will convert the 'phylo' object to an 'hclust' object
# which allows us to use the standard 'cutree' function.
message("Converting tree to hclust object...")
hclust_tree <- ape::as.hclust.phylo(mytree)

message("Cutting tree into k=2 groups...")
# This cuts the tree into exactly 2 groups (k = 2)
groups <- stats::cutree(hclust_tree, k = 2)

# --- 4. Format and Save Cluster File ---
# Convert the named vector into a proper data frame
# The 'groups' object looks like:
# cell-1 cell-2 cell-3
#      1      1      2
clone_assignments <- data.frame(
  cell_id = names(groups),
  clone_id = groups
)

# Rename '1' and '2' to 'A' and 'B' to match the screenshot
# (This is optional, but it's what you wanted)
message("Renaming clusters '1' and '2' to 'A' and 'B'...")
clone_assignments$clone_id <- recode(clone_assignments$clone_id,
                                     `1` = "A",
                                     `2` = "B")

# Save the final cluster file
write.csv(clone_assignments, output_cluster_file, row.names = FALSE, quote = FALSE)

message(paste("\nSuccess! Tree-based clone assignments saved to:", output_cluster_file))
message("The first few rows of your new file:")
print(head(clone_assignments))