library(ggplot2)
library(ggtree)
library(ape)      # For the read.tree function
library(here)
source(here("src/heatmap_plot.R"))
source(here("src/col_palettes.R"))

# Set your file path
treefile <- "tree.newick" # Or "results/latest/tree.newick"

# Read the tree data
mytree <- ape::read.tree(file = treefile)
cat("Tree has", length(mytree$tip.label), "tips\n")


if (normalize_tree == T) {
  tree <- format_tree(tree, branch_length)
}

tree_ggplot <- make_tree_ggplot(tree,
                                as.data.frame(clusters),
                                clone_pal = clone_pal,
                                ladderize = ladderize)
tree_plot_dat <- tree_ggplot$data

message("Creating tree...")
tree_hm <- make_corrupt_tree_heatmap(tree_ggplot, tree_width = tree_width)
ordered_cell_ids <- get_ordered_cell_ids(tree_plot_dat)