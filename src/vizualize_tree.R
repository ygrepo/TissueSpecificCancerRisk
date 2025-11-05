library(ggplot2)
library(ggtree)
library(ape)      # For the read.tree function


# Set your file path
newick_file <- "tree.newick" # Or "results/latest/tree.newick"

# Read the tree data
my_tree <- read.tree(newick_file)

# --- Basic Plot (Rectangular) ---
ggtree(my_tree) +
  geom_tiplab(size = 3, align = TRUE) + # Add tip labels (cell names)
  geom_treescale() +                     # Add a scale bar
  theme_tree2()

# --- Fancier Plot (Circular) ---
ggtree(my_tree, layout = "circular") +
  geom_tiplab(size = 3, align = TRUE, linesize = 0.2) +
  theme(legend.position = "none") # Hide legend if there is one
