#!/usr/bin/env Rscript

library(ape)
library(ggtree)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)

# --- 1. Define File Paths ---
tree_file <- here("trees/B216/tree.newick")
# This is the file you've been passing to --cnv_data
cnv_file <- here("data/allele_specific_cn/B2HET16-hscn.csv") 
output_pdf <- here("output/figures/B216_tree_with_deletions.pdf")

# --- 2. Load Data ---
message(paste("Loading tree:", tree_file))
if (!file.exists(tree_file)) stop("Tree file not found: ", tree_file)
tree <- read.tree(tree_file)

message(paste("Loading CNV data:", cnv_file))
if (!file.exists(cnv_file)) stop("CNV file not found: ", cnv_file)
B2.16_df <- read.csv(cnv_file)

# --- 3. Implement Your Deletion Logic ---

# Step 1: Convert deletions to 1, anything else to 0
message("Filtering for chr13/17 deletions...")
cnv_filtered <- B2.16_df %>%
  filter(chr %in% c("17", "13")) %>%
  mutate(binary = ifelse(state < 2, 1, 0)) %>%
  dplyr::select(cell_id, chr, binary) %>%
  group_by(cell_id, chr) %>%
  summarise(binary = max(binary), .groups = "drop") # 1 if any deletion, else 0

# Step 2: Pivot to wide format
message("Pivoting data to wide format...")
cnv_wide <- cnv_filtered %>%
  pivot_wider(
    names_from = chr,
    values_from = binary,
    names_prefix = "chr",
    values_fill = 0 # Ensure cells missing a chr get 0, not NA
  )

# Step 3: Remove cell label on tree
tree$tip.label <- sub("^cell_", "", tree$tip.label)

# Step 4: Make a data frame with tip labels and join
message("Joining deletion data to tree tips...")
tip_data <- data.frame(label = tree$tip.label) %>%
  left_join(cnv_wide, by = c("label" = "cell_id"))

# Handle any cells in the tree that weren't in the CNV file
tip_data <- tip_data %>% 
  mutate(
    chr13 = ifelse(is.na(chr13), 0, chr13),
    chr17 = ifelse(is.na(chr17), 0, chr17)
  )

# --- 4. Generate and Save Tree Plot ---
message("Generating ggtree plot...")
tree_plot <- ggtree(tree) %<+% tip_data +
  geom_tiplab(aes(color = case_when(
    chr17 == 1 & chr13 == 1 ~ "Both", # Added this case for clarity
    chr17 == 1 ~ "17q",
    chr13 == 1 ~ "13q",
    TRUE ~ "none"
  )), size = 2) +
  scale_color_manual(
    name = "Chromosomal Deletion",
    values = c("17q" = "red", "13q" = "blue", "Both" = "purple", "none" = "gray70"),
    labels = c("17q" = "17q deletion", "13q" = "13q deletion", "Both" = "13q & 17q deletion", "none" = "No Deletion")
  ) +
  theme_tree2() +
  ggtitle("B2.16: Phylogenetic Tree with chr17/13 deletions") +
  theme(legend.position = "right") # Added a legend

# --- 5. Save the Plot ---
# We use a tall height to make sure the tip labels are readable
message(paste("Saving plot to:", output_pdf))
ggsave(
  output_pdf,
  plot = tree_plot,
  width = 5,
  height = 8, # Start with a tall plot
  units = "in",
  limitsize = FALSE
)

message("Done.")

