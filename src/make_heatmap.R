#!/usr/bin/env Rscript

library(argparse)
library(ape)
library(glue)
library(dplyr)
library(grid)
library(ComplexHeatmap)
#library(signals)
library(magick)
library(tidyr)
library(ggtree)
# Suppress ComplexHeatmap messages
ht_opt$message = FALSE
library(here)
source(here("src/heatmap_plot.R"))
source(here("src/col_palettes.R"))


#' Generate a Heatmap Aligned with a Phylogenetic Tree
#' 
#' @param treefile Path to the tree file (Newick format)
#' @param clusters Optional data frame with cell_id column for subsetting
#' @param cnv_data Data frame containing CNV data with cell_id column (required only if clusters not provided)
#' @param output_file Optional path to save the output plot
#' @param plot_width Width of the plot (default: 89 * 0.039)
#' @param plot_height Height of the plot (default: 2)
#' @param chroms Character vector of chromosome labels
#' @param title Plot title (default: empty string)
#' @param tree_width Width of the tree portion (default: 1)
#' @param linkheight Height of the linking lines (default: 2)
#' @param clone_colors Named vector of colors for clones
#' @return List containing heatmap object and tree object
make_heatmap_tree <- function(treefile, 
                              clusters = NULL, 
                              use_umap_clusters = FALSE,
                              cnv_data = NULL,
                              output_file = NULL,
                              plot_width = 89 * 0.039,
                              plot_height = 2,
                              chroms = c(paste0(1:11), "13", "15", "17", "20", "X"),
                              title = "",
                              plot_tree = TRUE,
                              tree_width = 1,
                              linkheight = 2,
                              chr13_17_deletion = FALSE,
                              chr13_17_deletion_threshold = 0.5,
                              clone_colors = c("A" = "firebrick4", "B" = "deepskyblue4")) {
  
  # ---- Basic input sanity checks ----
  if (is.null(cnv_data)) {
    stop("cnv_data must be provided.")
  }
  if (!"cell_id" %in% colnames(cnv_data)) {
    stop("cnv_data must contain 'cell_id' column")
  }
  
  # ---- Tree vs. no-tree (UMAP) mode ----
  if (use_umap_clusters) {
    cat("Using UMAP+HDBSCAN clustering; no external treefile will be read.\n")
    mytree <- NULL
  } else {
    if (is.null(treefile) || !file.exists(treefile)) {
      stop("Tree file does not exist: ", treefile)
    }
    cat("Reading tree from:", treefile, "\n")
    mytree <- ape::read.tree(file = treefile)
    cat("Tree has", length(mytree$tip.label), "tips\n")
    
    # Remove "cell_" prefix from tree tip labels if present
    if (any(grepl("^cell_", mytree$tip.label))) {
      cat("Removing 'cell_' prefix from tree tip labels\n")
      mytree$tip.label <- gsub("^cell_", "", mytree$tip.label)
      cat("Updated tree tip labels (first 3):",
          paste(head(mytree$tip.label, 3), collapse = ", "), "\n")
    }
    
    # If clusters is given, restrict the tree
    if (!is.null(clusters)) {
      if (!"cell_id" %in% colnames(clusters)) {
        stop("clusters data frame must contain 'cell_id' column")
      }
      cat("Filtering tree to", nrow(clusters), "cells from clusters\n")
      mytree <- ape::keep.tip(mytree, clusters$cell_id)
    }
  }
  
  cat("CNV data has", nrow(cnv_data), "rows and",
      length(unique(cnv_data$cell_id)), "unique cells\n")
  
  # ---- Align CNV with tree, or use all cells in UMAP mode ----
  if (!is.null(mytree)) {
    overlap_cells <- intersect(mytree$tip.label, unique(cnv_data$cell_id))
    cat("Found", length(overlap_cells),
        "overlapping cells between tree and CNV data\n")
    
    if (length(overlap_cells) == 0) {
      cat("Tree tips (first 5):",
          paste(head(mytree$tip.label, 5), collapse = ", "), "\n")
      cat("CNV cell_ids (first 5):",
          paste(head(unique(cnv_data$cell_id), 5), collapse = ", "), "\n")
      stop("No overlapping cells found between tree and CNV data")
    }
    cat("Tip:", length(mytree$tip.label))               # tips in the tree
    cnaneuploid <- cnv_data %>%
      dplyr::filter(cell_id %in% mytree$tip.label)
    cat("cnaneuploid:", length(unique(cnaneuploid$cell_id)))     # rows in the heatmap
  } else {
    # No tree: use all CNV cells
    cnaneuploid <- cnv_data
  }
  
  cat("Filtered CNV data has", nrow(cnaneuploid), "rows for",
      length(unique(cnaneuploid$cell_id)), "cells\n")
  
  # ---- chr13/17 deletion annotations (works with or without tree) ----
  annotations_df <- NULL
  if (chr13_17_deletion) {
    message("Generating chr13/17 deletion annotations...")
    
    cnv_filtered_lumi <- cnaneuploid %>%
      mutate(del = ifelse(state < 2, 1, 0)) %>%
      group_by(cell_id, chr) %>%
      summarise(binary = mean(del) >= chr13_17_deletion_threshold,
                .groups = "drop") %>%
      mutate(binary = as.numeric(binary))
    
    cnv_wide <- cnv_filtered_lumi %>%
      dplyr::filter(chr %in% c("13", "17")) %>%
      tidyr::pivot_wider(
        names_from  = chr,
        values_from = binary,
        names_prefix = "chr"
      )
    
    # Base set of cells for annotations: tree tips if tree exists, else CNV cells
    if (!is.null(mytree)) {
      all_cells <- data.frame(cell_id = mytree$tip.label)
    } else {
      all_cells <- data.frame(cell_id = unique(cnaneuploid$cell_id))
    }
    
    cnv_wide <- all_cells %>%
      left_join(cnv_wide, by = "cell_id") %>%
      mutate(
        chr13 = tidyr::replace_na(chr13, 0),
        chr17 = tidyr::replace_na(chr17, 0)
      )
    
    cnv_wide_deletions <- cnv_wide %>%
      mutate(
        Deletion = dplyr::case_when(
          chr17 == 1 & chr13 == 1 ~ "Both",
          chr17 == 1             ~ "17q",
          chr13 == 1             ~ "13q",
          TRUE                   ~ "None"
        )
      ) %>%
      dplyr::select(cell_id, Deletion)
    
    deletion_summary <- table(cnv_wide_deletions$Deletion)
    message(
      "Deletion summary: ",
      paste(names(deletion_summary), "=", deletion_summary, collapse = ", ")
    )
    
    # Join clusters (if present) + deletions
    if (!is.null(clusters)) {
      message("Joining clone clusters and deletion data...")
      annotations_df <- all_cells %>%
        left_join(as.data.frame(clusters), by = "cell_id") %>%
        left_join(as.data.frame(cnv_wide_deletions), by = "cell_id")
      
      annotations_df$clone_id[is.na(annotations_df$clone_id)] <- "None"
      annotations_df$Deletion[is.na(annotations_df$Deletion)] <- "None"
    } else {
      message("Joining deletion data only...")
      annotations_df <- all_cells %>%
        left_join(as.data.frame(cnv_wide_deletions), by = "cell_id")
      annotations_df$Deletion[is.na(annotations_df$Deletion)] <- "None"
    }
    
    message("Final annotations_df has ", nrow(annotations_df), " rows")
    message("Final deletion summary: ",
            paste(table(annotations_df$Deletion), collapse = ", "))
  }
  
  # ---- Chromosome diagnostics and auto-adjust ----
  cat("=== CHROMOSOME DIAGNOSTIC ===\n")
  actual_chroms <- sort(unique(cnaneuploid$chr))
  expected_chroms <- chroms
  
  cat("Chromosomes in data:", paste(actual_chroms, collapse = ", "), "\n")
  cat("Expected chromosomes:", paste(expected_chroms, collapse = ", "), "\n")
  
  missing_in_data <- setdiff(expected_chroms, actual_chroms)
  extra_in_data   <- setdiff(actual_chroms, expected_chroms)
  
  if (length(missing_in_data) > 0) {
    cat("Missing from data:", paste(missing_in_data, collapse = ", "), "\n")
  }
  if (length(extra_in_data) > 0) {
    cat("Extra in data:", paste(extra_in_data, collapse = ", "), "\n")
  }
  
  if (length(missing_in_data) > 0 || length(extra_in_data) > 0) {
    cat("🔧 Auto-adjusting chromosome list to match data\n")
    chroms <- actual_chroms
    cat("Updated chroms:", paste(chroms, collapse = ", "), "\n")
  }
  cat("============================\n")
  
  cat("Chromosome column type:", class(cnaneuploid$chr), "\n")
  
  if (nrow(cnaneuploid) == 0) {
    stop("No CNV data remains after filtering to match tree cells")
  }
  
  if (!exists("plotHeatmap")) {
    stop("plotHeatmap function not found. Please ensure it's loaded.")
  }
  
  # ---- Plotting ----
  p <- plotHeatmap(
    cnaneuploid,
    tree = if (use_umap_clusters) NULL else mytree %>% ape::ladderize(),
    clusters = NULL,
    column_title = title,
    column_title_gp = gpar(fontsize = 10),
    linkheight = linkheight,
    chrlabels = chroms,
    plotfrequency = FALSE,
    frequency_height = 0.5,
    anno_width = 0.6,
    annofontsize = 8,
    labeladjust = -8,
    show_heatmap_legend = TRUE,
    show_legend = TRUE,
    show_clone_text = TRUE,
    show_library_label = FALSE,
    show_clone_label = TRUE,
    plottree = if (use_umap_clusters) FALSE else plot_tree,
    reorderclusters = TRUE,
    tree_width = tree_width,
    clone_pal = clone_colors,
    annotations = annotations_df,
    chr13_17_deletion = chr13_17_deletion,
    labels_rot = 0,
    extend_val = 0.15
  )
  
  pout <- grid.grabExpr(draw(p), width = plot_width, height = plot_height)
  
  if (!is.null(output_file)) {
    ext <- tools::file_ext(output_file)
    if (ext %in% c("pdf", "PDF")) {
      pdf(output_file, width = plot_width, height = plot_height)
      draw(p)
      dev.off()
    } else if (ext %in% c("png", "PNG")) {
      png(output_file, width = plot_width * 100,
          height = (plot_height + 1) * 100, res = 300)
      draw(p)
      dev.off()
    }
  }
  
  return(list(hm = pout, tree = mytree))
}

# Main function with argparse
main <- function() {
  # Set up for both interactive and command-line use
  if (interactive()) {
    # Mimic command-line input for testing
    argv <- c("-t", "tree.nwk",
              "-c", "clusters.csv",
              "-o", "output.pdf",
              "--width", "5",
              "--height", "3",
              "--title", "Test Heatmap")
  } else {
    # Get real command-line arguments
    argv <- commandArgs(trailingOnly = TRUE)
  }
  
  # Create argument parser
  parser <- ArgumentParser(description = "Generate a heatmap aligned with a phylogenetic tree")
  
  # Required arguments
  parser$add_argument("-t", "--treefile", type = "character", 
                      help = "Path to the tree file (Newick format)")
  
  # Optional arguments
  parser$add_argument("-c", "--clusters", type = "character", default = NULL,
                      help = "CSV file with cell_id column for subsetting")
  
  parser$add_argument("-d", "--cnv_data", type = "character", default = NULL,
                      help = "CSV file with CNV data (required only if clusters not provided)")
  
  parser$add_argument("-o", "--output", type = "character", default = NULL,
                      help = "Output file path (pdf or png)")
  
  parser$add_argument("--width", type = "double", default = 5,
                      help = "Plot width (default: 3.471)")
  
  parser$add_argument("--height", type = "double", default = 5,
                      help = "Plot height (default: 2)")
  
  parser$add_argument("--title", type = "character", default = "",
                      help = "Plot title (default: empty)")
  
  parser$add_argument("--tree_width", type = "double", default = 6,
                      help = "Width of the tree portion (default: 1)")
  
  parser$add_argument("--linkheight", type = "double", default = 8,
                      help = "Height of the linking lines (default: 4)")
  
  parser$add_argument("--chroms", type = "character", nargs = "*",
                      default = NULL,
                      help = "Chromosome labels (space-separated). If not provided, will auto-detect from data")
  
  parser$add_argument("--plot_tree", action = "store_true", default = TRUE,
                      help = "Plot the tree alongside the heatmap")
  
  parser$add_argument("--chr13_17_deletion", action = "store_true", default = FALSE,
                      help = "Include chr13/17 deletion annotations")
 
  parser$add_argument("--chr13_17_deletion_threshold", type = "double", default = .5,
                      help = "chr13/17 deletion threshold")
  
  parser$add_argument("--use_umap_clusters", action = "store_true",
                       help = "Ignore external tree and cluster cells by UMAP+HDBSCAN")
   
   
  # Parse arguments
  args <- parser$parse_args(argv)
  
  # Load data files
  clusters <- NULL
  cnv_data <- NULL
  
  if (!is.null(args$clusters)) {
    if (!file.exists(args$clusters)) {
      stop("Clusters file does not exist: ", args$clusters)
    }
    clusters <- read.csv(args$clusters)
  }
  
  if (!is.null(args$cnv_data)) {
    if (!file.exists(args$cnv_data)) {
      stop("CNV data file does not exist: ", args$cnv_data)
    }
    cnv_data <- read.csv(args$cnv_data)
  }
  
  if (args$use_umap_clusters) {
    # UMAP mode: require CNV, ignore treefile
    if (is.null(cnv_data)) {
      stop("When --use_umap_clusters is set, you must provide --cnv_data.")
    }
  } else {
    # Tree mode: require treefile and CNV
    if (is.null(args$treefile)) {
      stop("You must provide --treefile when not using --use_umap_clusters.")
    }
    if (is.null(cnv_data)) {
      stop("You must provide --cnv_data (copy-number bins) to plot the heatmap.")
    }
  }
  
  # Validate that either clusters or cnv_data is provided
  # if (is.null(clusters) && is.null(cnv_data)) {
  #   stop("Either --clusters or --cnv_data must be provided")
  # }
  
  # Auto-detect or use provided chromosomes
  if (is.null(args$chroms) && !is.null(cnv_data)) {
    # Auto-detect chromosomes from data
    chroms_to_use <- sort(unique(as.character(cnv_data$chr)))
    cat("Auto-detected chromosomes from data:", paste(chroms_to_use, collapse = ", "), "\n")
  } else if (!is.null(args$chroms)) {
    chroms_to_use <- args$chroms
  } else {
    # Fallback to default
    chroms_to_use <- c(paste0(1:11), "13", "15", "17", "20", "X")
  }
  
  # Set up clone colors
  clone_colors <- c("A" = "firebrick4", "B" = "deepskyblue4")
  
  
  # Call the function
  tryCatch({
    result <- make_heatmap_tree(
      treefile = args$treefile,
      clusters = clusters,
      use_umap_clusters = args$use_umap_clusters,
      cnv_data = cnv_data,
      output_file = args$output,
      plot_width = args$width,
      plot_height = args$height,
      chroms = chroms_to_use, 
      plot_tree = args$plot_tree,
      chr13_17_deletion = args$chr13_17_deletion,
      chr13_17_deletion_threshold = args$chr13_17_deletion_threshold,
      title = args$title,
      tree_width = args$tree_width,
      linkheight = args$linkheight,
      clone_colors = clone_colors
    )
    
    cat("Heatmap tree generated successfully!\n")
    if (!is.null(args$output)) {
      cat("Output saved to:", args$output, "\n")
    }
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
    quit(status = 1)
  })
}
# Run main function if script is executed directly
if (!interactive()) {
  main()
}