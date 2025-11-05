#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to display usage
show_usage <- function() {
  cat("Usage: Rscript cnv_processor.R <input_file> <output_file> [--seed <number>]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  input_file   Path to input CSV file with CNV data\n")
  cat("  output_file  Path for output CSV file\n")
  cat("  --seed       Optional random seed (default: 42)\n")
  cat("\n")
  cat("Example:\n")
  cat("  Rscript cnv_processor.R data/SA501.tbnc.cnv.csv data/output.csv\n")
  cat("  Rscript cnv_processor.R data/input.csv data/output.csv --seed 123\n")
  quit(status = 1)
}

# Check if help is requested or insufficient arguments
if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
  show_usage()
}

if (length(args) < 2) {
  cat("Error: Missing required arguments\n\n")
  show_usage()
}

# Parse arguments
input_file <- args[1]
output_file <- args[2]
seed_value <- 42  # default

# Parse optional seed argument
if ("--seed" %in% args) {
  seed_idx <- which(args == "--seed")
  if (seed_idx == length(args)) {
    cat("Error: --seed requires a value\n")
    quit(status = 1)
  }
  seed_value <- as.numeric(args[seed_idx + 1])
  if (is.na(seed_value)) {
    cat("Error: Invalid seed value\n")
    quit(status = 1)
  }
}

# Validate input file exists
if (!file.exists(input_file)) {
  cat("Error: Input file does not exist:", input_file, "\n")
  quit(status = 1)
}

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyverse)
  library(tidyr)
  library(stringr)
  library(readr)
  library(rlang)
  library(data.table)
})

# Clear environment and set seed
rm(list = ls()[!ls() %in% c("input_file", "output_file", "seed_value")])
set.seed(seed_value)
date <- Sys.Date()

cat("Processing CNV data...\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Seed:", seed_value, "\n")
cat("Date:", as.character(date), "\n\n")

# Read data
cat("Reading input data...\n")
tryCatch({
  data <- read.table(input_file, 
                     sep = ",",
                     header = TRUE,
                     stringsAsFactors = FALSE)
  cat("Successfully read", nrow(data), "rows and", ncol(data), "columns\n")
}, error = function(e) {
  cat("Error reading input file:", e$message, "\n")
  quit(status = 1)
})

#' Process Long CNV Data to a Sorted Wide Format
#'
#' Takes a long-format CNV data frame and transforms it into a wide format
#' where rows are sorted genomic loci and columns are cell IDs.
#'
#' @param long_cnv_data A data frame or tibble in long format, containing
#'   at least columns: cell_id, chr, start, end, state.
#'
#' @return A tibble in wide format. The first column is "loci"
#'   (e.g., "1_1000_2000"), and subsequent columns are cell_ids,
#'   with numeric CNV states as values. The rows are sorted
#'   numerically by chromosome (1-22, X, Y) and then by start position.
#'
format_cnv_wide <- function(long_cnv_data) {
  
  # --- 1. Pivot from Long to Wide ---
  cnv_wide <- long_cnv_data %>%
    transmute(
      cell_id = as.character(cell_id),
      chr     = as.character(chr),
      start   = as.integer(start),
      end     = as.integer(end),
      state   = as.numeric(state)
    ) %>%
    # normalize chr labels like "chr13" -> "13"
    mutate(chr = str_replace(chr, "^chr", "")) %>%
    # build locus id "CHR_START_END"
    mutate(loci = str_c(chr, start, end, sep = "_")) %>%
    select(loci, cell_id, state) %>%
    # handle accidental duplicates deterministically
    distinct(loci, cell_id, .keep_all = TRUE) %>%
    # wide: rows=loci, columns=cell_id, values=state
    pivot_wider(names_from = cell_id, values_from = state, values_fill = NA_real_)
  
  # --- 2. Sort Loci Properly (Genomic Order) ---
  
  # Use %||% (null-coalescing operator) to safely get loci.
  # This works because pivot_wider creates a tibble where 'loci' is a column.
  # rownames(cnv_wide) will be NULL, so it correctly uses cnv_wide$loci.
  all_loci <- cnv_wide$loci %||% rownames(cnv_wide)
  
  chr_part   <- str_replace(all_loci, "^([0-9XY]+)_.*$", "\\1")
  start_part <- as.integer(str_replace(all_loci,
                                       "^[0-9XY]+_([0-9]+)_.*$", "\\1"))
  
  # Map chr to sortable numbers (1-22, X=23, Y=24)
  chr_num <- case_when(
    chr_part %in% as.character(1:22) ~ as.integer(chr_part),
    chr_part == "X" ~ 23L,
    chr_part == "Y" ~ 24L,
    TRUE ~ 999L # Other/unknown chromosomes sort last
  )
  
  # Get the row order
  ord <- order(chr_num, start_part, na.last = TRUE)
  
  # Apply the sorting
  cnv_wide <- cnv_wide[ord, ]
  
  # --- 3. Final Formatting for Output ---
  
  # This block ensures 'loci' is the first column and there are no row names.
  # It's a robust way to prepare for writing to CSV.
  out <- tibble(loci = cnv_wide$loci %||% rownames(cnv_wide)) %>%
    bind_cols(as_tibble(cnv_wide %>% select(-loci), .name_repair = "minimal"))
  
  # Ensure all cell columns are numeric
  num_cols <- setdiff(names(out), "loci")
  out[num_cols] <- lapply(out[num_cols], as.numeric)
  
  # Return the final processed tibble
  return(out)
}

# Process the data
cat("Processing CNV data to wide format...\n")
tryCatch({
  processed_cnv_data <- format_cnv_wide(data)
  cat("Successfully processed data to", nrow(processed_cnv_data), "loci and", 
      ncol(processed_cnv_data)-1, "cells\n")
}, error = function(e) {
  cat("Error processing data:", e$message, "\n")
  quit(status = 1)
})

# Write output
cat("Writing output file...\n")
tryCatch({
  # Create output directory if it doesn't exist
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  write_csv(processed_cnv_data, output_file)
  cat("Successfully wrote output to:", output_file, "\n")
  
  # Show first few lines as verification
  cat("\nFirst 3 lines of output:\n")
  first_lines <- readLines(output_file, n = 3)
  cat(paste(first_lines, collapse = "\n"), "\n")
  
}, error = function(e) {
  cat("Error writing output file:", e$message, "\n")
  quit(status = 1)
})

cat("\nProcessing completed successfully!\n")