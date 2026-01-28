#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to display usage
show_usage <- function() {
  cat("Usage: Rscript read_cnv_data.R <input_file> <output_file> [--patient <id> | --sample_id <id>] [--seed <number>]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  input_file   Path to input CSV file with CNV data\n")
  cat("  output_file  Path for output CSV file\n")
  cat("\n")
  cat("Options:\n")
  cat("  --patient <id>     Filter to a single patient (matches column 'patient')\n")
  cat("  --sample_id <id>   Filter to a single sample_id (matches column 'sample_id')\n")
  cat("  --seed <number>    Set RNG seed (optional)\n")
  cat("\n")
  cat("Examples:\n")
  cat("  Rscript read_cnv_data.R data/SA501.tbnc.cnv.csv data/output.csv --patient SA501\n")
  cat("  Rscript read_cnv_data.R data/SA501.tbnc.cnv.csv data/output.csv --sample_id SA501\n")
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

# Parse required positional arguments
input_file <- args[1]
output_file <- args[2]

# --- Parse optional flags ---
get_flag_value <- function(args, flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  if (idx[1] == length(args)) {
    cat("Error: Missing value for", flag, "\n")
    quit(status = 1)
  }
  args[idx[1] + 1]
}

patient_id <- get_flag_value(args, "--patient")
sample_id  <- get_flag_value(args, "--sample_id")

seed_val <- get_flag_value(args, "--seed")
if (!is.null(seed_val)) {
  if (!grepl("^[0-9]+$", seed_val)) {
    cat("Error: --seed must be an integer\n")
    quit(status = 1)
  }
  set.seed(as.integer(seed_val))
}

# Prevent conflicting filters (optional but safer)
if (!is.null(patient_id) && !is.null(sample_id)) {
  cat("Error: Provide only one of --patient or --sample_id (not both)\n")
  quit(status = 1)
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

# Clear environment and keep args vars we need
rm(list = ls()[!ls() %in% c("input_file", "output_file", "patient_id", "sample_id", "seed_val")])
date <- Sys.Date()

cat("Processing CNV data...\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("Date:", as.character(date), "\n")
if (!is.null(patient_id)) cat("Filter patient:", patient_id, "\n")
if (!is.null(sample_id))  cat("Filter sample_id:", sample_id, "\n")
if (!is.null(seed_val))   cat("Seed:", seed_val, "\n")
cat("\n")

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

# --- Filter on patient or sample_id (NEW) ---
if (!is.null(patient_id)) {
  if (!"patient" %in% names(data)) {
    cat("Error: --patient provided but column 'patient' not found in input.\n")
    cat("Available columns:", paste(names(data), collapse = ", "), "\n")
    quit(status = 1)
  }
  before_n <- nrow(data)
  data <- data %>% filter(.data$patient == patient_id)
  after_n <- nrow(data)
  cat("Filtered by patient:", patient_id, "->", after_n, "rows (from", before_n, ")\n")
  if (after_n == 0) {
    cat("Error: No rows found for patient =", patient_id, "\n")
    quit(status = 1)
  }
}

if (!is.null(sample_id)) {
  if (!"sample_id" %in% names(data)) {
    cat("Error: --sample_id provided but column 'sample_id' not found in input.\n")
    cat("Available columns:", paste(names(data), collapse = ", "), "\n")
    quit(status = 1)
  }
  before_n <- nrow(data)
  data <- data %>% filter(.data$sample_id == sample_id)
  after_n <- nrow(data)
  cat("Filtered by sample_id:", sample_id, "->", after_n, "rows (from", before_n, ")\n")
  if (after_n == 0) {
    cat("Error: No rows found for sample_id =", sample_id, "\n")
    quit(status = 1)
  }
}

# ------------------------------------------------------------------------------
# Save filtered data
# ------------------------------------------------------------------------------

cat("\nSaving filtered data...\n")

# Ensure output directory exists
outdir <- dirname(output_file)
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# Write CSV (keep comma-separated, with header)
# Use readr::write_csv for speed and consistent quoting
tryCatch({
  readr::write_csv(data, output_file)
  cat("Saved", nrow(data), "rows and", ncol(data), "columns to:", output_file, "\n")
}, error = function(e) {
  cat("Error writing output file:", e$message, "\n")
  quit(status = 1)
})

# Optional: quick sanity summary
cat("\nSummary:\n")
if ("patient" %in% names(data)) {
  cat("  Unique patients:", length(unique(data$patient)), "\n")
}
if ("sample_id" %in% names(data)) {
  cat("  Unique sample_id:", length(unique(data$sample_id)), "\n")
}
cat("Done.\n")

