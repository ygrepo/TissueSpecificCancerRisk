library(data.table)
library(dplyr)
library(broom)
library(lme4)
library(matrixStats)
library(sandwich)
library(lmtest)
library(ggplot2)
library(dplyr)
library(here)

rm(list = ls())


# ==============================================================================
# 1. DATA LOADING & PREPROCESSING
# ==============================================================================

#' Load and standardize CNV data
#' @param filepath Path to the CSV file
#' @param col_map Named vector to rename columns to standard internal names (chr, start, end, state, sample_id, cell_id)
load_and_preprocess_cnv <- function(filepath, col_map = NULL) {
  
  raw_data <- fread(filepath, stringsAsFactors = FALSE)
  
  # Rename columns based on input map (e.g., convert "sample" to "sample_id")
  if (!is.null(col_map)) {
    setnames(raw_data, old = names(col_map), new = col_map, skip_absent = TRUE)
  }
  
  # Standardize Chromosome names (remove "chr" prefix)
  raw_data[, chr := gsub("^chr", "", as.character(chr))]
  
  # Create integer chromosome column for filtering autosomes (1-22)
  # Suppress warnings for X/Y becoming NA during integer conversion
  raw_data[, chr_int := suppressWarnings(as.integer(chr))]
  
  # Calculate segment length if not present (Crucial for your new file format)
  if (!"seg_len" %in% names(raw_data)) {
    if (all(c("end", "start") %in% names(raw_data))) {
      raw_data[, seg_len := as.numeric(end - start + 1)]
    } else if ("segment_size_bp" %in% names(raw_data)) {
      raw_data[, seg_len := as.numeric(segment_size_bp)]
    } else {
      stop("Cannot calculate segment length: missing start/end or segment_size_bp columns")
    }
  }
  
  return(raw_data)
}

# ==============================================================================
# 2. METRIC CALCULATIONS
# ==============================================================================

#' Calculate baseline ploidy per cell
calculate_baseline_ploidy <- function(dt) {
  # Use autosomes only (1-22)
  baseline <- dt[
    !is.na(chr_int) & chr_int >= 1 & chr_int <= 22,
    .(baseline_cn = weightedMedian(state, w = seg_len, na.rm = TRUE)),
    by = .(sample_id, cell_id)
  ]
  return(baseline)
}

#' Compute Genome-Wide Metrics
calc_genome_metrics <- function(dt, baseline_dt) {
  
  # Join baseline info
  setkey(baseline_dt, sample_id, cell_id)
  setkey(dt, sample_id, cell_id)
  dt_merged <- baseline_dt[dt]
  
  # Constants
  HIGH_CUTOFF <- 4
  
  # Calculate metrics on Autosomes
  metrics <- dt_merged[!is.na(chr_int) & chr_int >= 1 & chr_int <= 22, .(
    # (A) Ploidy / endopolyploidy proxies
    mean_cn     = weighted.mean(state, w = seg_len, na.rm = TRUE),
    frac_cn_ge4 = sum(seg_len[state >= HIGH_CUTOFF], na.rm = TRUE) / sum(seg_len, na.rm = TRUE),
    
    # (B) CNV burden relative to baseline
    frac_altered = sum(seg_len[abs(state - baseline_cn) >= 1], na.rm = TRUE) / sum(seg_len, na.rm = TRUE),
    cn_mad       = weightedMedian(abs(state - baseline_cn), w = seg_len),
    
    # (C) Directional burdens
    frac_loss = sum(seg_len[state <= (baseline_cn - 1)], na.rm = TRUE) / sum(seg_len, na.rm = TRUE),
    frac_gain = sum(seg_len[state >= (baseline_cn + 1)], na.rm = TRUE) / sum(seg_len, na.rm = TRUE)
  ), by = .(sample_id, cell_id, baseline_cn)] # Add other grouping vars here if they exist in input (e.g. patient)
  
  # Add WGD status
  metrics[, WGD := baseline_cn >= 3.5]
  
  return(metrics)
}

#' Determine specific chromosome deletion status (e.g., 13q, 17)
call_chr_deletion <- function(dt, metrics_df, target_chr, col_name_prefix, del_cutoff=1, thresh=0.5) {
  
  chr_stat <- dt[chr == as.character(target_chr),
                 .(frac_del = sum(seg_len[state <= del_cutoff], na.rm = TRUE) / sum(seg_len, na.rm = TRUE)),
                 by = .(sample_id, cell_id)]
  
  chr_stat[, is_deleted := (frac_del >= thresh)]
  
  # Merge into main metrics table
  setkey(metrics_df, sample_id, cell_id)
  setkey(chr_stat, sample_id, cell_id)
  
  metrics_df <- chr_stat[metrics_df]
  
  # Rename columns dynamically
  setnames(metrics_df, c("frac_del", "is_deleted"), 
           c(paste0("frac_", col_name_prefix, "_del"), paste0(col_name_prefix, "_del")))
  
  # Fill NAs (cells with no segments for that chromosome)
  del_col <- paste0(col_name_prefix, "_del")
  metrics_df[is.na(get(del_col)), (del_col) := FALSE]
  
  return(metrics_df)
}

# ==============================================================================
# 3. STATISTICAL ANALYSES
# ==============================================================================

#' Run Wilcoxon tests within samples
run_within_sample_tests <- function(metrics_df, group_col, metric_cols) {
  
  rank_biserial <- function(x, g) {
    x1 <- x[g]; x0 <- x[!g]
    if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
    wt <- suppressWarnings(wilcox.test(x1, x0, exact = FALSE))
    W <- unname(wt$statistic)
    n1 <- length(x1); n0 <- length(x0)
    U <- W - n1*(n1+1)/2
    2*U/(n1*n0) - 1
  }
  
  results <- rbindlist(lapply(split(metrics_df, by = "sample_id"), function(df) {
    df <- as.data.table(df)
    out <- lapply(metric_cols, function(m) {
      x <- df[[m]]
      g <- df[[group_col]]
      
      # Skip if not enough data points in both groups
      if (sum(g, na.rm=TRUE) < 3 || sum(!g, na.rm=TRUE) < 3) return(NULL)
      
      wt <- suppressWarnings(wilcox.test(x ~ g, exact = FALSE))
      data.table(
        sample_id = df$sample_id[1],
        metric = m,
        n_group_TRUE = sum(g, na.rm=TRUE),
        n_group_FALSE = sum(!g, na.rm=TRUE),
        med_TRUE = median(x[g], na.rm = TRUE),
        med_FALSE = median(x[!g], na.rm = TRUE),
        p_value = wt$p.value,
        rank_biserial = rank_biserial(x, g)
      )
    })
    rbindlist(out, fill = TRUE)
  }), fill = TRUE)
  
  if(nrow(results) > 0) {
    results[, q_value := p.adjust(p_value, method = "BH")]
  }
  return(results)
}

#' Check Specific Gene Loci (BRCA2 / RB1)
check_gene_loci <- function(dt, metrics_df) {
  # Add 13q status to the main data for filtering
  dt_merged <- merge(dt, metrics_df[, .(sample_id, cell_id, chr13_del)], by = c("sample_id", "cell_id"))
  
  # Define bins
  brca2_bins <- dt_merged[chr == "13" & start <= 32500000 & end >= 32300000]
  rb1_bins   <- dt_merged[chr == "13" & start <= 49000000 & end >= 48800000]
  
  print("BRCA2 State Table for 13q Del Cells:")
  print(brca2_bins[chr13_del == TRUE, table(state)])
  
  print("RB1 State Table for 13q Del Cells:")
  print(rb1_bins[chr13_del == TRUE, table(state)])
  
  # Statistical test for concordance
  gene_cn_summary <- dt_merged[
    chr == "13" & ((start <= 32500000 & end >= 32300000) | (start <= 49000000 & end >= 48800000)),
    .(gene_cn = median(state)),
    by = .(cell_id, chr13_del)
  ]
  
  if (length(unique(gene_cn_summary$chr13_del)) > 1) {
    print(wilcox.test(gene_cn ~ chr13_del, data = gene_cn_summary))
  }
}

# ==============================================================================
# 4. MAIN EXECUTION
# ==============================================================================

# Input File from your screenshot
#input_file <- here("data", "Screenshot_File_Name.csv") # UPDATE THIS PATH

input_file <- here("data","SA501.tbnc.cnv.csv")

# MAPPING: Map the column names in your screenshot to the internal names
# Screenshot format: "sample" -> "sample_id", "chr" -> "chr", etc.
# Note: "segment_size_bp" is missing in screenshot, so code calculates it via end-start
column_mapping <- c(
  "sample" = "sample_id", 
  "cell_id" = "cell_id",
  "chr" = "chr",
  "start" = "start",
  "end" = "end",
  "state" = "state"
)

# 1. Load
dt <- load_and_preprocess_cnv(input_file, col_map = column_mapping)

# 2. Baseline
cell_baseline <- calculate_baseline_ploidy(dt)

# 3. Metrics
cell_metrics <- calc_genome_metrics(dt, cell_baseline)

# 4. Call Deletions (13q and 17)
#call_chr_deletion <- function(dt, metrics_df, target_chr, col_name_prefix, del_cutoff=1, thresh=0.5) {
  
cell_metrics <- call_chr_deletion(dt, cell_metrics, 
                                  target_chr = 13, 
                                  col_name_prefix = "chr13",
                                  del_cutoff = 1,
                                  thresh = 0.5)
cell_metrics <- call_chr_deletion(dt, cell_metrics, 
                                  target_chr = 17, 
                                  col_name_prefix = "chr17",
                                  del_cutoff = 1,
                                  thresh = 0.5)

# 5. Run Statistics (Within Sample)
stats_results <- run_within_sample_tests(
  cell_metrics, 
  group_col = "chr13_del", 
  metric_cols = c("mean_cn", "frac_cn_ge4", "frac_altered", "cn_mad", "frac_loss", "frac_gain")
)

# View Results
print("Statistics Results for chr13_del:")
print(head(stats_results[order(q_value)]))

# 8. Save
fwrite(stats_results[order(q_value)], file = here("output/data", "chr13_del_within_sample_results.csv"))

# 5. Run Statistics (Within Sample)
stats_results <- run_within_sample_tests(
  cell_metrics, 
  group_col = "chr17_del", 
  metric_cols = c("mean_cn", "frac_cn_ge4", "frac_altered", "cn_mad", "frac_loss", "frac_gain")
)

# View Results
print("Statistics Results for chr17_del:")
print(head(stats_results[order(q_value)]))
fwrite(stats_results[order(q_value)], file = here("output/data", "chr17_del_within_sample_results.csv"))

# 6. Global Linear Model (Across Samples)
# Ensure you have enough samples/variance for this to work
if (length(unique(cell_metrics$chr13_del)) > 1 && length(unique(cell_metrics$WGD)) > 1) {
  m <- lm(frac_altered ~ chr13_del + WGD + chr17_del, data = cell_metrics)
  print(summary(m))
  print(coeftest(m, vcov = vcovHC(m, type = "HC3")))
}

# 7. Check specific genes
check_gene_loci(dt, cell_metrics)


# 1. Tidy the model results and get Confidence Intervals
# (Use the model object 'm' from your previous step)
plot_data <- tidy(m, conf.int = TRUE) %>%
  filter(term != "(Intercept)") %>%  # Remove intercept to focus on drivers
  mutate(
    # Rename terms for the slide (No more "TRUE" suffixes)
    term_label = case_when(
      term == "WGDTRUE" ~ "Whole Genome Doubling",
      term == "chr13_delTRUE" ~ "Chr13q Deletion",
      term == "chr17_delTRUE" ~ "Chr17q Deletion",
      TRUE ~ term
    )
  )

p <- ggplot(plot_data, aes(x = reorder(term_label, estimate), y = estimate)) +
  
  # The zero line
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  
  # The Confidence Intervals
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, size = 0.8) +
  
  # The Estimates
  geom_point(size = 4, color = "#2E86C1") +
  
  # Labels
  coord_flip() +
  labs(
    title = "Fraction Altered",
    #subtitle = "Whole Genome Doubling has the largest impact on genomic alteration",
    y = "Effect Size (Increase in Fraction Altered)",
    x = ""
  ) +
  
  # Theme Settings
  theme_bw(base_size = 14) +  # Sets the base font size for axes/ticks
  theme(
    panel.grid.minor = element_blank(),
    
    # --- FIX IS HERE ---
    plot.title = element_text(size = 14, face = "bold"),   
    plot.subtitle = element_text(size = 12, color = "gray30") # Smaller subtitle
  )

ggsave(here("output/figures", "genomic_instability_drivers.png"), plot = p, width = 7, height = 4, dpi = 300)


