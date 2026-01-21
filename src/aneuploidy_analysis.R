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
# DATA LOADING & PREPROCESSING
# ==============================================================================

#' Load and standardize CNV data
#' @param filepath Path to the CSV file
#' @param col_map Named vector to rename columns to standard internal names (chr, start, end, state, sample_id, cell_id)
load_and_preprocess_cnv <- function(filepath, col_map = NULL, drop_filtered = TRUE) {
  dt <- fread(filepath, stringsAsFactors = FALSE)
  
  # --------------------------
  # Optional user mapping
  # --------------------------
  if (!is.null(col_map)) {
    setnames(dt, old = names(col_map), new = col_map, skip_absent = TRUE)
  }
  
  # --------------------------
  # Auto-standardize schema
  # --------------------------
  # sample_id
  if (!"sample_id" %in% names(dt)) {
    if ("sample" %in% names(dt)) {
      setnames(dt, "sample", "sample_id")
    } else if ("Sample" %in% names(dt)) {
      setnames(dt, "Sample", "sample_id")
    }
  }
  
  # cell_id
  if (!"cell_id" %in% names(dt)) {
    if ("cell" %in% names(dt)) setnames(dt, "cell", "cell_id")
    if ("barcode" %in% names(dt)) setnames(dt, "barcode", "cell_id")
  }
  
  # state (copy number state)
  # Prefer "state" if present; otherwise use "copy"
  if (!"state" %in% names(dt)) {
    if ("copy" %in% names(dt)) {
      setnames(dt, "copy", "state")
    } else if ("cn" %in% names(dt)) {
      setnames(dt, "cn", "state")
    }
  }
  
  # start/end
  if (!"start" %in% names(dt) && "Start" %in% names(dt)) setnames(dt, "Start", "start")
  if (!"end"   %in% names(dt) && "End"   %in% names(dt)) setnames(dt, "End",   "end")
  
  # --------------------------
  # Basic validation
  # --------------------------
  req <- c("sample_id", "cell_id", "chr", "start", "end", "state")
  missing <- setdiff(req, names(dt))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # --------------------------
  # Optional: remove filtered segments
  # --------------------------
  if (drop_filtered && "filt" %in% names(dt)) {
    # "filt" often is TRUE for "remove", FALSE for "keep"
    dt <- dt[is.na(filt) | filt == FALSE | filt == 0]
  }
  
  # --------------------------
  # Type coercion + chr normalization
  # --------------------------
  dt[, chr := gsub("^chr", "", as.character(chr))]
  dt[, chr_int := suppressWarnings(as.integer(chr))]
  
  dt[, start := as.numeric(start)]
  dt[, end   := as.numeric(end)]
  dt[, state := as.numeric(state)]
  
  # --------------------------
  # Segment length
  # --------------------------
  if (!"seg_len" %in% names(dt)) {
    if ("segment_size_bp" %in% names(dt)) {
      dt[, seg_len := as.numeric(segment_size_bp)]
    } else {
      dt[, seg_len := as.numeric(end - start + 1)]
    }
  }
  
  # Remove rows with unusable segments
  dt <- dt[!is.na(start) & !is.na(end) & !is.na(state) & seg_len > 0]
  
  return(dt)
}

# ==============================================================================
# METRIC CALCULATIONS
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
# STATISTICAL ANALYSES
# ==============================================================================

#' Run Wilcoxon tests within samples
run_within_sample_tests <- function(metrics_df, group_col, metric_cols) {
  
  # --- Robust rank-biserial (Cliff-like), always in [-1, 1] ---
  # Positive => deleted group tends to have larger values than no_del group
  rank_biserial <- function(x, g) {
    keep <- !is.na(x) & !is.na(g)
    x <- x[keep]
    g <- as.logical(g[keep])
    
    x_del <- x[g]
    x_no  <- x[!g]
    
    n1 <- length(x_del)
    n0 <- length(x_no)
    if (n1 < 2 || n0 < 2) return(NA_real_)
    
    # Compute U for deleted group from ranks (handles ties via average ranks)
    r <- rank(x, ties.method = "average")
    W <- sum(r[g])                          # rank-sum for deleted group
    U <- W - n1 * (n1 + 1) / 2              # Mann–Whitney U for deleted group
    
    # Rank-biserial correlation
    # Range is exactly [-1, 1]
    rb <- (2 * U) / (n1 * n0) - 1
    return(as.numeric(rb))
  }
  
  
  results <- rbindlist(lapply(split(metrics_df, by = "sample_id"), function(df) {
    df <- as.data.table(df)
    out <- lapply(metric_cols, function(m) {
      x <- df[[m]]
      g <- df[[group_col]]
      
      keep <- !is.na(x) & !is.na(g)
      x <- x[keep]
      g <- as.logical(g[keep])
      
      n_del <- sum(g)
      n_no  <- sum(!g)
      if (n_del < 3 || n_no < 3) return(NULL)

      # Wilcoxon test (deleted vs no_del)
      wt <- suppressWarnings(wilcox.test(x[g], x[!g], exact = FALSE))
      
      data.table(
        sample_id       = df$sample_id[1],
        metric          = m,
        n_group_del     = n_del,
        n_group_no_del  = n_no,
        med_del         = median(x[g],  na.rm = TRUE),
        med_no_del      = median(x[!g], na.rm = TRUE),
        p_value         = wt$p.value,
        rank_biserial   = rank_biserial(x, g),
        effect_direction = {
          rb <- rank_biserial(x, g)
          if (is.na(rb)) NA_integer_
          else if (rb > 0) 1L
          else if (rb < 0) -1L
          else 0L
        }
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

pick_threshold <- function(frac_del, candidates = c(0.5, 0.25, 0.2, 0.1), min_true = 3) {
  for (t in candidates) {
    n_true <- sum(frac_del >= t, na.rm = TRUE)
    if (n_true >= min_true) return(t)
  }
  return(tail(candidates, 1))  # fallback to loosest threshold
}


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================


# MAPPING: Map the column names in your screenshot to the internal names
column_mapping <- c(
  "sample" = "sample_id", 
  "cell_id" = "cell_id",
  "chr" = "chr",
  "start" = "start",
  "end" = "end",
  "state" = "state"
)


# Load
#input_file <- here("data","SA501.tbnc.cnv.csv")
input_file <- here("data/allele_specific_cn","B2HET16-hscn.csv")
dt <- load_and_preprocess_cnv(input_file, col_map = column_mapping)
print(names(dt))
print(dt[1:3, .(sample_id, cell_id, chr, start, end, state, seg_len)])

# Baseline
cell_baseline <- calculate_baseline_ploidy(dt)

# Metrics
cell_metrics <- calc_genome_metrics(dt, cell_baseline)

# chr17_frac <- dt[chr == "17",
#                  .(frac_del = sum(seg_len[state <= 1], na.rm = TRUE) / sum(seg_len, na.rm = TRUE)),
#                  by = .(sample_id, cell_id)]
# 
# threshold <- pick_threshold(chr17_frac$frac_del)
# cat("Auto-picked chr17 threshold:", threshold, "\n")

threshold <- 0.25
# Call Deletions (13q and 17q)


cell_metrics <- call_chr_deletion(dt, cell_metrics, 
                                  target_chr = 13, 
                                  col_name_prefix = "chr13",
                                  del_cutoff = 1,
                                  thresh = threshold)

cell_metrics <- call_chr_deletion(dt, cell_metrics, 
                                  target_chr = 17, 
                                  col_name_prefix = "chr17",
                                  del_cutoff = 1,
                                  thresh = threshold)

# 5. Run Statistics (Within Sample)
stats_results <- run_within_sample_tests(
  cell_metrics, 
  group_col = "chr13_del", 
  metric_cols = c("mean_cn", "frac_cn_ge4", "frac_altered", "cn_mad", "frac_loss", "frac_gain")
)

# View Results
print("Statistics Results for chr13_del:")
print(stats_results)

# Save
#sample_id <- "SA501"
sample_id <- unique(cell_metrics$sample_id)

fwrite(
  stats_results,
  file = here("output/data", paste0(sample_id, "_chr13_del_within_sample_results.csv"))
)

# Run Statistics (Within Sample)
stats_results <- run_within_sample_tests(
  cell_metrics, 
  group_col = "chr17_del", 
  metric_cols = c("mean_cn", "frac_cn_ge4", "frac_altered", "cn_mad", "frac_loss", "frac_gain")
)

print("Statistics Results for chr17_del:")
print(stats_results)
fwrite(
  stats_results,
  file = here("output/data", paste0(sample_id, "_chr17_del_within_sample_results.csv"))
)

#  Global Linear Model (Across Samples)
if (length(unique(cell_metrics$chr13_del)) > 1 && length(unique(cell_metrics$WGD)) > 1) {
  m <- lm(frac_altered ~ chr13_del + WGD + chr17_del, data = cell_metrics)
  print(summary(m))
  print(coeftest(m, vcov = vcovHC(m, type = "HC3")))
}

# Check specific genes
check_gene_loci(dt, cell_metrics)

# Unadjusted median deltas (from your Wilcoxon stats output)
chr17_tbl <- fread(here("output/data", 
                        paste0(sample_id, "_chr17_del_within_sample_results.csv"))) %>%
  filter(metric == "frac_altered") %>%
  transmute(
    driver = "Chr17q Deletion",
    unadj_delta_median = med_del - med_no_del
  )

chr13_tbl <- fread(here("output/data", 
                        paste0(sample_id, "_chr13_del_within_sample_results.csv"))) %>%
  filter(metric == "frac_altered") %>%
  transmute(
    driver = "Chr13q Deletion",
    unadj_delta_median = med_del - med_no_del
  )

# Adjusted marginal effects from the multivariable model
avg_marginal_effect <- function(model, data, var) {
  d1 <- copy(as.data.table(data))
  d0 <- copy(as.data.table(data))
  
  d1[[var]] <- TRUE
  d0[[var]] <- FALSE
  
  mean(predict(model, newdata = d1) - predict(model, newdata = d0), na.rm = TRUE)
}

adj_tbl <- data.frame(
  driver = c("Chr17q Deletion", "Chr13q Deletion", "Whole Genome Doubling"),
  adjusted_effect = c(
    avg_marginal_effect(m, cell_metrics, "chr17_del"),
    avg_marginal_effect(m, cell_metrics, "chr13_del"),
    avg_marginal_effect(m, cell_metrics, "WGD")
  )
)

bridge_tbl <- bind_rows(chr17_tbl, chr13_tbl) %>%
  left_join(adj_tbl, by = "driver")

print(bridge_tbl)


ct <- coeftest(m, vcov = vcovHC(m, type = "HC3"))

plot_data <- data.frame(
  term = rownames(ct),
  estimate = ct[, 1],
  std.error = ct[, 2],
  p.value = ct[, 4]
) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    conf.low = estimate - 1.96 * std.error,
    conf.high = estimate + 1.96 * std.error,
    term_label = case_when(
      term == "WGDTRUE" ~ "Whole Genome Doubling",
      term == "chr13_delTRUE" ~ "Chr13q Deletion",
      term == "chr17_delTRUE" ~ "Chr17q Deletion",
      TRUE ~ term
    )
  )


p <- ggplot(plot_data, aes(x = reorder(term_label, -estimate), y = estimate)) +
  
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
    y = paste0(sample_id, ": Adjusted Effect on Fraction Altered"),
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

ggsave(here("output/figures", 
            paste0(sample_id, "_genomic_instability_drivers.png")), 
            plot = p, width = 7, height = 4, dpi = 300)


