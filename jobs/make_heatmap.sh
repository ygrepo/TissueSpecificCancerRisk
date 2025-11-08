#!/bin/bash

# --- LSF Job Options ---
#BSUB -P acc_DiseaseGeneCell
#BSUB -J make_heatmap           # Job name
#BSUB -n 4                      # Request 4 CPU cores
#BSUB -W 2:00                   # Walltime of 2 hours
#BSUB -M 32000                  # Request 32 GB of memory
#BSUB -o logs/heatmap.%J.out    # Standard output log file
#BSUB -e logs/heatmap.%J.err    # Standard error log file
# --- End LSF Options ---

echo "Running heatmap plotting on host $(hostname)..."

set -e

PROJECT_DIR="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/TissueSpecificCancerRisk"
echo "Changing to working directory: ${PROJECT_DIR}"
cd "${PROJECT_DIR}"

# Load required modules
module purge
module load R/4.3.0

# Create logs and output directories
mkdir -p logs plots

echo "Starting heatmap generation..."

# Run the heatmap script with command-line arguments
Rscript src/make_heatmap.R \
    --treefile "trees/SA501/tree.newick" \
    --cnv_data "data/SA501.tbnc.cnv.csv" \
    --output "data/figures/SA501_heatmap.pdf" \
    --width 20 \
    --height 12 \
    --title "SA501 Copy Number Heatmap" \
    --tree_width 6 \
    --linkheight 8

echo "Heatmap generation completed."
echo "Job finished."