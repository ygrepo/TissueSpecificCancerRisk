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
# 
# # Load required modules
module purge
module load R/4.4.3

# Create logs and output directories
mkdir -p logs output/figures

echo "Starting heatmap generation..."

# Run the heatmap script with command-line arguments
# Rscript src/make_heatmap.R \
#     --cnv_data "data/allele_specific_cn/B2HET16-hscn.csv" \
#     --use_umap_clusters \
#     --output "output/figures/B216_heatmap_chr13_17q_del_umap.pdf" \
#     --plot_tree \
#     --width 8 \
#     --height 6 \
#     --title "B216 Clonal Evolution chr13/17q Deletions" \
#     --tree_width 2.5 \
#     --linkheight 2

# SA501
Rscript src/make_heatmap.R \
    --treefile "trees/SA501/SA501.tree.2.newick" \
    --cnv_data "data/SA501.tbnc.cnv.csv" \
    --output "output/figures/SA501_heatmap_chr13_17q_del.pdf" \
    --chr13_17_deletion \
    --chr13_17_deletion_threshold 0.5 \
    --plot_tree \
    --width 8 \
    --height 6 \
    --title "SA501 TBC Clonal Evolution chr13/17q Deletions" \
    --tree_width 3 \
    --linkheight 1
     
# Rscript src/make_heatmap.R \
#     --treefile "trees/B218/tree.newick" \
#     --cnv_data "data/B2HET18-hscn.csv" \
#     --output "output/figures/B218_heatmap_chr13_17q_del.pdf" \
#     --chr13_17_deletion \
#     --chr13_17_deletion_threshold 0.25 \
#     --plot_tree \
#     --width 8 \
#     --height 6 \
#     --title "B218 Clonal Evolution chr13/17q Deletions" \
#     --tree_width 3 \
#     --linkheight 1

# Rscript src/make_heatmap.R \
#     --treefile "trees/B216/tree.newick" \
#     --cnv_data "data/allele_specific_cn/B2HET16-hscn.csv" \
#     --output "output/figures/B216_heatmap_chr13_17q_del_2.pdf" \
#     --chr13_17_deletion \
#    --chr13_17_deletion_threshold 0.25 \
#     --plot_tree \
#     --width 8 \
#     --height 6 \
#     --title "B216 Clonal Evolution chr13/17 Deletions" \
#     --tree_width 3 \
#     --linkheight 1

echo "Heatmap generation completed."
echo "Job finished."