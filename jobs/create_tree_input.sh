#!/bin/bash

# --- LSF Job Options ---
#BSUB -P acc_DiseaseGeneCell
#BSUB -J create_tree_input           # Job name
#BSUB -n 4                      # Request 4 CPU cores
#BSUB -W 2:00                   # Walltime of 2 hours
#BSUB -M 32000                  # Request 32 GB of memory
#BSUB -o logs/create_tree_input.%J.out    # Standard output log file
#BSUB -e logs/create_tree_input.%J.err    # Standard error log file
# --- End LSF Options ---

echo "Running on create_tree_input host $(hostname)..."

set -e

PROJECT_DIR="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/TissueSpecificCancerRisk"
echo "Changing to working directory: ${PROJECT_DIR}"
cd "${PROJECT_DIR}"
# 
# # Load required modules
module purge
module load R/4.4.3


echo "Starting..."

Rscript src/cnv_processor.R data/genotype.data.cnv.csv data/SA039.tbnc.cnv.csv --patient SA039

echo "Job finished."