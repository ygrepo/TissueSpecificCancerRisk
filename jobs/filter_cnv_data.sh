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
mkdir -p logs


echo "Starting..."
input_file="data/hgsc.cnv.csv"
output_file="data/SA1052.hgsc.cnv.csv"
patient_id="SA1052"
Rscript src/read_cnv_data.R "${input_file}" "${output_file}" --patient "${patient_id}"

echo "Job finished."