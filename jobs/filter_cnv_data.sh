#!/bin/bash

# --- LSF Job Options ---
#BSUB -P acc_DiseaseGeneCell
#BSUB -J read_cnv_data           # Job name
#BSUB -n 4                      # Request 4 CPU cores
#BSUB -W 2:00                   # Walltime of 2 hours
#BSUB -M 32000                  # Request 32 GB of memory
#BSUB -R "rusage[mem=32000]"
#BSUB -R "span[hosts=1]"
#BSUB -o logs/read_cnv_data.%J.out    # Standard output log file
#BSUB -e logs/read_cnv_data.%J.err    # Standard error log file
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
input_file="data/genotype.data.cnv.csv"
output_file="data/SA039.geno.cnv.csv"
patient_id="SA039"
Rscript src/read_cnv_data.R "${input_file}" "${output_file}" --patient "${patient_id}"

echo "Job finished."