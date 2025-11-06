#!/bin/bash

# --- LSF Job Options ---
#BSUB -P acc_DiseaseGeneCell
#BSUB -J sitka_filter      # Job name
#BSUB -n 1                    # Request 1 CPU core
#BSUB -W 2:00                 # Walltime of 2 hours
#BSUB -M 32000                # Request 32 GB of memory
#BSUB -o logs/filter.%J.out    # Standard output log file
#BSUB -e logs/filter.%J.err    # Standard error log file
# --- End LSF Options ---

echo "Running corrupt-filter on host $(hostname)..."

PROJECT_DIR="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/TissueSpecificCancerRisk"
echo "Changing to working directory: ${PROJECT_DIR}"
cd "${PROJECT_DIR}"

# Load the required Java module
module purge
module load java/1.8.0_151


# Add the sitkatree binaries to your PATH
export SITKA_BIN_PATH="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/sitkatree/sitka/build/install/nowellpack/bin"
export PATH=$SITKA_BIN_PATH:$PATH

# Run the command
FN="${PROJECT_DIR}/output.csv"

corrupt-filter --input $FN --lowerFraction 0.05
cp results/latest/filtered.csv ./

echo "Job finished."