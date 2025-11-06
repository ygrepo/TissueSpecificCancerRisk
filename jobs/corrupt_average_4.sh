#!/bin/bash

#BSUB -P acc_DiseaseGeneCell
#BSUB -J sitka_average
#BSUB -n 1
#BSUB -W 2:00
#BSUB -M 32000
#BSUB -o logs/average.%J.out
#BSUB -e logs/average.%J.err

echo "Running corrupt-average on host $(hostname)..."

PROJECT_DIR="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/TissueSpecificCancerRisk"
echo "Changing to working directory: ${PROJECT_DIR}"
cd "${PROJECT_DIR}"

# Load the required Java module
module purge
module load java/1.8.0_151

# Add the sitkatree binaries to your PATH
export SITKA_BIN_PATH="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/sitkatree/sitka/build/install/nowellpack/bin"
export PATH=$SITKA_BIN_PATH:$PATH

mkdir -p logs

# Run the main inference command
FN="${PROJECT_DIR}/phylo.csv"
echo "Running"

# Run the command
corrupt-average --csvFile "${FN}" --logisticTransform false
cp results/latest/average.csv ./

echo "Job finished."