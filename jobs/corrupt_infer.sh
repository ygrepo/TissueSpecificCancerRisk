#!/bin/bash

# --- LSF Job Options ---
#BSUB -P acc_DiseaseGeneCell
#BSUB -J sitka_infer           # Job name
#BSUB -n 10                   # Request 12 CPU cores
#BSUB -W 24:00                # Walltime of 24 hours
#BSUB -M 64000                # Request 64 GB of memory
#BSUB -o logs/infer.%J.out         # Standard output log file
#BSUB -e logs/infer.%J.err         # Standard error log file
# --- End LSF Options ---

echo "Running corrupt-infer on host $(hostname)..."

# Load the required Java module
module purge
module load java/1.8.0_151

# Add the sitkatree binaries to your PATH
export SITKA_BIN_PATH="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/sitkatree/sitka/build/install/nowellpack/bin"
export PATH=$SITKA_BIN_PATH:$PATH

# Run the main inference command
# We set --engine.nChains to 4 to match our '-n 4' CPU request
FILTERED_PATRH=/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/TissueSpecificCancerRisk/results/all/2025-11-05-10-32-39-xmCj8wy3.exec/filtered.csv
corrupt-infer-with-noisy-params \
    --model.globalParameterization true \
    --model.binaryMatrix $FILTERED_PATRH \
    --model.fprBound 0.1 \
    --model.fnrBound 0.5 \
    --engine PT \
    --engine.initialization FORWARD \
    --engine.nScans 10000 \
    --engine.nPassesPerScan 1 \
    --engine.nChains 10
 cp results/latest/samples/phylo.csv ./
echo "Job finished."
