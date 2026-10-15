#!/usr/bin/env bash

###############################################################################
# submit.sh
#
# Example SLURM submission script used for ABN analyses.
#
# Users will likely need to adapt:
#   - partition
#   - account
#   - R installation
#   - resource allocation
#
###############################################################################

set -euo pipefail

#SBATCH --job-name=abn_fit
#SBATCH --time=00-120:00:00
#SBATCH --constraint=rhel8
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
# Replace with local HPC settings:
#SBATCH --partition=your_partition
#SBATCH --account=your_account


# Project directory
export PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Optional local R installation
export R_BIN="${PROJECT_DIR}/R/bin"

# Run information
export START_STAMP="$(date +'%Y-%m-%d_%H%M%S')"
export RUN_NAME="abn_fit"

# Parallelisation
export OMP_NUM_THREADS=1

# MCMC settings
export MCMC_SEEDS="081514,123456,2025,3141"
export MCMC_SCHEME="100000,0,0"
export BURNIN_LENGTH=20000
export THINNING_STEPS=2

# Load HPC environment
source "${PROJECT_DIR}/hpc/load_env.sh"


# Output directory
mkdir -p "${PROJECT_DIR}/results"

# Log file
ABN_OUTPUT="${PROJECT_DIR}/results/${RUN_NAME}_${START_STAMP}.log"

echo "========================================"
echo "R analysis"
echo "========================================"
echo "Start time:    $(date)"
echo "Project:       ${PROJECT_DIR}"
echo "R:             ${R_BIN}/Rscript"
echo "R library:     ${R_LIBS_USER}"
echo "Run name:      ${RUN_NAME}"
echo "========================================"

{
    echo
    echo "--------------------------"
    echo "Running analysis: ${RUN_NAME}"
    echo "--------------------------"
    echo "Start time: $(date)"
    echo

    srun "${R_BIN}/Rscript" --vanilla \
        "${PROJECT_DIR}/hpc/scripts/01_RIA_me.R"

    echo
    echo "Analysis completed: $(date)"
} > "${ABN_OUTPUT}" 2>&1

echo "Finished: $(date)"
echo "Log: ${ABN_OUTPUT}"