#!/bin/bash
#SBATCH --job-name=genespace
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

REF="${1:?Usage: sbatch genespace.sh <REF_GENOME> <SPECIES>}"
SPECIES="${2:?Usage: sbatch genespace.sh <REF_GENOME> <SPECIES>}"

SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts/1c_genespace.wholegenome.r

WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/${SPECIES}

export WORKING_DIR="${WORKDIR}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

source myconda
mamba activate genespace

Rscript "${SCRIPT}" "${SPECIES}" "${WORKDIR}/" "${REF}"