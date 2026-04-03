#!/bin/bash
#SBATCH --job-name=gs_chicken_wg
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

PAIRFILE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv
SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/scripts/genespace.synteny_to_chicken.wholegenome.sexchr.r
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/plots_sexchr
SLURM_ARRAY_TASK_ID=20

# Grab the Nth line (array index) and parse "Species Chromosome"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PAIRFILE}")
SPECIES=$(awk '{print $1}' <<< "${line}")
SEXCHR=$(awk '{print $2}' <<< "${line}")
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/${SPECIES}

export WORKING_DIR="${WORKDIR}"

cd "${OUTDIR}"

source myconda
mamba activate genespace

Rscript "${SCRIPT}" "${SPECIES}" "${OUTDIR}/" ${SEXCHR}