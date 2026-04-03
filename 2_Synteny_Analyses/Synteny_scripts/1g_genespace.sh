#!/bin/bash
#SBATCH --job-name=genespace
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

source myconda
mamba activate genespace_py3.10

OUTDIR="${1%/}"

echo "Bash received OUTDIR: ${OUTDIR}"

SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts/1h_genespace.R

Rscript "${SCRIPT}" "${OUTDIR}/"