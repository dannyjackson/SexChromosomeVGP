#!/bin/bash
#SBATCH --job-name=genespace_liftover
#SBATCH --output=slurm_output/%x_%A.out
#SBATCH --error=slurm_output/%x_%A.err
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

REF="${1:?Usage: sbatch $0 <REF_GENOME>}"

PAIRFILE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv
SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts/1c_genespace.wholegenome.r
SEXCHRSCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts/1d_genespace.wholegenome.sexchr.r
BASE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace
INDIR=${BASE}/${REF}/sexshared_liftover
OUTDIR=${BASE}/${REF}/sexshared_plots

mkdir -p "${OUTDIR}" slurm_output

source myconda
mamba activate genespace_py3.10

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

shopt -s nullglob

for WORKDIR in "${INDIR}"/*; do
    [ -d "${WORKDIR}" ] || continue

    SPECIES=$(basename "${WORKDIR}")

    SEXCHR=$(awk -F',' -v s="${SPECIES}" '$1 == s {print $2; exit}' "${PAIRFILE}")

    if [ -z "${SEXCHR}" ]; then
        echo "Skipping ${SPECIES}: no sex chromosome found in ${PAIRFILE}" >&2
        continue
    fi

    echo "Running ${SPECIES} (sex chromosome: ${SEXCHR})"

    export WORKING_DIR="${WORKDIR}"

    cd "${OUTDIR}"

    Rscript "${SCRIPT}" "${SPECIES}" "${OUTDIR}/" "${REF}"

    r="${WORKDIR}/orthofinder"
    out="${r}.tar.gz"
    if [ -d "${r}" ]; then
        tar -C "$(dirname "${r}")" -cf - "$(basename "${r}")" | bgzip -c > "${out}" && rm -rf "${r}"
    fi

    Rscript "${SEXCHRSCRIPT}" "${SPECIES}" "${WORKDIR}/" "${SEXCHR}" "${REF}"
done