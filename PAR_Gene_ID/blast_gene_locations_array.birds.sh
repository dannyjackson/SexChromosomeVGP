#!/bin/bash
#SBATCH --job-name=gene_blast
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

module load blast

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis

BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"
SPECIES_LIST="species_for_blast_array.txt"

mkdir -p slurm_output blast_dbs blast_results species_results

sp=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_LIST")

if [ -z "$sp" ]; then
    echo "ERROR: No species found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi
blastn -query ARK2C.fa -db blast_dbs/Mergus_octosetaceus -out test.ARK2C.tsv \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
    -max_target_seqs 5 \
    -evalue 1e-20 \
    -num_threads 4
genes=(
    ACAA2
    ALPK2
    ARK2C
    ARK2N
    ATP8A1
    C18orf32
    CCDC68
    CPLX4
    CTIF
    DCC
    DYM
    EEF2
    ELAC1
    EPG5
    FECH
    GPN1
    GRP
    HAUS1
    HDHD2
    IER3IP1
    KATNAL2
    LAS2
    LIPG
    LMAN1
    LOXHD1
    LUZP2
    MALT1
    MAPK4
    ME2
    MECP2
    MEX3C
    MIR122
    MYO5B
    NARS1
    NEDD4
    ONECUT2
    P53
    PIAS2
    PIK3C3
    POLI
    PSTPIP2
    RAB27B
    RAX
    RIT2
    Rx2
    SEC11C
    SETBP1
    SIGLEC15
    SKA1
    SKOR2
    SLC14A2
    SMAD2
    SMAD4
    SMAD7
    ST8SIA3
    STARD6
    SYT4
    TCF4
    TSPAN36
    TXNL1
    WDR7
    ZBTB7C
    ZNF532
)

genome="${BASE}/${sp}/${sp}.fna"
OUT="species_results/${sp}.gene_locations.csv"

echo "Species,Gene,Hit_rank,Chrom,Start_pos,Stop_pos" > "$OUT"

if [ ! -f "$genome" ]; then
    echo "WARNING: genome not found for ${sp}: ${genome}" >&2

    for gene in "${genes[@]}"; do
        echo "${sp},${gene},NA,NA,NA,NA" >> "$OUT"
    done

    exit 0
fi

db="blast_dbs/${sp}"

if [ ! -f "${db}.nin" ] && [ ! -f "${db}.00.nin" ]; then
    echo "Making BLAST DB for ${sp}"

    makeblastdb \
        -in "$genome" \
        -dbtype nucl \
        -parse_seqids \
        -out "$db"
fi

for gene in "${genes[@]}"; do
    query="${gene}.fa"

    if [ ! -f "$query" ]; then
        echo "WARNING: query FASTA not found for ${gene}: ${query}" >&2
        echo "${sp},${gene},NA,NA,NA,NA" >> "$OUT"
        continue
    fi

    blast_out="blast_results/${sp}.${gene}.blast.tsv"

    echo "BLAST: ${sp} ${gene}"

    blastn \
        -query "$query" \
        -db "$db" \
        -out "$blast_out" \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
        -max_target_seqs 5 \
        -evalue 1e-20 \
        -num_threads "${SLURM_CPUS_PER_TASK}"

done

echo "Done for ${sp}. Results written to ${OUT}"