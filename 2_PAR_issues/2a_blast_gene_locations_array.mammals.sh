#!/bin/bash
#SBATCH --job-name=gene_blast
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

module load blast

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals/

BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"
SPECIES_LIST="species_for_blast_array.txt"

mkdir -p slurm_output blast_dbs blast_results species_results

sp=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SPECIES_LIST")

if [ -z "$sp" ]; then
    echo "ERROR: No species found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

genes=(
    AKAP17A
    ANOS1
    APOA1A
    ASMT
    ASMTL
    ASMTL-AS1
    BOSD2_1
    BOSD2_2
    BOSD2_3
    BOSD2_4
    CRLF2
    CSF2RA
    DHRSX
    GPR143
    GTPBP6
    GYG2
    HRG
    IL3RA
    LINC00102
    LYL1
    MALRD1
    CD99
    MIR3690
    MXRA5
    NLGN4X
    OBP
    P2RY8
    PLCXD1
    PNPLA4
    PPP2R3B
    PRKX
    PUDP
    SHOX
    SHROOM2
    SLC25A6
    STS
    TBL1X
    XG
    ZBED1
    ZNF665
)

genome="${BASE}/${sp}/${sp}.fna"

if [ ! -f "$genome" ]; then
    echo "WARNING: genome not found for ${sp}: ${genome}" >&2
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

echo "Done for ${sp}."
