# Use SCINKD to detect instances of PAR dropout in Cetaceans

# Load SCINKD
source myconda
mamba activate scinkd 

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/

mkdir cetaceans

cd cetaceans

# 0. Collate all the data needed to evaluate this 
## Species of interest
Eubalaena_glacialis
Eschrichtius_robustus
Balaenoptera_physalus
Megaptera_novaeangliae
Balaenoptera_musculus
Mesoplodon_densirostris
Mesoplodon_bidens
Inia_geoffrensis
Delphinus_delphis
Stenella_coeruleoalba
Grampus_griseus
Pseudorca_crassidens
Globicephala_melas

## Predicted PAR status
Eubalaena glacialis,PRESENT
Eschrichtius robustus,DROPOUT
Balaenoptera physalus,PRESENT
Megaptera novaeangliae,DROPOUT
Balaenoptera musculus,DROPOUT
Mesoplodon densirostris,DROPOUT
Mesoplodon bidens,PRESENT
Inia geoffrensis,PRESENT
Delphinus delphis,DROPOUT
Stenella coeruleoalba,DROPOUT
Grampus griseus,PRESENT
Pseudorca crassidens,PRESENT
Globicephala melas,DROPOUT

# Prep files for SCINKD
```

SPECIES_LIST="cetacean.species.txt"

BASE_ANALYSIS_DIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/cetaceans"
SEXCHR_DIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs"
GENOME_ROOT="/data/Wilson_Lab/data/VGP_genomes"

SEXCHROM_CSV="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/sexchrom_accessions.csv"

if [[ ! -f "$SPECIES_LIST" ]]; then
  echo "ERROR: species list not found: $SPECIES_LIST" >&2
  exit 1
fi
if [[ ! -f "$SEXCHROM_CSV" ]]; then
  echo "ERROR: mapping CSV not found: $SEXCHROM_CSV" >&2
  exit 1
fi

mkdir -p "$SEXCHR_DIR"

get_acc() {
  # usage: get_acc <species> <X|Y>
  local sp="$1"
  local chrom="$2"
  awk -F',' -v sp="$sp" -v c="$chrom" '
    $1==sp && $2==c { gsub(/\r/,"",$3); print $3; exit }
  ' "$SEXCHROM_CSV"
}

while IFS= read -r SPECIES; do
  [[ -z "${SPECIES:-}" ]] && continue
  [[ "${SPECIES:0:1}" == "#" ]] && continue

  GENOME_FA="${GENOME_ROOT}/${SPECIES}/${SPECIES}.fna"
  if [[ ! -f "$GENOME_FA" ]]; then
    echo "WARNING: genome not found for $SPECIES at $GENOME_FA; skipping." >&2
    continue
  fi

  X_ACC="$(get_acc "$SPECIES" "X")"
  Y_ACC="$(get_acc "$SPECIES" "Y")"

  if [[ -z "${X_ACC:-}" || -z "${Y_ACC:-}" ]]; then
    echo "WARNING: missing X/Y accession in $SEXCHROM_CSV for $SPECIES (X='${X_ACC:-}', Y='${Y_ACC:-}'); skipping." >&2
    continue
  fi

  echo "==> Processing: $SPECIES  (X=$X_ACC, Y=$Y_ACC)"

  SPECIES_DIR="${BASE_ANALYSIS_DIR}/${SPECIES}"
  mkdir -p "$SPECIES_DIR"
  cd "$SPECIES_DIR"

  # Ensure genome is indexed
  samtools faidx "$GENOME_FA"

  # Extract X and Y by accession/contig name from CSV
  samtools faidx "$GENOME_FA" "$X_ACC" > "${SEXCHR_DIR}/${SPECIES}.X.fa"
  bgzip -f "${SEXCHR_DIR}/${SPECIES}.X.fa"

  samtools faidx "$GENOME_FA" "$Y_ACC" > "${SEXCHR_DIR}/${SPECIES}.Y.fa"
  bgzip -f "${SEXCHR_DIR}/${SPECIES}.Y.fa"

  # Symlinks to hap FASTAs
  ln -sf "${SEXCHR_DIR}/${SPECIES}.X.fa.gz" "${SPECIES}.hap1.fasta.gz"
  ln -sf "${SEXCHR_DIR}/${SPECIES}.Y.fa.gz" "${SPECIES}.hap2.fasta.gz"

  # Index gz FASTAs
  samtools faidx "${SPECIES}.hap1.fasta.gz"
  samtools faidx "${SPECIES}.hap2.fasta.gz"

  # Write config
  cat > "${SPECIES_DIR}/config.${SPECIES}.json" <<EOF
{
  "prefix": "${SPECIES}"
}
EOF

done < "$SPECIES_LIST"

echo "Done."

```

# 1. Run SCINKD
Save the following as scinkd_submit.cetaceans.sh:
```
#!/bin/bash
#SBATCH -J runscinkd_cetaceans
#SBATCH -c 8
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -o runscinkd_cetaceans.%j.out
#SBATCH -e runscinkd_cetaceans.%j.err

# Parse command-line arguments
while getopts s: option; do
    case "${option}" in
        s) SPECIES=${OPTARG};;
        *) echo "Invalid option: -${OPTARG}" >&2; exit 1;;
    esac
done

if [ -z "${SPECIES}" ]; then
    echo "Error: No species name provided." >&2
    exit 1
fi

source myconda

mamba activate scinkd

#run SCINKD in greedy mode for quick testing

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/cetaceans/${SPECIES}

snakemake --use-conda -c "${SLURM_CPUS_PER_TASK:-1}" \
    -s /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/SCINKD.v2.1.0.GREEDY.snakefile \
    --configfile /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/cetaceans/${SPECIES}/config.${SPECIES}.json --printshellcmds --show-failed-logs
```
Run it as a slurm array:
```
SPECIES_LIST="cetacean.species.txt"

while IFS= read -r SPECIES; do
    [[ -z "$SPECIES" ]] && continue
    [[ "${SPECIES:0:1}" == "#" ]] && continue

    sbatch scinkd_submit.cetaceans.sh -s "$SPECIES"
done < "$SPECIES_LIST"
```
Plot it:
```

while IFS= read -r SPECIES; do
    [[ -z "$SPECIES" ]] && continue
    [[ "${SPECIES:0:1}" == "#" ]] && continue

OUTDIR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/cetaceans/$SPECIES

Rscript plot_kmer_density.r -s "$SPECIES" -o "$OUTDIR"
done < "$SPECIES_LIST"

```