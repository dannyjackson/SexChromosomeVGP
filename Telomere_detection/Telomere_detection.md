# VGP Genome Telomere Detection
This script infers if the X and Y chromosomes are bounded by telomeres in VGP Data Freeze 1 genomes.

## Outline of steps
1. Make a reference list of all VGP genome fasta files
2. Make a reference list of all sex chromosome accession numbers in these files
3. Subset out each sex chromosome from each fasta
4. Use the Telomere Identification toolKit (tidk) to annotate telomeres on each sex chromosome (https://github.com/tolkit/telomeric-identifier)
5. Visualize the data for each sex chromosome


## 1. Make a reference list of all VGP genome fasta files
```
cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/
wget https://raw.githubusercontent.com/VGP/vgp-phase1/main/VGPPhase1-freeze-1.0.tsv

ls /data/Wilson_Lab/data/VGP_genomes/*/*.fna | grep -Ev '\.(cds|rna)\.fna$' | wc -l

base=/data/Wilson_Lab/data/VGP_genomes

for d in "$base"/*/; do
   name=$(basename "$d")
   [[ -f "$d/$name.fna" ]] || echo "$name"
 done
ls $base/Myotis_mystacinus

# Missing Myotis_mystacinus
https://github.com/VGP/vgp-phase1/blob/main/VGPPhase1-freeze-1.0.tsv

ls /data/Wilson_Lab/data/VGP_genomes/*/*.fna | grep -Ev '\.(cds|rna)\.fna$' > /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/VGP_phase1_fastas.txt

sed 's#.*/##; s/\.fna$//' VGP_phase1_fastas.txt | sort -u > fasta_species.list
awk -F'\t' 'NR>1 {
  split($10,a," ");
  print a[1] "_" a[2]
}' VGPPhase1-freeze-1.0.tsv | sort -u > vgp_species.list

comm -23 vgp_species.list fasta_species.list | wc -l

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/
comm -23 vgp_species.2.list fasta_species.list | wc -l

# 8 missing species in fasta list
Aechmophorus_clarkii
Datnioides_microlepis
Myotis_mystacinus
Pelodiscus_sinensis
Pelomedusa_somalica
Rhinoderma_darwinii
Squalus_suckleyi
Tenrec_ecaudatus


# these are missing accession numbers:
Aechmophorus clarkii
Datnioides microlepis
Pelodiscus sinensis
Pelomedusa somalica
Rhinoderma darwinii
Squalus suckleyi
Tenrec ecaudatus

# so it is just the Myotis_mystacinus genome  
module load ncbi-datasets/18.13.0

cd /data/Wilson_Lab/data/VGP_genomes/Myotis_mystacinus

datasets download genome accession GCA_964094495.3

unzip ncbi_dataset.zip

mv ncbi_dataset/data/* .

mv GCA_964094495.3/GCA_964094495.3_mMyoMys1.hap1.2_genomic.fna Myotis_mystacinus.fna
```

## 2. Make a reference list of all sex chromosome accession numbers in these files

Index genomes:
```
BASE="/data/Wilson_Lab/data/VGP_genomes"
LIST="fna_missing_fai.txt"

# Create a list of all fna files that are missing an fai file

find "$BASE" -type f -path "$BASE/*/*.fna" ! -path "$BASE/*/*.*/" -print0 \
  | while IFS= read -r -d '' fna; do
      dir="$(dirname "$fna")"
      species="$(basename "$dir")"
      file="$(basename "$fna")"

      # Keep only SPECIES/SPECIES.fna exactly
      if [[ "$file" == "${species}.fna" ]]; then
        [[ -f "${fna}.fai" ]] || printf '%s\n' "$fna"
      fi
    done > "$LIST"

wc -l "$LIST"

# Use that list to submit a slurm array to index all genomes

mkdir -p logs
LIST="fna_missing_fai.txt"
N="$(wc -l < "$LIST")"

# Run up to 50 tasks at once (adjust %50 for your cluster)
sbatch --array=1-"$N"%50 faidx_array.sbatch "$LIST"
```

Use this code to check the status of the array -- it counts how many genome directories lack an fai file.

```
BASE="/data/Wilson_Lab/data/VGP_genomes"
LIST="fna_missing_fai.progress.txt"

find "$BASE" -type f -path "$BASE/*/*.fna" ! -path "$BASE/*/*.*/" -print0 \
  | while IFS= read -r -d '' fna; do
      dir="$(dirname "$fna")"
      species="$(basename "$dir")"
      file="$(basename "$fna")"

      # Keep only SPECIES/SPECIES.fna exactly
      if [[ "$file" == "${species}.fna" ]]; then
        [[ -f "${fna}.fai" ]] || printf '%s\n' "$fna"
      fi
    done > "$LIST"

wc -l "$LIST"

# after the script is finished, confirm that no errors were noted and then remove the log directory:

rm -r logs
```

Run this to create a TSV with columns "Species,Chromosome,Accession"
```
export TSV=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/VGP_freeze_hap1_combined_sexchroms_seq_reports.tsv
export GENOMEDIR=/data/Wilson_Lab/data/VGP_genomes

GITREP="~/dannys_githubs/SexChromosomeVGP/"

bash $GITREPO/make_sexchrom_accessions.sh $OUTDIR/sexchrom_accessions.csv
```
Some values were missing from this, due to mismatch between the VGP seq reports accession numbers and the version of fasta in our genome directory:
```
Meles_meles,X,NA
Meles_meles,Y,NA
Patagioenas_fasciata,W,NA
Patagioenas_fasciata,Z,NA
Tursiops_truncatus,Y,NA
Pongo_pygmaeus,X,NA
Pongo_pygmaeus,Y,NA
Aquila_chrysaetos,W,NA
Aquila_chrysaetos,Z,NA
Scyliorhinus_canicula,X,NA
Lemur_catta,X,NA
Lemur_catta,Y,NA
Pan_troglodytes,X,NA
Pan_troglodytes,Y,NA
Canis_lupus_familiaris,X,NA
Canis_lupus_baileyi,X,NA
Canis_lupus_baileyi,Y,NA
Myotis_mystacinus,X,NA
Myotis_mystacinus,Y,NA
Canis_lupus_orion,X,NA
Canis_lupus_orion,Y,NA
Pongo_abelii,X,NA
Pongo_abelii,Y,NA
```
Use these two lines of code to identify sex chromosome accession numbers for the missing species:
```
grep 'chromosome Z' $GENOMEDIR/Aquila_chrysaetos/Aquila_chrysaetos.fna
grep 'chromosome W' $GENOMEDIR/Aquila_chrysaetos/Aquila_chrysaetos.fna

grep 'chromosome X' $GENOMEDIR/Canis_lupus/Canis_lupus.fna
grep 'chromosome Y' $GENOMEDIR/Canis_lupus/Canis_lupus.fna
```
Run the following code to remove rows with NA values and insert updated ones. All *Canis lupus* genomes except for *baileyi* were dropped from this version of the freeze.
```
# Remove the placeholder rows
grep -v ',NA$' sexchrom_accessions.csv > tmp.csv

# Append the corrected entries
cat >> tmp.csv << 'EOF'
Meles_meles,X,NC_060087.1
Meles_meles,Y,NC_060088.1
Patagioenas_fasciata,W,NC_092559.1
Patagioenas_fasciata,Z,NC_092560.1
Tursiops_truncatus,Y,NW_022983135.1
Pongo_pygmaeus,X,NC_072396.2
Pongo_pygmaeus,Y,NC_072397.2
Aquila_chrysaetos,W,NC_054457.1
Aquila_chrysaetos,Z,NC_044030.1
Lemur_catta,X,NC_059155.1
Lemur_catta,Y,NC_059156.1
Pan_troglodytes,X,NC_072421.2
Pan_troglodytes,Y,NC_072422.2
Myotis_mystacinus,X,OZ075425.2
Myotis_mystacinus,Y,OZ075426.2
Pongo_abelii,X,NC_072008.2
Pongo_abelii,Y,NC_072009.2
Scyliorhinus_canicula,X,NC_052173.1
EOF

mv tmp.csv sexchrom_accessions.csv

echo 'Canis_lupus,X,NC_132876.1' >> sexchrom_accessions.csv
echo 'Canis_lupus,Y,NC_132877.1' >> sexchrom_accessions.csv
```
## 3. Subset out each sex chromosome from each fasta
```
SEXCHRFILE="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/sexchrom_accessions.csv"
OUTDIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs"

mkdir -p "${OUTDIR}"

# Skip header, read CSV
tail -n +2 "${SEXCHRFILE}" | while IFS=, read -r SPECIES SEXCHR ACCESSION; do
    FASTA="${GENOMEDIR}/${SPECIES}/${SPECIES}.fna"
    CHROM="${ACCESSION}"
    OUT="${OUTDIR}/${SPECIES}.${SEXCHR}.fa"

    # index if needed
    [[ -f "${FASTA}.fai" ]] || samtools faidx "${FASTA}"

    echo "Extracting ${SPECIES} ${SEXCHR}"
    samtools faidx "${FASTA}" "${CHROM}" > "${OUT}"
done
```

## 4. Use the Telomere Identification toolKit (tidk) to annotate telomeres on each sex chromosome (https://github.com/tolkit/telomeric-identifier)
```
FASTADIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs"
OUTDIR="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/Telomere_detection/sexchrs/tidk_results"

mkdir -p "${OUTDIR}"
    
for FA in "${FASTADIR}"/*.fa; do
    BASENAME=$(basename "${FA}" .fa)   # SPECIES.SEXCHR
    SPECIES="${BASENAME%%.*}"          # text before first dot

    echo "Running TIDK on ${FA}"
    tidk search "${FA}" \
        --string TTAGGG \
        --output "${BASENAME}" \
        --dir "${OUTDIR}"/"${SPECIES}"
done
```


# 

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/Telomere_detection

module load mamba_install
source myconda

conda install -c bioconda tidk

FASTA=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13/Xchr.fa

tidk search $FASTA --string TTAGGG --output CHM13_Xchr --dir Homo_sapiens

FASTA=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13/Xchr.fa

tidk search $FASTA --string TTAGGG --output CHM13_Xchr --dir Homo_sapiens


## 5. Visualize the data for each sex chromosome
```
#!/bin/bash

module load R

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/Telomere_detection/sexchrs/tidk_results

Rscript plot_telomere.faceted.r

```

sbatch plot_telomere.faceted.sh

9583517