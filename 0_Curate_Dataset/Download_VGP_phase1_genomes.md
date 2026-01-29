# Download all VGP Phase 1 genomes

### Using the list of accession numbers in the VGP_list_sex_chroms_curated file, download all genomes and annotation files from the VGP phase 1 datafreeze
This script uses the "Accession_num_main_haplotype" column in the reference datasheet, which is in format GCA_*. It preferrentially pulls RefSeq files by dropping GCA and using GCF first, but will fall back on GCA files if no GCF match is found.
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes

ACCESSION_LIST="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD.csv"
N=$(( $(wc -l < "$ACCESSION_LIST") - 1 ))
sbatch --array=1-"$N"%25 download_vgp_array.sh
```
### Merge log files
```
jobid=10389948  # your array job id
GENOME_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1"

# OK
{
  printf "Scientific.Name\tAccessionUsed\tZip\n"
  awk 'FNR>1 { print }' "${GENOME_DIR}/download_ok.${jobid}."*.tsv
} > "${GENOME_DIR}/download_ok.${jobid}.tsv"

# FAIL

{
  printf "Scientific.Name\tAccessionUsed\tZip\n"
  awk 'FNR>1 { print }' "${GENOME_DIR}/download_fail.${jobid}."*.tsv
} > "${GENOME_DIR}/download_fail.${jobid}.tsv"

wc -l download_ok.${jobid}.tsv
wc -l download_fail.${jobid}.tsv

# clear out individual jobid output files 

rm download_ok.${jobid}.*.tsv
rm download_fail.${jobid}.*.tsv
```
### Create summary file of genomes with sequence and annotation files
```
set -u

OUT="genome_file_summary.csv"

echo "Species,fna,gff,faa" > "$OUT"

for d in */; do
  # remove trailing slash
  species="${d%/}"

  # skip non-directories just in case
  [[ -d "$species" ]] || continue

  fna="N"
  gff="N"
  faa="N"

  # check for files recursively
  if find "$species" -type f \( -iname "*.fna" -o -iname "*.fa" -o -iname "*.fasta" \) -print -quit | grep -q .; then
    fna="Y"
  fi

  if find "$species" -type f \( -iname "*.gff" -o -iname "*.gff3" \) -print -quit | grep -q .; then
    gff="Y"
  fi

  if find "$species" -type f -iname "*.faa" -print -quit | grep -q .; then
    faa="Y"
  fi

  echo "${species},${fna},${gff},${faa}" >> "$OUT"
done

set +u

echo "Wrote $OUT"
```
### Check if all genomes downloaded an fna:
```
grep 'N,N,N' genome_file_summary.csv
```
All genomes that should have had a genome did download an fna.

These did not download an fna but they are not part of this datafreeze; they do not have published genomes.
```
Scientific.Name AccessionUsed   Zip
Aechmophorus_clarkii,N,N,N
Datnioides_microlepis,N,N,N
Pelomedusa_somalica,N,N,N
Rhinoderma_darwinii,N,N,N
Squalus_suckleyi,N,N,N
Tenrec_ecaudatus,N,N,N
```

Spot check these based on Slack from Simone:

Heads up these accession #s have been updated in the most recent data freeze version: "GCA_051911825.1" "GCA_964094495.3" "GCA_014108245.2" "GCA_014176215.2" "GCA_003287225.3" "GCA_014706295.2" "GCA_054100595.1" "GCA_964204865.2"



### Remove zip files:
This line outputs what will be deleted, so I can be sure it's not going to delete more than I want:
```
find . -mindepth 2 -maxdepth 2 -type f -name "*.zip"
```
This line deletes it:
```
find . -mindepth 2 -maxdepth 2 -type f -name "*.zip" -delete
```

# Create a list of GFF files
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/

mkdir -p reference_lists 

echo -e "Species\tAccession\tGFF" > /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/gff_file_list.tsv

find . -mindepth 5 -maxdepth 5 -type f -name "*.gff" | while read -r gff; do
  species=$(echo "$gff" | cut -d/ -f2)
  accession=$(echo "$gff" | cut -d/ -f5)
  fullpath=$(readlink -f "$gff")
  echo -e "${species}\t${accession}\t${fullpath}"
done >> /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/gff_file_list.tsv

mv ../genomes/*tsv .
mv ../genomes/*csv .
```
# Create a list of species lacking GFF files
```
awk -F',' '$3=="N"{print $1}' genome_file_summary.csv > species_requiring_lifted_gff.txt
```