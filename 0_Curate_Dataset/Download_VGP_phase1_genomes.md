# Download all VGP Phase 1 genomes

## Using the list of accession numbers in the VGP_list_sex_chroms_curated file, dentify latest GCF and GCA versions of all genomes
```
ACCESSION_LIST="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD.csv"

mamba activate ncbi_datasets   # or wherever datasets is installed
chmod +x make_latest_gca_gcf_csv.sh
./make_latest_gca_gcf_csv.sh $ACCESSION_LIST > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv

head $ACCESSION_LIST
head /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv
```
### Check which species are in the VGP list but didn't make it to the list of latest accessions
```
export ACCESSION_LIST=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD.csv
python - "$ACCESSION_LIST" <<'PY'
import csv, sys

file1 = sys.argv[1]
file2 = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv"

def read_set(path, col, transform=lambda s: s):
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        if col not in r.fieldnames:
            raise SystemExit(f"{path}: missing column '{col}'. Columns: {r.fieldnames}")
        out = set()
        for row in r:
            v = (row.get(col) or "").strip()
            if v:
                out.add(transform(v))
        return out

# File1 has spaces in Scientific.Name (e.g., "Squalus sukleyi")
s1 = read_set(file1, "Scientific.Name")

# File2 has underscores in Species (e.g., "Squalus_sukleyi") -> convert to spaces
s2 = read_set(file2, "Species", transform=lambda s: s.replace("_", " "))

missing = sorted(s1 - s2)
print("\n".join(missing))
print(f"\n# Missing: {len(missing)} (Scientific.Name in file1 not in file2)")
PY

```
These did not make it, which is expected -- they lack an accession in the datasheet and are not part of this freeze.
```
Aechmophorus clarkii
Datnioides microlepis
Pelodiscus sinensis
Pelomedusa somalica
Rhinoderma darwinii
Squalus suckleyi
Tenrec ecaudatus
```
### Check which species made it to the latest accessions but do not have a "latest accession"
```
awk -F',' 'NR==1{print; next} $2=="" && $3=="" {print $1 ",,"}' \
  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv
```
These are all missing the latest accession:
```
Species,Accession_GCF,Accession_GCA
Saccopteryx_bilineata,,
Saccopteryx_leptura,,
Artibeus_lituratus,,
Artibeus_intermedius,,
Hoplias_malabaricus,,
Hypanus_sabinus,,
Lepisosteus_oculatus,,
Acridotheres_tristis,,
Zootoca_vivipara,,
```
Save the following as the file "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/did_not_run.csv" 
```
Species,Accession_GCF,Accession_GCA
Saccopteryx_bilineata,GCF_036850995.1,
Saccopteryx_leptura,GCF_036850765.1,
Artibeus_lituratus,,GCA_038363095.4
Artibeus_intermedius,,GCA_038363145.1
Hoplias_malabaricus,GCF_029633855.1,
Hypanus_sabinus,GCF_030144855.1,
Lepisosteus_oculatus,GCF_040954835.1,
Acridotheres_tristis,,GCA_027559615.1
Zootoca_vivipara,GCF_963506605.1,
```
Then run the following:
```
chmod +x download_vgp_array.troubleshoot.sh
./download_vgp_array.troubleshoot.sh
```
## Download all genomes and annotation files from the most recent GCF (preferred) or, if no GCF exists, most recent GCA
```
N=$(( $(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv) - 1 ))

sbatch --array=1-"$N"%25 download_vgp_array.sh
```
### Some inevitably fail, due to pinging NCBI too many times. Run the following to catch it up:
```
head -n 1 /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/did_not_run.csv
for sp in */; do
  sp="${sp%/}"
  if ! find "$sp" -type f -name '*.fna' \
        ! -name 'cds_from_genomic.fna' ! -name 'rna.fna' \
        ! -path '*/cds_from_genomic*' \
        -print -quit | grep -q .; then
    grep "$sp" /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv >> /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/did_not_run.csv
  fi
done

chmod +x download_vgp_array.troubleshoot.sh
./download_vgp_array.troubleshoot.sh
```
#### Check progress and confirm output
```
find . -type f -name '*.fna' | grep -v 'cds_from_genomic' | grep -v 'rna.fna' | wc -l

find . -type f -name '*.gff' | wc -l

find . -type f -name '*.zip' | wc -l

ls -d */ | wc -l

```
## Download the Western Chorus Frog
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes

echo 'Species,Accession_GCF,Accession_GCA' > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/did_not_run.csv
echo 'Pseudacris_triseriata,,GCA_053478255.1' >> /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/did_not_run.csv
./download_vgp_array.troubleshoot.sh
```
## Create reference files about downloads (e.g. paths to fna files)
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

mv genome_file_summary.csv ../reference_lists/

```
### Use this for some quick data confirmations and cleaning
#### Confirm that all genome directories contain an fna:
```
grep 'N,N,N' ../reference_lists/genome_file_summary.csv
```

Spot check these based on Slack from Simone:

Heads up these accession #s have been updated in the most recent data freeze version: "GCA_051911825.1" "GCA_964094495.3" "GCA_014108245.2" "GCA_014176215.2" "GCA_003287225.3" "GCA_014706295.2" "GCA_054100595.1" "GCA_964204865.2"

grep 'GCA_964204865' /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species_latest_accessions.csv

#### Remove zip files:
This line outputs what will be deleted, so I can be sure it's not going to delete more than I want:
```
find . -mindepth 2 -maxdepth 2 -type f -name "*.zip"
```
This line deletes it (ideally none should still exist):
```
find . -mindepth 2 -maxdepth 2 -type f -name "*.zip" -delete
```

### Create a list of GFF files
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes

mkdir -p ../reference_lists 

echo -e "Species\tAccession\tGFF" > /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/gff_file_list.tsv

find . -mindepth 5 -maxdepth 5 -type f -name "*.gff" | while read -r gff; do
  species=$(echo "$gff" | cut -d/ -f2)
  accession=$(echo "$gff" | cut -d/ -f5)
  fullpath=$(readlink -f "$gff")
  echo -e "${species}\t${accession}\t${fullpath}"
done >> /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/gff_file_list.tsv
```
### Create a list of species lacking GFF files
```
cd ../reference_lists
awk -F',' '$3=="N"{print $1}' genome_file_summary.csv > species_requiring_lifted_gff.txt
```
#### Compare and make sure that all species are represented in the "Has GFF" vs "Needs GFF" lists
```

LIST1="genome_file_summary.csv"
LIST2="/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/gff_file_list.tsv"
LIST3="species_requiring_lifted_gff.txt"

# Extract species names
cut -d',' -f1 "$LIST1" | tail -n +2 | sort > list1.species
awk '{print $1}' "$LIST2" | sort > list2.species
sort "$LIST3" > list3.species

echo "=== Species in LIST1 missing from gff_file_list.tsv ==="
comm -23 list1.species list2.species > needs_gff.txt

echo
echo "=== Species in LIST1 missing from species_requiring_lifted_gff.txt ==="
comm -23 needs_gff.txt list3.species
comm -23 list3.species needs_gff.txt

# Cleanup
rm list1.species list2.species list3.species needs_gff.txt
```
### Create a list of FNA files
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes

out=/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/fna_file_list.tsv
echo -e "Species\tAccession\tFNA" > "$out"

find "${PWD}" -mindepth 5 -maxdepth 5 -type f -name "GC*_genomic.fna" \
| awk -F/ '
{
  for (i=1; i<=NF; i++) {
    if ($i == "genomes") {
      species = $(i+1)
      acc     = $(i+4)
      break
    }
  }

  fna = $0

  # parse acc like GCA_009914755.4 or GCF_009914755.4
  if (match(acc, /^(GC[AF])_([0-9]+)\.([0-9]+)$/, m)) {
    prefix = m[1]
    id  = m[2] + 0
    ver = m[3] + 0

    # Prefer GCF repository: rank 0 for GCF, 1 for GCA
    # (smaller rank sorts first)
    rank = (prefix=="GCF" ? 0 : 1)
  } else {
    id = 0; ver = 0; rank = 9
  }

  # key fields:
  # species, id, ver (desc), rank (GCF first), then acc, fna
  printf("%s\t%010d\t%06d\t%d\t%s\t%s\n",
         species, id, ver, rank, acc, fna)
}
' \
| sort -t$'\t' -k1,1 -k4,4n -k2,2n -k3,3nr \
| awk -F'\t' '!seen[$1]++ { print $1 "\t" $5 "\t" $6 }' \
>> "$out"

```
### Create a list of protein files
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes

mkdir -p ../reference_lists 

echo -e "Species\tAccession\tFAA" > /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/faa_file_list.tsv

find . -mindepth 5 -maxdepth 5 -type f -name "protein.faa" | while read -r faa; do
  species=$(echo "$faa" | cut -d/ -f2)
  accession=$(echo "$faa" | cut -d/ -f5)
  fullpath=$(readlink -f "$faa")
  echo -e "${species}\t${accession}\t${fullpath}"
done >> /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/faa_file_list.tsv
```
### Create a list of cds files
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes

mkdir -p ../reference_lists 

echo -e "Species\tAccession\tCDS" > /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/cds_file_list.tsv

find . -mindepth 5 -maxdepth 5 -type f -name "cds_from_genomic.fna" | while read -r cds; do
  species=$(echo "$cds" | cut -d/ -f2)
  accession=$(echo "$cds" | cut -d/ -f5)
  fullpath=$(readlink -f "$cds")
  echo -e "${species}\t${accession}\t${fullpath}"
done >> /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/cds_file_list.tsv
```
## Recreate data structure for Genespace using symlinks
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/

tail -n +2 /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/fna_file_list.tsv \
| while read -r species accession fna; do
    mkdir -p "${species}"
    ln -sf "${fna}" "${species}/${species}.fna"
done

tail -n +2 /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/gff_file_list.tsv \
| while read -r species accession gff; do
    mkdir -p "${species}"
    ln -sf "${gff}" "${species}/${species}.gff"
done

tail -n +2 /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/faa_file_list.tsv \
| while read -r species accession faa; do
    mkdir -p "${species}"
    ln -sf "${faa}" "${species}/${species}.faa"
done

tail -n +2 /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/cds_file_list.tsv \
| while read -r species accession cds; do
    mkdir -p "${species}"
    ln -sf "${cds}" "${species}/${species}.cds"
done

```
## Rename X chromosome in Narcine bancroftii
```
# Check how the X chromosome is noted in the fna file
grep 'chromosome 12' GCF_036971445.1_sNarBan1.hap1_genomic.fna

# test a replacement of it to confirm that it is replacing the text we expect and none of that we don't expect
sed 's/chromosome 12/chromosome X/g' GCF_036971445.1_sNarBan1.hap1_genomic.fna | grep 'chromosome 12'
sed 's/chromosome 12/chromosome X/g' GCF_036971445.1_sNarBan1.hap1_genomic.fna | grep 'chromosome X' 

# replace the text
sed -i 's/chromosome 12/chromosome X/g' GCF_036971445.1_sNarBan1.hap1_genomic.fna

# Check which other files annotate this chromosome as 12 instead of X
grep 'chromosome 12' *

# Replace it in matching files using the same script used for the fna file
sed -i 's/chromosome 12/chromosome X/g' genomic.gff
sed -i 's/chromosome 12/chromosome X/g' rna.fna

# Confirm that no instances of chromosome 12 remain in any file
grep 'chromosome 12' *
```
# Lift over annotations for western chrous frog
**Same family annotations:**
- Dendropsophus ebraccatus GCF_027789765.1
- Hyla sarda GCF_029499605.1
- Use Hyla sarda (Hifiasm, hi-C phasing); Dendropsophus was CLR

```
minimap2 -d Pseudacris_triseriata.fna.mmi Pseudacris_triseriata.fna
```

```
#!/bin/bash

cd /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Pseudacris_triseriata

source myconda
mamba activate lifton

FNA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Hyla_sarda/ncbi_dataset/data/GCF_029499605.1/GCF_029499605.1_aHylSar1.hap1_genomic.fna

GFF_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Hyla_sarda/ncbi_dataset/data/GCF_029499605.1/genomic.gff

FAA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Hyla_sarda/ncbi_dataset/data/GCF_029499605.1/protein.faa

GFF_QRY=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Pseudacris_triseriata/Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.gff

FNA_QRY=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Pseudacris_triseriata/Pseudacris_triseriata.fna

SC=0.5

lifton \
    -g ${GFF_REF} \
    -o ${GFF_QRY} \
    -P ${FAA_REF} \
    -t 16 \
    -sc ${SC} \
    ${FNA_QRY} \
    ${FNA_REF}
```

```
sbatch \
  -c 16 \
  -t 5:00:00 \
  --mem-per-cpu=24G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --output=slurm_output/lifton_anuran.%j \
  submit_lifton.sh 
  ```

# NOT DOING THIS ANYMORE
### Download alternate haplotypes 
```
cd /data/Wilson_Lab/data/VGP_genomes_phase1/Alternate_Haplotypes

ACCESSION_LIST="/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/VGP_Phase1_Alternate_Haplotype_Accessions.csv"
N=$(( $(wc -l < "$ACCESSION_LIST") - 1 ))
sbatch --array=1-"$N"%25 download_vgp_array.althap.sh
```
