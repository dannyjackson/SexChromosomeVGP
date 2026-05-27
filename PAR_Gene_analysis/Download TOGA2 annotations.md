# Download TOGA2 annotations
# Reference: Yury Malovichko et al. "Accurate, comprehensive gene annotation and ortholog identification across thousands of vertebrate genomes with TOGA2", in preparation

## Human reference
```
mkdir -p /data/Wilson_Lab/data/TOGA2_Hiller/Homo_sapiens_hg38

cd /data/Wilson_Lab/data/TOGA2_Hiller/Homo_sapiens_hg38
```
### Get list of all subdirs in TOGA database
```
base="https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/"

wget -q -O - "$base" \
  | grep -o 'href="[^"]*/"' \
  | sed 's/^href="//; s/"$//' \
  | grep -v '^\.\./$' \
  | sed 's#/$##; s#.*/##' \
  > subdirs.txt
```

### Submit slurm array
```
echo 'Pongo_pygmaeus' >> /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals/species_for_blast_array.txt
echo 'Symphalangus_syndactylus' >> /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals/species_for_blast_array.txt
echo 'Gorilla_gorilla' >> /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals/species_for_blast_array.txt
echo 'Pongo_abelii' >> /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals/species_for_blast_array.txt


while read -r species; do
    grep "${species}" subdirs.txt >> mammals_PAR.txt
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals/species_for_blast_array.txt

n=$(wc -l < mammals_PAR.txt)
sbatch --array=1-"$n"%20 download_query_annotation_array.sh
```
### Slurm array script:
```
#!/bin/bash
#SBATCH --job-name=toga_query_ann
#SBATCH --output=logs/toga_query_ann_%A_%a.out
#SBATCH --error=logs/toga_query_ann_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-51%20

set -euo pipefail

base="https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/"
subdirs_file="mammals_PAR.txt"

mkdir -p logs

subdir="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$subdirs_file")"

if [[ -z "${subdir}" ]]; then
  echo "No subdir found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: downloading query_annotation files from ${subdir}"

wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "${base}${subdir}/"
```

# Rename unannotated genes in gff, or use
```
wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Mustela_nivalis__Least_weasel__HLmusNiva2A__GCA_964662115.1/"
```


## Zebra finch reference
```
mkdir -p /data/Wilson_Lab/data/TOGA2_Hiller/zebrafinch_HLtaeGut5-GCF_003957565.2/

cd /data/Wilson_Lab/data/TOGA2_Hiller/zebrafinch_HLtaeGut5-GCF_003957565.2
```
### Get list of all subdirs in TOGA database
```
base="https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_zebrafinch_HLtaeGut5-GCF_003957565.2/"

wget -q -O - "$base" \
  | grep -o 'href="[^"]*/"' \
  | sed 's/^href="//; s/"$//' \
  | grep -v '^\.\./$' \
  | sed 's#/$##; s#.*/##' \
  > subdirs.txt
```

### Submit slurm array
```
cat /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/birds/species_for_blast_array.txt

while read -r species; do
    grep "${species}" subdirs.txt >> birds_PAR.txt
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/birds/species_for_blast_array.txt

n=$(wc -l < birds_PAR.txt)
sbatch --array=1-"$n"%20 download_query_annotation_array.sh
```
### Slurm array script:
```
#!/bin/bash
#SBATCH --job-name=toga_query_ann
#SBATCH --output=logs/toga_query_ann_%A_%a.out
#SBATCH --error=logs/toga_query_ann_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=1-51%20

set -euo pipefail

base="https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_zebrafinch_HLtaeGut5-GCF_003957565.2/"
subdirs_file="birds_PAR.txt"

cd /data/Wilson_Lab/data/TOGA2_Hiller/zebrafinch_HLtaeGut5-GCF_003957565.2

mkdir -p logs

subdir="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$subdirs_file")"

if [[ -z "${subdir}" ]]; then
  echo "No subdir found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: downloading query_annotation files from ${subdir}"

wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "${base}${subdir}/"
```
### Remove directories of unused accessions
```
rm -r Anas_platyrhynchos_platyrhynchos__common_mallard__HLanaPla2__GCA_002743455.1
rm -r Aythya_marila__Greater_scaup__HLaytMari1A__GCA_965140915.1
rm -r Calonectris_borealis__Corys_shearwater__HLcalBor1__GCA_013401115.1
rm -r Calonectris_borealis__Corys_shearwater__HLcalBor2
rm -r Colius_striatus__speckled_mousebird__colStr1__GCF_000690715.1
rm -r Coloeus_monedula__jackdaw__HLcolMon1__GCA_013407035.1
rm -r Columba_livia__rock_pigeon__HLcolLiv2__GCA_000337935.2
rm -r Lathamus_discolor__Swift_parrot__HLlatDis1B__GCA_037157625.1
rm -r Opisthocomus_hoazin__hoatzin__opiHoa1__GCF_000692075.1
rm -r Passer_domesticus__House_sparrow__HLpasDom1__GCA_001700915.1
rm -r Patagioenas_fasciata__Band-tailed_pigeon__HLpatFas2__GCA_036971685.1
rm -r Patagioenas_fasciata_monilis__band-tailed_pigeon__HLpatFas1__GCA_002029285.1
rm -r Poecile_atricapillus__Black-capped_chickadee__HLpoeAtr1__GCA_011421415.1
rm -r Rhynochetos_jubatus__kagu__HLrhyJub1__GCA_013398095.1
rm -r Zosterops_lateralis_melanops__silver-eye__HLzosLat1__GCA_001281735.1

wget Pseudopipra_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1/
wget Aythya_marila__Greater_scaup__HLaytMari1A__GCA_965140915.1/

base="https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_zebrafinch_HLtaeGut5-GCF_003957565.2/"
subdir="Pseudopipra_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1"
subdir="Aythya_marila__Greater_scaup__HLaytMari1A__GCA_965140915.1"

wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "${base}${subdir}/"

# Rename Dixiphia_pipra directory to match VGP nomenclature
mv Pseudopipra_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1 Dixiphia_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1
