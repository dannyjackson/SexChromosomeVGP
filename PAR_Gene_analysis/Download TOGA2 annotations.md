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
  "https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Mustela_nivalis_vulgaris__Least_weasel__HLmusNivaVul3A__GCA_057128415.1/"
  
```
# Get Molossus nigricans
wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Molossus_nigricans__Northern_black_mastiff_bat__HLmolNig2A__GCA_039880945.1/"

# Get Siamang
wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Symphalangus_syndactylus__siamang__HLsymSyn4__GCA_028878055.3/"

# Get genomes from round2 of HalfDeep

cd /data/Wilson_Lab/data/TOGA2_Hiller/Homo_sapiens_hg38

# Bos_taurus -- updated TOGA2 version just finished running
wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Bos_indicus_x_Bos_taurus__Hybrid_cattle__HLbosIndiX1__GCF_003369695.1/"


# Lycaon_pictus
wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Lycaon_pictus__African_hunting_dog__HLlycPict4A__GCA_040955705.1/"

# Panthera_onca 
wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_human_hg38/Panthera_onca__Jaguar__HLpanOnca3A__GCA_046562875.1/"



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

# Retrieve genomes that were missed 

base="https://genome.senckenberg.de/download/TOGA2/TOGA2/reference_zebrafinch_HLtaeGut5-GCF_003957565.2/"
subdir="Pseudopipra_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1"
subdir="Aythya_marila__Greater_scaup__HLaytMari1A__GCA_965140915.1"
subdir="Struthio_camelus__African_ostrich__HLstrCame2A__GCF_040807025.1"
subdir="Anas_platyrhynchos__Mallard__HLanaPlat4A__GCA_964188345.1"

wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "${base}${subdir}/"

# Rename Dixiphia_pipra and Bos_taurus directories to match VGP nomenclature
mv Pseudopipra_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1 Dixiphia_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1

mv Bos_indicus_x_Bos_taurus__Hybrid_cattle__HLbosIndiX1__GCF_003369695.1 Bos_taurus__Hybrid_cattle__HLbosIndiX1__GCF_003369695.1/
```
## Integrated reference
```
mkdir -p /data/Wilson_Lab/data/TOGA2_Hiller/zebrafinch_HLtaeGut5-GCF_003957565.2/

cd /data/Wilson_Lab/data/TOGA2_Hiller/zebrafinch_HLtaeGut5-GCF_003957565.2
```
### Get list of all subdirs in TOGA database
```
cd /data/Wilson_Lab/data/TOGA2_Hiller/Integrated/birds
base="https://genome.senckenberg.de/download/TOGA2/TOGA2integration/v1/Aves"

wget -q -O - "$base" \
  | grep -o 'href="[^"]*/"' \
  | sed 's/^href="//; s/"$//' \
  | grep -v '^\.\./$' \
  | sed 's#/$##; s#.*/##' \
  > subdirs.birds.txt


cd /data/Wilson_Lab/data/TOGA2_Hiller/Integrated/mammals
base="https://genome.senckenberg.de/download/TOGA2/TOGA2integration/v1/Mammalia"

wget -q -O - "$base" \
  | grep -o 'href="[^"]*/"' \
  | sed 's/^href="//; s/"$//' \
  | grep -v '^\.\./$' \
  | sed 's#/$##; s#.*/##' \
  > subdirs.mammals.txt
```

### Submit slurm array for birds
```
cd /data/Wilson_Lab/data/TOGA2_Hiller/Integrated/birds

while read -r species; do
    grep "${species}" subdirs.birds.txt >> birds_PAR.txt
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/birds/species_for_blast_array.txt

n=$(wc -l < birds_PAR.txt)
sbatch --array=1-"$n"%20 download_query_annotation_array.integrated.birds.sh
```
### Slurm array script for birds:
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

base="https://genome.senckenberg.de/download/TOGA2/TOGA2integration/v1/Aves/"
subdirs_file="birds_PAR.txt"

cd /data/Wilson_Lab/data/TOGA2_Hiller/Integrated/birds
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

### Submit slurm array for mammals
```
cd /data/Wilson_Lab/data/TOGA2_Hiller/Integrated/mammals

while read -r species; do
    grep "${species}" subdirs.mammals.txt >> mammals_PAR.txt
done < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals/species_for_blast_array.txt

n=$(wc -l < mammals_PAR.txt)
sbatch --array=1-"$n"%20 download_query_annotation_array.integrated.mammals.sh
```
### Slurm array script for mammals:
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

base="https://genome.senckenberg.de/download/TOGA2/TOGA2integration/v1/Mammalia/"
subdirs_file="mammals_PAR.txt"

cd /data/Wilson_Lab/data/TOGA2_Hiller/Integrated/mammals
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

# remove directories of non-used accessions

### Remove directories of unused mammalian accessions 
```
cd /data/Wilson_Lab/data/TOGA2_Hiller/Integrated/mammals/Mammalia

rm -r Balaenoptera_physalus__Fin_whale__HLbalPhy4__GCA_965194765.1
rm -r Balaenoptera_physalus__Fin_whale__HLbalPhy3__GCA_023338255.1
rm -r Balaenoptera_physalus__Fin_whale__HLbalPhy2
rm -r Callithrix_jacchus__white-tufted-ear_marmoset__HLcalJac5__GCF_011100555.1
rm -r Callithrix_jacchus__white-tufted-ear_marmoset__calJac4__GCF_009663435.1
rm -r Camelus_dromedarius__Arabian_camel__HLcamDro2__GCA_000803125.3
rm -r Camelus_dromedarius__Arabian_camel__HLcamDro3B__GCA_036321565.1
rm -r Capra_hircus__goat__HLcapHir2__GCA_001704415.1
rm -r Gorilla_gorilla_gorilla__western_lowland_gorilla__gorGor6__GCF_008122165.1
rm -r Grampus_griseus__Rissos_dolphin__HLgraGri1
rm -r Inia_geoffrensis__Boutu__HLiniGeo2B__GCA_036417475.1
rm -r Loxodonta_africana__African_savanna_elephant__HLloxAfr5B__GCA_030020305.1
rm -r Macaca_nemestrina__pig-tailed_macaque__macNem1__GCF_000956065.1
rm -r Marmota_flaviventris__Yellow-bellied_marmot__HLmarFlav2B__GCA_047512175.1
rm -r Marmota_flaviventris__yellow-bellied_marmot__HLmarFla1__GCA_003676075.2
rm -r Mesoplodon_bidens__Sowerbys_beaked_whale__HLmesBid1__GCA_004027085.1
rm -r Mesoplodon_bidens__Sowerbys_beaked_whale__HLmesBid2B__GCA_964148845.1
rm -r Myotis_nattereri__Natterers_bat__HLmyoNatt2B__GCA_964212025.2
rm -r Myotis_nattereri__Natterers_bat__HLmyoNat1B__GCA_964212025.1
rm -r Myotis_nattereri__Natterers_bat__HLmyoNat1A__GCA_964212035.1
rm -r Ovis_aries__sheep__HLoviAri6__GCF_016772045.2
rm -r Ovis_canadensis__bighorn_sheep__HLoviCan2__GCA_004026945.1
rm -r Pan_troglodytes__chimpanzee__panTro6__GCF_002880755.1
rm -r Pongo_abelii__Sumatran_orangutan__ponAbe3__GCF_002880775.1
rm -r Pongo_pygmaeus__Bornean_orangutan__HLponPyg1__GCA_023767775.1
rm -r Pseudorca_crassidens__false_killer_whale__HLpseCra1
rm -r Rhynchocyon_petersi__Black_and_rufous_elephant_shrew__HLrhyPet1B__GCA_043290065.1
rm -r Rhynchonycteris_naso__Proboscis_bat__HLrhyNas2B__GCA_037038555.1
rm -r Symphalangus_syndactylus__siamang__HLsymSyn3__GCF_028878055.1
rm -r Trichechus_inunguis__Amazon_manatee__HLtriInu1B__GCA_046562985.1
rm -r Urocitellus_parryii__Arctic_ground_squirrel__HLuroPar1__GCA_003426925.1
rm -r Urocitellus_parryii__Arctic_ground_squirrel__HLuroPar2B__GCA_045843765.1


base="https://genome.senckenberg.de/download/TOGA2/TOGA2integration/v1/Mammalia/"
subdir="Mustela_nivalis__Least_weasel__HLmusNiva2A__GCA_964662115.1"

wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "${base}${subdir}/"
```
# Download birds that didn't download in array

```
base="https://genome.senckenberg.de/download/TOGA2/TOGA2integration/v1/Aves/"
subdir="Pseudopipra_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1"
# subdir="Aythya_marila__Greater_scaup__HLaytMari1A__GCA_965140915.1"
# subdir="Struthio_camelus__African_ostrich__HLstrCame2A__GCF_040807025.1"
wget -c -r -l 1 -np -nH --cut-dirs=4 \
  --reject "index.html*" \
  --accept "query_annotation*" \
  "${base}${subdir}/"

# Rename Dixiphia_pipra directory to match VGP nomenclature
mv Pseudopipra_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1 Dixiphia_pipra__white-crowned_manakin__HLpsePip1A__GCF_036250125.1