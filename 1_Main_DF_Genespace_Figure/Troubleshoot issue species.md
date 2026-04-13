# Troubleshoot issue species

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

sbatch 1b_genespace.singleinstance.sh Gallus_gallus ${SPECIES} ${OUTDIR}

# Anurans
## Hyla sarda
Whichever chr is syntenic to Xenopus chr 4 is the SC
https://academic.oup.com/mbe/article/38/1/192/5882026?guestAccessKey=
When Sex Chromosomes Recombine Only in the Heterogametic Sex: Heterochiasmy and Heterogamety in Hyla Tree Frogs 
```
####################################################################################
## Synteny to Pseudacristis
####################################################################################
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Hyla_sarda
mkdir -p ${WORKDIR}
cd ${WORKDIR}

mkdir -p bed peptide

while read -r SPECIES; do
  cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed" bed/
  cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa" peptide/
done <<'EOF'
Hyla_sarda
EOF

cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/bed/Pseudacris_triseriata.bed bed/
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/peptide/Pseudacris_triseriata.fa peptide/


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}

####################################################################################
# Synteny to chicken
####################################################################################
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Hyla_sarda
mkdir -p ${WORKDIR}
cd ${WORKDIR}

mkdir -p bed peptide

while read -r SPECIES; do
  cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed" bed/
  cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa" peptide/
done <<'EOF'
Hyla_sarda
Gallus_gallus
EOF

# Replace chr 7 with X

awk 'BEGIN{OFS="\t"} {$1=($1=="7"?"X":$1); print}' \
bed/Hyla_sarda.bed \
> tmp && mv tmp bed/Hyla_sarda.bed


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}
```
# Sharks/Rays
## Scyliorhinus canicula
No ID'd sex chr -- ID with synteny in shark analysis; rename to X and rerun
```
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Scyliorhinus_canicula/
cd ${WORKDIR}

# Replace chr 28 with X

awk 'BEGIN{OFS="\t"} {$1=($1=="28"?"X":$1); print}' \
bed/Scyliorhinus_canicula.bed \
> tmp && mv tmp bed/Scyliorhinus_canicula.bed

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}
```
## Narcine bancroftii
```
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Narcine_bancroftii/
cd ${WORKDIR}

# Replace chr 28 with X

awk 'BEGIN{OFS="\t"} {$1=($1=="12"?"X":$1); print}' \
bed/Narcine_bancroftii.bed \
> tmp.bed && mv tmp.bed bed/Narcine_bancroftii.bed

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}
```
## Amblyraja radiata
No ID'd sex chr -- ID with synteny in shark analysis; rename to X and rerun
```
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Amblyraja_radiata/
cd ${WORKDIR}

mkdir OG_Naming
mv * OG_Naming
mv OG_Naming/bed .
mv OG_Naming/peptide .

# Replace chr X with Unnamed

awk 'BEGIN{OFS="\t"} {$1=($1=="X"?"Unnamed":$1); print}' \
bed/Amblyraja_radiata.bed \
> tmp.bed && mv tmp.bed bed/Amblyraja_radiata.bed

# Replace chr 46 with X
awk 'BEGIN{OFS="\t"} {$1=($1=="46"?"X":$1); print}' \
bed/Amblyraja_radiata.bed \
> tmp.bed && mv tmp.bed bed/Amblyraja_radiata.bed

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}
```
## Raja brachyura
```
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Raja_brachyura/
cd ${WORKDIR}

mkdir OG_Naming
mv * OG_Naming
mv OG_Naming/bed .
mv OG_Naming/peptide .

# Replace chr X with Unnamed

awk 'BEGIN{OFS="\t"} {$1=($1=="X"?"Unnamed":$1); print}' \
bed/Raja_brachyura.bed \
> tmp.bed && mv tmp.bed bed/Raja_brachyura.bed

# Replace chr 46 with X
awk 'BEGIN{OFS="\t"} {$1=($1=="40"?"X":$1); print}' \
bed/Raja_brachyura.bed \
> tmp.bed && mv tmp.bed bed/Raja_brachyura.bed

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}
```
# Mammals
## No annotated sex chr:
### Marsupial
```

WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Marsupial_NoSexAnn
mkdir -p ${WORKDIR}
cd ${WORKDIR}

mkdir -p bed peptide

while read -r SPECIES; do
  cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed" bed/
  cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa" peptide/
done <<'EOF'
Dasyurus_maculatus
Sarcophilus_harrisii
Sminthopsis_crassicaudata
EOF


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}
```
### Eutherian
```
source myconda
mamba activate genespace

WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Eutherians_NoSexAnn
mkdir -p ${WORKDIR}
cd ${WORKDIR}

mkdir -p bed peptide

while read -r SPECIES; do
  cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed" bed/
  cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa" peptide/
done <<'EOF'
Acomys_minous
Doryrhina_cyclops
Rhinolophus_yonghoiseni
Antrozous_pallidus
Eptesicus_fuscus
Myotis_daubentonii
EOF

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}


########################################################################
# Findings
########################################################################
Acomys minous -- chr 7 = X
Doryrhina cyclops -- chr 7 = X
Rhinolophus_yonghoiseni -- chr 8 = X
Antrozous pallidus -- chr 2 = X
Eptesicus fuscus -- chr 1 = X


########################################################################
# Acomys minous
########################################################################
SPECIES=Acomys_minous
X_CHR_NUM=7

while IFS=, read -r SPECIES X_CHR_NUM; do
    WORKDIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}"
    mkdir -p "$WORKDIR"
    cd "$WORKDIR" || exit 1
    
    mkdir -p bed peptide

    while read -r SP; do
        cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SP}.bed" bed/
        cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SP}.fa" peptide/
    done <<EOF
${SPECIES}
Gallus_gallus
EOF

    awk -v xchr="$X_CHR_NUM" 'BEGIN{OFS="\t"} {$1=($1==xchr?"X":$1); print}' \
        "bed/${SPECIES}.bed" > tmp && mv tmp "bed/${SPECIES}.bed"

    cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts || exit 1
    sbatch 1g_genespace.sh "$WORKDIR"
done <<EOF
Furcifer_pardalis,10
EOF
Doryrhina_cyclops,7
Rhinolophus_yonghoiseni,8
Antrozous_pallidus,2
Eptesicus_fuscus,1
EOF

Acomys_minous,7


```
# Squamates
## Furcifer pardalis
Run, rename chr 10 to Z and rerun again

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

sbatch 1b_genespace.singleinstance.sh Gallus_gallus Furcifer_pardalis /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/





# Run all of these from the beginning (Genome to bed to genespace)
## Mammals
Elephas maximus indicus # liftover annotations
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

sbatch 1b_genespace.singleinstance.sh Gallus_gallus Elephas_maximus_indicus /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared_liftover


Aulostomus maculatus # rename chr to sex chr
Lepisoteus oculatus # rename chr to sex chr
```
mamba activate genespace

SPECIES=Lepisosteus_oculatus
OUTDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/${SPECIES}/
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/${SPECIES}/ncbi_dataset/data/GC*/
transeq -sequence cds_from_genomic.fna -outseq cds_from_genomic.translated.cds
ln -sf /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Lepisosteus_oculatus/ncbi_dataset/data/GCF_040954835.1/cds_from_genomic.translated.cds ${OUTDIR}/${SPECIES}.translated.cds

SPECIES=Aulostomus_maculatus
OUTDIR=/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/${SPECIES}/
cd /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/${SPECIES}/ncbi_dataset/data/GC*/
transeq -sequence cds_from_genomic.fna -outseq cds_from_genomic.translated.cds
ln -sf /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Aulostomus_maculatus/ncbi_dataset/data/GCA_048301465.1/cds_from_genomic.translated.cds ${OUTDIR}/${SPECIES}.translated.cds

```
## Prepare all bed and peptide files
```

library(GENESPACE)

chrIdDictionary <- c(
  "NC_090696.1" = "1",
  "NC_090697.1" = "2",
  "NC_090698.1" = "3",
  "NC_090699.1" = "4",
  "NC_090700.1" = "5",
  "NC_090701.1" = "6",
  "NC_090702.1" = "7",
  "NC_090703.1" = "8",
  "NC_090704.1" = "9",
  "NC_090705.1" = "10",
  "NC_090706.1" = "11",
  "NC_090707.1" = "12",
  "NC_090708.1" = "13",
  "NC_090709.1" = "14",
  "NC_090710.1" = "15",
  "NC_090711.1" = "16",
  "NC_090712.1" = "17",
  "NC_090713.1" = "18",
  "NC_090714.1" = "19",
  "NC_090715.1" = "20",
  "NC_090716.1" = "21",
  "NC_090717.1" = "22",
  "NC_090718.1" = "23",
  "NC_090719.1" = "24",
  "NC_090720.1" = "25",
  "NC_090721.1" = "26",
  "NC_090722.1" = "27",
  "NC_090723.1" = "28",
  "NC_090724.1" = "29"
)

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/",
  genomeDirs    = "Lepisosteus_oculatus",
  genomeIDs     = "Lepisosteus_oculatus",
  gffString     = "gff",
  faString      = "translated.cds",
  presets       = "none",
  genespaceWd   = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/bed_chrfixed/",
  gffIdColumn      = "gene",
  headerSep        = " ",
  headerEntryIndex = 2,
  headerStripText  = "gene=|\\[|\\]",
  gffStripText     = "",
  chrIdDictionary  = chrIdDictionary
)


chrIdDictionary <- c(
  "CM107617.1" = "1",
  "CM107618.1" = "2",
  "CM107619.1" = "3",
  "CM107620.1" = "4",
  "CM107621.1" = "5",
  "CM107622.1" = "6",
  "CM107623.1" = "7",
  "CM107624.1" = "8",
  "CM107625.1" = "9",
  "CM107626.1" = "10",
  "CM107627.1" = "11",
  "CM107628.1" = "12",
  "CM107629.1" = "13",
  "CM107630.1" = "14",
  "CM107631.1" = "15",
  "CM107632.1" = "16",
  "CM107633.1" = "17",
  "CM107634.1" = "18",
  "CM107635.1" = "19",
  "CM107636.1" = "20",
  "CM107637.1" = "21",
  "CM107638.1" = "22",
  "CM107639.1" = "23",
  "CM107640.1" = "24"
)

parsedPaths2 <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/",
  genomeDirs    = "Aulostomus_maculatus",
  genomeIDs     = "Aulostomus_maculatus",
  gffString     = "gff",
  faString      = "translated.cds",
  presets       = "none",
  genespaceWd   = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/bed_chrfixed/",
  gffIdColumn        = "protein_id",
  headerSep          = " ",
  headerEntryIndex   = 1,
  headerStripText    = "^.*_cds_|_[0-9]+_[0-9]+$",
  gffStripText       = "",
  chrIdDictionary  = chrIdDictionary
)
```

## Remove scaffolds
```
grep -v 'NW' Lepisosteus_oculatus.bed > 
```
# Lepisosteus_oculatus
## Create requisite genespace directory structure and submit genespace
```
SPECIES=Lepisosteus_oculatus

GENESPACEDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/
cd ${GENESPACEDIR}/${SPECIES}
grep -v 'NW' /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/bed_chrfixed/bed/${SPECIES}.bed > bed/${SPECIES}.bed

# Replace chr 14 with X1
awk 'BEGIN{OFS="\t"} {$1=($1=="14"?"X1":$1); print}' \
bed/Lepisosteus_oculatus.bed \
> tmp.bed && mv tmp.bed bed/Lepisosteus_oculatus.bed

# Replace chr 18 with X2
awk 'BEGIN{OFS="\t"} {$1=($1=="18"?"X2":$1); print}' \
bed/Lepisosteus_oculatus.bed \
> tmp.bed && mv tmp.bed bed/Lepisosteus_oculatus.bed


# 2) Filter peptide FASTAs to match filtered BEDs (keep only genes still present)
#    - assumes peptide file name matches bed base name: SPECIES.bed <-> SPECIES.fa
  base="$(basename "$bed" .bed)"
pep_in="${INPEP}/${base}.fa"
pep_out="${OUTPEP}/${base}.fa"

# IDs to keep = col 4 from filtered BED
keep_ids="$(mktemp)"
awk -F'\t' '$4!=""{print $4}' bed/Lepisosteus_oculatus.bed | sort -u > "$keep_ids"

# Subset FASTA by header name (requires seqkit)
seqkit grep -n -f "$keep_ids" peptide/Lepisosteus_oculatus.fa > peptide/Lepisosteus_oculatus.fa.tmp

rm -f "$keep_ids"

mv peptide/Lepisosteus_oculatus.fa.tmp peptide/Lepisosteus_oculatus.fa

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

sbatch 1b_genespace.singleinstance.sh Gallus_gallus ${SPECIES} ${GENESPACEDIR}
```
# Aulostomus_maculatus
## Create requisite genespace directory structure and submit genespace
```
SPECIES=Aulostomus_maculatus

GENESPACEDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/
mkdir -p ${GENESPACEDIR}/${SPECIES}
cd ${GENESPACEDIR}/${SPECIES}
mkdir -p bed peptide
grep -v "JB" /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/bed_chrfixed/bed/${SPECIES}.bed > bed/${SPECIES}.bed

# Replace chr 14 with Z
awk 'BEGIN{OFS="\t"} {$1=($1=="14"?"Z":$1); print}' \
bed/Aulostomus_maculatus.bed \
> tmp.bed && mv tmp.bed bed/Aulostomus_maculatus.bed

# copy peptide file
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/usable_peptide/${SPECIES}.fa peptide/

# copy chicken files
cp ../Lepisosteus_oculatus/bed/Gallus_gallus.bed bed/
cp ../Lepisosteus_oculatus/peptide/Gallus_gallus.fa peptide/

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

sbatch 1b_genespace.singleinstance.sh Gallus_gallus ${SPECIES} ${GENESPACEDIR}
```
15979846
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared//Lepisosteus_oculatus
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared//plots//Lepisosteus_oculatus.Gallus_gallus.5x3.wholegenome.pdf



# Tiliqua_scincoides
SPECIES=Tiliqua_scincoides
WORKDIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1
mkdir -p bed peptide

cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed" bed/
cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa" peptide/


cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/Gallus_gallus.bed" bed/
cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/Gallus_gallus.fa" peptide/

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts || exit 1

sbatch 1g_genespace.sh "$WORKDIR"

while IFS=, read -r SPECIES X_CHR_NUM; do
    WORKDIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}"
    mkdir -p "$WORKDIR"
    cd "$WORKDIR" || exit 1
    
    mkdir -p bed peptide

    while read -r SP; do
        cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SP}.bed" bed/
        cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SP}.fa" peptide/
    done <<EOF
${SPECIES}
Gallus_gallus
EOF

    awk -v xchr="$X_CHR_NUM" 'BEGIN{OFS="\t"} {$1=($1==xchr?"X":$1); print}' \
        "bed/${SPECIES}.bed" > tmp && mv tmp "bed/${SPECIES}.bed"

    cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts || exit 1
    sbatch 1g_genespace.sh "$WORKDIR"
done <<EOF
Tiliqua_scincoides,7
EOF

# Ambystoma_mexicanum_x_Amblystoma_tigrinum

SPECIES=Ambystoma_mexicanum_x_Ambystoma_tigrinum
WORKDIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}"
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1
mkdir -p bed peptide

cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed" bed/
cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa" peptide/


cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/Gallus_gallus.bed" bed/
cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/Gallus_gallus.fa" peptide/

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts || exit 1

sbatch 1g_genespace.sh "$WORKDIR"

# syntenic with chicken 7, 9, and 24
# syntenic with newt LG736
while IFS=, read -r SPECIES X_CHR_NUM; do
    WORKDIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}"
    mkdir -p "$WORKDIR"
    cd "$WORKDIR" || exit 1
    
    mkdir -p bed peptide

    while read -r SP; do
        cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SP}.bed" bed/
        cp "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SP}.fa" peptide/
    done <<EOF
${SPECIES}
Gallus_gallus
EOF

    awk -v xchr="$X_CHR_NUM" 'BEGIN{OFS="\t"} {$1=($1==xchr?"X":$1); print}' \
        "bed/${SPECIES}.bed" > tmp && mv tmp "bed/${SPECIES}.bed"

    cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts || exit 1
    sbatch 1g_genespace.sh "$WORKDIR"
done <<EOF
Tiliqua_scincoides,7
EOF