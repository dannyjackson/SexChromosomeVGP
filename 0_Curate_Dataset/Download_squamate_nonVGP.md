# Download squamate genomes that are not in VGP
Xantusia_vigilis,GCA_051940855.1
Cyrtopodion_scabrum,GCA_054791115.1
Pogona_vitticeps,GCF_051106095.1
Paroedura_picta,GCF_049243985.1
Hemicordylus_capensis,GCF_027244095.1
Chamaeleo_calyptratus,GCA_043643385.1
Sphaerodactylus_townsendi,GCF_021028975.2
Varanus_acanthurus,GCA_050042745.1


## Download all genomes and annotation files from the most recent GCF (preferred) or, if no GCF exists, most recent GCA
```
source myconda
mamba activate ncbi_datasets

GENOME_DIR=/data/Wilson_Lab/data/VGP_genomes_phase1/squamate_nonVGP_genomes
mkdir -p "${GENOME_DIR}"
FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"

# Set this somewhere above if it is not already set
# FILES_TO_DOWNLOAD="genome,gff3,rna,cds,protein,seq-report"

while IFS=, read -r SPECIES_NAME ACCESSION; do
  OUTDIR="${GENOME_DIR}/${SPECIES_NAME}"
  mkdir -p "${OUTDIR}"

  ZIPFILE="${OUTDIR}/${ACCESSION}.zip"

  echo "=== ${SPECIES_NAME} -> ${ACCESSION}"
  echo "Downloading: ${FILES_TO_DOWNLOAD}"

  datasets download genome accession "${ACCESSION}" \
    --include "${FILES_TO_DOWNLOAD}" \
    --filename "${ZIPFILE}"

done <<'EOF'
Xantusia_vigilis,GCA_051940855.1
Cyrtopodion_scabrum,GCA_054791115.1
Pogona_vitticeps,GCF_051106095.1
Paroedura_picta,GCF_049243985.1
Hemicordylus_capensis,GCF_027244095.1
Chamaeleo_calyptratus,GCA_043643385.1
Sphaerodactylus_townsendi,GCF_021028975.2
Varanus_acanthurus,GCA_050042745.1
EOF


while IFS=, read -r SPECIES_NAME ACCESSION; do
  OUTDIR="${GENOME_DIR}/${SPECIES_NAME}"
  mkdir -p "${OUTDIR}"

  cd ${OUTDIR}
  ZIPFILE="${OUTDIR}/${ACCESSION}.zip"
  unzip "${ZIPFILE}"

done <<'EOF'
Xantusia_vigilis,GCA_051940855.1
Cyrtopodion_scabrum,GCA_054791115.1
Pogona_vitticeps,GCF_051106095.1
Paroedura_picta,GCF_049243985.1
Hemicordylus_capensis,GCF_027244095.1
Chamaeleo_calyptratus,GCA_043643385.1
Sphaerodactylus_townsendi,GCF_021028975.2
Varanus_acanthurus,GCA_050042745.1
EOF
```
All species:
# non VGP -- annotated
Pogona_vitticeps,GCF_051106095.1,Has_GFF
Paroedura_picta,GCF_049243985.1,Has_GFF
Hemicordylus_capensis,GCF_027244095.1,Has_GFF
Sphaerodactylus_townsendi,GCF_021028975.2,Has_GFF

# Liftover these
## non VGP -- unannotated
Xantusia_vigilis,GCA_051940855.1,Lacks_GFF
Cyrtopodion_scabrum,GCA_054791115.1,Lacks_GFF
Chamaeleo_calyptratus,GCA_043643385.1,Lacks_GFF
Varanus_acanthurus,GCA_050042745.1,Lacks_GFF

# Build liftover db
cd /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/scripts

## Run on one species to build gff_db
```
while IFS=, read -r SPECIES_NAME ACCESSION; do
  sbatch liftover.squamates.nonVGP.sh -s "$SPECIES_NAME" -a "$ACCESSION"
done <<'EOF'
Xantusia_vigilis,GCA_051940855.1
EOF
```
# perform liftover
```
while IFS=, read -r SPECIES_NAME ACCESSION; do
  sbatch liftover.squamates.nonVGP.useDB.sh -s "$SPECIES_NAME" -a "$ACCESSION"
done <<'EOF'
Cyrtopodion_scabrum,GCA_054791115.1
Chamaeleo_calyptratus,GCA_043643385.1
EOF

# Chamaeleo_calyptratus,GCA_043643385.1
# Varanus_acanthurus,GCA_050042745.1ls 

```
## in VGP

while IFS=, read -r SPECIES_NAME ACCESSION; do
  sbatch liftover.squamates.VGP.useDB.sh -s "$SPECIES_NAME" -a "$ACCESSION"
done <<'EOF'
Anniella_stebbinsi,GCA_051312515.2
Shinisaurus_crocodilurus,GCA_021292165.1
Furcifer_pardalis,GCA_030440675.1
Cyclura_pinguis,GCA_030412105.1
EOF

