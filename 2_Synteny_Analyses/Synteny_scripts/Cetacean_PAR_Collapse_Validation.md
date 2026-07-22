# A. Cetacean Y-to-X PAR inference plot

# B. Cetacean synteny PAR validation plot

## 1 Prepare bed and peptide for species with missing files (Inia geoffrensis)

### Write R script for parsing gff

```
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
repo <- args[2]
stopifnot(!is.na(outdir), nzchar(outdir))

SPECIES <- sort(basename(list.dirs(repo, full.names = TRUE, recursive = FALSE)))

library(GENESPACE)

parsedPaths2 <- parse_annotations(
  rawGenomeRepo = repo,
  genomeDirs    = SPECIES,
  genomeIDs     = SPECIES,
  gffString     = "gff",
  faString      = "translated_cds",
  presets       = "none",
  genespaceWd   = outdir,
  gffIdColumn      = "orig_transcript_id",
  headerSep        = " ",
  headerEntryIndex = 1,

  headerStripText = "^rna-|_[0-9]+$",
  gffStripText     = ".*\\|"
)
```
### Parse GFF
```

mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans

source myconda
mamba activate genespace

# translate nucleic acid seq 

transeq -sequence Inia_geoffrensis.cds -outseq Inia_geoffrensis.translated_cds

Rscript prepare_genespace.R /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/Inia_working/ /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/Inia_working/

```
# 2. Filter genomes to X and Y chromosomes
```

cat > species.txt <<'EOF'
Tursiops_truncatus
Pseudorca_crassidens
Lagenorhynchus_acutus
Globicephala_melas
Grampus_griseus
Inia_geoffrensis
Delphinus_delphis
Balaenoptera_musculus
Stenella_coeruleoalba
Eschrichtius_robustus
Megaptera_novaeangliae
Balaenoptera_physalus
Mesoplodon_densirostris
Mesoplodon_bidens
Eubalaena_glacialis
Monodon_monoceros
EOF


mkdir -p genespace/bed genespace/peptide 

for SPECIES in \
    Tursiops_truncatus \
    Pseudorca_crassidens \
    Globicephala_melas \
    Delphinus_delphis \
    Balaenoptera_musculus \
    Eschrichtius_robustus \
    Megaptera_novaeangliae \
    Mesoplodon_densirostris \
    Eubalaena_glacialis
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/usable_bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/genespace/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/usable_peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/genespace/peptide/
done


for SPECIES in \
    Lagenorhynchus_acutus \
    Grampus_griseus \
    Inia_geoffrensis \
    Stenella_coeruleoalba \
    Balaenoptera_physalus \
    Mesoplodon_bidens 
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/bed_chrfixed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/genespace/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/genespace/peptide/
done

for SPECIES in \
    Inia_geoffrensis
do
  cp Inia_working/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/genespace/bed/

  cp Inia_working/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/genespace/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Cetaceans/genespace


```

