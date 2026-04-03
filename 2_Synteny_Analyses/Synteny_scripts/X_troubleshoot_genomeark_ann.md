# TROUBLESHOOT GENESPACE

# Does haplotype have X/Z?
## No:
```
Anniella stebbinsi # Download hap 2
Doryrhina cyclops # no sex chrs
```
## Genome file is empty:
```
Artibeus lituratus # GCA_038363095.4
```
## Maybe it's bc of sex chrs
```
Hoplias malabaricus
```
# Is it an issue with the annotation?
## Needs sex chrs annotated:
```
Panthera onca 
```
# Not present in usable bed
## NCBI annot
```
Peromyscus maniculatus
Lemur catta
```
## Genomeark annot
```
Monodon monocero
Lycocorax pyrrhopterus
Willisornis vidua
Larus fuscus
Anas platyrhynchos # missing annotation?
Echiichthys vipera
Inia geoffrensis
Corynorhinus townsendii
Myotis emarginatus
Myotis mystacinus
Pipistrellus pygmaeus
Vespertilio murinus
Phascolarctos cinereus
Tupaia tana
Podarcis siculus
```
# Present in usable bed
## Error found in slurm output:
```
Chlamydotis macqueenii #  ERROR: diamond makedb failed
Girardinichthys multiradiatus # Error in rbindlist
Amblyraja radiata # Error in rbindlist
Tursiops truncatus # Error in rbindlist
```
## No error:
```
Ammospiza maritima
Mergus octosetaceus
```

# Troubleshoot:
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/

## Echiichthys_vipera
Notes: Genespace ran but it broke on MCScanX I think? No clear reason why, I reran genespace single threaded and it still broke.

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/Echiichthys_vipera

ls /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Echiichthys_vipera/ncbi_dataset/data/GCA_963691815.1/
ls /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/fEchVip8/

cp /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/fEchVip8/fEchVip8.gff Echiichthys_vipera.gff
cp /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/fEchVip8/fEchVip8.cds Echiichthys_vipera.translated.cds

/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/ /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/

source myconda
mamba activate genespace_py3.10

library(GENESPACE)

SPECIES <- sort(basename(list.dirs(repo, full.names = TRUE, recursive = FALSE)))

parsedPaths2 <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting",
  genomeDirs    = "Echiichthys_vipera",
  genomeIDs     = "Echiichthys_vipera",
  gffString     = "gff",
  faString      = "cds",
  presets       = "none",
  genespaceWd   = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting",
  gffIdColumn        = "protein_id",
  headerSep          = " ",
  headerEntryIndex   = 1,
  headerStripText    = ".*\\|",   # strip everything up to last '|'
  gffStripText       = ".*\\|"    # same normalization on protein_id values
)

```
library(GENESPACE)
library(ggplot2)
library(data.table)

wd="/vf/users/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting"

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/"
)

out <- run_genespace(gpar, overwrite = FALSE)

gpar <- init_genespace(
  wd = wd,
  nCores = 16,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/"
)

```
```
#!/bin/bash
#SBATCH --job-name=genespace_Echiichthys_vipera
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

source myconda
mamba activate genespace_py3.10

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/scripts

Rscript genespace.Echiichthys_vipera.R
```
cp ../sexshared/Echiichthys_vipera/bed/Gallus_gallus.bed bed/
cp ../sexshared/Echiichthys_vipera/peptide/Gallus_gallus.fa peptide/