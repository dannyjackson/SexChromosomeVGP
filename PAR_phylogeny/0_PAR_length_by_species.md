# Vertebrate PAR lengths testing
This script creates plots demonstrating the length of pseudoautosomal regions (PARs) across vertebrate clades using VGP data from the 2025 data freeze.

*Created: December 3rd, 2025* <br>
*Last edited: December 3rd, 2025*
```
cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/
mkdir analyses referencelists datafiles
cd datafiles

curl -L -o roadies_v1.1.16b.nwk \
  https://raw.githubusercontent.com/VGP/vgp-trees/main/phase-1/roadies_v1.1.16b.nwk

curl -L -o annotations.tsv \
  https://raw.githubusercontent.com/VGP/vgp-trees/main/phase-1/annotations.tsv
```

I also manually copied over a file PAR_annotations.zip from google drive:
https://drive.google.com/file/d/1TPe9rXQc82eJWCov5zCq8CTgf1iRJYjO/view?usp=sharing

### Plot the trees using specific clades 
```

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_phylogeny

module load R

export TREE=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/roadies_v1.1.16b.nwk
export ANN=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/annotations.tsv
export BEDDIR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/PAR_annotations
export SEXCHRLEN=VGP_freeze_hap1_combined_sexchroms_seq_reports.tsv

Rscript 0a_ggplotTree_PARs.R Birds
Rscript 0a_ggplotTree_PARs.R Mammals
Rscript 0a_ggplotTree_PARs.R Fishes
Rscript 0a_ggplotTree_PARs.R Amphibians
Rscript 0a_ggplotTree_PARs.R Reptiles

Add * for we should look at it more
Add X for nothing here at all



```
All of these have complex sex chromosome systems. Sort out how to plot them such that there are two bars for these (at half the width).

Warning messages:
1: In get_par_segment_lengths(sci) :
Multiple BED files for Pseudorca crassidens; using first: Pseudorca_crassidens_fastatofasta.aln.id98_5.len10k.refqry.bed
2: In get_par_segment_lengths(sci) :
Multiple BED files for Artibeus lituratus; using first: Artibeus_lituratus_Y2toX.aln.id98_5.len10k.refqry.bed
3: In get_par_segment_lengths(sci) :
Multiple BED files for Artibeus intermedius; using first: Artibeus_intermedius_Y1toX.aln.id98_5.len10k.refqry.bed
4: In get_par_segment_lengths(sci) :
Multiple BED files for Rhynchocyon petersi; using first: Rhynchocyon_petersi_YchromtoXchrom.aln.id98_5.len10k.refqry.bed
Wrote tree + PAR bar plot + small legend to: subtree_mammals_tree_with_PARbars_byOrder.pdf 

```
Rscript ggplotTree_PARs.R birds
```

Issues to clean up:
2. The genomes with complex sex chromosomes (i.e. X1X2) are skipping the latter and just using the former.
3. Overlap between plots is imperfect and hiding tip labels.
4. Add a version with ratio of PAR to entire sex chromosome (probably entire X/Z?)