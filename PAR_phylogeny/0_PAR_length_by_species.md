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


Why isn't ostrich in it?
Are the Xs technical or biological?
Spot check neosex chromosomes -- at minimum star the Neosex
  See Simone's annotation of the
Rather than PAR ratio, add size of the Y
Violin plot across species of PAR ratios in major groupings
Is the Y annotated for each of these? 
Add length of Y and add color of sequencing technology
Add whether there is a telomere
Standardize 
Potentially add length of non-PAR X/Z as stacked gray
Important to 


## Investigating issues, December 15th
# Which species have more than one bed file and why?
``` ls | awk -F '_' '{print $1,$2}' ```

These have more than one bed file:
```
Artibeus intermedius # Y1toX and Y2toX
Artibeus lituratus # Y1toX and Y2toX
Canis lupus # two subspecies -- baileyi and orion
Hoplias malabaricus # YtoX1 and YtoX2
Hypanus sabinus # YtoX1 and YtoX2
Pongo abelii # fastatofasta.aln and YtoX
Pongo pygmaeus # fastatofasta.aln and YtoX
Pseudorca crassidens # fastatofasta.aln and YtoX
Rhynchocyon petersi # YchromtoXchrom and YtoX
```
The fastatofasta.aln and YchromotoXchrom files are duplicates

Sort each PAR file
Collect the value that is the distance between the end of one region and the next
Make a histogram of these distances
Identify a cutoff from the histogram that can be used to separate PARs on each end of the chromosome
Write it up as both a methods section and a protocol on github

Synteny Pipeline
1. Go through list of neo sex column 
1a. Identify list(s) of genomes that we care about
2. Liftover annotations to unannotated genomes of interest
3. Engage Adrian 