# Infer PAR boundaries from sequence alignment

## 1. Align Y/W to X/Z chromosome to identify putative PAR regions

### A. Align the sex chromosomes and output bed files with regions that matched
Simone performed the alignment -- add in her code here.

### Join contiguous regions at either end of the sex chromosomes
The above code output a collection of bedfiles, found here:
```
export BEDDIR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/align_PAR_99thr/align_PAR_inference
```
There should be only one bed per genome. To confirm, I ran the following:

``` 
ls $BEDDIR | awk -F '_' '{print $1,$2}' 
```

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
The fastatofasta.aln and YchromotoXchrom files are duplicates of each other.

#### Create a list that excludes duplicates:
```
printf "%s\n" "$BEDDIR"/* | grep -v 'YchromtoXchrom' | grep -v 'fastatofasta' > /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/alignment_beds.txt
```

Check to confirm that all bed files only contain one chromsome each:
```
while read -r bed; do
    # Count unique chromosomes (first column)
    num_chr=$(awk '!/^#/ && NF {print $1}' "$bed" | sort -u | wc -l)

    if [ "$num_chr" -gt 1 ]; then
        echo "$bed"
    fi
done < /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/alignment_beds.txt
```
Sort all bed files, and merge overlapping regions
```
chmod +x 1x_sort_merge_regions.sh
./1x_sort_merge_regions.sh
```
Fix errors in bed file names. Fringilla is misnammed, two copies of the Tammar wallaby exist (Notamacropus and Macropus), and the PAR_job_list file was carried over into the original bed list and needs to not exist here.

one misnamed bed file and remove a file made in error:
```
mv sorted_beds/Fringilla_coelebs_WtoZ.aln.id99.len10k.refqry.bed sorted_beds/Tamandua_tetradactyla_YtoX.aln.id99.len10k.refqry.bed
rm sorted_beds/Notamacropus_eugenii_YtoX.aln.id99.len10k.refqry.bed

```

Compute the distance between each region of the bed files across all bed files:
```
chmod +x 1x_compute_bed_spacers.sh
./1x_compute_bed_spacers.sh

```
## Infer if chromosome confidently has only one PAR. If all regions of interest are on only one side of the midpoint of the chromosome, infer that these are all one par.
Create a reference data table that links file names to SEXCHRLEN file, which will have the following structure:

| Accession | Chromosome_name | Seq_length | Species | Filename |
| --- | --- | --- | --- | --- |
| Assembly accession ID (e.g., GCA_*) | Name of the chromosome (e.g., W, Z, X, etc.) | Length of the chromosome sequence | Species name associated with the chromosome | BED filename from the `sorted_beds` directory |


```
chmod +x 1x_combine_sex_chr_len_beds.sh
./1x_combine_sex_chr_len_beds.sh

```
This file (bed_sexchr_info.tsv) has the following structure:
```
Accession       GenBank_accession       RefSeq_accession        Seq_length      Species Filename
```
### Confirm that either the Genbank or RefSex accession matches the chromosome name in the bed file. 
```
Rscript 1x_confirm_chromosome_match.R
```
Genomes that did not match:
| Species | Chrom in bed | GenBank in metadata | Cause of error | Solution |
| --- | --- | --- | --- | --- |
| Fringilla   | CM043353.1 | OY740727.1 | Mislabelled Tamandua tetradactyla | rename to Tamandua, Simone is realigning |
| Meles meles |  NC_060087.1 | OV277448.1 |  matches to another haplotype.  Chr length of NC_060087: 132,955,750; of PV277448.1 132,955,750 | Nothing -- this one is fine |
| Pan_troglodytes | NC_072421.2 | CM054457.2 | Is refseq chrom, just missing from metadata file | Nothing - is fine |
| Patagioenas_fasciata | NC_089087.1 | CM104797.1 |  Matches with an old version of the genome https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036971685.1/. chr len metadata: 85,896,199; chr len current version: 85,896,199; chr len previous version: 85,910,070 | ? |
| Pongo_abelii | NC_072008.2 | CM054702.2 | Is refseq chrom, just missing from metadata file | Nothing - is fine |
| Pongo_pygmaeus | NC_072396.2 | CM054653.2 | Is refseq chrom, just missing from metadata file | Nothing - is fine |
| Sus_scrofa | CM060728.1 | CM107810.1 | CM107810 is in T2T_pig on ncbi with corresponding GCA and is 137,020,761 bp ; CM060728 is a different genome on NCBI and is 135,116,946 bp (https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_030704935.2/) | ? |

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/
mv Fringilla_coelebs_YtoX.aln.paf Tamandua_tetradactyla_YtoX.aln.paf
mv Notamacropus_eugenii_YtoX.aln.paf Macropus_eugenii_YtoX.aln.paf
mkdir -p junk_continuous_percentID
mv Bos_taurus_YtoX.aln.paf junk_continuous_percentID/
mv Dama_dama_YtoX.aln.paf junk_continuous_percentID/
mv Lycaon_pictus_YtoX.aln.paf junk_continuous_percentID/

mv Pongo_abelii_fastatofasta.aln.paf junk_continuous_percentID/
mv Pongo_pygmaeus_fastatofasta.aln.paf junk_continuous_percentID/
mv Pseudorca_crassidens_fastatofasta.aln.paf junk_continuous_percentID
mv Rhynchocyon_petersi_YchromtoXchrom.aln.paf junk_continuous_percentID/

mv Monodelphis_domestica_YtoX.aln.paf junk_continuous_percentID/
grep -v 'NW_' junk_continuous_percentID/Monodelphis_domestica_YtoX.aln.paf > Monodelphis_domestica_YtoX.aln.paf

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/continuous_percentID
mv Fringilla_coelebs_YtoX.refqry.csv Tamandua_tetradactyla_YtoX.aln.refqry.csv
mv Pongo_abelii_fastatofasta.aln.refqry.csv ../junk_continuous_percentID/
mv Pongo_pygmaeus_fastatofasta.aln.refqry.csv ../junk_continuous_percentID/
mv Pseudorca_crassidens_fastatofasta.aln.refqry.csv ../junk_continuous_percentID

mv Monodelphis_domestica_YtoX.aln.refqry.csv ../junk_continuous_percentID/
grep -v 'NW_' ../junk_continuous_percentID/Monodelphis_domestica_YtoX.aln.refqry.csv > Monodelphis_domestica_YtoX.aln.refqry.csv

# Next steps (don't do this)
 - Identify all bed files with regions on only one half of the chromosome. 
    - If the regions are closer to start of the chromosome, and then output a new bed file to dir Dir/one_PAR that is just chr start to highest end pos in region file. 
    - If the regions are closer to the end of the chromosome, output a new bed file to Dir/one_PAR that is lowest start pos in region file to chromosome end
 - Then for all others, test if the largest spacer is >0.5 of the chromosome. 
    - Make a bed file in Dir/two_PAR that is start of the chromosome to start of largest spacer (one line of bed) and end of largest spacer to end of chromosome (second line of bed)

If there are ones that remain... messy.


# Next steps 
 - Choose genomes with high quality PAR inferences using chromosome_regions.png
 - Compare to TE approach
 - Assume that Y and W are poorly generated
 - PAR collapse onto the X or Z
 - Do human, canis lupis (not bad one), maybe apes too

 Half-deep:
 1. Byung June Penn with Kateryna
    Ran half deep on every chromosome of every genome if 


Rscript X_continuous_identity_plotting.r
Rscript X_continuous_identity_plotting.thresholds.r
Rscript X_continuous_identity_plotting.tech.r
Rscript X_AlignmentTech.r





















# Junk scripts (for now)
Make a histogram in R:
```
module load R/4.5.0
R

# Read CSV
df <- read.csv("bed_spacers.csv", header = FALSE,
               col.names = c("file", "positions", "distance"))

# Convert distance to numeric if needed
df$distance <- as.numeric(df$distance)

# Output file
png("bed_spacers_hist.png", width = 1200, height = 900)

# Histogram
hist(df$distance,
     breaks = 100,
     main = "Histogram of BED Region Spacers",
     xlab = "Distance Between Regions",
     col = "steelblue",
     border = "white")

dev.off()

```
This histogram is pretty useless because the data is heavily skewed toward 0. I think a summary stat approach would work better.

Make a file ```bed_spacer_stats.csv``` that has the following structure:
| file | spacer_count | max_spacer_size | second_max | diff_largest_second | average |
| --- | --- | --- | --- | --- | --- |
| file name | number of spacers | maximum spacer size | second max spacer size | $3 - $4 | average spacer size
```
echo "file,spacer_count,max_spacer_size,second_max,third_max,diff_largest_second,diff_second_third,delta_diff,average" \
    > bed_spacer_stats.csv

awk -F',' -v OFS=',' '
{
    file = $1
    dist = $3 + 0

    count[file]++
    sum[file] += dist

    # track top three distances
    if (dist > max1[file]) {
        max3[file] = max2[file]
        max2[file] = max1[file]
        max1[file] = dist
    } else if (dist > max2[file]) {
        max3[file] = max2[file]
        max2[file] = dist
    } else if (dist > max3[file]) {
        max3[file] = dist
    }
}
END {
    for (f in count) {
        c = count[f]
        regions = c + 1

        m1 = max1[f]
        if (c > 1) {
            m2 = max2[f]
        } else {
            m2 = m1
        }
        if (c > 2) {
            m3 = max3[f]
        } else {
            m3 = m2
        }

        d12 = m1 - m2                 # max - 2nd max
        d23 = m2 - m3                 # 2nd max - 3rd max
        delta = d12 - d23             # (max-2nd) - (2nd-3rd)
        avg = sum[f] / c

        printf "%s,%d,%d,%d,%d,%d,%d,%d,%.6f\n", \
               f, regions, m1, m2, m3, d12, d23, delta, avg
    }
}
' bed_spacers.csv >> bed_spacer_stats.csv


```
Make histograms of the maximum spacer size and the difference between the max and 2nd max

```
R

# Read stats file
df <- read.csv("bed_spacer_stats.csv", header = TRUE)

# -------------------------
# Histogram 1: max_spacer_size
# -------------------------

png("max_spacer_size_hist.png", width = 1200, height = 900)
hist(df$max_spacer_size,
     breaks = 100,
     main = "Histogram of Max Spacer Size",
     xlab = "Max Spacer Size",
     col = "steelblue",
     border = "white")
dev.off()

# -------------------------
# Histogram 2: diff_largest_second
# -------------------------

png("diff_largest_second_hist.png", width = 1200, height = 900)
hist(df$diff_largest_second,
     breaks = 500,
     main = "Histogram of Difference Between Largest & Second Largest Spacer",
     xlab = "Difference (max - second max)",
     col = "darkorange",
     border = "white")
dev.off()

# -------------------------
# Histogram 3: delta_diff
# -------------------------

png("delta_diff.png", width = 1200, height = 900)
hist(df$delta_diff,
     breaks = 500,
     main = "Histogram of Delta Difference",
     xlab = "delta difference (max - 2nd max) - (2nd max - 3rd max)",
     col = "darkorange",
     border = "white")
dev.off()

```
Infer natural breaks
```
library(classInt)

vals <- df$diff_largest_second

# Find natural breaks into 2 or 3 groups
classIntervals(vals, n = 25, style = "jenks")$brks




Using $SEXCHRLEN, compute 

```
cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference

export SEXCHRLEN=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_phylogeny/VGP_freeze_hap1_combined_sexchroms_seq_reports.tsv


```
export TREE=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/roadies_v1.1.16b.nwk
export ANN=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/annotations.tsv
