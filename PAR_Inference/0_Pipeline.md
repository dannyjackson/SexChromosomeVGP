# Working Pipeline

## Goal:
Summarize the state of the data with respect to PAR boundaries on sex chromosomes, infer PAR boundaries where absent (if possible), and produce a protocol for improving PAR annotations.


## Outline:
1. Download the csv containing information that associates sex chromosome information and subset it to relevant columns
2. For all genomes with both X and Y or Z and W assembled, estimate PAR boundaries using sequence alignment.
3. For these same genomes, estimate PAR boundaries using depth (halfdeep).
4. Compute metrics that compare across the three methods for inferring PAR boundaries.
5. Visualize the data.

## Step 1: Download the csv containing information that associates sex chromosome information and subset it to relevant columns
### Download the CSV:
December 15th I downloaded the main tab of the google sheet "VGP_list_sex_chroms_curated" as a csv, and then put it on biowulf. It can be found here on biowulf:
```/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/VGP_OrdinalList_Phase1Freeze_v1.1_sex_chroms_seqCov_HalfDeep_SCINKD_May8.25_DwnldDec16.25.csv```

It is stored in this GitHub repository in Datafiles/VGP_OrdinalList_Phase1Freeze_v1.1_sex_chroms_seqCov_HalfDeep_SCINKD_May8.25_DwnldDec16.25.csv

### Subset the datafile to only the columns relevant to this analysis:
I organized my thoughts regarding which data to keep from the "VGP_OrdinalList_Phase1Freeze_v1" csv file 

```
cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/

export VGP=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/VGP_OrdinalList_Phase1Freeze_v1.1_sex_chroms_seqCov_HalfDeep_SCINKD_May8.25_DwnldDec16.25.csv

module load R/4.5.0

Rscript 0a_Filter_Datafile.r
```
