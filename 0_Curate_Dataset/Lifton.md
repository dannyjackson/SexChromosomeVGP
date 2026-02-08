##  Subset split_gff_by_lineage.py

python3 split_gff_by_lineage.py \
  --ref "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD.csv" \
  --gff "/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/species_requiring_lifted_gff.txt" \
  --outdir ./gff_by_lineage


# environment
```
source myconda

conda env create -f /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/yamls/lifton.yml
conda activate lifton
```

```
Afrotheria Trichechus_inunguis # HiFi, Hifiasm Hi-C phasing; 695 contigs
Laurasiatheria Ovis_canadensis # HiFi, ONT, Hi-C, Verkko; 44 contigs
Marsupials Macropus eugenii # HiFi, Hifiasm solo; 829 contigs
Monotremes Ornithorhynchus anatinus # CLR II, FALCON unzip; 834 contigs
Supraprimates Microtus pennsylvanicus # HiFi, Hifiasm Hi-C phasing; 241 contigs
Xenarthra Dasypus novemcinctus # HiFi, Hifiasm Hi-C phasing; 1,310 contigs
```
list of just species to analyze
```

Afrotheria Trichechus_inunguis
Laurasiatheria Ovis_canadensis
Marsupials Macropus_eugenii
Monotremes Ornithorhynchus_anatinus
Supraprimates Microtus_pennsylvanicus
Xenarthra Dasypus_novemcinctus
```
# run lifton
```
#!/usr/bin/env bash
#SBATCH --job-name=vgp_dl
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/vgp_dl_%A_%a.out
#SBATCH --error=slurm_logs/vgp_dl_%A_%a.err

source myconda

conda activate lifton

# input tables
REFERENCE_FILE="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD.csv"

GFF_LIST="/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/gff_file_list.tsv"

FNA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Homo_sapiens/ncbi_dataset/data/GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna

GFF_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Homo_sapiens/ncbi_dataset/data/GCF_009914755.1/genomic.gff

FAA_REF=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Homo_sapiens/ncbi_dataset/data/GCF_009914755.1/protein.faa

GFF_QRY=/data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/Afrotheria/Trichechus_inunguis/lifted.Homo_sapiens.0_5.gff

FNA_QRY=/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Trichechus_inunguis/ncbi_dataset/data/GCA_046562895.1/GCA_046562895.1_mTriInu1_haplotype_2_genomic.fna

SC=0.5

lifton \
    -g ${GFF_REF} \
    -o ${GFF_QRY} \
    -P ${FAA_REF} \
    -t 2 \
    -sc ${SC} \
    ${FNA_QRY} \
    /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Homo_sapiens/ncbi_dataset/data/GCA_009914755.4/GCA_009914755.4_T2T-CHM13v2.0_genomic.fna

mv lifton_output/Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.gff .


grep -v '.1_' lifton_output/Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.gff > Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.no_dups.gff

wc -l lifton_output/Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.gff 
wc -l Pseudacris_triseriata.lifton.REF_Hyla_sarda.SC_0_5.no_dups.gff