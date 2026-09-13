# Is the PAR misassembled in the house finch?
# download the bam
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus
```
# Download bam script
```
#!/bin/bash
#SBATCH --job-name=download_bam
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.err


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/bamfile

wget -c https://sra-pub-src-2.s3.amazonaws.com/SRR25647252/m54306U_210616_143938.subreads.bam.1 \
  -O SRR25647252.subreads.bam
```
# download fastas script
```
#!/bin/bash
#SBATCH --job-name=sra_fasta
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-35
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.err

module load sratoolkit/3.3.0

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus || exit 1

mkdir -p logs

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

mkdir -p "${run}_out"

fasterq-dump "$run" -O "${run}_out"

```
# submit
```
sbatch download_fasta_array.sh
```
# accessions.txt
```
SRR25647249
SRR25647250
SRR25647251
```
# Index Haemorhous_mexicanus genome
```
#!/bin/bash
#SBATCH --job-name=index_genome
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/index.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/index.err

module load bwa
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Haemorhous_mexicanus/ncbi_dataset/data/GCF_027477595.1/GCF_027477595.1_bHaeMex1.pri_genomic.fna"

bwa index $REF_GENOME
```
# Fast QC
```
#!/bin/bash
#SBATCH --job-name=raw_QC
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=4:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-35
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.err

module load fastqc

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/ || exit 1

mkdir -p trimmed_fastas raw_QC_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

while read -r run; do
    fastqc \
    ${run}_out/${run}_1.fastq ${run}_out/${run}_2.fastq \
    --outdir raw_QC_files/
    echo $run" raw QC done" >> trim_and_QC_log.txt
done < accessions.txt
```
# align
```
#!/bin/bash
#SBATCH --job-name=align
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-4
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.err

module load bwa
module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Haemorhous_mexicanus/ncbi_dataset/data/GCF_027477595.1/GCF_027477595.1_bHaeMex1.pri_genomic.fna"
RUN="${run}_out/${run}.fastq"

if [ ! -f "$RUN" ] ; then
  echo "Missing FASTQ files for $run"
  exit 1
fi

bwa mem -t "${SLURM_CPUS_PER_TASK}" "$REF_GENOME" "$RUN" | \
  samtools view -@ "${SLURM_CPUS_PER_TASK}" -o "bam_files/${run}.bam" -S
```
# sort
```
#!/bin/bash
#SBATCH --job-name=sortbams
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-95
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/%x_%A_%a.err

module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Haemorhous_mexicanus/ncbi_dataset/data/GCF_027477595.1/GCF_027477595.1_bHaeMex1.pri_genomic.fna"

mkdir -p sorted_bam_files

samtools sort "bam_files/${run}.bam" -o "sorted_bam_files/${run}.sorted.bam"

samtools index "sorted_bam_files/${run}.sorted.bam"
```
# Call variants
```
#!/bin/bash
#SBATCH --job-name=callvariants
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/callvariants.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/logs/callvariants.err

module load bcftools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus || exit 1

bamdir="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Haemorhous_mexicanus/sorted_bam_files/"

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Myuchelys_georgesi/ncbi_dataset/data/GCA_040894355.2/GCA_040894355.2_rMyuGeo1.pri_genomic.fna"

bcftools mpileup -Ou -f "$REF_GENOME" -a FORMAT/AD,DP,INFO/AD,SP "$bamdir"*.bam | bcftools call -mv -V indels > Myuchelys_georgesi_snps_multiallelic.vcf
```
# Plot region
samtools index bamfile/SRR25647252.subreads.bam

samplot plot \
  -n SRR25647252 \
  -b bamfile/SRR25647252.subreads.bam \
  -c NC_082345.1 \
  -s 76532527 \
  -e 76907864 \
  -o SRR25647252.PAR.png


  76654594	76657426
  76659645	76661564

samplot plot \
  -n SRR25647249 \
  -b sorted_bam_files/SRR25647249.sorted.bam \
  -c NC_082345.1 \
  -s 76533000 \
  -e 76537000 \
  -o SRR25647249.chr5_PAR.junction.png \
  -t DEL
