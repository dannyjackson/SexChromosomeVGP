# Is the PAR misassembled in Willisornis vidua?

mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua
```
#!/bin/bash
#SBATCH --job-name=sra_fasta
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-1
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/%x_%A_%a.err

module load sratoolkit/3.3.0

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua || exit 1

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
SRR34768764
```
# Index Willisornis_vidua genome
```
#!/bin/bash
#SBATCH --job-name=index_genome
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --gres=lscratch:500
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/index.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/index.err

module load bwa
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Willisornis_vidua/ncbi_dataset/data/GCA_045364795.1/GCA_045364795.1_bWilVid1_haplotype_1_genomic.fna"

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
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/%x_%A_%a.err

module load fastqc

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/ || exit 1

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
#SBATCH --array=1
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_viduas/logs/%x_%A_%a.err

module load bwa
module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Willisornis_vidua/ncbi_dataset/data/GCA_045364795.1/GCA_045364795.1_bWilVid1_haplotype_1_genomic.fna"
RUN="${run}_out/${run}.fastq"

if [ ! -f "$RUN" ] ; then
  echo "Missing FASTQ files for $run"
  exit 1
fi

bwa mem -t "${SLURM_CPUS_PER_TASK}" "$REF_GENOME" "$RUN" | \
  samtools view -@ "${SLURM_CPUS_PER_TASK}" -o "bam_files/${run}.bam" -S
```
# Sort bams
```
#!/bin/bash
#SBATCH --job-name=sortbams
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:500
#SBATCH --array=1-95
#SBATCH --output=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/%x_%A_%a.out
#SBATCH --error=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua/logs/%x_%A_%a.err

module load samtools

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_issues/Willisornis_vidua || exit 1

mkdir -p logs bam_files

run=$(sed -n "${SLURM_ARRAY_TASK_ID}p" accessions.txt)

if [ -z "$run" ]; then
  echo "No accession found for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: ${run}"

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Willisornis_vidua/ncbi_dataset/data/GCA_045364795.1/GCA_045364795.1_bWilVid1_haplotype_1_genomic.fna"

mkdir -p sorted_bam_files

samtools sort "bam_files/${run}.bam" -o "sorted_bam_files/${run}.sorted.bam"

samtools index "sorted_bam_files/${run}.sorted.bam"
```
# Plot region
```
samplot plot \
  -n SRR34768764 \
  -b sorted_bam_files/SRR34768764.sorted.bam \
  -c CM098464.1 \
  -s 464609 \
  -e 970715 \
  -o SRR34768764.PAR.png


samplot plot \
  -n SRR25647252 \
  -b bamfile/SRR25647252.subreads.sorted.bam \
  -c CM098464.1 \
  -s 464609 \
  -e 970715 \
  -c chrB \
  -s 1 \
  -e 10000 \
  -t TRA \
  -o chrA_chrB_edge.png