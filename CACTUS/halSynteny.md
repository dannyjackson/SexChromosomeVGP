halSynteny: a fast, easy-to-use conserved synteny block construction method for multiple whole-genome alignments
https://doi.org/10.1093/gigascience/giaa047
https://s3.amazonaws.com/genomeark/index.html?prefix=downstream_analyses/genome_alignments/cactus/577way/


#!/bin/bash
#SBATCH --job-name=halSynteny_test
#SBATCH --output=slurm_output/halSynteny_test.out
#SBATCH --error=slurm_output/halSynteny_test.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --mail-type=ALL

module load cactus/3.0.0

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/cactus

HAL_FILE=/data/Wilson_Lab/data/VGP_genomes_phase1/cactus_alignments/vgp-577way-v1.hal

halValidate ${HAL_FILE}

halStats ${HAL_FILE}

halSynteny ${HAL_FILE} test.psl --queryGenome GCF_016700215.2 --targetGenome GCA_009914755.4