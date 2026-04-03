Get CACTUS alignment files
cd 
# vgp-577way-v1.hal

#!/bin/bash
#SBATCH --job-name=CACTUS_hal
#SBATCH --output=slurm_output/CACTUS_hal.out
#SBATCH --error=slurm_output/CACTUS_hal.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --mail-type=ALL

cd /data/Wilson_Lab/data/VGP_genomes_phase1/cactus_alignments
URL="https://s3.amazonaws.com/genomeark/downstream_analyses/genome_alignments/cactus/577way/vgp-577way-v1.hal"

aria2c -x 16 "$URL"

#!/bin/bash
#SBATCH --job-name=CACTUS_maf
#SBATCH --output=slurm_output/CACTUS_maf.out
#SBATCH --error=slurm_output/CACTUS_maf.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --mail-type=ALL

cd /data/Wilson_Lab/data/VGP_genomes_phase1/cactus_alignments
URL="https://s3.amazonaws.com/genomeark/downstream_analyses/genome_alignments/cactus/577way/vgp-577way-v1-hs1.maf.gz"

aria2c -x 16 "$URL"

wget https://s3.amazonaws.com/genomeark/downstream_analyses/genome_alignments/cactus/577way/vgp-577way-v1-hs1.maf.gz.tai