# Can SCINKD detect instances of PAR dropout?

## Logic:
If the pseudoautosomal region of a sex chromosome system was only assembled on one of the two systems, we would expect to see high hapmer density on the PAR region compared to the rest of that chromosome.

To test this, I will:

0. Install SCINKD
1. Run SCINKD on the human X and Y chromosomes
2. Run SCINKD on the X and PAR masked Y

# 0. Install SCINKD
git clone https://github.com/DrPintoThe2nd/SCINKD.git

source myconda

mamba create -n scinkd meryl=1.4.1 snakemake=7.32.4 pigz r r-dplyr r-ggplot2 samtools --yes

mamba activate scinkd 

cat > make_scinkd_env.sbatch <<'EOF'
#!/bin/bash
#SBATCH -J make_scinkd_env
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 01:00:00
#SBATCH -o make_scinkd_env.%j.out
#SBATCH -e make_scinkd_env.%j.err

source ~/.bashrc
conda activate base

mamba create -n scinkd -c conda-forge -c bioconda \
  meryl=1.4.1 snakemake=7.32.4 pigz r r-dplyr r-ggplot2 samtools --yes
EOF

sbatch make_scinkd_env.sbatch
mamba activate scinkd

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/

cp -r ~/programs/SCINKD/ .

cd SCINKD


## 1. Run SCINKD on the human X and Y chromosomes

### Test run
wget --content-disposition https://ndownloader.figshare.com/files/49948980
wget --content-disposition https://ndownloader.figshare.com/files/49948983

ln -s Anniella_stebbinsi_HiFi_2024.asm.hic.hap1.p_ctg.FINAL.Genbank.fasta.gz Anniella_stebbinsi_HiFi_2024.asm.hic.hap1.fasta.gz
ln -s Anniella_stebbinsi_HiFi_2024.asm.hic.hap2.p_ctg.FINAL.Genbank.fasta.gz Anniella_stebbinsi_HiFi_2024.asm.hic.hap2.fasta.gz

cat <<'EOF' > SCINKD/config.json
{
	"prefix": "Anniella_stebbinsi_HiFi_2024.asm.hic"
}
EOF

cp ~/programs/SCINKD/config.json SCINKD/

snakemake \
  --use-conda \
  -np \
  -s SCINKD/SCINKD.v2.1.0.FULL.snakefile \
  --configfile SCINKD/config.json

snakemake --use-conda -c 24 -s SCINKD/SCINKD.v2.1.0.GREEDY.snakefile --configfile SCINKD/config.json       #run SCINKD in greedy mode for quick testing

### Actual run
```
cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/Homo_sapiens/raw

bgzip /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs/Homo_sapiens.X.fa

bgzip /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs/Homo_sapiens.Y.fa


ln -s \
/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs/Homo_sapiens.X.fa.gz \
Homo_sapiens.hap1.fasta.gz
ln -s \
/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs/Homo_sapiens.Y.fa.gz \
Homo_sapiens.hap2.fasta.gz

samtools faidx Homo_sapiens.hap1.fasta.gz
samtools faidx Homo_sapiens.hap2.fasta.gz
```

Then, edit the SCINKD/config.json file to the name of the species of interest:
```
printf '%s\n' '{' '  "prefix": "Homo_sapiens"' '}' > /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/config.json
```
Run SCINKD
```
#dry-run to test inputs
snakemake --use-conda   -np \
    -s /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/SCINKD.v2.1.0.FULL.snakefile \
    --configfile /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/config.json
```

Run the following as a slurm script:

```
#!/bin/bash
#SBATCH -J runscinkd_raw
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 2:00:00
#SBATCH -o runscinkd_raw.%j.out
#SBATCH -e runscinkd_raw.%j.err

set -euo pipefail

source myconda

mamba activate scinkd

#run SCINKD in greedy mode for quick testing

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/Homo_sapiens/raw

snakemake --use-conda -c 24 \
    -s /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/SCINKD.v2.1.0.GREEDY.snakefile \
    --configfile /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/config.json --printshellcmds --show-failed-logs
```
```
sbatch scinkd_submit.sh
```
Plot it:
```
Rscript plot_kmer_density_raw.r
```

## 2. Run SCINKD on the X and PAR masked Y
```
cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/Homo_sapiens/masked_PAR

wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/masked_DJ_rDNA_PHR_5S_wi_rCRS/chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa

bgzip chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa

samtools faidx chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa.gz

samtools faidx chm13v2.0_masked_DJ_5S_rDNA_PHR_PAR_wi_rCRS.fa.gz chrY > Homo_sapiens.masked.Y.fa

bgzip Homo_sapiens.masked.Y.fa

ln -s \
/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/sex_chrs/Homo_sapiens.X.fa.gz \
Homo_sapiens.masked.hap1.fasta.gz
ln -s \
/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/Homo_sapiens/masked_PAR/Homo_sapiens.masked.Y.fa.gz \
Homo_sapiens.masked.hap2.fasta.gz

samtools faidx Homo_sapiens.masked.hap1.fasta.gz
samtools faidx Homo_sapiens.masked.hap2.fasta.gz
```

Revise config file:
```
printf '%s\n' '{' '  "prefix": "Homo_sapiens.masked"' '}' > /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/config.masked.json
```

Run the following using slurm:
```
#!/bin/bash
#SBATCH -J runscinkd_masked
#SBATCH -c 2
#SBATCH --mem=16G
#SBATCH -t 1:00:00
#SBATCH -o runscinkd_masked.%j.out
#SBATCH -e runscinkd_masked.%j.err

set -euo pipefail

source myconda

mamba activate scinkd

#run SCINKD in greedy mode for quick testing

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/Homo_sapiens/masked_PAR

snakemake --use-conda -c 24 \
    -s /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/SCINKD.v2.1.0.GREEDY.snakefile \
    --configfile /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/SCINKD_PAR_dropout/SCINKD/config.masked.json --printshellcmds --show-failed-logs
```

```
sbatch scinkd_submit.masked.sh
```
Plot it:
```
Rscript plot_kmer_density_masked.r
```
