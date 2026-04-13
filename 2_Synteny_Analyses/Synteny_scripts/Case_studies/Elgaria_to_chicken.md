# Elgaria to chicken
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Elgaria
mkdir raw_genome
cd raw_genome

source myconda
mamba activate ncbi_datasets

FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"


datasets download genome accession GCF_023053635.1 \
    --include "${FILES_TO_DOWNLOAD}" 

unzip ncbi_dataset.zip

cd ncbi_dataset/data/GCF_023053635.1/

mamba activate genespace

transeq -sequence cds_from_genomic.fna -outseq cds_from_genomic.translated.cds

```
## Prepare all bed and peptide files
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Elgaria

mkdir -p Elgaria_multicarinata_webbii
cd Elgaria_multicarinata_webbii

SPECIES=Elgaria_multicarinata_webbii
ln -sf /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Elgaria/raw_genome/ncbi_dataset/data/GCF_023053635.1/cds_from_genomic.translated.cds ${SPECIES}.translated.cds
ln -sf /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Elgaria/raw_genome/ncbi_dataset/data/GCF_023053635.1/genomic.gff ${SPECIES}.gff
ln -sf /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Elgaria/raw_genome/ncbi_dataset/data/GCF_023053635.1/GCF_023053635.1_rElgMul1.1.pri_genomic.fna ${SPECIES}.fa

SCRIPTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Elgaria/

Rscript ${SCRIPTDIR}/0b_prepare_genespace.R "$OUTDIR" "$OUTDIR"

```
# Create requisite genespace directory structure
```
GENESPACEDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/case_studies/Elgaria/Genespace/
mkdir -p ${GENESPACEDIR}/${SPECIES}
mv bed ${GENESPACEDIR}/${SPECIES}
mv peptide ${GENESPACEDIR}/${SPECIES}

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/bed/Gallus_gallus.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/peptide/Gallus_gallus.fa"

cd ${GENESPACEDIR}

# Symlink references (handy to keep alongside each species)
ln -sf "$REF_BED" "${SPECIES}/bed/Gallus_gallus.bed"
ln -sf "$REF_PEP" "${SPECIES}/peptide/Gallus_gallus.fa"
```
# Submit genespace
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts

sbatch 1b_genespace.singleinstance.sh Gallus_gallus ${SPECIES} ${GENESPACEDIR}
