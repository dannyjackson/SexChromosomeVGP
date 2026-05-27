# PAR annotation

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_annotation

# Troubleshoot with X of marmoset

# Filter protein.faa to just proteins mapped to the X in the gff
## Filter GFF to X chromosome
```
module load agat/1.4.0

base="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Callithrix_jacchus/ncbi_dataset/data/GCF_049354715.1"
chr="NC_133524.1"

gff="${base}/genomic.gff"
faa="${base}/protein.faa"

mkdir -p XZ_chrs

agat_sp_filter_record_by_coordinates.pl \
  --gff "$gff" \
  --tsv <(echo -e "${chr}\t1\t999999999") \
  -o XZ_chrs/Callithrix_jacchus.Xchr.gff
```
## Extract protein IDs from X chr gff
```
agat_sp_extract_attributes.pl \
  --gff XZ_chrs/Callithrix_jacchus.Xchr.gff/NC_133524.1_1_999999999.gff3 \
  --att protein_id \
  -t cds \
  -o XZ_chrs/Callithrix_jacchus.Xchr.txt
```
## Filter protein fasta
```
seqkit grep -f XZ_chrs/Callithrix_jacchus.Xchr_protein_id.txt "$faa" \
  > XZ_chrs/Callithrix_jacchus.Xchr.protein.faa
```
# Run emapper

base="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Callithrix_jacchus/ncbi_dataset/data/GCF_049354715.1"

emapper.py -m diamond --itype proteins \
  -i XZ_chrs/Callithrix_jacchus.Xchr.protein.faa \
  -o Callithrix_jacchus.eggnog \
  --output_dir eggnog_annotations

# rename all genes
base="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Callithrix_jacchus/ncbi_dataset/data/GCF_049354715.1"

emapper.py -m diamond --itype proteins \
  -i "${base}/protein.faa" \
  -o Callithrix_jacchus.eggnog \
  --output_dir eggnog_annotations



######################################################################
# prior attempts
######################################################################


/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_annotation
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv

mkdir -p XZ_chrs
mkdir -p YW_chrs

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

species=Callithrix_jacchus

echo 'NC_133524.1' > ${species}.xchr.txt
samtools faidx ${base}/${species}/${species}.fna -r ${species}.xchr.txt > XZ_chrs/${species}.Xchr.fasta
echo 'NC_133525.1' > ${species}.ychr.txt
samtools faidx ${base}/${species}/${species}.fna -r ${species}.ychr.txt > YW_chrs/${species}.Ychr.fasta

module load eggnog-mapper/2.1.12

mkdir -p eggnog_annotations

emapper.py -m diamond --itype genome --genepred prodigal \
-i XZ_chrs/${species}.Xchr.fasta -o test --output_dir eggnog_annotations

