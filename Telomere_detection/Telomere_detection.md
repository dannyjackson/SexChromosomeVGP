# Use the Telomere Identification toolKit (tidk) to annotate telomeres on the VGP genomes
# https://github.com/tolkit/telomeric-identifier

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/Telomere_detection

module load mamba_install
source myconda

conda install -c bioconda tidk

FASTA=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13/Xchr.fa

tidk search $FASTA --string TTAGGG --output CHM13_Xchr --dir Homo_sapiens


FASTA=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13/Xchr.fa

tidk search $FASTA --string TTAGGG --output CHM13_Xchr --dir Homo_sapiens

# Make a reference list of all VGP fasta files
cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/reference_lists/
wget https://raw.githubusercontent.com/VGP/vgp-phase1/main/VGPPhase1-freeze-1.0.tsv

First, make sure that my search criteria
ls /data/Wilson_Lab/data/VGP_genomes/*/*.fna | grep -Ev '\.(cds|rna)\.fna$' | wc -l

base=/data/Wilson_Lab/data/VGP_genomes

for d in "$base"/*/; do
   name=$(basename "$d")
   [[ -f "$d/$name.fna" ]] || echo "$name"
 done
ls $base/Myotis_mystacinus

Missing Myotis_mystacinus
https://github.com/VGP/vgp-phase1/blob/main/VGPPhase1-freeze-1.0.tsv

ls /data/Wilson_Lab/data/VGP_genomes/*/*.fna | grep -Ev '\.(cds|rna)\.fna$' > /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/reference_lists/VGP_phase1_fastas.txt

sed 's#.*/##; s/\.fna$//' VGP_phase1_fastas.txt | sort -u > fasta_species.list
awk -F'\t' 'NR>1 {
  split($10,a," ");
  print a[1] "_" a[2]
}' VGPPhase1-freeze-1.0.tsv | sort -u > vgp_species.list

comm -23 vgp_species.list fasta_species.list | wc -l
