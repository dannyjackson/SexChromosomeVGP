# VGP Sex Chromosome Synteny Analysis

awk -F',' 'NR>1{print $1}' /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv | sort -u | wc -l

sbatch /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts/filter_to_sex_chrs.all.sh
