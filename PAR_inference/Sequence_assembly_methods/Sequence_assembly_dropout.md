# Influence of sequencing / assembly method on PAR dropout
```
cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/Sequence_Assembly_PAR_dropout

source myconda

conda env create -f /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/yamls/python_env.yml

conda activate python_env

```

## 1. Create a datatable of all sequencing and assembly methods for each species
```
chmod +x summarize_methods.py
./summarize_methods.py \
  --root /data/Wilson_Lab/data/VGP_genomes_phase1/genomes \
  --out vgp_phase1_methods.tsv

python3 classify_seq_method.py vgp_phase1_methods.tsv vgp_phase1_methods.classified.tsv

./summarize_methods.py \
  --root /data/Wilson_Lab/data/VGP_genomes_phase1/Alternate_Haplotypes \
  --out vgp_phase1_methods.althaps.tsv

python3 classify_seq_method.py vgp_phase1_methods.althaps.tsv vgp_phase1_methods.althaps.classified.tsv

chmod +x merge_classifications.py
./merge_classifications.py vgp_phase1_methods.althaps.classified.tsv \
  vgp_phase1_methods.classified.tsv \
  -o vgp_phase1_methods.althaps.merged.tsv
```

Okay finally giving this up for the weekend. Note for myself on post-shutdown-day: 345 of the genomes are unclassified in sequencing/assembly method. 112 of those are Darwin Tree of Life Project genomes. Pulled up the methods of one of them which has a lot more information than they are putting in the NCBI Assembly metadata field: https://wellcomeopenresearch.org/articles/9-399/v1


datasets summary genome accession --inputfile VGP_Phase1_Accessions.mainhap.txt --as-json-lines | dataformat tsv genome --fields organism-name,assminfo-name,accession,annotinfo-name,annotinfo-release-date,annotinfo-pipeline,assminfo-paired-assm-accession,assminfo-paired-assm-name | awk 'BEGIN{FS="\t";OFS="\t"}{if($6=="NCBI EGAPx"){print $1,$2,$3,$5}else{if($7~/GCF/){print $1,$2,$7,$8}}}' > NCBI_list_of_annotated_genomes.tsv
