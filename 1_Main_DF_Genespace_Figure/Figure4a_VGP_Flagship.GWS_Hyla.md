# Prepare bed and peptide files
```
WORKDIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Synteny_Figure"
mkdir -p ${WORKDIR}/bed ${WORKDIR}/peptide

for SPECIES in \
    Gallus_gallus \
    Homo_sapiens \
    Carcharodon_carcharias \
    Podarcis_raffonei \
    Gasterosteus_aculeatus \
    Hyla_sarda
do
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed ${WORKDIR}/bed
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa ${WORKDIR}/peptide
done

cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Hyla_sarda/bed/Hyla_sarda.bed ${WORKDIR}/bed
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Hyla_sarda/peptide/Hyla_sarda.fa ${WORKDIR}/peptide

# replace chr 7 with X in Hyla sarda

awk 'BEGIN{OFS="\t"} {$1=($1=="7"?"X":$1); print}' \
bed/Hyla_sarda.bed \
> tmp && mv tmp bed/Hyla_sarda.bed


```
# Run genespace
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch genespace.sh ${WORKDIR}

Rscript plot.VGP_MainDF_Figure.r