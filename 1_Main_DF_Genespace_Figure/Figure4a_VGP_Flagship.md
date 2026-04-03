
Narcine_bancroftii
Gasterosteus_aculeatus
Pseudacris_triseriata
Podarcis_raffonei
Gallus_gallus
Homo_sapiens

# Prepare bed and peptide files
See Prepare_Pseudacris_tiseriata.md for details on preparing that genome.
```
WORKDIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure"
mkdir -p ${WORKDIR}/bed ${WORKDIR}/peptide

for SPECIES in \
    Gallus_gallus \
    Homo_sapiens \
    Narcine_bancroftii \
    Podarcis_raffonei \
    Gasterosteus_aculeatus
do
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed ${WORKDIR}/bed
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa ${WORKDIR}/peptide
done

# rename chr 12 in Narcine to X
awk 'BEGIN{OFS="\t"} {$1=($1==12 ? "X" : $1); print}' bed/Narcine_bancroftii.bed > bed/Narcine_bancroftii.tmp 

awk '{print $1}' bed/Narcine_bancroftii.tmp | sort -u

mv bed/Narcine_bancroftii.tmp bed/Narcine_bancroftii.bed
```
# Run genespace
## try using mamba env genespace instead of genespace_py3.10
sbatch 1g_genespace.sh ${WORKDIR}
# broke on Narcine bancroftii and Gasterosteus aculeatus
# Everything else works in a pairwise comparison (human chicken, Narcine v Podarcis)
