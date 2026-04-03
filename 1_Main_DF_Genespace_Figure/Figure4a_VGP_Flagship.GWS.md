Carcharodon_carcharias
Gasterosteus_aculeatus
Pseudacris_triseriata
Podarcis_raffonei
Gallus_gallus
Homo_sapiens

# Prepare bed and peptide files
See Prepare_Pseudacris_tiseriata.md for details on preparing that genome.
```
WORKDIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure_GWS"
mkdir -p ${WORKDIR}/bed ${WORKDIR}/peptide

for SPECIES in \
    Gallus_gallus \
    Homo_sapiens \
    Carcharodon_carcharias \
    Podarcis_raffonei \
    Gasterosteus_aculeatus
do
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_bed/${SPECIES}.bed ${WORKDIR}/bed
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/sexshared/usable_peptide/${SPECIES}.fa ${WORKDIR}/peptide
done

cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/bed/Pseudacris_triseriata.bed ${WORKDIR}/bed
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/peptide/Pseudacris_triseriata.fa ${WORKDIR}/peptide

```
# Run genespace
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/scripts
sbatch 1g_genespace.sh ${WORKDIR}

Rscript plot.VGP_MainDF_Figure.GWS.r