# All Independent SC
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Independent_SC
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Independent_SC/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Independent_SC/peptide

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Independent_SC

for SPECIES in \
    Heptranchias_perlo \
    Coilia_mystus \
    Hoplias_malabaricus \
    Argentina_silus \
    Aulostomus_maculatus \
    Girardinichthys_multiradiatus \
    Echiichthys_vipera \
    Gasterosteus_aculeatus \
    Hyla_sarda \
    Ornithorhynchus_anatinus \
    Homo_sapiens \
    Dibamus_smithi \
    Tiliqua_scincoides \
    Shinisaurus_crocodilurus \
    Anniella_stebbinsi \
    Furcifer_pardalis \
    Anolis_sagrei \
    Cyclura_pinguis \
    Erythrolamprus_reginae \
    Podarcis_bocagei \
    Gallus_gallus
do

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Independent_SC/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Independent_SC/peptide/
done

cd ../scripts

sbatch IndepSC.genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Independent_SC/

IndepSC.genespace.R
```