# Chondricthyes
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Chondrichthyes
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Chondrichthyes/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Chondrichthyes/peptide

for SPECIES in \
  Raja_brachyura \
  Amblyraja_radiata \
  Narcine_bancroftii \
  Pristis_pectinata \
  Hypanus_sabinus \
  Mobula_birostris \
  Heptranchias_perlo \
  Pristiophorus_japonicus \
  Heterodontus_francisci \
  Stegostoma_tigrinum \
  Hemiscyllium_ocellatum \
  Mustelus_asterias \
  Scyliorhinus_canicula \
  Cetorhinus_maximus \
  Carcharodon_carcharias
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Chondrichthyes/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Chondrichthyes/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Chondrichthyes/

```
# Shrews
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Shrews
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Shrews/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Shrews/peptide

for SPECIES in \
  Suncus_etruscus \
  Sorex_araneus
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Shrews/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Shrews/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Shrews/
```
# Bat neo sex
```
CLADE="Bat_Neo"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Desmodus_rotundus \
  Phyllostomus_discolor \
  Glossophaga_mutica \
  Artibeus_intermedius \
  Artibeus_lituratus
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Snakes
```
CLADE="Snakes"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Vipera_latastei \
  Vipera_berus \
  Vipera_ursinii \
  Erythrolamprus_reginae \
  Thamnophis_elegans \
  Natrix_helvetica
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Squamates
```
CLADE="Squamates"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Rhineura_floridana \
  Lacerta_agilis \
  Zootoca_vivipara \
  Podarcis_siculus \
  Podarcis_muralis \
  Podarcis_liolepis \
  Podarcis_bocagei \
  Podarcis_vaucheri \
  Podarcis_pityusensis \
  Podarcis_tiliguerta \
  Podarcis_raffonei \
  Podarcis_filfolensis \
  Podarcis_melisellensis \
  Podarcis_gaigeae \
  Podarcis_cretensis \
  Podarcis_erhardii
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Strisores
```
CLADE="Strisores"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Caprimulgus_europaeus \
  Nyctibius_grandis \
  Podargus_strigoides \
  Aegotheles_albertisi \
  Hemiprocne_comata \
  Apus_apus \
  Phaethornis_superciliosus \
  Heliangelus_exortis \
  Calypte_anna
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Pigeons
```
CLADE="Pigeons"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Pterocles_gutturalis \
  Caloenas_nicobarica \
  Patagioenas_fasciata \
  Columba_livia \
  Nesoenas_mayeri \
  Streptopelia_turtur \
  Streptopelia_decaocto
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Spoonbills
```
CLADE="Spoonbills"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Pelecanus_crispus \
  Theristicus_caerulescens \
  Platalea_leucorodia \
  Morus_bassanus
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# manually run chondrichthyes and snakes
```
CLADE="Chondrichthyes"
./1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/

CLADE="Snakes"
./1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Squamates_broad
```
CLADE="Squamates"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Dibamus_smithi \
  Anolis_sagrei \
  Cyclura_pinguis \
  Vipera_latastei
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Hawks
```
CLADE="Hawks"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Athene_noctua \
  Strix_aluco \
  Sarcoramphus_papa \
  Gypaetus_barbatus \
  Aquila_chrysaetos \
  Morphnus_guianensis \
  Harpia_harpyja \
  Accipiter_gentilis \
  Buteo_buteo \
  Haliaeetus_albicilla
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Coraciimorphs
```
CLADE="Coraciimorphs"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Leptosomus_discolor \
  Pogoniulus_pusillus \
  Dryobates_pubescens \
  Colius_striatus
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Parrots
```
CLADE="Parrots"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Strigops_habroptilus \
  Amazona_ochrocephala \
  Guaruba_guaruba \
  Ara_ararauna \
  Psittacula_echo \
  Lathamus_discolor \
  Melopsittacus_undulatus \
  Acanthisitta_chloris
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Sylviida
```
CLADE="Sylviida"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Poecile_atricapillus \
  Hirundo_rustica \
  Zosterops_lateralis \
  Sylvia_atricapilla \
  Sylvia_borin
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Muscicapoids
```
CLADE="Muscicapoids"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Acridotheres_tristis \
  Cinclus_cinclus \
  Catharus_ustulatus \
  Erithacus_rubecula \
  Oenanthe_melanoleuca
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Monotremes
```
CLADE="Monotremes"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Tachyglossus_aculeatus \
  Ornithorhynchus_anatinus
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Marsupials
```
CLADE="Marsupials"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Monodelphis_domestica \
  Dromiciops_gliroides \
  Sminthopsis_crassicaudata \
  Sarcophilus_harrisii \
  Dasyurus_maculatus \
  Macrotis_lagotis \
  Notoryctes_typhlops \
  Phascolarctos_cinereus \
  Trichosurus_vulpecula \
  Macropus_eugenii
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Mammals_broad
```
CLADE="Mammals_broad"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Ornithorhynchus_anatinus \
  Monodelphis_domestica \
  Macrotis_lagotis \
  Homo_sapiens
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```
# Eutherian_fusion
```

CLADE="Eutherian_fusion"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Monodelphis_domestica \
  Dromiciops_gliroides \
  Sminthopsis_crassicaudata \
  Sarcophilus_harrisii \
  Trichosurus_vulpecula \
  Macropus_eugenii \
  Homo_sapiens
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/




```
CLADE="Ray_finned"
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide

for SPECIES in \
  Lepisosteus_oculatus \
  Coilia_mystus
do
  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/bed/${SPECIES}.bed \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/bed/

  cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/${SPECIES}/peptide/${SPECIES}.fa \
     /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/peptide/
done

sbatch 1g_genespace.sh /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/${CLADE}/
```