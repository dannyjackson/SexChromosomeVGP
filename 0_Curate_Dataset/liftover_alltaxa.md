# Add liftover species to genespace-to-chicken dataset
cd /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/scripts

# Bats
```
while IFS=, read -r species refgenome; do
  sbatch liftover.sh -s "$species" -r "$refgenome"
done <<'EOF'
Rousettus_aegyptiacus,Megaptera_novaeangliae
Doryrhina_cyclops,Hipposideros_larvatus
Rhinolophus_yonghoiseni,Rhinolophus_ferrumequinum
Artibeus_lituratus,Glossophaga_mutica
Molossus_molossus,Molossus_nigricans
Myotis_mystacinus,Myotis_daubentonii
EOF
```
# FILES_TO_DOWNLOAD="genome,gff3,rna,cds,protein,seq-report"

while IFS=, read -r SPECIES_NAME ACCESSION; do
  echo "=== ${SPECIES_NAME} -> ${ACCESSION}"
  echo "Downloading: ${FILES_TO_DOWNLOAD}"

  datasets download genome accession "${ACCESSION}" \
    --include "${FILES_TO_DOWNLOAD}" \

done <<'EOF'
Artibeus_lituratus,GCA_038363095.4
EOF

# use the gff_db for these
```
Myotis_emarginatus,Myotis_daubentonii

while IFS=, read -r species refgenome; do
  sbatch liftover.sh -s "$species" -r "$refgenome"
done <<'EOF'
Myotis_myotis,Myotis_daubentonii
Antrozous_pallidus,Myotis_daubentonii
Corynorhinus_townsendii,Myotis_daubentonii
Eptesicus_fuscus,Myotis_daubentonii
Vespertilio_murinus,Myotis_daubentonii
Nyctalus_leisleri,Myotis_daubentonii
Pipistrellus_nathusii,Myotis_daubentonii
Pipistrellus_kuhlii,Myotis_daubentonii
Pipistrellus_pygmaeus,Myotis_daubentonii
EOF
```
# Random
ls /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Pseudorca*/ncbi_dataset/data/GC*/*gff

```
while IFS=, read -r species refgenome; do
  sbatch liftover.sh -s "$species" -r "$refgenome"
done <<'EOF'
Larus_fuscus,Larus_michahellis
Anas_platyrhynchos,Anas_acuta
Willisornis_vidua,Dixiphia_pipra
Inia_geoffrensis,Pseudorca_crassidens
Monodon_monocero,Phocoena_phocoena
EOF
```


Podarcis_siculus  Cannot find the following peptide files in /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/sexshared/Podarcis_siculus/peptide:

cd /data/Wilson_Lab/data/VGP_genomes_phase1/lifted_gffs/scripts

# Everything else
ls /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/analyses/lacks_gff/ > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/reference_lists/lacks_gff.txt

while IFS=, read -r species refgenome; do
  sbatch liftover.sh -s "$species" -r "$refgenome"
done <<'EOF'
Equus_caballus,Diceros_bicornis_minor
Podarcis_siculus,Podarcis_raffonei
EOF

Elephas_maximus_indicus,Loxodonta_africana
Capra_hircus,Ovis_canadensis
Sus_scrofa,Pseudorca_crassidens
Notoryctes_typhlops,Macropus_eugenii
Mus_musculus,Apodemus_sylvaticus
Peromyscus_maniculatus,Onychomys_torridus
Lycocorax_pyrrhopterus,Corvus_hawaiiensis
Echiichthys_vipera,Gasterosteus_aculeatus
Tupaia_tana,Homo_sapiens
Chlamydotis_macqueenii,Columba_livia
EOF

# must redo:
while IFS=, read -r species refgenome; do
  sbatch liftover.sh -s "$species" -r "$refgenome"
done <<'EOF'
Notoryctes_typhlops,Macropus_eugenii
Tupaia_tana,Homo_sapiens
EOF


## Duplicates
Phascolarctos_cinereus,Macropus_eugenii
Hoplias_malabaricus,Gasterosteus_aculeatus
Echiichthys_vipera,Gasterosteus_aculeatus
Girardinichthys_multiradiatus,Gasterosteus_aculeatus
Lemur_catta,Homo_sapiens

while IFS=, read -r species refgenome; do
  sbatch liftover.sh -s "$species" -r "$refgenome"
done <<'EOF'
Phascolarctos_cinereus,Macropus_eugenii
Hoplias_malabaricus,Gasterosteus_aculeatus
Echiichthys_vipera,Gasterosteus_aculeatus
Girardinichthys_multiradiatus,Gasterosteus_aculeatus
Lemur_catta,Homo_sapiens
EOF

Amblyraja radiata
Tursiops truncatus
Panthera onca