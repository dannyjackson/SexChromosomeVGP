# Simone's code for PAR inference via alignment

srun ~/minimap2-2.28_x64-linux/minimap2 -x asm5 -c --eqx -k 23 -t 3 ${species_name}.${small_chr}chrom.fasta ${species_name}.${big_chr}chrom.fasta \
 > /scratch/sgable3/WGA_Benchmarking/minimap2/PAR_analysis/${species_name}_${small_chr}to${big_chr}.aln.paf

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/

base="$in"
infile="${base}.paf"

outbed="${base}.id98_5.len10k.refqry.bed"

# Filter and write BED:
# Output BED columns:
# ref_id  ref_start  ref_end  query_id  query_start  query_end  identity_percent  ref_span_len
awk -v thr=0.985 -v minlen=10000 -v OFS='\t' '
  function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
  function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
  {
    # PAF fields:
    # 1 qname 2 qlen 3 qstart 4 qend 5 strand 6 tname 7 tlen 8 tstart 9 tend 10 nmatch 11 alnlen 12 mapq
    tspan = $9 - $8
    de=""; for(i=13;i<=NF;i++) if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }
    id = ident_from_de(de); if (id < 0) id = ident_from_core($10,$11)

    if (id >= thr && tspan >= minlen)
      printf "%s\t%d\t%d\t%s\t%d\t%d\t%.4f\t%d\n", $6, $8, $9, $1, $3, $4, id*100, tspan
  }
' "${base}.paf" > "PAR_annotations/$outbed"

# output continuous percent identity
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2
mkdir -p continuous_percentID

for infile in *.paf; do
  base=$(basename "$infile" .paf)
  outbed="continuous_percentID/${base}.refqry.csv"

  # Write header
  printf "chrom_qry\tlen_qry\tbp_start_qry\tbp_end_qry\tpercent_identity_qry\tchrom_ref\tlen_ref\tbp_start_ref\tbp_end_ref\tpercent_identity_ref\n" \
    > "$outbed"

  # Process file
  awk -v OFS='\t' '
    function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
    function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
    {
      # Extract de:f tag if present
      de=""
      for(i=13;i<=NF;i++)
        if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }

      # Compute identity
      id = ident_from_de(de)
      if (id < 0) id = ident_from_core($10,$11)

      pid = (id >= 0 ? id * 100 : "NA")

      printf "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\n", \
             $1, $2, $3, $4, pid, \
             $6, $7, $8, $9, pid
    }
  ' "$infile" >> "$outbed"

  echo "Wrote $outbed"
done



# plot continuous percent ID

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

# Identify all species missing from Simone's PAF files that are heterogametous

## 1. ID all heterogametous species
## 2. Compare counts of heterogametous genomes and unique PAF files (216)
## 3. Intersect the two and ID unique in [A] and unique in [B]

List of heterogametic genomes
```
Tachyglossus_aculeatus
Ornithorhynchus_anatinus
Ambystoma_mexicanum
Anniella_stebbinsi
Mustelus_asterias
Lycaon_pictus
Dama_dama
Argentina_silus
Girardinichthys_multiradiatus
Gasterosteus_aculeatus
Echiichthys_vipera
Elephas_maximus_indicus
Loxodonta_africana
Heterohyrax_brucei
Trichechus_inunguis
Rhynchocyon_petersi
Hemiscyllium_ocellatum
Carcharodon_carcharias
Pristiophorus_japonicus
Narcine_bancroftii
Tursiops_truncatus
Pseudorca_crassidens
Lagenorhynchus_acutus
Globicephala_melas
Grampus_griseus
Inia_geoffrensis
Delphinus_delphis
Balaenoptera_musculus
Stenella_coeruleoalba
Eschrichtius_robustus
Megaptera_novaeangliae
Balaenoptera_physalus
Mesoplodon_densirostris
Mesoplodon_bidens
Eubalaena_glacialis
Monodon_monoceros
Ovis_canadensis
Ovis_aries
Capra_hircus
Sus_scrofa
Camelus_dromedarius
Equus_caballus
Phyllostomus_discolor
Desmodus_rotundus
Molossus_molossus
Molossus_nigricans
Tadarida_brasiliensis
Myotis_daubentonii
Myotis_mystacinus
Myotis_nattereri
Pipistrellus_pygmaeus
Pipistrellus_hanaki
Eptesicus_nilssonii
Rhynchonycteris_naso
Corynorhinus_townsendii
Vespertilio_murinus
Nyctalus_leisleri
Miniopterus_schreibersii
Rhinolophus_trifoliatus
Rhinolophus_perniger_lanosus
Rousettus_aegyptiacus
Aselliscus_stoliczkanus
Lutra_lutra
Mustela_erminea
Mustela_lutreola
Mustela_nivalis_vulgaris
Meles_meles
Canis_lupus_baileyi
Zalophus_californianus
Lynx_canadensis
Neofelis_nebulosa
Panthera_onca
Manis_pentadactyla
Suncus_etruscus
Martes_martes
Trichosurus_vulpecula
Macropus_eugenii
Monodelphis_domestica
Anolis_sagrei
Cyclura_pinguis
Dibamus_smithi
Homo_sapiens
Pan_troglodytes
Pan_paniscus
Gorilla_gorilla_gorilla
Pongo_abelii
Pongo_pygmaeus
Symphalangus_syndactylus
Macaca_nemestrina
Callithrix_jacchus
Saimiri_boliviensis
Nycticebus_coucang
Lemur_catta
Cynocephalus_volans
Erethizon_dorsatum
Muscardinus_avellanarius
Sciurus_vulgaris
Sciurus_carolinensis
Urocitellus_parryii
Callospermophilus_lateralis
Marmota_flaviventris
Apodemus_sylvaticus
Arvicanthis_niloticus
Rattus_norvegicus
Microtus_pennsylvanicus
Chionomys_nivalis
Peromyscus_maniculatus_sonoriensis
Jaculus_jaculus
Thomomys_bottae
Ochotona_princeps
Tupaia_tana
Macaca_fascicularis
Tamandua_tetradactyla
Choloepus_didactylus
Dasypus_novemcinctus
Artibeus_lituratus
Artibeus_intermedius
Hoplias_malabaricus
Hypanus_sabinus
Calypte_anna
Heliangelus_exortis
Apus_apus
Hemiprocne_comata
Caprimulgus_europaeus
Aegotheles_albertisi
Podargus_strigoides
Nyctibius_grandis
Streptopelia_turtur
Nesoenas_mayeri
Patagioenas_fasciata
Columba_livia
Phoenicopterus_ruber
Taeniopygia_guttata
Vidua_chalybeata
Haemorhous_mexicanus
Melospiza_georgiana
Melospiza_melodia_melodia
Ammospiza_nelsoni
Passerculus_sandwichensis
Zonotrichia_albicollis
Passer_domesticus
Sylvia_atricapilla
Hirundo_rustica
Zosterops_lateralis
Geothlypis_trichas
Erithacus_rubecula
Catharus_ustulatus
Poecile_atricapillus
Agelaius_phoeniceus
Lycocorax_pyrrhopterus_obiensis
Corvus_hawaiiensis
Corvus_moneduloides
Cyanocitta_cristata
Coloeus_monedula
Chiroxiphia_lanceolata
Dixiphia_pipra
Willisornis_vidua
Acanthisitta_chloris
Strigops_habroptilus
Melopsittacus_undulatus
Amazona_ochrocephala
Lathamus_discolor
Guaruba_guaruba
Falco_peregrinus
Falco_rusticolus
Falco_biarmicus
Falco_cherrug
Falco_naumanni
Cariama_cristata
Dryobates_pubescens
Pogoniulus_pusillus
Merops_nubicus
Bucorvus_abyssinicus
Trogon_surrucura
Colius_striatus
Strix_aluco
Aquila_chrysaetos_chrysaetos
Haliaeetus_albicilla
Harpia_harpyja
Morphnus_guianensis
Accipiter_gentilis
Gypaetus_barbatus
Sarcoramphus_papa
Pelecanus_crispus
Platalea_leucorodia
Gulosus_aristotelis
Phalacrocorax_carbo
Ciconia_maguari
Calonectris_borealis
Spheniscus_humboldti
Phaethon_aethereus
Rhynochetos_jubatus
Pluvialis_apricaria
Numenius_arquata
Sterna_hirundo
Rissa_tridactyla
Larus_michahellis
Larus_argentatus
Balearica_regulorum_gibbericeps
Porphyrio_hochstetteri
Opisthocomus_hoazin
Lagopus_muta
Aythya_fuligula
Aythya_ferina
Aythya_marila
Netta_rufina
Anas_platyrhynchos
Cygnus_olor
Anas_acuta
Anser_cygnoides
Anser_brachyrhynchus
Anser_anser
Cygnus_columbianus
Mergus_octosetaceus
Chlamydotis_macqueenii
Tauraco_erythrolophus
Cuculus_canorus
Struthio_camelus_australis
Eudromia_elegans
Dromaius_novaehollandiae
Vipera_ursinii
Vipera_berus
Thamnophis_elegans
Natrix_helvetica
Lacerta_agilis
Podarcis_raffonei
Podarcis_siculus
Podarcis_muralis
Podarcis_bocagei
Podarcis_gaigeae
Podarcis_pityusensis
Podarcis_erhardii
Podarcis_melisellensis
Podarcis_filfolensis
Podarcis_liolepis
Podarcis_vaucheri
Podarcis_tiliguerta
Gallus_gallus
Zootoca_vivipara
```
List of putatively heterogametic genomes (require curation)
```
Rhinolophus_yonghoiseni
Eptesicus_fuscus
Antrozous_pallidus
Rhinopoma_microphyllum
Doryrhina_cyclops
Dasyurus_maculatus
Saccopteryx_leptura
Saccopteryx_bilineata
Sarcophilus_harrisii
Rhinolophus_affinis
Mops_condylurus
Megaderma_spasma
Stegostoma_tigrinum
Scyliorhinus_canicula
Sternotherus_odoratus
Shinisaurus_crocodilurus
Acomys_minous
```
# Species in paf list
ls minimap2/*paf \
  | sed 's|.*/||; s/_WtoZ\.aln\.paf$//' \
  | sed 's|.*/||; s/_YtoX\.aln\.paf$//' \
  | sed 's|.*/||; s/_Y1toX\.aln\.paf$//' \
  | sed 's|.*/||; s/_Y2toX\.aln\.paf$//' \
  | sed 's|.*/||; s/_YtoX1\.aln\.paf$//' \
  | sed 's|.*/||; s/_YtoX2\.aln\.paf$//' \
  | sort -u > paf_species.txt

# Species in heterogametic list
sort -u heterogametic.txt > heterogametic_species.txt

# Heterogametic species missing from paf list
comm -23 heterogametic_species.txt paf_species.txt > heterogametic_missing_in_paf.txt

# PAF species missing from heterogametic list
comm -13 heterogametic_species.txt paf_species.txt > paf_missing_in_heterogametic.txt

wc -l heterogametic_missing_in_paf.txt
wc -l paf_missing_in_heterogametic.txt
```
# Remove species that were actually run but this is missing them due to nomenclature issues
```
Amblystoma_mexicanum
Gorilla_gorilla_gorilla
Balearica_regulorum_gibbericeps
Melospiza_melodia_melodia
Rhinolophus_perniger_lanosus
```
# Rename species in file that include subspecies in name
```
# Aquila_chrysaetos_chrysaetos
sed -i 's/chrysaetos_chrysaetos/chrysaetos/g' heterogametic_missing_in_paf.txt
# Lycocorax_pyrrhopterus_obiensis
sed -i 's/_obiensis//g' heterogametic_missing_in_paf.txt
# Peromyscus_maniculatus_sonoriensis
sed -i 's/_sonoriensis//g' heterogametic_missing_in_paf.txt
```
# Make paf for species that were missed

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"


## Gorilla gorilla
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles

species=Gorilla_gorilla

echo 'NC_073247.2' > Gorilla.xchr.txt
samtools faidx ${base}/${species}/${species}.fna -r Gorilla.xchr.txt > XZ_chrs/${species}.Xchr.fasta
echo 'NC_073248.2' > Gorilla.ychr.txt
samtools faidx ${base}/${species}/${species}.fna -r Gorilla.ychr.txt > YW_chrs/${species}.Ychr.fasta

srun ~/minimap2-2.28_x64-linux/minimap2 

minimap2 -x asm5 -c -eqx -k 23 -t 3 YW_chrs/${species}.Ychr.fasta XZ_chrs/${species}.Xchr.fasta \
 > minimap2/${species}_YtoX.aln.paf


species=Gorilla_gorilla

sbatch --cpus-per-task=3 --mem=8G --time=12:00:00 --job-name="${species}_YtoX" \
  --wrap="minimap2 -x asm5 -c -eqx -k 23 -t 3 \
  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/YW_chrs/${species}.Ychr.fasta \
  /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/${XZ_chrs}/${species}.Xchr.fasta \
  > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/${species}_YtoX.aln.paf"

# Output continuous percent id file 
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2

infile=Gorilla_gorilla_YtoX.aln.paf
infile=Pongo_abelii_YtoX.aln.paf

base=$(basename "$infile" .paf)
outbed="continuous_percentID/${base}.refqry.csv"

# Write header
printf "chrom_qry\tlen_qry\tbp_start_qry\tbp_end_qry\tpercent_identity_qry\tchrom_ref\tlen_ref\tbp_start_ref\tbp_end_ref\tpercent_identity_ref\n" \
  > "$outbed"

# Process file
awk -v OFS='\t' '
  function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
  function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
  {
    # Extract de:f tag if present
    de=""
    for(i=13;i<=NF;i++)
      if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }

    # Compute identity
    id = ident_from_de(de)
    if (id < 0) id = ident_from_core($10,$11)

    pid = (id >= 0 ? id * 100 : "NA")

    printf "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\n", \
            $1, $2, $3, $4, pid, \
            $6, $7, $8, $9, pid
  }
' "$infile" >> "$outbed"

echo "Wrote $outbed"
```
# Zosterops_lateralis


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/"
species=Zosterops_lateralis

echo 'OZ246513.1' > Zosterops_lateralis.zchr.txt
samtools faidx ${base}/${species}/${species}.fna -r ${species}.zchr.txt > XZ_chrs/${species}.Zchr.fasta
echo 'OZ246514.1' > Zosterops_lateralis.wchr.txt
samtools faidx ${base}/${species}/${species}.fna -r ${species}.wchr.txt > YW_chrs/${species}.Wchr.fasta


sbatch --cpus-per-task=3 --mem=8G --time=12:00:00 --job-name="${species}_WtoZ" \
  --wrap="module load minimap2; minimap2 -x asm5 -c -eqx -k 23 -t 3 /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/YW_chrs/${species}.Wchr.fasta /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/XZ_chrs/${species}.Zchr.fasta > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/${species}_WtoZ.aln.paf"


# Output continuous percent id file 
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2

infile=Zosterops_lateralis_WtoZ.aln.paf


base=$(basename "$infile" .paf)
outbed="continuous_percentID/${base}.refqry.csv"

# Write header
printf "chrom_qry\tlen_qry\tbp_start_qry\tbp_end_qry\tpercent_identity_qry\tchrom_ref\tlen_ref\tbp_start_ref\tbp_end_ref\tpercent_identity_ref\n" \
  > "$outbed"

# Process file
awk -v OFS='\t' '
  function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
  function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
  {
    # Extract de:f tag if present
    de=""
    for(i=13;i<=NF;i++)
      if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }

    # Compute identity
    id = ident_from_de(de)
    if (id < 0) id = ident_from_core($10,$11)

    pid = (id >= 0 ? id * 100 : "NA")

    printf "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\n", \
            $1, $2, $3, $4, pid, \
            $6, $7, $8, $9, pid
  }
' "$infile" >> "$outbed"

echo "Wrote $outbed"


# Rhynochetos_jubatus


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/"
species=Rhynochetos_jubatus

echo 'CM050690.1' > Rhynochetos_jubatus.zchr.txt
samtools faidx ${base}/${species}/${species}.fna -r ${species}.zchr.txt > XZ_chrs/${species}.Zchr.fasta
echo 'CM050689.1' > Rhynochetos_jubatus.wchr.txt
samtools faidx ${base}/${species}/${species}.fna -r ${species}.wchr.txt > YW_chrs/${species}.Wchr.fasta


sbatch --cpus-per-task=3 --mem=8G --time=12:00:00 --job-name="${species}_WtoZ" \
  --wrap="module load minimap2; minimap2 -x asm5 -c -eqx -k 23 -t 3 /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/YW_chrs/${species}.Wchr.fasta /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/XZ_chrs/${species}.Zchr.fasta > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/${species}_WtoZ.aln.paf"


# Output continuous percent id file 
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2

infile=Rhynochetos_jubatus_WtoZ.aln.paf

base=$(basename "$infile" .paf)
outbed="continuous_percentID/${base}.refqry.csv"

# Write header
printf "chrom_qry\tlen_qry\tbp_start_qry\tbp_end_qry\tpercent_identity_qry\tchrom_ref\tlen_ref\tbp_start_ref\tbp_end_ref\tpercent_identity_ref\n" \
  > "$outbed"

# Process file
awk -v OFS='\t' '
  function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
  function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
  {
    # Extract de:f tag if present
    de=""
    for(i=13;i<=NF;i++)
      if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }

    # Compute identity
    id = ident_from_de(de)
    if (id < 0) id = ident_from_core($10,$11)

    pid = (id >= 0 ? id * 100 : "NA")

    printf "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\n", \
            $1, $2, $3, $4, pid, \
            $6, $7, $8, $9, pid
  }
' "$infile" >> "$outbed"

echo "Wrote $outbed"


# Run as a loop
```
########################################################################
# Resolved species
########################################################################
# Corynorhinus_townsendii
echo 'Corynorhinus_townsendii,X,CM133721.1' >> ${sexchr_csv}
echo 'Corynorhinus_townsendii,Y,CM133722.1' >> ${sexchr_csv}
# Lycocorax_pyrrhopterus
echo 'Lycocorax_pyrrhopterus,Z,CM133637.1' >> ${sexchr_csv}
echo 'Lycocorax_pyrrhopterus,W,CM133636.1' >> ${sexchr_csv}
# Macaca_fascicularis
echo 'Macaca_fascicularis,Y,NC_132903.1' >> ${sexchr_csv}
# Molossus_molossus
echo 'Molossus_molossus,X,CM138289.1' >> ${sexchr_csv}
echo 'Molossus_molossus,Y,CM138290.1' >> ${sexchr_csv}
# Peromyscus_maniculatus
echo 'Peromyscus_maniculatus,X,CM138667.1' >> ${sexchr_csv}
echo 'Peromyscus_maniculatus,Y,CM138668.1' >> ${sexchr_csv}
# Rousettus_aegyptiacus
echo 'Rousettus_aegyptiacus,X,CM133704.1' >> ${sexchr_csv}
echo 'Rousettus_aegyptiacus,Y,CM133705.1' >> ${sexchr_csv}
# Tupaia_tana
echo 'Tupaia_tana,X,CM133136.1' >> ${sexchr_csv}
echo 'Tupaia_tana,Y,CM133137.1' >> ${sexchr_csv}

# Elephas_maximus_indicus
sed -i 's/NC_064846/CM123151/g' ${sexchr_csv}
sed -i 's/NC_064847/CM123152/g' ${sexchr_csv}
# Lemur_catta
sed -i 's/NC_059155/CM036499/g' ${sexchr_csv}
sed -i 's/NC_059156/CM036500/g' ${sexchr_csv}
# Pelecanus_crispus
sed -i 's/CM060033/NC_134676/g' ${sexchr_csv}
sed -i 's/CM060032/NC_134675/g' ${sexchr_csv}
# Tursiops_truncatus
sed -i 's/NW_022983135/NC_133428/g' ${sexchr_csv}

# Monodon_monoceros
mv ${base}/Monodon_monocero/ ${base}/Monodon_monoceros
mv ${base}/Monodon_monoceros/Monodon_monocero.fna ${base}/Monodon_monoceros/Monodon_monoceros.fna

########################################################################
# Unsolved species
########################################################################

species=Tupaia_tana
grep '>' ${base}/${species}/${species}.fna | grep 'chromosome'
grep ${species} ${sexchr_csv}

# Anniella_stebbinsi -- missing genome
# Dama_dama  -- only has X
# Lycaon_pictus -- only has X
# Mustelus_asterias -- only has X
# Narcine_bancroftii -- only has X (Y in alt hap)
# Ornithorhynchus_anatinus -- Neo sex
# Tachyglossus_aculeatus -- Neo sex
# Thamnophis_elegans -- only has Z
# Trichosurus_vulpecula -- only has X
########################################################################

species=Tupaia_tana
grep '>' ${base}/${species}/${species}.fna | grep 'chromosome'
grep ${species} ${sexchr_csv}


# Make paf for species that were missed

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"
workdir="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles"
sexchr_csv="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"
missing_file="heterogametic_missing_in_paf.txt"

cd "${workdir}" || exit 1

mkdir -p XZ_chrs YW_chrs minimap2 logs

while read -r species; do
    [[ -z "${species}" ]] && continue

    echo "Processing ${species}"

    # Get accessions
    x_acc=$(awk -F',' -v sp="${species}" '$1 == sp && $2 == "X" {print $3}' "${sexchr_csv}")
    y_acc=$(awk -F',' -v sp="${species}" '$1 == sp && $2 == "Y" {print $3}' "${sexchr_csv}")
    z_acc=$(awk -F',' -v sp="${species}" '$1 == sp && $2 == "Z" {print $3}' "${sexchr_csv}")
    w_acc=$(awk -F',' -v sp="${species}" '$1 == sp && $2 == "W" {print $3}' "${sexchr_csv}")

    if [[ -n "${x_acc}" && -n "${y_acc}" ]]; then
        ref_chr="X"
        query_chr="Y"
        ref_acc="${x_acc}"
        query_acc="${y_acc}"
        paf_name="${species}_YtoX.aln.paf"

    elif [[ -n "${z_acc}" && -n "${w_acc}" ]]; then
        ref_chr="Z"
        query_chr="W"
        ref_acc="${z_acc}"
        query_acc="${w_acc}"
        paf_name="${species}_WtoZ.aln.paf"

    else
        echo "WARNING: Could not find complete X/Y or Z/W pair for ${species}; skipping"
        continue
    fi

    # Write accession files
    echo "${ref_acc}" > "${species}.${ref_chr}chr.txt"
    echo "${query_acc}" > "${species}.${query_chr}chr.txt"

    # Extract chromosomes
    samtools faidx \
        "${base}/${species}/${species}.fna" \
        -r "${species}.${ref_chr}chr.txt" \
        > "XZ_chrs/${species}.${ref_chr}chr.fasta"

    samtools faidx \
        "${base}/${species}/${species}.fna" \
        -r "${species}.${query_chr}chr.txt" \
        > "YW_chrs/${species}.${query_chr}chr.fasta"
done < "${missing_file}"


while read -r species; do
    # Submit minimap2 job: heterogametic chromosome as query, homogametic chromosome as target
    sbatch \
        --cpus-per-task=3 \
        --mem=8G \
        --time=12:00:00 \
        --job-name="${species}_${query_chr}to${ref_chr}" \
        --output="logs/${species}_${query_chr}to${ref_chr}.%j.out" \
        --error="logs/${species}_${query_chr}to${ref_chr}.%j.err" \
        --wrap="minimap2 -x asm5 -c -eqx -k 23 -t 3 \
        ${workdir}/XZ_chrs/${species}.${ref_chr}chr.fasta \
        ${workdir}/YW_chrs/${species}.${query_chr}chr.fasta \
        > ${workdir}/minimap2/${paf_name}"

done < "${missing_file}"
```