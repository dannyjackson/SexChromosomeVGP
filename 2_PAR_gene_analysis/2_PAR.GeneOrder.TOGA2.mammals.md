# Identify all genes found in any avian PAR, curate fastas for blast analysis, then identify genes within PARs of all genomes
## 0. Quantify gappiness of each PAR
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/TOGA2_annotations/mammals

TOGA_DIR="/data/Wilson_Lab/data/TOGA2_Hiller/Homo_sapiens_hg38"

PAF_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/continuous_percentID"

rm -r Callithrix_jacchus__white-tufted-ear_marmoset__HLcalJac5__GCF_011100555.1
rm -r Callithrix_jacchus__white-tufted-ear_marmoset__calJac4__GCF_009663435.1/
rm -r Pan_troglodytes__chimpanzee__panTro6__GCF_002880755.1/
rm -r Loxodonta_africana__African_savanna_elephant__HLloxAfr5B__GCA_030020305.1/
rm -r Urocitellus_parryii__Arctic_ground_squirrel__HLuroPar1__GCA_003426925.1/
rm -r Urocitellus_parryii__Arctic_ground_squirrel__HLuroPar2B__GCA_045843765.1/
rm -r Pongo_pygmaeus__Bornean_orangutan__HLponPyg1__GCA_023767775.1/
rm -r Symphalangus_syndactylus__siamang__HLsymSyn4__GCA_028878055.3/
rm -r Gorilla_gorilla_gorilla__western_lowland_gorilla__gorGor6__GCF_008122165.1/
rm -r Pongo_abelii__Sumatran_orangutan__ponAbe3__GCF_002880775.1/
rm -r Myotis_nattereri__Natterers_bat__HLmyoNat1A__GCA_964212035.1
rm -r Myotis_nattereri__Natterers_bat__HLmyoNat1B__GCA_964212025.1
rm -r Myotis_nattereri__Natterers_bat__HLmyoNatt2B__GCA_964212025.2
rm -r Rhynchonycteris_naso__Proboscis_bat__HLrhyNas2B__GCA_037038555.1

cat > PAR.species_chr_region.gaps.txt <<'EOF'
Balaenoptera_physalus,OZ239531:0-7257636
Bos_taurus,NC_040105:0-6843483
Callithrix_jacchus,CM111807:0-1757279
Camelus_dromedarius,NC_087472:108540129-114202744
Capra_hircus,CP168640:0-7202443
Eubalaena_glacialis,NC_083736:0-7069966
Gorilla_gorilla,NC_073247:0-12039748
Grampus_griseus,OZ206318:0-7148355
Homo_sapiens,NC_060947:0-2394410
Inia_geoffrensis,CM070920:0-7142434
Loxodonta_africana,CM057446:0-10314587
Lycaon_pictus,CM082710:0-6575129
Macaca_nemestrina,NC_092145:158114029-159757195 
Manis_pentadactyla,NC_080038:0-5004179
Marmota_flaviventris,NC_092518:0-10556221 
Meles_meles,OV277448:0-6399546
Mesoplodon_bidens,OZ073217:135117924-142816029
Molossus_nigricans,CM078089:0-4370896
Mustela_nivalis_vulgaris,CM169857:131948904-138477583
Myotis_nattereri,OZ125678:0-2687171
Ovis_aries,CP162266:0-7021881
Ovis_canadensis,NC_091727:0-7072606
Pan_paniscus,CM055495:0-2524164 
Pan_troglodytes,CM054457:0-3170188 
Panthera_onca,CM102116:122788602-130846311
Pongo_abelii,NC_072008:0-2382235
Pongo_pygmaeus,CM054653:0-2356740
Pseudorca_crassidens,NC_090317:127639906-136059651 
Rhynchocyon_petersi,CM091802:0-20182068
Rhynchonycteris_naso,CM073052:138579092-142670400
Symphalangus_syndactylus,NC_072447:0-16994066
Trichechus_inunguis,CM102173:0-9624845
Urocitellus_parryii,CM099876:122998411-131433435
Canis_lupus_baileyi,NC_132876:118611595-125236116
Nyctalus_leisleri,OZ183628:0-1665649
Dasypus_novemcinctus,NC_080704:0-9679476
Artibeus_intermedius,CM076326:147359214-151570710
Artibeus_lituratus,CM076392:0-4918541
Myotis_mystacinus,OZ075425:125458419-127971689
Corynorhinus_townsendii,CM133721.1:0-1829836
Miniopterus_schreibersii,OZ071095:103813303-106815783
EOF
```
# Identify gaps in the assembly indicative of a Hi-C joining
```
chmod +x 2a_find_PAR_gaps.sh 
./2a_find_PAR_gaps.sh

Rscript 2b_plot_PAR_gaps.R

# Make a file containing just the gaps in the sex chromosomes,
out="PAR_gap_results/SexChr_species_chr_gaps.bed" 
: > "$out"

while IFS=',' read -r species region; do
    # Skip blank lines
    [ -z "$species" ] && continue
    echo $species
    # Trim whitespace from region
    region=$(printf '%s' "$region" | tr -d '[:space:]')

    # Extract chromosome from chr:start-end
    chr=${region%%:*}

    # Match the corresponding species gap file
    bed="PAR_gap_results/per_species/${species}.all_N.bed"

    if [ -f "$bed" ]; then
        awk -v species="$species" -v chr="$chr" '
            $1 == chr {
                print species, $0
            }
        ' OFS='\t' "$bed" >> "$out"
    else
        echo "Warning: missing file for $species: $bed" >&2
    fi
done < PAR.species_chr_region.gaps.txt

# Condense gaps to contiguous sequences

module load bedtools

PAR_FILE="PAR.species_chr_region.gaps.txt"
indir="PAR_gap_results"
infile="$indir/SexChr_species_chr_gaps.bed"

out_merged="$indir/PAR_only_species_chr_gaps.merged.bed"
summary="$indir/PAR_only_species_chr_gaps.merged.summary.tsv"

mkdir -p "$indir"/per_species "$indir"/logs

: > "$out_merged"

echo -e "species\tchrom\tPAR_start\tPAR_end\tPAR_len\tn_gaps\ttotal_gap_bp\tmax_gap_bp" > "$summary"

echo "Merging contiguous gaps from $infile..."

while IFS=, read -r species region_raw; do
    [[ -z "${species:-}" ]] && continue

    # Trim whitespace and repair accidental quote typo, if present
    region="$(printf '%s' "$region_raw" | tr -d '[:space:]' | sed 's/"/:/g')"

    chrom="${region%%:*}"
    coords="${region#*:}"
    start="${coords%-*}"
    end="${coords#*-}"

    par_len=$((end - start))

    tmp_raw="$indir/per_species/${species}.PAR_only_chr.gaps.raw.bed"
    tmp_merged="$indir/per_species/${species}.PAR_only_chr.gaps.merged.bed"

    echo "  $species  $chrom:$start-$end"

    # Pull this species/chromosome from the combined file,
    # trim to the PAR interval, then merge adjacent/overlapping gap bases.
    #
    # Input columns are expected to be:
    # species chrom gap_start gap_end name score strand
    awk -v sp="$species" -v c="$chrom" -v s="$start" -v e="$end" '
        BEGIN { OFS="\t" }
        $1 == sp && $2 == c && $4 > s && $3 < e {
            a = ($3 < s ? s : $3)
            b = ($4 > e ? e : $4)
            if (b > a) print $2, a, b
        }
    ' "$infile" \
        | sort -k1,1 -k2,2n \
        | bedtools merge -d 1 -i - \
        > "$tmp_merged"

    # Add species name and relative PAR coordinates.
    # Output columns:
    # species chrom gap_start gap_end gap_len PAR_start PAR_end rel_gap_start rel_gap_end
    awk -v sp="$species" -v ps="$start" -v pe="$end" '
        BEGIN { OFS="\t" }
        {
            print sp, $1, $2, $3, $3-$2, ps, pe, $2-ps, $3-ps
        }
    ' "$tmp_merged" >> "$out_merged"

    n_gaps=$(wc -l < "$tmp_merged" | tr -d ' ')
    total_gap_bp=$(awk '{sum += $3-$2} END {print sum+0}' "$tmp_merged")
    max_gap_bp=$(awk 'BEGIN {max=0} {len=$3-$2; if (len>max) max=len} END {print max+0}' "$tmp_merged")

    echo -e "${species}\t${chrom}\t${start}\t${end}\t${par_len}\t${n_gaps}\t${total_gap_bp}\t${max_gap_bp}" >> "$summary"

done < "$PAR_FILE"

echo "Done."
echo "Merged gap BED-like table: $out_merged"
echo "Summary table:            $summary"

```
## PAR boundaries
### Identify PAR boundaries in bat species
#### Artibeus_lituratus
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/TOGA2_annotations/mammals/Bats

awk 'BEGIN{OFS="\t"} $11>=10000 && ($10/$11)>0.90 {print $1,$3,$4}' \
"Artibeus_lituratus_GCA_038363095.4_YtoX.aln.paf" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> Artibeus_lituratus.X_gt90_10kb_merged.bed
```
#### Artibeus_intermedius
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/TOGA2_annotations/mammals/Bats

awk 'BEGIN{OFS="\t"} $11>=10000 && ($10/$11)>0.90 {print $1,$3,$4}' \
"Artibeus_intermedius_Y1toX.aln.paf" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> Artibeus_intermedius.X_gt90_10kb_merged.bed
```
#### Corynorhinus_townsendii
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/TOGA2_annotations/mammals/Bats

awk 'BEGIN{OFS="\t"} $11>=10000 && ($10/$11)>0.90 {print $1,$3,$4}' \
"Corynorhinus_townsendii_GCA_026230055.2_hap1_Y_to_GCA_026230055.2_hap1_X.aln.paf" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> Corynorhinus_townsendii.X_gt90_10kb_merged.bed
```
#### Miniopterus_schreibersii
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/TOGA2_annotations/mammals/Bats

awk 'BEGIN{OFS="\t"} $11>=10000 && ($10/$11)>0.90 {print $1,$3,$4}' \
"Corynorhinus_townsendii_GCA_026230055.2_hap1_Y_to_GCA_026230055.2_hap1_X.aln.paf" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> Corynorhinus_townsendii.X_gt90_10kb_merged.bed
```
#### Nyctalus_leisleri
```
awk 'BEGIN{OFS="\t"} $11>=10000 && ($10/$11)>0.90 {print $1,$3,$4}' \
"Nyctalus_leisleri_YtoX.aln.paf" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> Nyctalus_leisleri.X_gt90_10kb_merged.bed
```
#### Myotis_mystacinus 
```
awk 'BEGIN{OFS="\t"} $11>=10000 && ($10/$11)>0.90 {print $1,$3,$4}' \
"Myotis_mystacinus_GCA_964094495.3_hap1_Y_to_GCA_964094495.3_hap1_X.aln.paf" \
| sort -k1,1 -k2,2n \
| bedtools merge -i - \
> Myotis_mystacinus.X_gt90_10kb_merged.bed
```
# Identify PAR genes in the TOGA annotations
```
cat > PAR.species_chr_region.txt <<'EOF'
Artibeus_lituratus,CM076392:0-4918541
Balaenoptera_physalus,OZ239531:0-7257636
Bos_taurus,NC_040105:0-6843483
Callithrix_jacchus,CM111807:0-1757279
Camelus_dromedarius,NC_087472:108540129-114202744
Capra_hircus,CP168640:0-7202443
Eubalaena_glacialis,NC_083736:0-7069966
Gorilla_gorilla,NC_073247:0-12039748
Grampus_griseus,OZ206318:0-7148355
Homo_sapiens,NC_060947:0-2394410
Inia_geoffrensis,CM070920:0-7142434
Loxodonta_africana,CM057446:0-10314587
Lycaon_pictus,CM082710:0-6575129
Macaca_nemestrina,NC_092145:158114029-159757195
Manis_pentadactyla,NC_080038:0-5004179
Marmota_flaviventris,NC_092518:0-10556221
Meles_meles,OV277448:0-6399546
Mesoplodon_bidens,OZ073217:135117924-142816029
Molossus_nigricans,CM078089:0-4370896
Mustela_nivalis_vulgaris,CM169857:131948904-138477583
Myotis_mystacinus,OZ075425:125458419-127971689
Myotis_nattereri,OZ125678:0-2687171
Ovis_aries,CP162266:0-7021881
Ovis_canadensis,NC_091727:0-7072606
Pan_paniscus,CM055495:0-2524164
Pan_troglodytes,CM054457:0-3170188
Panthera_onca,CM102116:122788602-130846311
Pongo_abelii,NC_072008:0-2382235
Pongo_pygmaeus,CM054653:0-2356740
Pseudorca_crassidens,NC_090317:127639906-136059651
Rhynchocyon_petersi,CM091802:0-20182068
Rhynchonycteris_naso,CM073052:138579092-142670400
Symphalangus_syndactylus,NC_072447:0-16994066
Trichechus_inunguis,CM102173:0-9624845
Urocitellus_parryii,CM099876:122998411-131433435
Artibeus_intermedius,CM076326:147359214-151570710
Miniopterus_schreibersii,OZ071095:103813303-106815783
Canis_lupus_baileyi,NC_132876:118611595-125236116
Nyctalus_leisleri,OZ183628:0-1665649
Dasypus_novemcinctus,NC_080704:0-9679476
Corynorhinus_townsendii,CM133721.1:0-1829836
EOF


chmod +x 2c_find_par_genes_toga.mammals.awk 
./2c_find_par_genes_toga.mammals.awk PAR.species_chr_region.txt > mammals_PAR_genes.tsv


awk -F',' '{print $1}' PAR.species_chr_region.txt > species.txt

while read -r species; do
  echo $species
  grep $species mammals_PAR_genes.tsv | wc -l
done < species.txt

# Redo the one(s) that didn't work due to naming issues

grep Mustela_nivalis PAR.species_chr_region.txt | sed 's/_vulgaris//g' > PAR.species_chr_region.redo.txt
./2c_find_par_genes_toga.mammals.awk PAR.species_chr_region.redo.txt > mammals_PAR_genes.redo.tsv


# Redo the ones that didn't work because of chr naming in TOGA2 files
cat > species.chrX.txt <<'EOF'
Pan_paniscus
Pan_troglodytes
Pongo_pygmaeus
Gorilla_gorilla
Pongo_abelii
Symphalangus_syndactylus
EOF


rm PAR.species_chr_region.chrX.txt
while read -r species; do
  grep $species PAR.species_chr_region.txt >> PAR.species_chr_region.chrX.txt
done < species.chrX.txt

chmod +x 2d_find_par_genes_toga.chrX.awk
./2d_find_par_genes_toga.chrX.awk PAR.species_chr_region.chrX.txt > mammals_PAR_genes.chrX.tsv

# Add human
gff="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Homo_sapiens/ncbi_dataset/data/GCF_009914755.1/genomic.gff"
out="Human_PAR_genes.tsv"

awk -F'\t' '
BEGIN {
    OFS = "\t"

    species = "Homo_sapiens"
    togadir = "Homo_sapiens__human__T2T-CHM13v2.0__GCF_009914755.1"

    par_chrom = "NC_060947.1"
    par_start = 1
    par_end = 2394410

    print "Species", "TOGADir", "Chromosome", "StartPos", "StopPos", "GeneName", "InPAR"
}

function get_attr(attrs, key,    n, a, i, kv) {
    n = split(attrs, a, ";")
    for (i = 1; i <= n; i++) {
        split(a[i], kv, "=")
        if (kv[1] == key) {
            return kv[2]
        }
    }
    return ""
}

function par_status(start, stop, par_start, par_end) {
    if (start >= par_start && stop <= par_end) {
        return "Y"
    }

    if (stop < par_start || start > par_end) {
        return "N"
    }

    return "Edge"
}

/^#/ { next }

$1 == par_chrom && $3 == "gene" {
    attrs = $9

    gene_biotype = get_attr(attrs, "gene_biotype")
    gbkey = get_attr(attrs, "gbkey")

    if (gbkey != "Gene") {
        next
    }

    if (gene_biotype != "protein_coding") {
        next
    }

    start = $4 + 0
    stop = $5 + 0

    gene = get_attr(attrs, "gene")
    if (gene == "") gene = get_attr(attrs, "Name")
    if (gene == "") gene = get_attr(attrs, "gene_name")
    if (gene == "") gene = get_attr(attrs, "ID")

    gsub(/%20/, " ", gene)

    chrom = $1
    sub(/\..*/, "", chrom)

    inpar = par_status(start, stop, par_start, par_end)

    print species, togadir, chrom, start, stop, gene, inpar
}
' "$gff" > "$out"


# Combine all PAR gene files

sort -u mammals_PAR_genes.chrX.tsv > mammals_PAR_genes.chrX.tsv.tmp
mv mammals_PAR_genes.chrX.tsv.tmp mammals_PAR_genes.chrX.tsv

cat mammals_PAR_genes.tsv | grep -v Desmodus > mammals_PAR_genes.all.tsv 
tail -n +2 mammals_PAR_genes.redo.tsv >> mammals_PAR_genes.all.tsv 
tail -n +2 mammals_PAR_genes.chrX.tsv >> mammals_PAR_genes.all.tsv 
tail -n +2 Human_PAR_genes.tsv >> mammals_PAR_genes.all.tsv 
sed -i 's/Mustela_nivalis/Mustela_nivalis_vulgaris/g' mammals_PAR_genes.all.tsv 


# check that all species ran
while read -r species; do
  echo $species
  grep $species mammals_PAR_genes.all.tsv | wc -l
done < species.txt

```
# Prepare data for plotting
```
# Clean data to standardize genes named for X and Y locations in human
sed -i 's/NLGN4Y/NLGN4/g' mammals_PAR_genes.all.tsv
sed -i 's/NLGN4X/NLGN4/g' mammals_PAR_genes.all.tsv
sed -i 's/TBL1Y/TBL1/g' mammals_PAR_genes.all.tsv
sed -i 's/TBL1X/TBL1/g' mammals_PAR_genes.all.tsv
sed -i 's/DDX3Y/DDX3/g' mammals_PAR_genes.all.tsv
sed -i 's/DDX3X/DDX3/g' mammals_PAR_genes.all.tsv
sed -i 's/USP9Y/USP9/g' mammals_PAR_genes.all.tsv
sed -i 's/USP9X/USP9/g' mammals_PAR_genes.all.tsv
sed -i 's/AMELY/AMEL/g' mammals_PAR_genes.all.tsv
sed -i 's/AMELX/AMEL/g' mammals_PAR_genes.all.tsv

{ head -n 1 mammals_PAR_genes.all.tsv; tail -n +2 mammals_PAR_genes.all.tsv | sort -u; } > mammals_PAR_genes.all.tsv.tmp
mv mammals_PAR_genes.all.tsv.tmp mammals_PAR_genes.all.tsv

wc -l mammals_PAR_genes.all.tsv
grep -v Desmodus mammals_PAR_genes.all.tsv > mammals_PAR_genes.all.tsv.tmp
mv mammals_PAR_genes.all.tsv.tmp mammals_PAR_genes.all.tsv

# Modify the file to combine major blocks of gene family expansions into a single point
{
  head -n 1 mammals_PAR_genes.all.tsv
  tail -n +2 mammals_PAR_genes.all.tsv | sort -t $'\t' -k1,1 -k2,2 -k3,3 -k4,4n
} | awk 'BEGIN{FS=OFS="\t"}
NR==1 { print; next }

function flush_block() {
    if (in_block) {
        split(block_row, f, OFS)
        f[4] = block_start
        f[5] = block_stop
        f[6] = block_name
        print f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]
        in_block = 0
    }
}

{
    key = $1 OFS $2 OFS $3

    is_array = 0
    name = ""

    if ($6 ~ /^ZNF/) {
        is_array = 1
        name = "ZF_array"
    } else if ($6 ~ /^VC/) {
        is_array = 1
        name = "VC_array"
    } else if ($6 ~ /^SET/) {
        is_array = 1
        name = "SET_array"
    }

    if (is_array) {
        if (!in_block || key != block_key || name != block_name) {
            flush_block()
            block_row = $0
            block_start = $4
            block_stop = $5
            block_key = key
            block_name = name
            in_block = 1
        } else {
            if ($4 < block_start) block_start = $4
            if ($5 > block_stop) block_stop = $5
        }
    } else {
        flush_block()
        print
    }
}

END {
    flush_block()
}' > mammals_PAR_genes.all.ZNF_arrays.tsv

cp mammals_PAR_genes.all.ZNF_arrays.tsv mammals_PAR_genes.all.ZNF_arrays.tsv.save

grep -v 'LOC' mammals_PAR_genes.all.ZNF_arrays.tsv | grep -v 'LINC' > mammals_PAR_genes.all.ZNF_arrays.tsv.noLOC
mv mammals_PAR_genes.all.ZNF_arrays.tsv.noLOC mammals_PAR_genes.all.ZNF_arrays.tsv
```
## set XG as the origin gene
```
awk -v GENE="XG" -f 2e_orient_to_gene_midpoint.awk \
  mammals_PAR_genes.all.ZNF_arrays.tsv \
  mammals_PAR_genes.all.ZNF_arrays.tsv \
  > mammals_PAR_genes.XG_oriented.tsv

# change 2nd shroom2 to shroom2-like
awk 'BEGIN{OFS="\t"} $1=="Mesoplodon_bidens" && $6=="SHROOM2" {n++; if (n==2) $6="SHROOM2-like"} {print}' \
  mammals_PAR_genes.all.ZNF_arrays.tsv > mammals_PAR_genes.all.ZNF_arrays.edited.tsv

cp mammals_PAR_genes.all.ZNF_arrays.tsv mammals_PAR_genes.all.ZNF_arrays.tsv.unedited
mv mammals_PAR_genes.all.ZNF_arrays.edited.tsv mammals_PAR_genes.all.ZNF_arrays.tsv
```
### Note if telomere is ID'd on the end of the PAR
```
telomere_file="species.telomeres.txt"

cat > "$telomere_file" <<'EOF'
Bos_taurus,L,NO
Lycaon_pictus,L,YES
Panthera_onca,L,YES
Eubalaena_glacialis,L,YES
Macaca_nemestrina,R,YES
Callithrix_jacchus,L,YES
Mesoplodon_bidens,R,NO
Inia_geoffrensis,L,YES
Camelus_dromedarius,R,YES
Ovis_canadensis,L,YES
Manis_pentadactyla,L,YES
Pan_paniscus,L,YES
Marmota_flaviventris,L,YES
Balaenoptera_physalus,L,YES
Pan_troglodytes,L,YES
Rhynchonycteris_naso,R,YES
Ovis_aries,L,YES
Trichechus_inunguis,L,YES
Pseudorca_crassidens,R,YES
Myotis_nattereri,L,NO
Rhynchocyon_petersi,L,YES
Grampus_griseus,L,YES
Mustela_nivalis_vulgaris,R,YES
Capra_hircus,L,YES
Loxodonta_africana,L,YES
Urocitellus_parryii,R,NO
Meles_meles,L,YES
Homo_sapiens,L,YES
Pongo_pygmaeus,L,YES
Symphalangus_syndactylus,L,YES
Gorilla_gorilla,L,YES
Pongo_abelii,L,YES
EOF
```
# Create the plot
```
Rscript 2f_PAR_GeneOrder_Mammals.August2026.R
```