# Identify all genes found in any avian PAR, curate fastas for blast analysis, then identify genes within PARs of all genomes
## 0. Quantify gappiness of each PAR
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/TOGA2_annotations/mammals

TOGA_DIR="/data/Wilson_Lab/data/TOGA2_Hiller/Homo_sapiens_hg38"

## Excluding:
Microtus_pennsylvanicus,CHROM:0-116367
Ochotona_princeps,CHROM:107724012-107801096 



# Confirm that accessions in sexchrfile still match correctly (genomes have been updated)
SPECIES=Urocitellus_parryii
grep '>' ${BASE}/${SPECIES}/${SPECIES}.fna | grep 'chromosome'

# Revise those that need revision
SEXCHR_FILE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"

grep 'Urocitellus_parryii' $SEXCHR_FILE
sed -i 's/CM099876/NC_135547/g' $SEXCHR_FILE
sed -i 's/CM099877/NC_135548/g' $SEXCHR_FILE

Capra_hircus # no matches at all

# These species need accession numbers changed from RefSeq to GenBank
Callithrix_jacchus          GCA_049354715.1    GCF_049354715.1     Yes:GCA-GCF
Pan_paniscus                GCA_029289425.3    GCF_029289425.2     Yes:GCA-GCF
Pan_troglodytes             GCA_028858775.2    GCF_028858775.2     Yes:GCA-GCF
Loxodonta_africana          GCA_030014295.1    GCF_030014295.1     Yes:GCA-GCF
Urocitellus_parryii         GCA_045843805.1    GCF_045843805.1     Yes:GCA-GCF
Meles_meles                 GCA_922984935.1    GCF_922984935.1     Yes:GCA-GCF
Pongo_pygmaeus              GCA_028885625.2    GCF_028885625.2     Yes:GCA-GCF
Symphalangus_syndactylus    GCF_028878055.3    GCA_028878055.3
Gorilla_gorilla             GCA_029281585.3    GCF_029281585.2    Yes:GCA-GCF
Pongo_abelii                GCF_028885655.2    GCA_028885655.2    Yes:GCA-GCF

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
Eubalaena_glacialis,NC_083736.1:0-7069966
Macaca_nemestrina,NC_092145.1:158149255-159757195 
Callithrix_jacchus,CM111807.1.1:0-1757279
Mesoplodon_bidens,OZ073217.1:135117924-142816029
Inia_geoffrensis,CM070920.1:0-7142434
Camelus_dromedarius,NC_087472.1:108540129-114202744
Ovis_canadensis,NC_091727.1:0-7072606
Manis_pentadactyla,NC_080038.1:0-4881534
Pan_paniscus,CM055495.2:0-2524164                      
Marmota_flaviventris,NC_092518.1:521224-10347709        
Balaenoptera_physalus,OZ239531.1:0-7257636
Pan_troglodytes,CM054457.2:0-3170188          
Rhynchonycteris_naso,CM073052.1:138579092-142670400
Ovis_aries,CP162266.1:0-7021881
Trichechus_inunguis,CM102173.1:0-9558467
Pseudorca_crassidens,NC_090317.1:127639906-136059651 
Myotis_nattereri,OZ125678.2:0-1449976
Rhynchocyon_petersi,CM091802.1:0-20182068
Grampus_griseus,OZ206318.1:0-7148355
Mustela_nivalis_vulgaris,OZ211688.1:127035994-138477583
Capra_hircus,CP168640.1:0-7040093
Loxodonta_africana,CM057446.1:0-10314587
Urocitellus_parryii,CM099876.1:122998411-131433435 
Meles_meles,OV277448.1:0-6387864
Homo_sapiens,NC_060947.1:0-2394410
Pongo_pygmaeus,CM054653.2.2:0-2356740
Symphalangus_syndactylus,CM054533.2:0-16994066
Gorilla_gorilla,CM054702.2:0-12039748
Pongo_abelii,NC_072008.2:0-2382235
EOF


chmod +x find_PAR_gaps.sh 
./find_PAR_gaps.sh

Rscript plot_PAR_gaps.R

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
### Species without gaps
Loxodonta africana has four gaps each of 13 base pairs... this feels a little suspicious to me? Since assemblers can put in any random number of bp for a gap? Excluding it for now.
```
cat > species_par.csv <<'EOF'
Pseudorca_crassidens,127639906-136059651 
Inia_geoffrensis,0-7142434
Ovis_canadensis,0-7072606
Capra_hircus,0-7040093
Ovis_aries,0-7021881
Meles_meles,0-6387864
Pan_troglodytes,0-3170188          
Pan_paniscus,0-2524164                      
Homo_sapiens,0-2394410
Callithrix_jacchus,0-1757279
Pongo_pygmaeus,NC_072396.2:0-2356740
Symphalangus_syndactylus,NC_072447.2:0-16994066
Gorilla_gorilla,NC_073247.2:0-12039748
Pongo_abelii,NC_072008.2:0-2382235
EOF
```
## 1. Curate annotated genes found in the PARs
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/GFFs/mammals


cat > PAR.species_chr_region.txt <<'EOF'
Eubalaena_glacialis,NC_083736:0-7069966
Macaca_nemestrina,NC_092145:158149255-159757195 
Callithrix_jacchus,CM111807:0-1757279
Mesoplodon_bidens,OZ073217:135117924-142816029
Inia_geoffrensis,CM070920:0-7142434
Camelus_dromedarius,NC_087472:108540129-114202744
Ovis_canadensis,NC_091727:0-7072606
Manis_pentadactyla,NC_080038:0-4881534
Pan_paniscus,CM055495:0-2524164                      
Marmota_flaviventris,NC_092518:521224-10347709        
Balaenoptera_physalus,OZ239531:0-7257636
Pan_troglodytes,CM054457:0-3170188          
Rhynchonycteris_naso,CM073052.1:138579092-142670400
Ovis_aries,CP162266:0-7021881
Trichechus_inunguis,CM102173:0-9558467
Pseudorca_crassidens,NC_090317:127639906-136059651 
Myotis_nattereri,OZ125678:0-1449976
Rhynchocyon_petersi,CM091802:0-20182068
Grampus_griseus,OZ206318:0-7148355
Mustela_nivalis_vulgaris,OZ211688:127035994-138477583
Capra_hircus,CP168640:0-7040093
Loxodonta_africana,CM057446:0-10314587
Urocitellus_parryii,CM099876:122998411-131433435 
Meles_meles,OV277448:0-6387864
Homo_sapiens,NC_060947:0-2394410
Pongo_pygmaeus,CM054653:0-2356740
Symphalangus_syndactylus,CM054533:0-16994066
Gorilla_gorilla,CM054702:0-12039748
Pongo_abelii,NC_072008:0-2382235
EOF

TOGA_DIR="/data/Wilson_Lab/data /TOGA2_Hiller/Homo_sapiens_hg38" 

chmod +x find_par_genes_toga.awk 
./find_par_genes_toga.awk PAR.species_chr_region.txt > mammals_PAR_genes.tsv



cat > species.txt <<'EOF'
Eubalaena_glacialis
Macaca_nemestrina
Callithrix_jacchus
Mesoplodon_bidensfvg
Inia_geoffrensis
Camelus_dromedarius
Ovis_canadensis
Manis_pentadactyla
Pan_paniscus
Marmota_flaviventris
Balaenoptera_physalus
Pan_troglodytes
Rhynchonycteris_naso
Ovis_aries
Trichechus_inunguis
Pseudorca_crassidens
Myotis_nattereri
Rhynchocyon_petersi
Grampus_griseus
Mustela_nivalis_vulgaris
Capra_hircus
Loxodonta_africana
Urocitellus_parryii
Meles_meles
Homo_sapiens
Pongo_pygmaeus
Symphalangus_syndactylus
Gorilla_gorilla
Pongo_abelii
EOF

while read -r species; do
  echo $species
  grep $species mammals_PAR_genes.tsv | wc -l
done < species.txt

# Redo the one(s) that didn't work due to naming issues

grep Mustela_nivalis PAR.species_chr_region.txt | sed 's/_vulgaris//g' > PAR.species_chr_region.redo.txt
./find_par_genes_toga.awk PAR.species_chr_region.redo.txt > mammals_PAR_genes.redo.tsv


grep Rhynchonycteris_naso PAR.species_chr_region.txt > PAR.species_chr_region.redo.txt
./find_par_genes_toga.awk PAR.species_chr_region.redo.txt >> mammals_PAR_genes.redo.tsv


# Redo the ones that didn't work because of chr naming in TOGA2 files
cat > species.chrX.txt <<'EOF'
Pan_paniscus
Pan_troglodytes
Pongo_pygmaeus
Symphalangus_syndactylus
Gorilla_gorilla
Pongo_abelii
EOF


while read -r species; do
  grep $species PAR.species_chr_region.txt >> PAR.species_chr_region.chrX.txt
done < species.chrX.txt

chmod +x find_par_genes_toga.chrX.awk
./find_par_genes_toga.chrX.awk PAR.species_chr_region.chrX.txt > mammals_PAR_genes.chrX.tsv

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


grep -v 'Myotis' mammals_PAR_genes.tsv > mammals_PAR_genes.tmp
grep 'GCA_964212035.2' mammals_PAR_genes.tsv >> mammals_PAR_genes.tmp
mv mammals_PAR_genes.tmp mammals_PAR_genes.tsv

sort -u mammals_PAR_genes.chrX.tsv > mammals_PAR_genes.chrX.tsv.tmp
mv mammals_PAR_genes.chrX.tsv.tmp mammals_PAR_genes.chrX.tsv

cat mammals_PAR_genes.tsv > mammals_PAR_genes.all.tsv 
tail -n +2 mammals_PAR_genes.redo.tsv >> mammals_PAR_genes.all.tsv 
tail -n +2 mammals_PAR_genes.chrX.tsv >> mammals_PAR_genes.all.tsv 
tail -n +2 Human_PAR_genes.tsv >> mammals_PAR_genes.all.tsv 

sed -i 's/Mustela_nivalis/Mustela_nivalis_vulgaris/g' mammals_PAR_genes.all.tsv 

```
# Plot
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

Rscript plot_combined.allgenelabels.mammals.R

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

Rscript plot_combined.mammals.R
Rscript plot_combined.allgenelabels.mammals.GeneArrays.R

```

cp mammals_PAR_genes.all.ZNF_arrays.tsv mammals_PAR_genes.all.ZNF_arrays.tsv.save

grep -v 'LOC' mammals_PAR_genes.all.ZNF_arrays.tsv | grep -v 'LINC' > mammals_PAR_genes.all.ZNF_arrays.tsv.noLOC
mv mammals_PAR_genes.all.ZNF_arrays.tsv.noLOC mammals_PAR_genes.all.ZNF_arrays.tsv
Rscript plot_combined.mammals.R
Rscript plot_combined.allgenelabels.mammals.GeneArrays.R

```
# set a specific gene as the origin
```

awk -v GENE="XG" -f orient_to_gene_midpoint.awk \
  mammals_PAR_genes.all.ZNF_arrays.tsv \
  mammals_PAR_genes.all.ZNF_arrays.tsv \
  > mammals_PAR_genes.XG_oriented.tsv



Rscript plot_combined.mammals.geneoriented.R 




# Compute LOC freq

```
library(dplyr)
library(readr)

# Input file from your AWK script
infile <- "mammals_PAR_genes.all.tsv"

df <- read_tsv(infile, show_col_types = FALSE) %>%
  filter(Species != "Homo_sapiens")
  

df2 <- df %>%
  mutate(
    Region = if_else(InPAR == "N", "nonPAR", "PAR"),
    Is_LOC = grepl("^LOC", GeneName)  )

loc_freq <- df2 %>%
  group_by(Region) %>%
  summarise(
    TotalGenes = n(),
    LOCGenes = sum(Is_LOC, na.rm = TRUE),
    NonLOCGenes = sum(!Is_LOC, na.rm = TRUE),
    LOCFrequency = LOCGenes / TotalGenes,
    LOCPercent = 100 * LOCFrequency,
    .groups = "drop"
  )

print(loc_freq, width = Inf)

# In all mammals except humans, PAR%LOC is 3.88%; nonPAR%LOC is 0.446%

df <- read_tsv(infile, show_col_types = FALSE) %>%
  filter(Species == "Homo_sapiens")
  

df2 <- df %>%
  mutate(
    Region = if_else(InPAR == "N", "nonPAR", "PAR"),
    Is_LOC = grepl("^LOC", GeneName)  )

loc_freq <- df2 %>%
  group_by(Region) %>%
  summarise(
    TotalGenes = n(),
    LOCGenes = sum(Is_LOC, na.rm = TRUE),
    NonLOCGenes = sum(!Is_LOC, na.rm = TRUE),
    LOCFrequency = LOCGenes / TotalGenes,
    LOCPercent = 100 * LOCFrequency,
    .groups = "drop"
  )

print(loc_freq, width = Inf)

# In just humans, PAR%LOC is 23.3%; nonPAR%LOC is 17.7%
















##########################################################################################
# Previous code
##########################################################################################

# Filter X genes to just genes in the PAR and the adjacent non-PAR genes
```
#!/usr/bin/env Rscript

library(dplyr)
library(readr)

infile <- "x_par_genes.with_inferred_names.tsv"
outfile <- "x_par_genes.with_inferred_names.PAR_plus_20_adjacent.tsv"

df <- read_tsv(infile, show_col_types = FALSE)

filtered <- df %>%
  group_by(Species) %>%
  arrange(StartPos, StopPos, .by_group = TRUE) %>%
  mutate(
    gene_index = row_number(),
    par_index = if_else(PAR_status == "Y", gene_index, NA_integer_)
  ) %>%
  mutate(
    first_par_index = min(par_index, na.rm = TRUE),
    last_par_index  = max(par_index, na.rm = TRUE),
    keep_gene = gene_index >= first_par_index - 20 &
                gene_index <= last_par_index + 20,
    RegionWindowStatus = case_when(
      PAR_status == "Y" ~ "PAR",
      gene_index < first_par_index ~ "upstream_20_adjacent",
      gene_index > last_par_index ~ "downstream_20_adjacent",
      TRUE ~ "within_PAR_window"
    )
  ) %>%
  filter(keep_gene) %>%
  ungroup() %>%
  select(-gene_index, -par_index, -first_par_index, -last_par_index, -keep_gene)

write_tsv(filtered, outfile)

print(filtered, n = Inf, width = Inf)
```














### Note if telomere is ID'd on the end of the PAR
```
telomere_file="species.telomeres.txt"

cat > "$telomere_file" <<'EOF'
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
## 3. Plot the data
### Plot a PCA of PAR status by gene and species
```
module load R

R

library(tidyverse)
library(ggrepel)

infile <- "gene_locations_by_species.with_chr_label.all.In_PAR.X_only.gapless_species.csv"

df <- read_csv(infile, show_col_types = FALSE)

# Collapse multiple Z hits per Species/Gene.
# Priority:
#   Y/Edge > N
z_states <- df %>%
  mutate(
    state = case_when(
      In_PAR %in% c("Y", "Edge") ~ 2,
      In_PAR == "N" ~ 1,
      In_PAR == "U" ~ 0,
      is.na(In_PAR) ~ 0,
      TRUE ~ 0
    )
  ) %>%
  group_by(Species, Gene) %>%
  summarise(
    state = max(state),
    .groups = "drop"
  )

# Make the full Species x Gene matrix.
# Missing Species/Gene combinations are interpreted as not found on Z.
mat <- z_states %>%
  pivot_wider(
    names_from = Gene,
    values_from = state,
    values_fill = 0
  )

species <- mat$Species

pca_mat <- mat %>%
  select(-Species) %>%
  as.data.frame()

rownames(pca_mat) <- species

# Remove genes with no variation across species
pca_mat <- pca_mat[, apply(pca_mat, 2, var) > 0]

pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("Species")

var_explained <- summary(pca)$importance[2, ] * 100

ggplot(pca_df, aes(PC1, PC2, label = Species)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.1,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_bw()

ggsave("In_PAR_PCA_from_X_only.pdf", width = 7, height = 6)

write_csv(
  mat,
  "species_by_gene.In_PAR_numeric_matrix.from_X_only.csv"
)
```
### Plot a PCA of PAR status by gene
```
library(ggrepel)

# mat is Species x Gene:
# first column = Species
# remaining columns = numeric states 0/1/2

gene_pca_mat <- mat %>%
  column_to_rownames("Species") %>%
  t() %>%
  as.data.frame()

# Remove species columns with no variation across genes
gene_pca_mat <- gene_pca_mat[, apply(gene_pca_mat, 2, var) > 0]

# Remove genes with no variation across species after transpose, if any
gene_pca_mat <- gene_pca_mat[apply(gene_pca_mat, 1, var) > 0, ]

gene_pca <- prcomp(gene_pca_mat, center = TRUE, scale. = TRUE)

gene_pca_df <- as.data.frame(gene_pca$x) %>%
  rownames_to_column("Gene")

gene_var_explained <- summary(gene_pca)$importance[2, ] * 100

ggplot(gene_pca_df, aes(PC1, PC2)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  xlab(paste0("PC1 (", round(gene_var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(gene_var_explained[2], 1), "%)")) +
  theme_bw()

ggsave("In_PAR_PCA_genes_labeled.pdf", width = 10, height = 8)

write_csv(gene_pca_df, "In_PAR_PCA_gene_coordinates.csv")
```
# Make combined plot
```
Rscript ../PAR_Combined_Plot.R mammals
```
# Permissions modification
chmod -R g+rwx /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR*
chmod -R g+rwx /data/Wilson_Lab/data/VGP_genomes_phase1/
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles