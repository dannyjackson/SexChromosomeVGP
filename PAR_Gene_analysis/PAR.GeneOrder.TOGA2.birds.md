# Identify all genes found in any avian PAR, curate fastas for blast analysis, then identify genes within PARs of all genomes
## 0. Quantify gappiness of each PAR
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/TOGA2_annotations/birds

TOGA_DIR="/data/Wilson_Lab/data/TOGA2_Hiller/zebrafinch_HLtaeGut5-GCF_003957565.2"


# Confirm that accessions in sexchrfile still match correctly (genomes have been updated)
SEXCHR_FILE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"
BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/"

SPECIES=Rhynochetos_jubatus
grep $SPECIES $SEXCHR_FILE
grep '>' ${BASE}/${SPECIES}/${SPECIES}.fna | grep 'chromosome'



# Revise those that need revision

grep 'Opisthocomus_hoazin' $SEXCHR_FILE
sed -i 's/CM061422/NC_134454/g' $SEXCHR_FILE
sed -i 's/CM061421/NC_134453/g' $SEXCHR_FILE
grep 'Opisthocomus_hoazin' $SEXCHR_FILE

# Excluding this species, which had no surviving regions after filtering for PID and length: Rhynochetos_jubatus

cat > PAR.species_chr_region.gaps.txt <<'EOF'
Aegotheles_albertisi,CM078494.1:0-3571166
Anas_platyrhynchos,OZ076978.1:0-1777881
Aythya_ferina,OZ124217.1:0-1863434
Aythya_marila,OZ223658.1:0-2273918
Calonectris_borealis,NC_134352.1:77280424-87715099
Colius_striatus,NC_084790.1:0-1337305
Coloeus_monedula,OZ238506.1:84644553-85378027
Columba_livia,NC_088642.1:0-1995161
Cyanocitta_cristata,CM100569.1:0-807097
Cygnus_columbianus,OZ223797.1:0-2524550
Dixiphia_pipra,NC_087581.1:0-503193
Falco_naumanni,NC_054080.1:0-813326
Heliangelus_exortis,NC_092454.1:0-854851
Larus_argentatus,OZ207420.1:0-4531783
Lathamus_discolor,NC_088909.1:111967751-112479095
Mergus_octosetaceus,CM072318.1:0-2607306
Morphnus_guianensis,CM098430.1:93319090-102737265
Numenius_arquata,NC_133616.1:0-3300171
Opisthocomus_hoazin,NC_134454.1:88690378-91338688
Passer_domesticus,NC_087512.1:0-667073
Patagioenas_fasciata,NC_092560.1:83825531-85896199
Platalea_leucorodia,OZ238966.1:0-3253226
Poecile_atricapillus,NC_081289.1:121777972-146584261
Phaethon_aethereus,OZ196914.1:0-10085609
Rissa_tridactyla,NC_071497.1:83212817-88208671
Sarcoramphus_papa,CM075626.1:87011595-96420541
Strix_aluco,NC_133971.1:94105843-97043613
Struthio_camelus_australis,NC_090982.1:35190081-89934196
Taeniopygia_guttata,NC_133063.1:0-2824462
Zosterops_lateralis,OZ246513.1:0-517530
EOF


# chmod +x find_PAR_gaps.sh 
# ./find_PAR_gaps.sh

# Rscript plot_PAR_gaps.R

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
```
## 1. Curate annotated genes found in the PARs

Calonectris_borealis # 	swap to OZ122106.1 from NC_134352.1
Numenius_arquata # swap to OZ067243.1 from NC_133616.1
Opisthocomus_hoazin # swap to CM061422.1 from NC_134454.1
Strix_aluco # swap to CM062916.1 from NC_133971.1
```




# Excluding Poecile_atricapillus: no clear PAR boundary, 2MB gap between first and second contiguous aligned seq
# Excluding Colius_striatus,NC_084790.1:0-1337305: inferred PAR is neo-sex
# Excluding Falco_naumanni,NC_054080.1:0-813326: Can't confirm PAR side... dig more into later


cat > PAR.species_chr_region.txt <<'EOF'
Aegotheles_albertisi,CM078494.1:0-3571166
Anas_platyrhynchos,OZ076978.1:0-1777881
Aythya_ferina,OZ124217.1:0-1863434
Aythya_marila,OZ223658.1:0-2273918
Calonectris_borealis,OZ122106.1:77280424-87715099
Coloeus_monedula,OZ238506.1:84644553-85378027
Columba_livia,NC_088642.1:0-1995161
Cyanocitta_cristata,CM100569.1:0-807097
Cygnus_columbianus,OZ223797.1:0-2524550
Dixiphia_pipra,NC_087581.1:0-503193
Heliangelus_exortis,NC_092454.1:0-854851
Larus_argentatus,OZ207420.1:0-4531783
Lathamus_discolor,NC_088909.1:111967751-112479095
Mergus_octosetaceus,CM072318.1:0-2607306
Morphnus_guianensis,CM098430.1:93319090-102737265
Numenius_arquata,OZ067243.1:0-3300171
Opisthocomus_hoazin,CM061422.1:88690378-91338688
Passer_domesticus,NC_087512.1:0-667073
Patagioenas_fasciata,NC_092560.1:83825531-85896199
Platalea_leucorodia,OZ238966.1:0-3253226
Phaethon_aethereus,OZ196914.1:0-10085609
Rissa_tridactyla,NC_071497.1:83212817-88208671
Sarcoramphus_papa,CM075626.1:87011595-96420541
Strix_aluco,CM062916.1:94105843-97043613
Struthio_camelus,NC_090982.1:35190081-89934196
Taeniopygia_guttata,NC_133063.1:0-2824462
Zosterops_lateralis,OZ246513.1:0-517530
EOF

export TOGA_DIR="/data/Wilson_Lab/data/TOGA2_Hiller/zebrafinch_HLtaeGut5-GCF_003957565.2" 

chmod +x find_par_genes_toga.birds.awk 
./find_par_genes_toga.birds.awk PAR.species_chr_region.txt > birds_PAR_genes.tsv



cat > species.txt <<'EOF'
Aegotheles_albertisi
Anas_platyrhynchos
Aythya_ferina
Aythya_marila
Calonectris_borealis
Colius_striatus
Coloeus_monedula
Columba_livia
Cyanocitta_cristata
Cygnus_columbianus
Dixiphia_pipra
Falco_naumanni
Heliangelus_exortis
Larus_argentatus
Lathamus_discolor
Mergus_octosetaceus
Morphnus_guianensis
Numenius_arquata
Opisthocomus_hoazin
Passer_domesticus
Patagioenas_fasciata
Phaethon_aethereus
Platalea_leucorodia
Rissa_tridactyla
Sarcoramphus_papa
Strix_aluco
Struthio_camelus
Taeniopygia_guttata
EOF

while read -r species; do
  echo $species
  grep $species birds_PAR_genes.tsv | wc -l
done < species.txt


## Anas_platyrhynchos -- TOGA2 has different assembly


cp birds_PAR_genes.tsv birds_PAR_genes.all.tsv 

sed -i 's/Struthio_camelus/Struthio_camelus_australis/g' birds_PAR_genes.all.tsv 

```
# Plot
```

telomere_file="species.telomeres.txt"

cat > "$telomere_file" <<'EOF'
Aegotheles_albertisi,L,NO
Anas_platyrhynchos,L,YES
Aythya_ferina,L,YES
Aythya_marila,L,YES
Calonectris_borealis,R,NO
Coloeus_monedula,R,YES
Columba_livia,L,YES
Cyanocitta_cristata,L,YES
Cygnus_columbianus,L,NO
Dixiphia_pipra,L,NO
Heliangelus_exortis,L,NO
Larus_argentatus,L,NO
Lathamus_discolor,R,NO
Mergus_octosetaceus,L,YES
Morphnus_guianensis,R,NO
Numenius_arquata,L,NO
Opisthocomus_hoazin,R,NO
Passer_domesticus,L,NO
Patagioenas_fasciata,R,YES
Phaethon_aethereus,L,NO
Platalea_leucorodia,L,NO
Rissa_tridactyla,R,YES
Sarcoramphus_papa,R,NO
Strix_aluco,R,NO
Struthio_camelus_australis,R,YES
Taeniopygia_guttata,L,YES
Zosterops_lateralis,R,NO
EOF

# Clean data to standardize genes named for X and Y locations in human
## sed -i 's/NLGN4Y/NLGN4/g' mammals_PAR_genes.all.tsv

Rscript plot_combined.allgenelabels.birds.R

Rscript plot_combined.birds.R



# Modify the file to combine major blocks of gene family expansions into a single point
{
  head -n 1 birds_PAR_genes.all.tsv
  tail -n +2 birds_PAR_genes.all.tsv | sort -t $'\t' -k1,1 -k2,2 -k3,3 -k4,4n
} | awk 'BEGIN{FS=OFS="\t"}
NR==1 { print; next }

function print_saved_rows(    i) {
    for (i = 1; i <= n_rows; i++) {
        print rows[i]
    }
}

function flush_block(    f) {
    if (n_rows == 0) return

    if (n_rows >= 3) {
        split(rows[1], f, OFS)
        f[4] = block_start
        f[5] = block_stop
        f[6] = block_name
        print f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8]
    } else {
        print_saved_rows()
    }

    n_rows = 0
    in_block = 0
}

{
    key = $1 OFS $2 OFS $3

    is_array = 0
    name = ""

    if ($6 ~ /^LOC/) {
        is_array = 1
        name = "LOC_array"
    }

    if (is_array) {
        if (!in_block || key != block_key || name != block_name) {
            flush_block()

            n_rows = 1
            rows[n_rows] = $0
            block_start = $4
            block_stop = $5
            block_key = key
            block_name = name
            in_block = 1
        } else {
            n_rows++
            rows[n_rows] = $0

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
}' > birds_PAR_genes.all.LOC_arrays.tsv


Rscript plot_combined.birds.arrays.R

```
```
# set a specific gene as the origin

cp birds_PAR_genes.all.tsv birds_PAR_genes.all.tsv.save

grep -v 'LOC' birds_PAR_genes.all.tsv | grep -v 'LINC' > birds_PAR_genes.all.tsv.noLOC
mv birds_PAR_genes.all.tsv.noLOC birds_PAR_genes.all.tsv


Rscript plot_combined.birds.geneoriented.R 




