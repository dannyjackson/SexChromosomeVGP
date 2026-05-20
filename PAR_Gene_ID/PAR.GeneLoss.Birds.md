# Identify all genes found in any avian PAR, curate fastas for blast analysis, then identify genes within PARs of all genomes
## 0. Quantify gappiness of each PAR
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/birds


cat > PAR.species_chr_region.txt <<'EOF'
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
EOF

chmod +x find_PAR_gaps.sh 
./find_PAR_gaps.sh

Rscript plot_PAR_gaps.R

# Colius striatus has a PAR assembled in the middle of the W. Check if there are gaps around it.
cat > PAR.species_chr_region.txt <<'EOF'
Colius_striatus,NC_084790.1:87870922-88807031
EOF

## No gaps within it... check if there are gaps around it.
cat > PAR.species_chr_region.txt <<'EOF'
Colius_striatus,NC_084790.1:87000922-89007031
EOF
# Gap found 793,927 bp before the first PAR gene: Colius_striatus NC_084790.1     87075675        87076995        1320    87000922        89007031        74753   7607

chmod +x find_PAR_gaps.sh 
./find_PAR_gaps.sh

# Passer domesticus may have had the PAR assembled inverted
cat > PAR.species_chr_region.txt <<'EOF'
Passer_domesticus,NC_087512.1:0-5000073
EOF
./find_PAR_gaps.sh

# First gap after the PAR is just shy of 2.7 Mb
Passer_domesticus       NC_087512.1     2695374 2695574 200     0       5000073 2695374 2695574

# Identify nearest gap outside of PARs
awk -F'[,\t :.-]+' '
NR==FNR {
    chr[$1]=$2
    start[$1]=$3
    end[$1]=$4
    edge[$1]=($3==0 ? $4 : $3)
    side[$1]=($3==0 ? "right" : "left")
    next
}

FNR==1 {
    sp=FILENAME
    sub(".*/","",sp)
    sub(/\.all_N\.bed$/,"",sp)
}

sp in chr && $1==chr[sp] {
    if (side[sp]=="right" && $2 >= edge[sp]) {
        d = $2 - edge[sp]
    } else if (side[sp]=="left" && $3 <= edge[sp]) {
        d = edge[sp] - $3
    } else {
        next
    }

    if (!(sp in best) || d < best[sp]) {
        best[sp]=d
        gap_s[sp]=$2
        gap_e[sp]=$3
    }
}

END {
    print "species","chr","PAR_edge","gap_start","gap_end","distance_bp"
    for (sp in chr) {
        if (sp in best)
            print sp, chr[sp], edge[sp], gap_s[sp], gap_e[sp], best[sp]
        else
            print sp, chr[sp], edge[sp], "NA", "NA", "NA"
    }
}
' OFS='\t' PAR.species_chr_region.txt PAR_gap_results/per_species/*.all_N.bed \
> PAR.nearest_gap_distance.tsv


# Do our gapless species have a PAR-NAR junction gap?

# <1Mb from PAR edge
Mergus_octosetaceus   0
Coloeus_monedula      51046
Cyanocitta_cristata   90937
Lathamus_discolor     136241

# >1Mb
Morphnus_guianensis   1332323
Passer_domesticus     2028301
Strix_aluco           2145690
Rissa_tridactyla      2616784
Larus_argentatus      3087714
Taeniopygia_guttata   7349429
Patagioenas_fasciata  64699483

# What genes are immediately proximate to the PAR in each species?
## gapless
Mergus_octosetaceus 0-2607306           dymeclin-like,CTIF,SMAD7,ZBTB7C,SMAD2,SKOR2,IER3IP1,HDHD2
Coloeus_monedula    84644553-85378027   TCF4,SMAD4,ELAC1,ME2,MAPK4,SKA1,TETRASPANIN-36-LIKE,CFAP53
Cyanocitta_cristata,0-807097            TCF4,MEX3C,SMAD4,ELAC1,ME2,MAPK4,SKA1,TETRASPANIN-36-LIKE,CFAP53
Lathamus_discolor,111967751-112479095   E3 ubiquitin-protein ligase NEDD4-like,phospholipid-transporting ATPase IC-like,
Morphnus_guianensis,93319090-102737265
Passer_domesticus,0-667073
Strix_aluco,94105843-97043613
Rissa_tridactyla,83212817-88208671
Larus_argentatus,0-4531783
Taeniopygia_guttata,0-2824462
Patagioenas_fasciata,83825531-85896199

## not gapless
Aegotheles_albertisi,0-3571166
Anas_platyrhynchos,0-1777881
Aythya_ferina,0-1863434
Aythya_marila,0-2273918
Calonectris_borealis,77280424-87715099
Colius_striatus,0-1337305
Columba_livia,0-1995161
Cygnus_columbianus,0-2524550
Dixiphia_pipra,0-503193
Falco_naumanni,0-813326
Heliangelus_exortis,0-854851
Numenius_arquata,0-3300171
Opisthocomus_hoazin,88690378-91338688
Phaethon_aethereus,0-10085609
Sarcoramphus_papa,87011595-96420541
Platalea_leucorodia,0-3253226
Poecile_atricapillus,121777972-146584261
Struthio_camelus_australis,35190081-89934196

```
## 1. Curate fastas of genes found in any avian PAR
```

species,PAR,Z_size,W_size
Aegotheles_albertisi,0-3571166,84215446,24988134

chmod +x extract_z_par_genes.sh
./extract_z_par_genes.sh species_par.csv > z_par_genes.tsv
Colius_striatus
Morphnus_guianensis
Platalea_leucorodia
Poecile_atricapillus

# Genes from z_par_genes.tsv found in 1+ species:
cat > genes.txt <<'EOF'
ACAA2
ALPK2
ARK2C
ARK2N
ATP8A1
C18orf32
CCDC68
CFAP53
CPLX4
CTIF
C18orf54
DCC
DYM
EEF2
ELAC1
EPG5
FECH
GPN1
GRP
HAUS1
HDHD2
IER3IP1
KATNAL2
LAS2
LIPG
LMAN1
LOXHD1
LUZP2
MALT1
MAPK4
ME2
MECP2
MEX3C
MIR122
MYO5B
NA
NARS1
NEDD4
ONECUT2
PIAS2
PIK3C3
POLI
PSTPIP2
RAB27B
RAX
RIT2
Rx2
SEC11C
SETBP1
SIGLEC15
SKA1
SKOR2
SLC14A2
SMAD2
SMAD4
SMAD7
SNORD58
ST8SIA3
STARD6
SYT4
TCF4
TSPAN36
TXNL1
WDR7
ZBTB7C
ZNF532
EOF

GFF="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Taeniopygia_guttata/Taeniopygia_guttata.gff"

awk -F'\t' '
    NR==FNR {
        genes[$1]
        order[++n] = $1
        next
    }

    $0 ~ /^#/ { next }

    $3 == "gene" {
        gene_name = ""

        if (match($9, /gene=([^;]+)/, m)) {
            gene_name = m[1]
        } else if (match($9, /Name=([^;]+)/, m)) {
            gene_name = m[1]
        } else if (match($9, /gene_name=([^;]+)/, m)) {
            gene_name = m[1]
        }

        if (gene_name in genes) {
            region = gene_name ";" $1 ":" $4 "-" $5

            # Always prefer NC_133063.1 if present
            if (!(gene_name in best) || $1 == "NC_133063.1") {
                best[gene_name] = region
                found[gene_name] = 1
            }
        }
    }

    END {
        for (i = 1; i <= n; i++) {
            gene = order[i]

            if (gene in found) {
                print best[gene]
            } else {
                print gene ";NOT_FOUND"
            }
        }
    }
' genes.txt "$GFF" >> genes.regions.Taeniopygia_guttata.txt

# repeat using Colius for any "NOT FOUND" genes

cat > genes.Colius_striatus.txt <<'EOF'
ALPK2
C18orf32
CCDC68
CZH18orf54
DCC
FECH
LAS2
LMAN1
MALT1
MEX3C
NA
ONECUT2
POLI
RAB27B
Rx2
SLC14A2
SNORD58
ST8SIA3
STARD6
TCF4
TSPAN36
TXNL1
WDR7
ZNF532
EOF

GFF="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Colius_striatus/Colius_striatus.gff"

awk -F'\t' '
    NR==FNR {
        genes[$1]
        order[++n] = $1
        next
    }

    $0 ~ /^#/ { next }

    $3 == "gene" {
        gene_name = ""

        if (match($9, /gene=([^;]+)/, m)) {
            gene_name = m[1]
        } else if (match($9, /Name=([^;]+)/, m)) {
            gene_name = m[1]
        } else if (match($9, /gene_name=([^;]+)/, m)) {
            gene_name = m[1]
        }

        if (gene_name in genes) {
            region = gene_name ";" $1 ":" $4 "-" $5

            # Always prefer NC_133063.1 if present
            if (!(gene_name in best) || $1 == "NC_133063.1") {
                best[gene_name] = region
                found[gene_name] = 1
            }
        }
    }

    END {
        for (i = 1; i <= n; i++) {
            gene = order[i]

            if (gene in found) {
                print best[gene]
            } else {
                print gene ";NOT_FOUND"
            }
        }
    }
' genes.Colius_striatus.txt "$GFF" >> genes.regions.Colius_striatus.txt


C18orf32;NOT_FOUND Sarcoramphus_papa WM294_016952
CZH18orf54;NOT_FOUND Patagioenas_fasciata
LAS2;NOT_FOUND Sarcoramphus_papa WM294_016714
NA;NOT_FOUND 
Rx2;NOT_FOUND Sarcoramphus_papa WM294_016892
SLC14A2;NOT_FOUND WM294_016706 Sarcoramphus_papa
TSPAN36;NOT_FOUND  WM294_016439	tetraspanin-36-like Sarcoramphus_papa

# Sarcoramphus_papa

cat > genes.Sarcoramphus_papa.txt <<'EOF'
WM294_016952
WM294_016714
WM294_016892
WM294_016706
WM294_016439
EOF


GFF="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Sarcoramphus_papa/Sarcoramphus_papa.gff"

awk -F'\t' '
    NR==FNR {
        genes[$1]
        order[++n] = $1
        next
    }

    $0 ~ /^#/ { next }

    $3 == "gene" {
        gene_name = ""

        if (match($9, /gene=([^;]+)/, m)) {
            gene_name = m[1]
        } else if (match($9, /Name=([^;]+)/, m)) {
            gene_name = m[1]
        } else if (match($9, /gene_name=([^;]+)/, m)) {
            gene_name = m[1]
        }

        if (gene_name in genes) {
            region = gene_name ";" $1 ":" $4 "-" $5

            # Always prefer NC_133063.1 if present
            if (!(gene_name in best) || $1 == "NC_133063.1") {
                best[gene_name] = region
                found[gene_name] = 1
            }
        }
    }

    END {
        for (i = 1; i <= n; i++) {
            gene = order[i]

            if (gene in found) {
                print best[gene]
            } else {
                print gene ";NOT_FOUND"
            }
        }
    }
' genes.Sarcoramphus_papa.txt "$GFF" >> genes.regions.Sarcoramphus_papa.txt

```
Pull fastas from reference genomes (Taeniopygia_guttata, Colius_striatus, Sarcoramphus_papa)
```
grep -v NOT_FOUND genes.regions.Taeniopygia_guttata.txt > genes.regions.Taeniopygia_guttata.noNA.txt

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Taeniopygia_guttata/Taeniopygia_guttata.fna"

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

# Index reference if needed
if [ ! -f "${REF_GENOME}.fai" ]; then
    samtools faidx "$REF_GENOME"
fi

while IFS=';' read -r gene region; do
    # Skip empty or malformed lines
    if [ -z "${gene}" ] || [ -z "${region}" ]; then
        echo "Skipping malformed line: gene='${gene}' region='${region}'" >&2
        continue
    fi

    out="${gene}.fa"

    echo "Extracting ${gene}: ${region} -> ${out}"

    samtools faidx "$REF_GENOME" "$region" \
        | awk -v gene="$gene" -v region="$region" '
            NR == 1 { print ">" gene "|" region; next }
            { print }
        ' > "$out"

done < genes.regions.Taeniopygia_guttata.noNA.txt




grep -v NOT_FOUND genes.regions.Colius_striatus.txt > genes.regions.Colius_striatus.noNA.txt

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Colius_striatus/Colius_striatus.fna"

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

# Index reference if needed
if [ ! -f "${REF_GENOME}.fai" ]; then
    samtools faidx "$REF_GENOME"
fi

while IFS=';' read -r gene region; do
    # Skip empty or malformed lines
    if [ -z "${gene}" ] || [ -z "${region}" ]; then
        echo "Skipping malformed line: gene='${gene}' region='${region}'" >&2
        continue
    fi

    out="${gene}.fa"

    echo "Extracting ${gene}: ${region} -> ${out}"

    samtools faidx "$REF_GENOME" "$region" \
        | awk -v gene="$gene" -v region="$region" '
            NR == 1 { print ">" gene "|" region; next }
            { print }
        ' > "$out"

done < genes.regions.Colius_striatus.noNA.txt







grep -v NOT_FOUND genes.regions.Sarcoramphus_papa.txt > genes.regions.Sarcoramphus_papa.noNA.txt

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Sarcoramphus_papa/Sarcoramphus_papa.fna"

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

# Index reference if needed
if [ ! -f "${REF_GENOME}.fai" ]; then
    samtools faidx "$REF_GENOME"
fi

while IFS=';' read -r gene region; do
    # Skip empty or malformed lines
    if [ -z "${gene}" ] || [ -z "${region}" ]; then
        echo "Skipping malformed line: gene='${gene}' region='${region}'" >&2
        continue
    fi

    out="${gene}.fa"

    echo "Extracting ${gene}: ${region} -> ${out}"

    samtools faidx "$REF_GENOME" "$region" \
        | awk -v gene="$gene" -v region="$region" '
            NR == 1 { print ">" gene "|" region; next }
            { print }
        ' > "$out"

done < genes.regions.Sarcoramphus_papa.noNA.txt



mv WM294_016952.fa C18orf32.fa
mv WM294_016714.fa LAS2.fa
mv WM294_016892.fa Rx2.fa
mv WM294_016706.fa SLC14A2.fa
mv WM294_016439.fa TSPAN36.fa
```
## 2. Blast all avian genomes for the fastas of putatively conserved PAR genes
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/birds

cat > species_for_blast_array.txt <<'EOF'
Mergus_octosetaceus
Aegotheles_albertisi
Sarcoramphus_papa
Larus_argentatus
Rissa_tridactyla
Numenius_arquata
Columba_livia
Patagioenas_fasciata
Rhynochetos_jubatus
Falco_naumanni
Anas_platyrhynchos
Aythya_ferina
Aythya_marila
Cygnus_columbianus
Opisthocomus_hoazin
Coloeus_monedula
Cyanocitta_cristata
Taeniopygia_guttata
Passer_domesticus
Dixiphia_pipra
Phaethon_aethereus
Calonectris_borealis
Strix_aluco
Struthio_camelus_australis
Heliangelus_exortis
Morphnus_guianensis
Colius_striatus
Poecile_atricapillus
Zosterops_lateralis
Platalea_leucorodia
Lathamus_discolor
EOF

N=$(wc -l < species_for_blast_array.txt)
sbatch --array=1-${N} blast_gene_locations_array.sh
```
### Filter blast outputs
```
mkdir -p blast_results_Zchr

for f in blast_results/*.blast.tsv; do
    base=$(basename "$f")
    species=$(echo "$base" | cut -d'.' -f1)
    out="blast_results_Zchr/$base"

    count=$(
        awk -v species="$species" '
            BEGIN {
                FS = "[,\t ]+"
            }

            NR == FNR {
                if ($1 == species && $2 == "Z") {
                    xacc[$3] = 1
                }
                next
            }

            {
                subj = $2

                # Handles:
                #   ref|NC_060947.1|
                #   emb|OZ239531.1|
                #   gb|ABC123.1|
                #   dbj|XYZ123.1|
                #   NC_060947.1
                if (subj ~ /^[A-Za-z_][A-Za-z0-9_]*\|[^|]+\|?$/) {
                    split(subj, a, "|")
                    acc = a[2]
                } else {
                    acc = subj
                }

                if (acc in xacc) {
                    print
                }
            }
        ' "$SEXCHR_FILE" "$f" | tee "$out" | wc -l
    )

    printf "%s\t%s\n" "$base" "$count"
done
```
### Combine blast outputs
```
chmod +x Blast_Match_inPAR.birds.sh
./Blast_Match_inPAR.birds.sh

```
### Add chromosome label to output
For each blast match, this column will have a value of Z, W, or NA. A value of NA here indicates that the chromosome with the PAR genes is not a sex chromosome. Some manual curation of this column will be necessary because I made the sexchrom_accessions file a while ago.

```
SEXCHR_FILE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"

LABELED_OUT="gene_locations_by_species.with_chr_label.all.csv"


awk -F',' '
    BEGIN {
        OFS = ","
    }

    # First file: sex chromosome accession map
    NR == FNR {
        if (FNR == 1) next

        species = $1
        chr_label = $2
        accession = $3

        gsub(/\r/, "", accession)
        gsub(/\r/, "", chr_label)

        key = species SUBSEP accession
        label[key] = chr_label

        next
    }

    # Second file: BLAST output
    FNR == 1 {
        print $0, "chr_label"
        next
    }

    {
        species = $1
        chrom = $4

        gsub(/\r/, "", chrom)

        key = species SUBSEP chrom

        if (key in label) {
            chr_label = label[key]
        } else {
            chr_label = "NA"
        }

        print $0, chr_label
    }
' "$SEXCHR_FILE" "$OUT_ALL" > "$LABELED_OUT"

echo "Labeled results written to ${LABELED_OUT}"
```
### Infer gene overlap
```
# Add neo-sex species to species PAR:
cp ../PAR_Gene_ID/birds/species_par.csv .

echo 'Morphnus_guianensis,93319090-102737265' >> species_par.csv
echo 'Platalea_leucorodia,0-3253226' >> species_par.csv
echo 'Poecile_atricapillus,121777972-146584261' >> species_par.csv
echo 'Struthio_camelus_australis,35190081-89934196' >> species_par.csv
```
### 
```
awk -F',' '
BEGIN {
    OFS = FS
}

# First file: species_par.csv
# Format: Species,start-stop
# These coordinates refer to the Z chromosome for that species.
NR == FNR {
    split($2, b, "-")
    rstart = b[1] + 0
    rstop  = b[2] + 0

    species = $1
    n[species]++
    start[species, n[species]] = rstart
    stop[species, n[species]]  = rstop

    next
}

# Second file: gene_locations_by_species.with_chr_label.all.csv
FNR == 1 {
    print $0, "In_PAR"
    next
}

{
    species = $1
    chr_label = $7
    gstart = $5 + 0
    gstop = $6 + 0

    # Only evaluate Z chromosomes.
    # Anything else, including W or NA, is unknown for this PAR test.
    if (chr_label != "Z") {
        status = "U"
    } else if (!(species in n)) {
        status = "U"
    } else {
        status = "N"

        for (i = 1; i <= n[species]; i++) {
            # Any overlap between gene interval and species Z PAR interval
            if (gstop >= start[species, i] && gstart <= stop[species, i]) {
                # Fully contained within the PAR interval
                if (gstart >= start[species, i] && gstop <= stop[species, i]) {
                    status = "Y"
                } else {
                    status = "Edge"
                }
                break
            }
        }
    }

    print $0, status
}
' species_par.csv \
  gene_locations_by_species.with_chr_label.all.csv \
  > gene_locations_by_species.with_chr_label.all.In_PAR.csv

# Filter to just Z matches
awk -F',' 'BEGIN { OFS = FS } NR == 1 || $(NF-1) == "Z"' \
  gene_locations_by_species.with_chr_label.all.In_PAR.csv \
  > gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.csv
```

# Filter to gapless species (plus the ostrich), also drop Colius_striatus ("PAR" appears to be assembled in the middle of the W... could be a loss? Could also be an assembly error -- a )
```

input="gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.csv"
output="gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.gapless_species.csv"
species_list="gapless_species.txt"

cat > "$species_list" <<'EOF'
Lathamus_discolor
Passer_domesticus
Coloeus_monedula
Cyanocitta_cristata
Patagioenas_fasciata
Mergus_octosetaceus
Taeniopygia_guttata
Strix_aluco
Larus_argentatus
Rissa_tridactyla
Morphnus_guianensis
EOF

# Keep header, then rows matching one of the gapless species.
{
    head -n 1 "$input"
    tail -n +2 "$input" | grep -Ff "$species_list"
} > "$output"

echo "Wrote: $output"
echo "Species retained:"
cut -d',' -f1 "$output" | tail -n +2 | sort -u
```
### Note if telomere is ID'd on the end of the PAR
```
telomere_file="gapless_species.telomeres.txt"

cat > "$telomere_file" <<'EOF'
Lathamus_discolor,R,NO
Passer_domesticus,L,NO
Coloeus_monedula,R,YES
Cyanocitta_cristata,L,YES
Patagioenas_fasciata,R,YES
Mergus_octosetaceus,L,YES
Taeniopygia_guttata,L,YES
Strix_aluco,R,YES
Larus_argentatus,L,YES
Rissa_tridactyla,R,YES
Morphnus_guianensis,R,NO
EOF
```

## 3. Plot the data
### Plot a PCA of PAR status by gene and species
```
module load R

R

library(tidyverse)
library(ggrepel)

infile <- "gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.gapless_species.csv"

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

ggsave("In_PAR_PCA_from_Z_only.pdf", width = 7, height = 6)

write_csv(
  mat,
  "species_by_gene.In_PAR_numeric_matrix.from_Z_only.csv"
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
Rscript ../PAR_Combined_Plot.R birds
Rscript ../PAR_nonPAR_Combined_Plot.R birds
Rscript ../PAR_nonPAR_Combined_Plot.minimaltext.R birds

Mergus_octosetaceus   0
Coloeus_monedula      51046
Cyanocitta_cristata   90937
Lathamus_discolor     136241

```
### Repeat with the ostrich included as an "outgroup" species
```
input="gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.csv"
output="gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.gapless_species.csv"
species_list="gapless_species.txt"

cat > "$species_list" <<'EOF'
Lathamus_discolor
Passer_domesticus
Coloeus_monedula
Cyanocitta_cristata
Patagioenas_fasciata
Mergus_octosetaceus
Opisthocomus_hoazin
Taeniopygia_guttata
Strix_aluco
Numenius_arquata
Larus_argentatus
Rissa_tridactyla
Morphnus_guianensis
Calonectris_borealis
Struthio_camelus_australis
EOF

# Keep header, then rows matching one of the gapless species.
{
    head -n 1 "$input"
    tail -n +2 "$input" | grep -Ff "$species_list"
} > "$output"

echo "Wrote: $output"
echo "Species retained:"
cut -d',' -f1 "$output" | tail -n +2 | sort -u
```
## 3. Plot the data
### Plot a PCA of PAR status by gene and species
```
module load R

R

library(tidyverse)
library(ggrepel)

infile <- "gene_locations_by_species.with_chr_label.all.In_PAR.Z_only.gapless_species.csv"

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

ggsave("In_PAR_PCA_from_Z_only.ostrich.pdf", width = 7, height = 6)

write_csv(
  mat,
  "species_by_gene.In_PAR_numeric_matrix.from_Z_only.ostrich.csv"
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

ggsave("In_PAR_PCA_genes_labeled.ostrich.pdf", width = 10, height = 8)

write_csv(gene_pca_df, "In_PAR_PCA_gene_coordinates.ostrich.csv")
```
# Make combined plot
```
Rscript PAR_Combined_Plot.birds.ostrich.R
```
# Permissions modification
chmod -R g+rwx /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR*
chmod -R g+rwx /data/Wilson_Lab/data/VGP_genomes_phase1/
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles