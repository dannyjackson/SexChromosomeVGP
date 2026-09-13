# Identify all genes found in any avian PAR, curate fastas for blast analysis, then identify genes within PARs of all genomes
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

```
## 1. Curate fastas of genes found in any avian PAR
```

chmod +x 1a_extract_z_par_genes.sh
./1a_extract_z_par_genes.sh PAR.species_chr_region.txt > z_par_genes.tsv

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
```
## 2. Blast all avian genomes for the fastas of putatively conserved PAR genes
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/birds

awk -F',' '{print $1}' PAR.species_chr_region.txt > species_for_blast_array.txt

N=$(wc -l < species_for_blast_array.txt)
sbatch --array=1-${N} 1b_blast_gene_locations_array.birds.sh
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
chmod +x 1c_Blast_Match_inPAR.birds.sh
./1c_Blast_Match_inPAR.birds.sh

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
### Infer gene overlap with PAR
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