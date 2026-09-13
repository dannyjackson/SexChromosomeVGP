BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

SEXCHR_FILE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"
PAR_FILE="PAR.species_chr_region.txt"

species=(
    Eubalaena_glacialis
    Macaca_nemestrina
    Callithrix_jacchus
    Mesoplodon_bidens
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
)

genes=(
    AKAP17A
    ANOS1
    APOA1A
    ASMT
    ASMTL
    ASMTL-AS1
    BOSD2_1
    BOSD2_2
    BOSD2_3
    BOSD2_4
    CRLF2
    CSF2RA
    DHRSX
    GPR143
    GTPBP6
    GYG2
    HRG
    IL3RA
    LINC00102
    LYL1
    MALRD1
    CD99
    MIR3690
    MXRA5
    NLGN4X
    OBP
    P2RY8
    PLCXD1
    PNPLA4
    PPP2R3B
    PRKX
    PUDP
    SHOX
    SHROOM2
    SLC25A6
    STS
    TBL1X
    XG
    ZBED1
    ZNF665
)

OUT_ALL="gene_locations_by_species.all.csv"
OUT_XREGION="gene_locations_by_species.X_region.csv"

echo "Species,Gene,Hit_rank,Chrom,Start_pos,Stop_pos" > "$OUT_ALL"
echo "Species,Gene,Hit_rank,Chrom,Start_pos,Stop_pos,Bitscore,Aln_len,X_region_start,X_region_stop" > "$OUT_XREGION"

for sp in "${species[@]}"; do
    genome="${BASE}/${sp}/${sp}.fna"

    if [ ! -f "$genome" ]; then
        echo "WARNING: genome not found for ${sp}: ${genome}" >&2
        continue
    fi

    for gene in "${genes[@]}"; do
        query="${gene}.fa"

        if [ ! -f "$query" ]; then
            echo "WARNING: query FASTA not found for ${gene}: ${query}" >&2
            continue
        fi

        blast_out="blast_results_Xchr/${sp}.${gene}.blast.tsv"

        echo "BLAST: ${sp} ${gene}"

        if [ ! -s "$blast_out" ]; then
            echo "${sp},${gene},NA,NA,NA,NA" >> "$OUT_ALL"
            echo "${sp},${gene},NA,NA,NA,NA,NA,NA,NA,NA" >> "$OUT_XREGION"
            continue
        fi

        ############################################################
        # Original output:
        # Pick best hit overall by highest bitscore, then length
        ############################################################

        awk -v sp="$sp" -v gene="$gene" '
            BEGIN {
                FS = OFS = "\t"
            }
            {
                chrom = $2
                gsub(/\r/, "", chrom)

                if (chrom ~ /\|/) {
                    split(chrom, a, "|")
                    chrom = a[2]
                }

                s = $9
                e = $10
                bits = $12
                len = $4

                if (s <= e) {
                    start = s
                    stop = e
                } else {
                    start = e
                    stop = s
                }

                print bits, len, chrom, start, stop
            }
        ' "$blast_out" \
        | sort -k1,1gr -k2,2gr \
        | head -n 1 \
        | awk -v sp="$sp" -v gene="$gene" '
            BEGIN {
                OFS = ","
                rank = 0
            }
            {
                rank++
                print sp, gene, rank, $3, $4, $5
            }
        ' >> "$OUT_ALL"

        ############################################################
        # New output:
        # Pick best hit on that species X chromosome within PAR region
        ############################################################

        best_x_region=$(
    awk \
        -v sp="$sp" \
        -v gene="$gene" \
        -v sexchr_file="$SEXCHR_FILE" \
        -v par_file="$PAR_FILE" '
        BEGIN {
            FS = OFS = "\t"

            while ((getline line < par_file) > 0) {
                gsub(/\r/, "", line)
                gsub(/[[:space:]]+$/, "", line)

                if (line == "" || line ~ /^#/) {
                    continue
                }

                split(line, f, ",")
                species_name = f[1]
                region = f[2]

                gsub(/^[[:space:]]+|[[:space:]]+$/, "", species_name)
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", region)

                split(region, r, "-")

                par_start[species_name] = r[1] + 0
                par_stop[species_name] = r[2] + 0
            }
            close(par_file)

            while ((getline line < sexchr_file) > 0) {
                gsub(/\r/, "", line)

                if (line == "" || line ~ /^Species,Chromosome,Accession/) {
                    continue
                }

                split(line, f, ",")

                species_name = f[1]
                chrom_name = f[2]
                accession = f[3]

                gsub(/^[[:space:]]+|[[:space:]]+$/, "", species_name)
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", chrom_name)
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", accession)

                if (chrom_name == "X") {
                    x_accession[species_name] = accession
                }
            }
            close(sexchr_file)

            if (!(sp in par_start) || !(sp in x_accession)) {
                missing_info = 1
            }
        }

        missing_info == 1 {
            next
        }

        {
            chrom = $2
            gsub(/\r/, "", chrom)

            if (chrom ~ /\|/) {
                split(chrom, a, "|")
                chrom = a[2]
            }

            s = $9
            e = $10
            bits = $12
            len = $4

            if (s <= e) {
                start = s
                stop = e
            } else {
                start = e
                stop = s
            }

            if (chrom == x_accession[sp] && stop >= par_start[sp] && start <= par_stop[sp]) {
                print bits, len, chrom, start, stop, par_start[sp], par_stop[sp]
            }
        }
    ' "$blast_out" \
    | sort -k1,1gr -k2,2gr \
    | head -n 1
)

        if [ -z "$best_x_region" ]; then
            echo "${sp},${gene},NA,NA,NA,NA,NA,NA,NA,NA" >> "$OUT_XREGION"
        else
            echo "$best_x_region" \
            | awk -v sp="$sp" -v gene="$gene" '
                BEGIN {
                    OFS = ","
                }
                {
                    rank = 1
                    bits = $1
                    len = $2
                    chrom = $3
                    start = $4
                    stop = $5
                    region_start = $6
                    region_stop = $7

                    print sp, gene, rank, chrom, start, stop, bits, len, region_start, region_stop
                }
            ' >> "$OUT_XREGION"
        fi

    done
done

echo "Done."
echo "Overall best hits written to ${OUT_ALL}"
echo "Best X-region hits written to ${OUT_XREGION}"