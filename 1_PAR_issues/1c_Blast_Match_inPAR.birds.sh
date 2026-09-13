BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

SEXCHR_FILE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"
PAR_FILE="PAR.species_chr_region.txt"


species=(
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
)

genes=(
    ACAA2
    ALPK2
    ARK2C
    ARK2N
    ATP8A1
    C18orf32
    CCDC68
    CPLX4
    CTIF
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
    NARS1
    NEDD4
    ONECUT2
    CFAP53
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
    ST8SIA3
    STARD6
    SYT4
    TCF4
    TSPAN36
    TXNL1
    WDR7
    ZBTB7C
    ZNF532
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

                if (chrom_name == "Z") {
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