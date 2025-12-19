#!/bin/bash

export SEXCHRLEN=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_phylogeny/VGP_freeze_hap1_combined_sexchroms_seq_reports.tsv

sed -i 's/Guaruba_guaruba/Guaruba_guarouba/g' $SEXCHRLEN

OUTFILE="bed_sexchr_info.tsv"
ERRFILE="bed_sexchr_info_errors.tsv"

# Headers
printf "Accession\tGenBank_accession\tRefSeq_accession\tSeq_length\tSpecies\tFilename\n" > "$OUTFILE"
printf "Filename\tSpecies_key\tChromosome\tError\tDetail\n" > "$ERRFILE"

for bed in sorted_beds/*; do
    [ -f "$bed" ] || continue
    base=$(basename "$bed")

    # Strip at first dot: Genus_species[_subspecies]_comparison...
    core=${base%%.*}

    # Split core into tokens
    IFS='_' read -r genus species third fourth rest <<< "$core"
    unset IFS

    if [[ -z "$genus" || -z "$species" ]]; then
        printf "%s\t\t\tMalformed filename (no genus/species)\t\n" "$base" >> "$ERRFILE"
        continue
    fi

    has_sub=0
    subspecies=""
    comparison_code=""

    if [[ -n "$third" ]]; then
        firstchar=${third:0:1}
        if [[ "$firstchar" =~ [a-z] ]]; then
            # third token is subspecies
            subspecies="$third"
            has_sub=1
            species_key="${genus}_${species}_${subspecies}"
            comparison_code="$fourth"
        else
            # third token is comparison code
            has_sub=0
            species_key="${genus}_${species}"
            comparison_code="$third"
        fi
    else
        has_sub=0
        species_key="${genus}_${species}"
        comparison_code=""
    fi

    # Determine chromosome name from comparison code
    chr=""

    if [[ "$comparison_code" == *"WtoZ"* ]]; then
        chr="Z"
    elif [[ "$comparison_code" == *"X1"* ]]; then
        chr="X1"
    elif [[ "$comparison_code" == *"X2"* ]]; then
        chr="X2"
    elif [[ "$comparison_code" == *"Y1toX"* ]]; then
        # treat Y1toX as standard X
        chr="X"
    elif [[ "$comparison_code" == *"Y2toX"* ]]; then
        # treat Y2toX as standard X
        chr="X"
    elif [[ "$comparison_code" == *"YtoX"* ]]; then
        chr="X"
    fi

    if [[ -z "$chr" ]]; then
        printf "%s\t%s\t\tUnable to determine chromosome from comparison code '%s'\t\n" \
            "$base" "$species_key" "$comparison_code" >> "$ERRFILE"
        continue
    fi

    # Search SEXCHRLEN
    match=$(awk -v FS="\t" -v OFS="\t" \
        -v spkey="$species_key" -v has_sub="$has_sub" -v ch="$chr" -v fname="$base" '
        NR == 1 { next }

        {
            full = $16  # Species field
            n = split(full, a, "_")

            # Genus_species key
            if (n >= 2)
                key2 = a[1] "_" a[2]
            else
                key2 = full

            # Genus_species_subspecies key (if present)
            if (n >= 3)
                key3 = key2 "_" a[3]
            else
                key3 = key2

            # Match logic:
            # - If filename has subspecies -> match on key3
            # - If filename has no subspecies -> match on key2 (any subspecies)
            if ( (has_sub == 1 && key3 == spkey && $4 == ch) ||
                 (has_sub == 0 && key2 == spkey && $4 == ch) ) {
                count++
                acc[count]  = $1   # Accession
                gb[count]   = $7   # GenBank_accession
                ref[count]  = $10  # RefSeq_accession
                len[count]  = $12  # Seq_length
                spec[count] = full # full species (with subspecies)
            }
        }

        END {
            if (count == 0) {
                print "NOTFOUND"
                exit
            }

            if (count == 1) {
                print acc[1], gb[1], ref[1], len[1], spec[1], fname
                print "FOUND_ONE"
                exit
            }

            # Multiple rows matched: check distinct full species (subspecies)
            uniq = 0
            for (i = 1; i <= count; i++) {
                exists = 0
                for (j = 1; j <= uniq; j++) {
                    if (spec[i] == uspec[j]) { exists = 1; break }
                }
                if (!exists) {
                    uniq++
                    uspec[uniq] = spec[i]
                }
            }

            if (uniq > 1) {
                s = uspec[1]
                for (k = 2; k <= uniq; k++) s = s ";" uspec[k]
                print "AMBIG\t" s
            } else {
                # multiple rows but same subspecies: OK, take first
                print acc[1], gb[1], ref[1], len[1], spec[1], fname
                print "FOUND_MULTI_SAME"
            }
        }
    ' "$SEXCHRLEN")

    status=$(echo "$match" | tail -n 1 | cut -f1)

    if [[ "$status" == "FOUND_ONE" || "$status" == "FOUND_MULTI_SAME" ]]; then
        echo "$match" | head -n -1 >> "$OUTFILE"
    elif [[ "$status" == "NOTFOUND" ]]; then
        printf "%s\t%s\t%s\tNo matching entry in SEXCHRLEN\t\n" \
            "$base" "$species_key" "$chr" >> "$ERRFILE"
    elif [[ "$status" == "AMBIG" ]]; then
        detail=$(echo "$match" | tail -n 1 | cut -f2-)
        printf "%s\t%s\t%s\tMultiple subspecies matched\t%s\n" \
            "$base" "$species_key" "$chr" "$detail" >> "$ERRFILE"
    fi

done