#!/bin/bash

export SEXCHRLEN=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_phylogeny/VGP_freeze_hap1_combined_sexchroms_seq_reports.tsv
PAFDIR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2

# Fix known spelling mismatch
sed -i 's/Guaruba_guaruba/Guaruba_guarouba/g' "$SEXCHRLEN"

OUTFILE="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/paf_sexchr_info.tsv"
ERRFILE="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/paf_sexchr_info_errors.tsv"

printf "Accession\tXZ_Accession\tXZ_Seq_length\tYW_Accession\tYW_Seq_length\tSpecies\tFilename\n" > "$OUTFILE"
printf "Filename\tSpecies_key\tComparison\tError\tDetail\n" > "$ERRFILE"

# Get assembly accession + canonical Species label from SEXCHRLEN (handles subspp ambiguity)
find_species_and_assembly() {
    local species_key="$1"
    local has_sub="$2"

    awk -v FS="\t" -v OFS="\t" -v spkey="$species_key" -v has_sub="$has_sub" '
        NR==1 { next }
        {
            full = $16
            n = split(full, a, "_")
            key2 = (n>=2 ? a[1]"_"a[2] : full)
            key3 = (n>=3 ? key2"_"a[3] : key2)

            if ((has_sub==1 && key3==spkey) || (has_sub==0 && key2==spkey)) {
                count++
                acc[count]  = $1
                spec[count] = full
            }
        }
        END {
            if (count==0) { print "NOTFOUND"; exit }

            # uniq species labels
            uniq=0
            for (i=1;i<=count;i++){
                exists=0
                for (j=1;j<=uniq;j++) if (spec[i]==uspec[j]) {exists=1; break}
                if (!exists) { uniq++; uspec[uniq]=spec[i]; uacc[uniq]=acc[i] }
            }

            if (uniq>1) {
                s=uspec[1]
                for (k=2;k<=uniq;k++) s=s";"uspec[k]
                print "AMBIG\t"s
            } else {
                print uacc[1], uspec[1]
                print "FOUND"
            }
        }
    ' "$SEXCHRLEN"
}

# Get expected length for a specific chromosome (X/Z or Y/W) from SEXCHRLEN
find_expected_len() {
    local species_key="$1"
    local has_sub="$2"
    local chr="$3"

    awk -v FS="\t" -v OFS="\t" -v spkey="$species_key" -v has_sub="$has_sub" -v ch="$chr" '
        NR==1 { next }
        {
            full = $16
            n = split(full, a, "_")
            key2 = (n>=2 ? a[1]"_"a[2] : full)
            key3 = (n>=3 ? key2"_"a[3] : key2)

            if ( ((has_sub==1 && key3==spkey) || (has_sub==0 && key2==spkey)) && $4==ch ) {
                count++
                len[count]  = $12
                spec[count] = full
            }
        }
        END {
            if (count==0) { print "NOTFOUND"; exit }

            uniq=0
            for (i=1;i<=count;i++){
                exists=0
                for (j=1;j<=uniq;j++) if (spec[i]==uspec[j]) {exists=1; break}
                if (!exists) { uniq++; uspec[uniq]=spec[i]; ulen[uniq]=len[i] }
            }

            if (uniq>1) {
                s=uspec[1]
                for (k=2;k<=uniq;k++) s=s";"uspec[k]
                print "AMBIG\t"s
            } else {
                print ulen[1]
                print "FOUND"
            }
        }
    ' "$SEXCHRLEN"
}

# Parse the first alignment line of a PAF: qname qlen tname tlen
parse_paf_first() {
    local paf="$1"
    awk -v FS="\t" '
        NF>=12 && $0 !~ /^#/ {
            print $1 "\t" $2 "\t" $6 "\t" $7
            exit
        }
    ' "$paf"
}

abs() { awk -v x="$1" 'BEGIN{print (x<0?-x:x)}'; }

for paf in "$PAFDIR"/*.paf; do
    [ -f "$paf" ] || continue
    base=$(basename "$paf")

    # Strip at first dot: Genus_species[_subspecies]_comparison...
    core=${base%%.*}

    # Split core into tokens
    IFS='_' read -r genus species third fourth rest <<< "$core"
    unset IFS

    if [[ -z "${genus:-}" || -z "${species:-}" ]]; then
        printf "%s\t\t\tMalformed filename (no genus/species)\t\n" "$base" >> "$ERRFILE"
        continue
    fi

    has_sub=0
    subspecies=""
    comparison_code=""

    if [[ -n "${third:-}" ]]; then
        firstchar=${third:0:1}
        if [[ "$firstchar" =~ [a-z] ]]; then
            subspecies="$third"
            has_sub=1
            species_key="${genus}_${species}_${subspecies}"
            comparison_code="${fourth:-}"
        else
            has_sub=0
            species_key="${genus}_${species}"
            comparison_code="$third"
        fi
    else
        has_sub=0
        species_key="${genus}_${species}"
        comparison_code=""
    fi

        # Determine which pair this is (XZ vs YW)
    chr_xz=""
    chr_yw=""

    if [[ "$comparison_code" == *"WtoZ"* ]]; then
        chr_xz="Z"
        chr_yw="W"

    # neo-Y cases: Y1/Y2 mapping to standard X
    elif [[ "$comparison_code" == *"Y1toX"* ]]; then
        chr_xz="X"
        chr_yw="Y1"
    elif [[ "$comparison_code" == *"Y2toX"* ]]; then
        chr_xz="X"
        chr_yw="Y2"

    # neo-X cases: Y to X1/X2
    elif [[ "$comparison_code" == *"YtoX1"* ]]; then
        chr_xz="X1"
        chr_yw="Y"
    elif [[ "$comparison_code" == *"YtoX2"* ]]; then
        chr_xz="X2"
        chr_yw="Y"

    # standard XY
    elif [[ "$comparison_code" == *"YtoX"* ]]; then
        chr_xz="X"
        chr_yw="Y"

    # if the filename comparison code is "X1..." or "X2..." (keep your original behavior)
    elif [[ "$comparison_code" == *"X1"* ]]; then
        chr_xz="X1"
        chr_yw="Y"
    elif [[ "$comparison_code" == *"X2"* ]]; then
        chr_xz="X2"
        chr_yw="Y"
    fi


    if [[ -z "$chr_xz" || -z "$chr_yw" ]]; then
        printf "%s\t%s\t%s\tUnable to determine XZ/YW chromosomes\t\n" \
            "$base" "$species_key" "$comparison_code" >> "$ERRFILE"
        continue
    fi

    # Assembly + Species label
    spasm=$(find_species_and_assembly "$species_key" "$has_sub")
    sp_status=$(echo "$spasm" | tail -n 1 | cut -f1)

    if [[ "$sp_status" == "NOTFOUND" ]]; then
        printf "%s\t%s\t%s\tNo matching species in SEXCHRLEN\t\n" \
            "$base" "$species_key" "$comparison_code" >> "$ERRFILE"
        continue
    elif [[ "$sp_status" == "AMBIG" ]]; then
        detail=$(echo "$spasm" | tail -n 1 | cut -f2-)
        printf "%s\t%s\t%s\tMultiple subspecies matched (species)\t%s\n" \
            "$base" "$species_key" "$comparison_code" "$detail" >> "$ERRFILE"
        continue
    fi

    asm_line=$(echo "$spasm" | head -n -1)
    accession=$(echo "$asm_line" | cut -f1)
    full_species=$(echo "$asm_line" | cut -f2-)

    # Expected lengths from SEXCHRLEN for disambiguating PAF direction
    exz=$(find_expected_len "$species_key" "$has_sub" "$chr_xz")
    exz_status=$(echo "$exz" | tail -n 1 | cut -f1)
    eyw=$(find_expected_len "$species_key" "$has_sub" "$chr_yw")
    eyw_status=$(echo "$eyw" | tail -n 1 | cut -f1)

    if [[ "$exz_status" != "FOUND" ]]; then
        detail=$(echo "$exz" | tail -n 1 | cut -f2-)
        printf "%s\t%s\t%s\tExpected length not resolved for %s\t%s\n" \
            "$base" "$species_key" "$comparison_code" "$chr_xz" "$detail" >> "$ERRFILE"
        continue
    fi
    if [[ "$eyw_status" != "FOUND" ]]; then
        detail=$(echo "$eyw" | tail -n 1 | cut -f2-)
        printf "%s\t%s\t%s\tExpected length not resolved for %s\t%s\n" \
            "$base" "$species_key" "$comparison_code" "$chr_yw" "$detail" >> "$ERRFILE"
        continue
    fi

    exz_val=$(echo "$exz" | head -n -1 | tail -n 1)
    eyw_val=$(echo "$eyw" | head -n -1 | tail -n 1)

    # Parse PAF first alignment
    paf_fields=$(parse_paf_first "$paf")
    if [[ -z "$paf_fields" ]]; then
        printf "%s\t%s\t%s\tEmpty/invalid PAF\t\n" \
            "$base" "$species_key" "$comparison_code" >> "$ERRFILE"
        continue
    fi

    qname=$(echo "$paf_fields" | cut -f1)
    qlen=$(echo "$paf_fields" | cut -f2)
    tname=$(echo "$paf_fields" | cut -f3)
    tlen=$(echo "$paf_fields" | cut -f4)

    # Decide which side is XZ vs YW by closest expected lengths
    # Option A: query=YW, target=XZ
    dA=$(awk -v q="$qlen" -v t="$tlen" -v ey="$eyw_val" -v ex="$exz_val" \
        'BEGIN{da=(q-ey); if(da<0) da=-da; db=(t-ex); if(db<0) db=-db; print da+db}')
    # Option B: query=XZ, target=YW
    dB=$(awk -v q="$qlen" -v t="$tlen" -v ex="$exz_val" -v ey="$eyw_val" \
        'BEGIN{da=(q-ex); if(da<0) da=-da; db=(t-ey); if(db<0) db=-db; print da+db}')

    if awk -v a="$dA" -v b="$dB" 'BEGIN{exit !(a<=b)}'; then
        # A wins: query=YW, target=XZ
        xz_acc="$tname"; xz_len="$tlen"
        yw_acc="$qname"; yw_len="$qlen"
    else
        # B wins: query=XZ, target=YW
        xz_acc="$qname"; xz_len="$qlen"
        yw_acc="$tname"; yw_len="$tlen"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$accession" "$xz_acc" "$xz_len" "$yw_acc" "$yw_len" "$full_species" "$base" >> "$OUTFILE"

done

# Manually add Panthera onca since it is not in the reference sheet
printf "NA\tCM102116.1\t130846311\tCM102117.1\t30836324\tPanthera_onca\tPanthera_onca_YtoX.aln.paf\n" >> "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/paf_sexchr_info.tsv"

# also manually add Bos taurus since it is weird in the reference sheet
printf "NA\tNC_040105.1\t146092946\tCM011803.1\t15658480\tBos_taurus\tBos_taurus_YtoX.aln.paf\n" >> "/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/paf_sexchr_info.tsv"

# Manually swap Macropus eugenii between accession convention
sed -i 's/NC_092879.1/CM051820.1/g' /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/paf_sexchr_info.tsv
sed -i 's/NC_092880.1/CM051821.1/g' /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/paf_sexchr_info.tsv

mv /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/junk_continuous_percentID/Bos_taurus_YtoX.aln.paf /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/