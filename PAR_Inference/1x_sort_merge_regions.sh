#!/bin/bash

LIST="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/alignment_beds.txt"
OUTDIR="sorted_beds"

mkdir -p "$OUTDIR"

while read -r bed; do
    [ -z "$bed" ] && continue
    [ ! -f "$bed" ] && continue

    base=$(basename "$bed")

    awk 'BEGIN{OFS="\t"} {print $4, $5, $6}' "$bed" \
        | sort -k1,1 -k2,2n -k3,3n \
        | awk '
            BEGIN {OFS="\t"}
            NR==1 {
                chr=$1; start=$2; end=$3;
                next
            }
            {
                if ($1 == chr && $2 <= end) {
                    # overlapping or contiguous: extend region
                    if ($3 > end) end = $3
                } else {
                    # emit previous
                    print chr, start, end
                    chr=$1; start=$2; end=$3
                }
            }
            END {
                # print final interval
                print chr, start, end
            }
        ' > "$OUTDIR/$base"

done < "$LIST"