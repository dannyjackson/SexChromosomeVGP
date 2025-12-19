#!/bin/bash

BEDDIR="sorted_beds"
OUTFILE="bed_spacers.csv"

# truncate / create output file
: > "$OUTFILE"

for bed in "$BEDDIR"/*; do
    [ -f "$bed" ] || continue

    awk -v fname="$(basename "$bed")" 'BEGIN{OFS=","}
        NR == 1 {
            prev_chr = $1
            prev_start = $2
            prev_end = $3
            next
        }
        {
            chr = $1
            start = $2
            end = $3

            # only compute distance within same chromosome
            if (chr == prev_chr) {
                dist = start - prev_end
                printf "%s,%s_%s,%d\n", fname, prev_end, start, dist
            }

            prev_chr = chr
            prev_start = start
            prev_end = end
        }
    ' "$bed" >> "$OUTFILE"

done