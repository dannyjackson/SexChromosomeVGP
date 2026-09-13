#!/usr/bin/env awk -f

# Usage:
# awk -v GENE="XG" -f orient_to_gene_midpoint.awk mammals_PAR_genes.all.ZNF_arrays.tsv > mammals_PAR_genes.XG_oriented.tsv

BEGIN {
    FS = OFS = "\t"
}

# First pass: record midpoint of anchor gene per Species + TOGADir + Chromosome
FNR == NR {
    if (FNR == 1) next

    key = $1 OFS $2 OFS $3

    if ($6 == GENE) {
        anchor_mid[key] = ($4 + $5) / 2
        anchor_start[key] = $4
        anchor_stop[key]  = $5
    }

    next
}

# Second pass: print genes relative to anchor midpoint
FNR == 1 {
    print $0, "AnchorGene", "AnchorMidpoint", "RelStart", "RelStop", "RelMidpoint"
    next
}

{
    key = $1 OFS $2 OFS $3

    # Only print rows from chromosome/scaffold records that contain the anchor gene
    if (!(key in anchor_mid)) next

    rel_start = $4 - anchor_mid[key]
    rel_stop  = $5 - anchor_mid[key]
    rel_mid   = (($4 + $5) / 2) - anchor_mid[key]

    print $0, GENE, anchor_mid[key], rel_start, rel_stop, rel_mid
}
