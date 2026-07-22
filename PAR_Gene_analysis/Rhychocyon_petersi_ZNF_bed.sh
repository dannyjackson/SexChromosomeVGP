# Code to create bed files for zinc finger genes in the Rhynchocyon petersi PAR

zcat query_annotation.gtf.gz |
awk -F '\t' '
$1 == "CM091802" && $0 ~ /ZNF/ {
    if (match($9, /gene_id "[^"]+"/)) {
        gene_id = substr($9, RSTART + 9, RLENGTH - 10)
    }

    if (match($9, /#[^#]*ZNF[^#]*#/)) {
        gene_name = substr($9, RSTART + 1, RLENGTH - 2)
        names[gene_id] = gene_name
    }
}

$1 == "CM091802" && $3 == "gene" {
    if (match($9, /gene_id "[^"]+"/)) {
        gene_id = substr($9, RSTART + 9, RLENGTH - 10)
        chrom[gene_id] = $1
        start[gene_id] = $4 - 1
        end[gene_id] = $5
    }
}

END {
    for (gene_id in names) {
        if (gene_id in chrom)
            print names[gene_id], chrom[gene_id], start[gene_id], end[gene_id]
    }
}' OFS='\t' > ZNF_genes.with_gene_name.bed


awk -F '\t' '{print $2,$3,$4}' OFS='\t' ZNF_genes.with_gene_name.bed > ZNF_genes.bed

awk -F '\t' '
$1 == "CM091802" && $3 < 20182068 {print $0}' ZNF_genes.bed > ZNF_genes.PAR.bed
