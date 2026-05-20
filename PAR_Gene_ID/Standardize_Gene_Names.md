# Build a database of gene descriptions and gene names
From the sex chromosomes of all VGP GFFs, build a database of gene descriptions and gene names. 
```
#!/usr/bin/env bash

base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"
accessions="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv"
outdir="sex_chr_gene_dict"

mkdir -p "$outdir"

tail -n +2 "$accessions" | cut -d, -f1 | sort -u | while read -r species
do
    gff="${base}/${species}/${species}.gff"
    outfile="${outdir}/${species}.tsv"

    if [ ! -f "$gff" ]; then
        echo "Warning: missing GFF for ${species}: ${gff}" >&2
        continue
    fi

    echo "Processing ${species}" >&2

    awk \
      -v species="$species" \
      -v accessions="$accessions" \
      -f gff_sex_chr_gene_dict.awk \
      "$gff" | \
    awk 'BEGIN{FS=OFS="\t"} NR==1 || !seen[$0]++' \
      > "$outfile"
done


```
# Combine across species
Combine each unique gene description and gene name pair across all species into a single file.

Output instances of multiple mapping to a separate file.

```
#!/usr/bin/env bash

input_dir="sex_chr_gene_dict"
database_file="gene_name_description_database.tsv"
conflicts_file="gene_name_description_database.conflicts.tsv"

tmp_all="$(mktemp)"
tmp_counts="$(mktemp)"
tmp_single="$(mktemp)"
tmp_conflicts="$(mktemp)"

awk '
BEGIN { FS=OFS="\t" }

FNR == 1 { next }

NF >= 2 && $1 != "" && $2 != "" {
    description = $1
    gene_name = $2

    # Trim leading/trailing whitespace
    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", description)
    gsub(/^[ \t\r\n]+|[ \t\r\n]+$/, "", gene_name)

    # Skip generic/predicted locus names
    if (gene_name ~ /^LOC[0-9]+$/) next
    if (gene_name ~ /^[A-Z0-9]+_[0-9]+$/) next

    # Collapse case-only gene-name differences:
    # CNP, cnp, Cnp -> CNP
    gene_name = toupper(gene_name)

    print description, gene_name
}
' "$input_dir"/*.tsv \
| sort -u > "$tmp_all"

# Count how many distinct gene names each description maps to after filtering/case normalization.
awk '
BEGIN { FS=OFS="\t" }
{
    desc = $1
    gene = $2
    seen[desc SUBSEP gene] = 1
}
END {
    for (k in seen) {
        split(k, a, SUBSEP)
        desc = a[1]
        count[desc]++
    }

    for (desc in count) {
        print desc, count[desc]
    }
}
' "$tmp_all" > "$tmp_counts"

# Descriptions with exactly one gene name.
awk '
BEGIN { FS=OFS="\t" }

NR == FNR {
    if ($2 == 1) single[$1] = 1
    next
}

$1 in single {
    print $1, $2
}
' "$tmp_counts" "$tmp_all" \
| sort -k1,1 > "$tmp_single"

# Descriptions with multiple gene-name matches.
awk '
BEGIN { FS=OFS="\t" }

NR == FNR {
    if ($2 > 1) conflict[$1] = 1
    next
}

$1 in conflict {
    print $1, $2
}
' "$tmp_counts" "$tmp_all" \
| sort -k1,1 -k2,2 > "$tmp_conflicts"

{
    printf "Description\tGeneName\n"
    cat "$tmp_single"
} > "$database_file"

{
    printf "Description\tGeneName\n"
    cat "$tmp_conflicts"
} > "$conflicts_file"

rm -f "$tmp_all" "$tmp_counts" "$tmp_single" "$tmp_conflicts"

echo "Wrote: $database_file"
echo "Wrote: $conflicts_file"
```
# Rename X genes using this database
```
#!/usr/bin/env bash

db_file="gene_name_description_database.tsv"
input_file="x_par_genes.tsv"
output_file="x_par_genes.with_inferred.tsv"

awk '
BEGIN {
    FS = OFS = "\t"
}

# -----------------------------
# First file: gene_name_description_database.tsv
# -----------------------------
ARGIND == 1 {
    if (FNR == 1) next

    desc = clean_description($1)
    gene = $2

    if (desc != "" && gene != "") {
        db[desc] = gene
    }

    next
}

# -----------------------------
# Second file: x_par_genes.tsv
# -----------------------------
ARGIND == 2 {
    if (FNR == 1) {
        print $0, "InferredGeneName"
        next
    }

    species = $1
    chromosome = $2
    start_pos = $3
    stop_pos = $4
    gene_name = $5
    gene_description = $6
    par_status = $7

    clean_desc = clean_description(gene_description)

    inferred_gene_name = ""

    if (clean_desc in db) {
        inferred_gene_name = db[clean_desc]
    }

    print $0, inferred_gene_name
}

# -----------------------------
# Cleaning function used on both files
# -----------------------------
function clean_description(s) {
    gsub(/%2C/, ",", s)
    gsub(/%2c/, ",", s)

    s = tolower(s)

    # Remove annotation suffixes/prefixes
    gsub(/-like/, "", s)
    gsub(/^putative /, "", s)
    gsub(/ putative /, " ", s)

    # Normalize common wording differences
    gsub(/anchor protein/, "anchoring protein", s)

    # Normalize punctuation/spaces
    gsub(/[ \t]+/, " ", s)
    gsub(/^[ \t\r\n]+/, "", s)
    gsub(/[ \t\r\n]+$/, "", s)

    return s
}
' "$db_file" "$input_file" > "$output_file"

echo "Wrote: $output_file"
```
# output unique gene descriptions with no matched gene names
```
awk -F'\t' '
BEGIN { OFS="\t" }

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "GeneDescription") desc_col = i
        if ($i == "InferredGeneName") inferred_col = i
    }

    print "GeneDescription"
    next
}

inferred_col && desc_col && $inferred_col == "" {
    desc = $desc_col

    gsub(/%2C/, ",", desc)
    gsub(/%2c/, ",", desc)
    gsub(/^[ \t\r\n]+/, "", desc)
    gsub(/[ \t\r\n]+$/, "", desc)
    gsub(/[ \t]+/, " ", desc)

    if (desc != "") seen[desc] = 1
}

END {
    for (desc in seen) print desc
}
' x_par_genes.with_inferred.tsv | sort > unmatched_gene_descriptions.tsv
```
# Identify errors
Output instances of conflict between GeneName and InferredGeneName for all instances without an LOC in GeneName
```
awk -F'\t' '
BEGIN { OFS="\t" }

NR == 1 {
    print $0
    for (i = 1; i <= NF; i++) {
        if ($i == "GeneName") gene_col = i
        if ($i == "InferredGeneName") inferred_col = i
    }
    next
}

gene_col && inferred_col {
    gene = $gene_col
    inferred = $inferred_col

    # Skip LOC-style gene names
    if (gene ~ /^LOC[0-9]+$/) next

    # Skip blank inferred names
    if (inferred == "") next

    # Compare case-insensitively
    if (toupper(gene) != toupper(inferred)) {
        print $0
    }
}
' x_par_genes.with_inferred.tsv > gene_name_inferred_conflicts.tsv

```
# Correct errors
The only error observed is that all genes named "APOO" are now named "APOOL"
```
awk -F'\t' '
BEGIN { OFS = "\t" }

NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "GeneName") gene_col = i
        if ($i == "InferredGeneName") inferred_col = i
    }
    print
    next
}

{
    if (gene_col && inferred_col && toupper($gene_col) == "APOO" && toupper($inferred_col) == "APOOL") {
        $inferred_col = "APOO"
    }

    print
}
' x_par_genes.with_inferred.tsv > x_par_genes.with_inferred.apoo_fixed.tsv

# Replace original with revised if it looks good
grep 'APOO' x_par_genes.with_inferred.apoo_fixed.tsv
mv x_par_genes.with_inferred.apoo_fixed.tsv x_par_genes.with_inferred.tsv 

```