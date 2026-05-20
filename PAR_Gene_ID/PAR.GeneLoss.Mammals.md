# Identify all genes found in any avian PAR, curate fastas for blast analysis, then identify genes within PARs of all genomes
## 0. Quantify gappiness of each PAR
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_analysis/mammals

## Excluding:
Microtus_pennsylvanicus,CHROM:0-116367
Ochotona_princeps,CHROM:107724012-107801096 



# Confirm that accessions in sexchrfile still match correctly (genomes have been updated)
SPECIES=Urocitellus_parryii
grep '>' ${BASE}/${SPECIES}/${SPECIES}.fna | grep 'chromosome'

# Revise those that need revision
SEXCHR_FILE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"

grep 'Urocitellus_parryii' $SEXCHR_FILE
sed -i 's/CM099876/NC_135547/g' $SEXCHR_FILE
sed -i 's/CM099877/NC_135548/g' $SEXCHR_FILE

Capra_hircus # no matches at all


cat > PAR.species_chr_region.txt <<'EOF'
Eubalaena_glacialis,NC_083736.1:0-7069966
Macaca_nemestrina,NC_092145.1:158149255-159757195 
Callithrix_jacchus,NC_133525.1.1:0-1757279
Mesoplodon_bidens,OZ073217.1:135117924-142816029
Inia_geoffrensis,CM070920.1:0-7142434
Camelus_dromedarius,NC_087472.1:108540129-114202744
Ovis_canadensis,NC_091727.1:0-7072606
Manis_pentadactyla,NC_080038.1:0-4881534
Pan_paniscus,NC_073272.2:0-2524164                      
Marmota_flaviventris,NC_092518.1:521224-10347709        
Balaenoptera_physalus,OZ239531.1:0-7257636
Pan_troglodytes,NC_072421.2:0-3170188          
Rhynchonycteris_naso,CM073052.1:0-10832599
Ovis_aries,CP162266.1:0-7021881
Trichechus_inunguis,CM102173.1:0-9558467
Pseudorca_crassidens,NC_090317.1:127639906-136059651 
Myotis_nattereri,OZ125678.2:0-1449976
Rhynchocyon_petersi,CM091802.1:0-20182068
Grampus_griseus,OZ206318.1:0-7148355
Mustela_nivalis_vulgaris,OZ211688.1:127035994-138477583
Capra_hircus,CP168640.1:0-7040093
Loxodonta_africana,NC_087369.1:0-10314587
Urocitellus_parryii,NC_135547.1:122998411-131433435 
Meles_meles,NC_060087.1:0-6387864
Homo_sapiens,NC_060947.1:0-2394410
EOF

chmod +x find_PAR_gaps.sh 
./find_PAR_gaps.sh

Rscript plot_PAR_gaps.R
```
### Species without gaps
Loxodonta africana has four gaps each of 13 base pairs... this feels a little suspicious to me? Since assemblers can put in any random number of bp for a gap? Excluding it for now.
```
cat > species_par.csv <<'EOF'
Pseudorca_crassidens,127639906-136059651 
Inia_geoffrensis,0-7142434
Ovis_canadensis,0-7072606
Capra_hircus,0-7040093
Ovis_aries,0-7021881
Meles_meles,0-6387864
Pan_troglodytes,0-3170188          
Pan_paniscus,0-2524164                      
Homo_sapiens,0-2394410
Callithrix_jacchus,0-1757279
EOF
```
## 1. Curate fastas of genes found in any avian PAR
```

chmod +x extract_z_par_genes.sh
./extract_x_par_genes.sh species_par.csv > x_par_genes.tsv
```
# Genes from x_par_genes.tsv that are not being used in blast
```
serine/arginine repetitive matrix protein 1-like
small nucleolar RNA U109
TRNAC-ACA
TRNAC-GCA
TRNAE-UUC
TRNAW-CCA
U6 spliceosomal RNA
SNORA70
SNORD82

# too many blast matches 
    ARSD
    ARSF
    ARSH
    ARSL
```
# Genes from x_par_genes.tsv found in 1+ species:
cat > genes.txt <<'EOF'
AKAP17A
ANOS1
APOA1A
ARSD
ARSF
ARSH
ARSL
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
NA
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
EOF

GFF="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Homo_sapiens/Homo_sapiens.gff"

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
' genes.txt "$GFF" >> genes.regions.Homo_sapiens.txt

# repeat using Ovis canadensis for any "NOT FOUND" genes
```
grep 'NOT_FOUND' genes.regions.Homo_sapiens.txt

cat > genes.Ovis_canadensis.txt <<'EOF'
LOC138930114
LOC138930123
LOC138930124
LOC138930131
LOC138930134
LOC138930121
EOF

GFF="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Ovis_canadensis/Ovis_canadensis.gff"

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

            # Always prefer NC_091727.1 if present
            if (!(gene_name in best) || $1 == "NC_091727.1") {
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
' genes.Ovis_canadensis.txt "$GFF" >> genes.regions.Ovis_canadensis.txt
```
# replace LOC with interpretable gene names
```
sed -i 's/LOC138930123/BOSD2_1/g' genes.regions.Ovis_canadensis.txt
sed -i 's/LOC138930124/BOSD2_2/g' genes.regions.Ovis_canadensis.txt
sed -i 's/LOC138930131/BOSD2_3/g' genes.regions.Ovis_canadensis.txt
sed -i 's/LOC138930134/BOSD2_4/g' genes.regions.Ovis_canadensis.txt
sed -i 's/LOC138930114/APOA1A/g' genes.regions.Ovis_canadensis.txt
sed -i 's/LOC138930121/OBP/g' genes.regions.Ovis_canadensis.txt
```
# Pull fastas from reference genomes (Homo_sapiens, Ovis_canadensis)
```
## Homo_sapiens
grep -v NOT_FOUND genes.regions.Homo_sapiens.txt > genes.regions.Homo_sapiens.noNA.txt

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Homo_sapiens/Homo_sapiens.fna"

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

done < genes.regions.Homo_sapiens.noNA.txt
```
## Ovis_canadensis
```
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Ovis_canadensis/Ovis_canadensis.fna"

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

done < genes.regions.Ovis_canadensis.txt

```
## 2. Blast all avian genomes for the fastas of putatively conserved PAR genes
```
awk -v FS=',' '{print $1}' PAR.species_chr_region.txt > species_for_blast_array.txt

N=$(wc -l < species_for_blast_array.txt)
sbatch --array=1-${N} blast_gene_locations_array.mammals.sh
```
### Filter blast outputs
```


mkdir -p blast_results_Xchr

for f in blast_results/*.blast.tsv; do
    base=$(basename "$f")
    species=$(echo "$base" | cut -d'.' -f1)
    out="blast_results_Xchr/$base"

    count=$(
        awk -v species="$species" '
            BEGIN {
                FS = "[,\t ]+"
            }

            NR == FNR {
                if ($1 == species && $2 == "X") {
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
chmod +x Blast_Match_inPAR.mammals.sh
./Blast_Match_inPAR.mammals.sh
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
### Infer gene overlap
```
awk -F',' '
BEGIN {
    OFS = FS
}

# First file: species_par.csv
# Format: Species,start-stop
# These coordinates refer to the X chromosome for that species.
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

    # Only evaluate X chromosomes.
    # Anything else, including W or NA, is unknown for this PAR test.
    if (chr_label != "X") {
        status = "U"
    } else if (!(species in n)) {
        status = "U"
    } else {
        status = "N"

        for (i = 1; i <= n[species]; i++) {
            # Any overlap between gene interval and species X PAR interval
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

# Filter to just X matches
awk -F',' 'BEGIN { OFS = FS } NR == 1 || $(NF-1) == "X"' \
  gene_locations_by_species.with_chr_label.all.In_PAR.csv \
  > gene_locations_by_species.with_chr_label.all.In_PAR.X_only.csv


# Filter to gapless species 
awk -v FS=',' '{print $1}' species_par.csv > gapless_species.txt

input="gene_locations_by_species.with_chr_label.all.In_PAR.X_only.csv"
output="gene_locations_by_species.with_chr_label.all.In_PAR.X_only.gapless_species.csv"
species_list="gapless_species.txt"

# Keep header, then rows matching one of the gapless species.
{
    head -n 1 "$input"
    tail -n +2 "$input" | grep -Ff "$species_list"
} > "$output"

echo "Wrote: $output"
echo "Species retained:"
cut -d',' -f1 "$output" | tail -n +2 | sort -u
```
### Note if telomere is ID'd on the end of the PAR
```
telomere_file="gapless_species.telomeres.txt"

cat > "$telomere_file" <<'EOF'
Pseudorca_crassidens,R,YES
Inia_geoffrensis,L,YES
Ovis_canadensis,L,YES
Capra_hircus,L,YES
Ovis_aries,L,YES
Meles_meles,L,YES
Pan_troglodytes,L,
Pan_paniscus,L,YES
Homo_sapiens,L,YES
Callithrix_jacchus,L,YES
EOF
```
## 3. Plot the data
### Plot a PCA of PAR status by gene and species
```
module load R

R

library(tidyverse)
library(ggrepel)

infile <- "gene_locations_by_species.with_chr_label.all.In_PAR.X_only.gapless_species.csv"

df <- read_csv(infile, show_col_types = FALSE)

# Collapse multiple Z hits per Species/Gene.
# Priority:
#   Y/Edge > N
z_states <- df %>%
  mutate(
    state = case_when(
      In_PAR %in% c("Y", "Edge") ~ 2,
      In_PAR == "N" ~ 1,
      In_PAR == "U" ~ 0,
      is.na(In_PAR) ~ 0,
      TRUE ~ 0
    )
  ) %>%
  group_by(Species, Gene) %>%
  summarise(
    state = max(state),
    .groups = "drop"
  )

# Make the full Species x Gene matrix.
# Missing Species/Gene combinations are interpreted as not found on Z.
mat <- z_states %>%
  pivot_wider(
    names_from = Gene,
    values_from = state,
    values_fill = 0
  )

species <- mat$Species

pca_mat <- mat %>%
  select(-Species) %>%
  as.data.frame()

rownames(pca_mat) <- species

# Remove genes with no variation across species
pca_mat <- pca_mat[, apply(pca_mat, 2, var) > 0]

pca <- prcomp(pca_mat, center = TRUE, scale. = TRUE)

pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("Species")

var_explained <- summary(pca)$importance[2, ] * 100

ggplot(pca_df, aes(PC1, PC2, label = Species)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.1,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  xlab(paste0("PC1 (", round(var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(var_explained[2], 1), "%)")) +
  theme_bw()

ggsave("In_PAR_PCA_from_X_only.pdf", width = 7, height = 6)

write_csv(
  mat,
  "species_by_gene.In_PAR_numeric_matrix.from_X_only.csv"
)
```
### Plot a PCA of PAR status by gene
```
library(ggrepel)

# mat is Species x Gene:
# first column = Species
# remaining columns = numeric states 0/1/2

gene_pca_mat <- mat %>%
  column_to_rownames("Species") %>%
  t() %>%
  as.data.frame()

# Remove species columns with no variation across genes
gene_pca_mat <- gene_pca_mat[, apply(gene_pca_mat, 2, var) > 0]

# Remove genes with no variation across species after transpose, if any
gene_pca_mat <- gene_pca_mat[apply(gene_pca_mat, 1, var) > 0, ]

gene_pca <- prcomp(gene_pca_mat, center = TRUE, scale. = TRUE)

gene_pca_df <- as.data.frame(gene_pca$x) %>%
  rownames_to_column("Gene")

gene_var_explained <- summary(gene_pca)$importance[2, ] * 100

ggplot(gene_pca_df, aes(PC1, PC2)) +
  geom_point(size = 3) +
  geom_text_repel(
    aes(label = Gene),
    size = 3,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.alpha = 0.5
  ) +
  xlab(paste0("PC1 (", round(gene_var_explained[1], 1), "%)")) +
  ylab(paste0("PC2 (", round(gene_var_explained[2], 1), "%)")) +
  theme_bw()

ggsave("In_PAR_PCA_genes_labeled.pdf", width = 10, height = 8)

write_csv(gene_pca_df, "In_PAR_PCA_gene_coordinates.csv")
```
# Make combined plot
```
Rscript ../PAR_Combined_Plot.R mammals
```
# Permissions modification
chmod -R g+rwx /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR*
chmod -R g+rwx /data/Wilson_Lab/data/VGP_genomes_phase1/
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists
chmod -R g+rwx /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles