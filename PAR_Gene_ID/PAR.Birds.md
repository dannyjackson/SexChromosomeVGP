# ID genomes with well-assembled PARs
This script follows the curation of genomes that show no evidence of PAR dropout in HalfDeep plots. Other forms of PAR misassembly can be missed by the HalfDeep assessment. We need to check if the ancestral PAR is found on the end of the Z and W chromosomes.

To do this, we will first identify genes conserved in the PARs across multiple genomes that pass the HalfDeep filter. We do not know what proportion of these genomes have misassemblies, so this analysis uses cutoffs that are inferred from the data itself, not predetermined.

## Set up the environment
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds
sinteractive --mem=20g
module load bcftools vcftools samtools R
```
# Plot regions of sequence alignment of sex chromosomes
mv /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/align_PAR_985thr/PAR_annotations/Mustela_nivalis_vulgaris_YtoX.aln.id98_5.len10k.refqry.bed /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/align_PAR_985thr/old/Mustela_nivalis_vulgaris_YtoX.aln.id98_5.len10k.refqry.bed 

mv /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/continuous_percentID/Mustela_nivalis_vulgaris_YtoX.aln.refqry.csv /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/old/Mustela_nivalis_vulgaris_YtoX.aln.refqry.csv

# Process Molossus
paf="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/PAFs_22JUNE2026/filtered_98.5_10kb_Molossus_nigricans_GCA_039880945.1_hap1_GCA_039880925.1_hap2_YtoX.aln.paf"

out="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/PAFs_22JUNE2026/filtered_98.5_10kb_Molossus_nigricans_GCA_039880945.1_hap1_GCA_039880925.1_hap2_YtoX.aln.id98_5.len10k.refqry.bed"

awk 'BEGIN{
  OFS="\t"
  print "chrom_qry","len_qry","bp_start_qry","bp_end_qry","percent_identity_qry","chrom_ref","len_ref","bp_start_ref","bp_end_ref","percent_identity_ref"
}
{
  pid = ($11 > 0 ? 100 * $10 / $11 : 0)

  print $1, $2, $3, $4, sprintf("%.4f", pid), $6, $7, $8, $9, sprintf("%.4f", pid)
}' "$paf" > "$out"

cp $out /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/align_PAR_985thr/PAR_annotations/Molossus_nigricans_YtoX.aln.id98_5.len10k.refqry.bed 

cp $out /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/continuous_percentID/Molossus_nigricans_YtoX.aln.refqry.csv


# Process Mustela
paf="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/PAFs_22JUNE2026/filtered_98.5_10kb_Mustela_nivalis_vulgaris_GCA_057128415.1_hap1_GCA_057128425.1_hap2_YtoX.aln.paf"

out="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/PAFs_22JUNE2026/filtered_98.5_10kb_Mustela_nivalis_vulgaris_GCA_057128415.1_hap1_GCA_057128425.1_hap2_YtoX.aln.id98_5.len10k.refqry.csv"

awk 'BEGIN{
  OFS="\t"
  print "chrom_qry","len_qry","bp_start_qry","bp_end_qry","percent_identity_qry","chrom_ref","len_ref","bp_start_ref","bp_end_ref","percent_identity_ref"
}
{
  pid = ($11 > 0 ? 100 * $10 / $11 : 0)

  print $1, $2, $3, $4, sprintf("%.4f", pid), $6, $7, $8, $9, sprintf("%.4f", pid)
}' "$paf" > "$out"


cp $out /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/align_PAR_985thr/PAR_annotations/Mustela_nivalis_vulgaris_YtoX.aln.id98_5.len10k.refqry.bed 

cp $out /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/continuous_percentID/Mustela_nivalis_vulgaris_YtoX.aln.refqry.csv

head -n 1 /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/continuous_percentID/Vipera_berus_WtoZ.aln.refqry.csv > \
    /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/continuous_percentID/Mustela_nivalis_vulgaris_YtoX.aln.refqry.csv
cat /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/align_PAR_985thr/PAR_annotations/Mustela_nivalis_vulgaris_YtoX.aln.id98_5.len10k.refqry.bed >> \
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/continuous_percentID/Mustela_nivalis_vulgaris_YtoX.aln.refqry.csv


# Plot with R script

Rscript PAR_Plotting_AllSpecies.R
# Use the output csv of surviving regions to identify PAR boundaries in the species with half-deep inferred no-PAR-dropout

# neo sex
```
Colius_striatus
Morphnus_guianensis
Platalea_leucorodia
Poecile_atricapillus
```
# NA
```
Anser_cygnoides,NA
Guaruba_guaruba,NA
Haemorhous_mexicanus,NA
Rhynochetos_jubatus,NA
Struthio_camelus_australis,NA
Willisornis_vidua,NA
Zonotrichia_albicollis,NA
Zosterops_lateralis,NA
```
# no neo sex

```
Aegotheles_albertisi,0-3571166
Anas_platyrhynchos,0-1777881
Aythya_ferina,0-1863434
Aythya_marila,0-2273918
Calonectris_borealis,77280424-87715099
Coloeus_monedula,84644553-85378027
Columba_livia,0-1995161
Cyanocitta_cristata,0-807097
Cygnus_columbianus,0-2524550
Dixiphia_pipra,0-503193
Falco_naumanni,0-813326
Heliangelus_exortis,0-854851
Larus_argentatus,0-4531783
Lathamus_discolor,111967751-112479095
Mergus_octosetaceus,0-2607306
Numenius_arquata,0-3300171
Opisthocomus_hoazin,88690378-91338688
Passer_domesticus,0-667073
Patagioenas_fasciata,83825531-85896199
Phaethon_aethereus,0-10085609
Rissa_tridactyla,83212817-88208671
Sarcoramphus_papa,87011595-96420541
Strix_aluco,94105843-97043613
Taeniopygia_guttata,0-2824462
```
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv

chmod +x extract_z_par_genes.sh
./extract_z_par_genes.sh species_par.csv > z_par_genes.tsv

# Manually match genes to same ID / description
# Top 15 candidate PAR genes:
Appear in 9-12 species with confident PARs
```
ALPK2
MALT1
ZNF532
GRP
LMAN1
SEC11C
ST8SIA3
CPLX4
FECH
NARS1
NEDD4
ONECUT2
Rx2
TXNL1
WDR7
```
# Ones from flycatcher
```
ALPK2
MALT1
ZNF532
SEC11C
GRP
# RAX
CPLX4
LMAN1
NEDD4
# ATP8B1
NARS # NARS1 in flycatcher
FECH
ONECUT2
ST8SIA3
WDR7
TXNL1
```
# ID which of the species with PARs have these
```
base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

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
    Anser_cygnoides
    Aythya_ferina
    Aythya_marila
    Cygnus_columbianus
    Opisthocomus_hoazin
    Coloeus_monedula
    Cyanocitta_cristata
    Taeniopygia_guttata
    Haemorhous_mexicanus
    Zonotrichia_albicollis
    Passer_domesticus
    Dixiphia_pipra
    Willisornis_vidua
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
    Guaruba_guaruba
    Lathamus_discolor
)

genes=(
  ALPK2
  MALT1
  ZNF532
  GRP
  LMAN1
  SEC11C
  ST8SIA3
  CPLX4
  FECH
  NARS1
  NEDD4
  ONECUT2
  Rx2
  TXNL1
  WDR7
)

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds

echo "species,gene,chromosome" > gene_chromosomes.csv

for sp in "${species[@]}"; do
  gff="${base}/${sp}/${sp}.gff"

  if [[ ! -s "$gff" ]]; then
    echo "WARNING: missing or empty GFF: $gff" >&2
    continue
  fi

  for gene in "${genes[@]}"; do
    awk -v sp="$sp" -v gene="$gene" '
      BEGIN { FS = "\t" }
      $0 ~ ("(^|[;[:space:]])gene=" gene "([;[:space:]]|$)") ||
      $0 ~ ("(^|[;[:space:]])Name=" gene "([;[:space:]]|$)") ||
      $0 ~ ("ID=gene-" gene "([;[:space:]]|$)") {
        print sp "," gene "," $1
      }
    ' "$gff" | sort -u
  done
done >> gene_chromosomes.csv



cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/gffs

FILES_TO_DOWNLOAD="gff3,rna,cds,protein,seq-report"

while IFS=, read -r SPECIES ACCESSION; do
  [[ -z "${SPECIES}" || -z "${ACCESSION}" ]] && continue

  echo "Downloading ${SPECIES} ${ACCESSION}"

  mkdir -p "${SPECIES}"
  cd "${SPECIES}"

  datasets download genome accession "${ACCESSION}" \
    --include "${FILES_TO_DOWNLOAD}" \
    --filename "${ACCESSION}.zip"

  cd ..
done <<'EOF'
Larus_argentatus,GCF_964417175.1
Anas_platyrhynchos,GCF_047663525.1
Aythya_ferina,GCA_964211825.1
Aythya_marila,GCF_965140915.1
Cygnus_columbianus,GCF_965151615.1
Coloeus_monedula,GCF_965178545.1
Willisornis_vidua,GCA_045364795.1
Phaethon_aethereus,GCF_964289735.1
Zosterops_lateralis,GCA_965231245.1
Platalea_leucorodia,GCF_965183815.1
EOF



while IFS=, read -r SPECIES ACCESSION; do
  [[ -z "${SPECIES}" || -z "${ACCESSION}" ]] && continue

  echo "Downloading ${SPECIES} ${ACCESSION}"

  mkdir -p "${SPECIES}"
  cd "${SPECIES}"

  unzip ${ACCESSION}.zip
  cd ..
done <<'EOF'
Larus_argentatus,GCF_964417175.1
Anas_platyrhynchos,GCF_047663525.1
Aythya_ferina,GCA_964211825.1
Aythya_marila,GCF_965140915.1
Cygnus_columbianus,GCF_965151615.1
Coloeus_monedula,GCF_965178545.1
Willisornis_vidua,GCA_045364795.1
Phaethon_aethereus,GCF_964289735.1
Zosterops_lateralis,GCA_965231245.1
Platalea_leucorodia,GCF_965183815.1
EOF


while IFS=, read -r SPECIES ACCESSION; do
  [[ -z "${SPECIES}" || -z "${ACCESSION}" ]] && continue
  ln -sf /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/gffs/${SPECIES}/ncbi_dataset/data/${ACCESSION}/genomic.gff symlinks/${SPECIES}.gff
done <<'EOF'
Larus_argentatus,GCF_964417175.1
Anas_platyrhynchos,GCF_047663525.1
Aythya_ferina,GCA_964211825.1
Aythya_marila,GCF_965140915.1a
Cygnus_columbianus,GCF_965151615.1
Coloeus_monedula,GCF_965178545.1
Willisornis_vidua,GCA_045364795.1
Phaethon_aethereus,GCF_964289735.1
Zosterops_lateralis,GCA_965231245.1
Platalea_leucorodia,GCF_965183815.1
EOF
```
# Manually curate a few
```
ln -sf /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/bAytFer1/bAytFer1.gff symlinks/Aythya_ferina.gff
ln -sf /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/bAytMar2/bAytMar2.gff symlinks/Aythya_marila.gff
ln -sf /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/bWilVid1/bWilVid1.gff symlinks/Willisornis_vidua.gff
ln -sf /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/bZosLat1/bZosLat1.gff symlinks/Zosterops_lateralis.gff
```

# Redo 
```
base="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/gffs/symlinks"

species=(
    Larus_argentatus
    Anas_platyrhynchos
    Aythya_ferina
    Aythya_marila
    Cygnus_columbianus
    Coloeus_monedula
    Willisornis_vidua
    Phaethon_aethereus
    Zosterops_lateralis
    Platalea_leucorodia
)

genes=(
  ALPK2
  MALT1
  ZNF532
  GRP
  LMAN1
  SEC11C
  ST8SIA3
  CPLX4
  FECH
  NARS1
  NEDD4
  ONECUT2
  Rx2
  TXNL1
  WDR7
)

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds


for sp in "${species[@]}"; do
  gff="${base}/${sp}.gff"

  if [[ ! -s "$gff" ]]; then
    echo "WARNING: missing or empty GFF: $gff" >&2
    continue
  fi

  for gene in "${genes[@]}"; do
    awk -v sp="$sp" -v gene="$gene" '
      BEGIN { FS = "\t" }
      $0 ~ ("(^|[;[:space:]])gene=" gene "([;[:space:]]|$)") ||
      $0 ~ ("(^|[;[:space:]])Name=" gene "([;[:space:]]|$)") ||
      $0 ~ ("ID=gene-" gene "([;[:space:]]|$)") {
        print sp "," gene "," $1
      }
    ' "$gff" | sort -u
  done
done >> gene_chromosomes.csv
```
# Blast PAR genes
## Curate fastas of each PAR gene from the T2T goose
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/T2T_goose

ALPK2;NC_089912.1:61791880-61881527
MALT1;NC_089912.1:61881580-61920160
ZNF532;NC_089912.1:61931500-61977720
GRP;NC_089912.1:62010049-62015796
LMAN1;NC_089912.1:62056284-62076307
SEC11C;NC_089912.1:61991217-61997428
ST8SIA3;NC_089912.1:61492255-61498740
CPLX4;NC_089912.1:62034152-62056048
FECH;NC_089912.1:61547325-61567295
NARS1;NC_089912.1:61569543-61582963
NEDD4;NC_089912.1:61649818-61804441
ONECUT2;NC_089912.1:61510801-61537550
RAX;NC_089912.1:62023452-62026108
TXNL1;NC_089912.1:61313084-61328909
WDR7;NC_089912.1:61330349-61444936

REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Anser_cygnoides/ncbi_dataset/data/GCF_040182565.1/GCF_040182565.1_Taihu_goose_T2T_genome_genomic.fna"


base="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/gffs/symlinks"


base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

species=(
    Mergus_octosetaceus
    Aegotheles_albertisi
    Sarcoramphus_papa
    Larus_argentatus
    Rissa_tridactyla
    Numenius_arquata
    Columba_livia
    Patagioenas_fasciata
)


cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/T2T_goose
module load samtools_1.23


REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Anser_cygnoides/ncbi_dataset/data/GCF_040182565.1/GCF_040182565.1_Taihu_goose_T2T_genome_genomic.fna"

# Index reference if needed
if [ ! -f "${REF_GENOME}.fai" ]; then
    samtools faidx "$REF_GENOME"
fi

cat > gene_regions.txt <<'EOF'
ALPK2;NC_089912.1:61791880-61881527
MALT1;NC_089912.1:61881580-61920160
ZNF532;NC_089912.1:61931500-61977720
GRP;NC_089912.1:62010049-62015796
LMAN1;NC_089912.1:62056284-62076307
SEC11C;NC_089912.1:61991217-61997428
ST8SIA3;NC_089912.1:61492255-61498740
CPLX4;NC_089912.1:62034152-62056048
FECH;NC_089912.1:61547325-61567295
NARS1;NC_089912.1:61569543-61582963
NEDD4;NC_089912.1:61649818-61804441
ONECUT2;NC_089912.1:61510801-61537550
RAX;NC_089912.1:62023452-62026108
TXNL1;NC_089912.1:61313084-61328909
WDR7;NC_089912.1:61330349-61444936
EOF

while IFS=';' read -r gene region; do
    # Skip empty or malformed lines
    if [ -z "${gene}" ] || [ -z "${region}" ]; then
        echo "Skipping malformed line: gene='${gene}' region='${region}'" >&2
        continue
    fi

    out="${gene}.T2Tgoose.fa"

    echo "Extracting ${gene}: ${region} -> ${out}"

    samtools faidx "$REF_GENOME" "$region" \
        | awk -v gene="$gene" -v region="$region" '
            NR == 1 { print ">" gene "|" region; next }
            { print }
        ' > "$out"

done < gene_regions.txt

```


# Run blast
```
module load blast

BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

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
    Anser_cygnoides
    Aythya_ferina
    Aythya_marila
    Cygnus_columbianus
    Opisthocomus_hoazin
    Coloeus_monedula
    Cyanocitta_cristata
    Taeniopygia_guttata
    Haemorhous_mexicanus
    Zonotrichia_albicollis
    Passer_domesticus
    Dixiphia_pipra
    Willisornis_vidua
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
    Guaruba_guaruba
    Lathamus_discolor
)

genes=(
    ALPK2
    MALT1
    ZNF532
    GRP
    LMAN1
    SEC11C
    ST8SIA3
    CPLX4
    FECH
    NARS1
    NEDD4
    ONECUT2
    RAX
    TXNL1
    WDR7
)

OUT="gene_locations_by_species.csv"

echo "Species,Gene,Hit_rank,Chrom,Start_pos,Stop_pos" > "$OUT"

mkdir -p blast_dbs blast_results

for sp in "${species[@]}"; do
    genome="${BASE}/${sp}/${sp}.fna"

    if [ ! -f "$genome" ]; then
        echo "WARNING: genome not found for ${sp}: ${genome}" >&2
        continue
    fi

    db="blast_dbs/${sp}"

    if [ ! -f "${db}.nin" ] && [ ! -f "${db}.00.nin" ]; then
        echo "Making BLAST DB for ${sp}"
        makeblastdb \
            -in "$genome" \
            -dbtype nucl \
            -parse_seqids \
            -out "$db"
    fi

    for gene in "${genes[@]}"; do
        query="${gene}.T2Tgoose.fa"

        if [ ! -f "$query" ]; then
            echo "WARNING: query FASTA not found for ${gene}: ${query}" >&2
            continue
        fi

        blast_out="blast_results/${sp}.${gene}.blast.tsv"

        echo "BLAST: ${sp} ${gene}"

        blastn \
            -query "$query" \
            -db "$db" \
            -out "$blast_out" \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" \
            -max_target_seqs 5 \
            -evalue 1e-20 \
            -num_threads 4

        if [ ! -s "$blast_out" ]; then
            echo "${sp},${gene},NA,NA,NA,NA" >> "$OUT"
            continue
        fi

        # Pick best hit by highest bitscore, then longest alignment.
        # Convert reverse-strand hits so Start_pos < Stop_pos.
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
        | sort -k1,1nr -k2,2nr \
        | head -n 2 \
        | awk -v sp="$sp" -v gene="$gene" '
            BEGIN {
                OFS = ","
                rank = 0
            }
            {
                rank++
                print sp, gene, rank, $3, $4, $5
            }
        ' >> "$OUT"

    done
done

echo "Done. Results written to ${OUT}"

```
# add chromosome label
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv

SEXCHR_FILE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"

LABELED_OUT="gene_locations_by_species.with_chr_label.csv"

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
' "$SEXCHR_FILE" "$OUT" > "$LABELED_OUT"

echo "Labeled results written to ${LABELED_OUT}"

NC_089873.1;NC_089874.1;NC_089875.1;NC_089876.1;NC_089883.1;NC_089892.1;NC_089899.1


grep NW_027428339.1 /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Zonotrichia_albicollis/Zonotrichia_albicollis.fna
# >NW_027428339.1 Zonotrichia albicollis isolate bZonAlb1 chromosome Z unlocalized genomic scaffold, bZonAlb1.hap1 SUPER_Z_unloc_1, whole genome shotgun sequence