# ID genomes with well-assembled PARs
This script follows the curation of genomes that show no evidence of PAR dropout in HalfDeep plots. HalfDeep identified patterns associated with PAR assembly on only one of the two sex chromosomes. However, other forms of PAR misassembly would be missed by the HalfDeep assessment. We need to check if the ancestral PAR is found on the end of both of the Z and W chromosomes.

To do this, we will first identify genes conserved in the PARs across multiple genomes that pass the HalfDeep filter. We do not know what proportion of these genomes have misassemblies, so this analysis uses cutoffs that are inferred from the data itself, not predetermined.

## 1. Set up the environment
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds

sinteractive --mem=20g

module load bcftools vcftools samtools R blast
```
## 2. Plot regions of sequence alignment of sex chromosomes (a method for inferring the boundaries of the PAR)
This step produces files in the directory ```/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_inference/alignment/continuous_percentID```. 

I used the output file surviving_regions.qry.thr98p5.len10000.ALLFILES.csv for all downstream analyses. The file name with "ref" instead of "qry" lists PAR boundaries in the W/Y chromosome. 

The thr98p5 section indicates that I used a threshold of 98.5% identity as the cutoff between alignments for potential PAR identification, and the len10000 section indicates that I used a cutoff of 10kb size for any matching region between the alignments. 

I tested various other thresholds including 50% ID and 0kb and found that the PAR boundaries do not change, but small regions throughout the chromosomes match and create noise. 
```
Rscript PAR_Plotting_AllSpecies.R
```
## 3. Curate a list of genes found in the PAR regions across all genomes that may have a PAR assembled on both sex chromosomes
Use the output csv of surviving regions to identify PAR boundaries in the species with half-deep inferred no-PAR-dropout. Then, using only genomes without a neo-sex fusion, identify genes within the PAR for each species.

First, subset the list of species that survived the HalfDeep filter to just those that do not have a neo sex chromosome AND that have some evidence of a PAR in the alignment file.
```
awk -F'[,\t]' 'NR > 1 {
  if (!($1 in min) || $2 < min[$1]) min[$1] = $2
  if (!($1 in max) || $3 > max[$1]) max[$1] = $3
}
END {
  for (s in min) print s "," min[s] "-" max[s]
}' /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_inference/alignment/continuous_percentID/surviving_regions.qry.thr98p5.len50000.ALLFILES.csv > ../InferredParBoundaries.csv
```
### Subset the above list to just the species of interest:
```
awk -F',' '
ARGIND == 1 {
  keep[$1] = 1
  next
}
ARGIND == 2 && $1 in keep
' - ../InferredParBoundaries.csv <<'EOF' > InferredParBoundaries.Birds.SurvivedHalfDeep.csv
Aegotheles_albertisi
Anas_platyrhynchos
Anser_cygnoides
Aythya_ferina
Aythya_marila
Calonectris_borealis
Colius_striatus
Coloeus_monedula
Columba_livia
Cyanocitta_cristata
Cygnus_columbianus
Dixiphia_pipra
Falco_naumanni
Guaruba_guaruba
Haemorhous_mexicanus
Heliangelus_exortis
Larus_argentatus
Lathamus_discolor
Mergus_octosetaceus
Morphnus_guianensis
Numenius_arquata
Opisthocomus_hoazin
Passer_domesticus
Patagioenas_fasciata
Phaethon_aethereus
Platalea_leucorodia
Poecile_atricapillus
Rhynochetos_jubatus
Rissa_tridactyla
Sarcoramphus_papa
Strix_aluco
Struthio_camelus_australis
Taeniopygia_guttata
Willisornis_vidua
Zonotrichia_albicollis
Zosterops_lateralis
EOF
```
### Species with evidence of a neo sex chromosome (exclude)
```
Colius_striatus
Morphnus_guianensis
Platalea_leucorodia
Poecile_atricapillus
```
### Species with no PAR identified by alignment; exclude
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
### Final species list (save as species_par.csv)
PARs can only exist on the end of a chromosome, and so I manually edited each boundary to extend to the end of the chromosome. If the inferred PAR is at the start of the chromosome, I edited the minimum position to be 0. If it is at the end of the chromosome, I edited the maximum position to match the total number of base pairs on the Z chromosome. I did this manually by looking up the genome on NCBI and finding the length of the Z.

I saved this file as ```species_par.csv```:
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
## 4. Extract all genes within the PAR regions 
Extract genes within the inferred PAR for these species using the GFF for that species (found in /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv)
```
chmod +x extract_z_par_genes.sh
./extract_z_par_genes.sh species_par.csv > z_par_genes.tsv
```
## 5. Manually match genes to same ID / description
Gene names in these GFFs are not consistent across species. The above script outputs both gene name and gene description. I downloaded the CSV to my local disk and then used the gene descriptions to recode genes to similar name. 

I googled the description and find the gene name, using sites that I trust like Gene Cards or UniProt. Afterwards, I put the entire list of GeneName, GeneDescription, and GeneName_Revised columns into ChatGPT and asked it to double check my conversions. ChatGPT caught a few that I'd accidentally labelled using the protein name instead of the gene name, but it also gave incorrect gene names for some as well so be sure to thoroughly check Chat's work.

For example, I added a new column to z_par_genes.tsv named "GeneName_Revised" and recoded these variously named genes to ACAA2.
```
| Species | Chromosome | StartPos | StopPos | GeneName | GeneDescription | GeneName_Revised |
|---|---|---:|---:|---|---|---|
| Aegotheles_albertisi | CM078494.1 | 2445494 | 2458490 | AAHN32_012953 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Rissa_tridactyla | NC_071497.1 | 84859140 | 84877649 | LOC128902953 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Columba_livia | NC_088642.1 | 1921455 | 1937181 | LOC135577390 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Patagioenas_fasciata | NC_092560.1 | 83911857 | 83925983 | LOC136114726 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Mergus_octosetaceus | CM072318.1 | 2435566 | 2447305 | V3H86_014483 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Sarcoramphus_papa | CM075626.1 | 93109717 | 93131945 | WM294_016710 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
```
## 6. Curate a list of putatively conserved PAR genes (top 15):
These genes each appear in 9-12 species with confident PARs. This cutoff is arbitrary, but I figured 15 conserved genes would be sufficient to identify the PAR region in a genome if it is there.
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
## 7. Compare the list to the extant literature
The avian PAR is much less well understood than the mammalian PAR(s). I found a study that identified PAR genes in the collared flycatcher (Smeds et al 2014) and I compared my list to theirs.

Citation:

Smeds, L., Kawakami, T., Burri, R., Bolivar, P., Husby, A., Qvarnström, A., ... & Ellegren, H. (2014). Genomic identification and characterization of the pseudoautosomal region in highly differentiated avian sex chromosomes. Nature communications, 5(1), 5448.

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
## 8. Use the list of PAR genes to identify the list of the ancestral avian PAR in our short list of species

First, curate fastas of each PAR gene from the T2T zebra finch genome. To do this, I need to generate a file of the chromosome and bp boundaries for each gene. I created this manually by looking up the genes by their descriptions on the annotations on NCBI.

```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/PAR_fastas

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/PAR_fastas

# Note that CPLX4 comes from the W not the Z, because it is not annotated on the Z.

cat > gene_regions.txt <<'EOF'
ALPK2;NC_133063.1:2717251-2744194
MALT1;NC_133063.1:2697857-2716479
ZNF532;NC_133063.1:2652873-2686892
GRP;NC_133063.1:2634400-2639278
LMAN1;NC_133063.1:2596856-2621284
SEC11C;NC_133063.1:2643018-2647109
ST8SIA3;NC_133063.1:2933267-2940220
CPLX4;NC_133064.1:2498281-2506420
FECH;NC_133063.1:2897203-2908590
NARS1;NC_133063.1:2888371-2895562
NEDD4;NC_133063.1:2767383-2853173
ONECUT2;NC_133063.1:2921246-2930515
RAX;NC_133063.1:2627901-2629552
WDR7;NC_133063.1:2963348-3005220
EOF

# pull TXNL1 from Colius striatus, no matches in the zebra finch
## Colius striatus was chosen only because I found this gene in this genome previously during a manual search

echo 'TXNL1;NC_084790.1:87870999-87891810' > TXNL1_regions.txt
```
Next, use these gene boundaries to create a fasta file for each gene of interest from the respective source genome (zebra finch for most, Colius striatus for TXNL1)
```
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Taeniopygia_guttata/ncbi_dataset/data/GCF_048771995.1/GCF_048771995.1_bTaeGut7.mat_genomic.fna"

base="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/gffs/symlinks"

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

done < gene_regions.txt
```
Repeat with TXNL1 from Colius striatus
```
REF_GENOME="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Colius_striatus/ncbi_dataset/data/GCF_028858725.1/GCF_028858725.1_bColStr4.1.hap1_genomic.fna"

base="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_Gene_ID/birds/gffs/symlinks"

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

done < TXNL1_regions.txt

```
Run blast to search each of the bird genomes that survived the HalfDeep filtering step for the putative PAR genes.
```
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
        query="${gene}.fa"

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
# Add a chromosome label for each blast match. 
For each blast match, this column will have a value of Z, W, or NA. A value of NA here indicates that the chromosome with the PAR genes is not a sex chromosome. Some manual curation of this column will be necessary because I made the sexchrom_accessions file a while ago.

```
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
```
The output was analyzed manually in excel to identify potential PAR misassemblies.









######################################################
# Draft code

## I used this code in step 8 to attempt to find PAR genes in the gffs. The blast approach is much more effective.

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