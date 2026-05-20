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
Reduce the output to PAR regions, integrating across multiple blocks of inferred PARs.
```
awk -F'[,\t]' 'NR > 1 {
  if (!($1 in min) || $2 < min[$1]) min[$1] = $2
  if (!($1 in max) || $3 > max[$1]) max[$1] = $3
}
END {
  for (s in min) print s "," min[s] "-" max[s]
}' /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_inference/alignment/continuous_percentID/surviving_regions.qry.thr98p5.len10000.ALLFILES.csv > ../InferredParBoundaries.csv
```
Subset the above list to just the species of interest:
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
Manually inspect the output for reasonable PAR inferences. Confirm that the PAR boundaries align with the expectations based on the png output of the R script above.

### Species with evidence of a neo sex chromosome (exclude)
```
Colius_striatus
Guaruba_guaruba
Lathamus_discolor
Morphnus_guianensis
Platalea_leucorodia
Poecile_atricapillus
```
### Species with no PAR identified by alignment; exclude
```
Anser_cygnoides,NA
Haemorhous_mexicanus,NA
Rhynochetos_jubatus,NA
Willisornis_vidua,NA
Zonotrichia_albicollis,NA
Zosterops_lateralis,NA
```
### Final species list (save as species_par.csv)
PARs can only exist on the end of a chromosome, and so I manually edited each boundary to extend to the end of the chromosome. If the inferred PAR is at the start of the chromosome, I edited the minimum position to be 0. If it is at the end of the chromosome, I edited the maximum position to match the total number of base pairs on the Z chromosome. I did this manually by looking up the genome on NCBI and finding the length of the Z. I also removed Struthio_camelus_australis from this step, as it is known to have a substantially larger PAR than other birds.

I saved this file as ```species_par.csv```:
```
species,PAR,Z_size,W_size
Aegotheles_albertisi,0-3571166,84215446,24988134
Anas_platyrhynchos,0-1777881,85226422,19630300
Aythya_ferina,0-1863434,85915276,19376790
Aythya_marila,0-2273918,89030289,19204845
Calonectris_borealis,77280424-87715099,87715099,56726463
Coloeus_monedula,84644553-85378027,85378027,34140743
Columba_livia,0-1995161,84824678,23867768
Cyanocitta_cristata,0-807097,84313221,37504850
Cygnus_columbianus,0-2524550,88185588,22579564
Dixiphia_pipra,0-503193,75011098,66133916
Falco_naumanni,0-813326,86597978,29619203
Heliangelus_exortis,0-854851,73898144,14480854
Larus_argentatus,0-4531783,88668981,35824208
Mergus_octosetaceus,0-2607306,87928280,16130779
Numenius_arquata,0-3300171,85602926,35829655
Opisthocomus_hoazin,88690378-91338688,91338688,50397613
Passer_domesticus,0-667073,84808660,47418805
Patagioenas_fasciata,83825531-85896199,85896199,70564437
Phaethon_aethereus,0-10085609,90167876,43195556
Rissa_tridactyla,83212817-88208671,88208671,31882723
Sarcoramphus_papa,87011595-96420541,96420541,49671510
Strix_aluco,94105843-97043613,97043613,49279228
Taeniopygia_guttata,0-2824462,82218670,33777440
```
## 4. Extract all genes within the PAR regions 
Extract genes within the inferred PAR for these species using the GFF for that species (found in /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv)
```
chmod +x extract_z_par_genes.sh
./extract_z_par_genes.sh species_par.csv > z_par_genes.tsv

# Used this in mammals
./extract_x_par_genes.sh species_par.csv > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/nora/x_par_genes.tsv

```
## 5. Manually match genes to same ID / description
Gene names in these GFFs are not consistent across species. The above script outputs both gene name and gene description. I downloaded the CSV to my local disk and then used the gene descriptions to recode genes to similar name. 

I googled the description and find the gene name, using sites that I trust like Gene Cards or UniProt. Afterwards, I put the entire list of GeneName, GeneDescription, and GeneName_Revised columns into ChatGPT and asked it to double check my conversions. ChatGPT caught a few that I'd accidentally labelled using the protein name instead of the gene name, but it also gave incorrect gene names for some as well so be sure to thoroughly check Chat's work.

For example, I added a new column to z_par_genes.tsv named "GeneName_Revised" and recoded these variously named genes to ACAA2.

| Species | Chromosome | StartPos | StopPos | GeneName | GeneDescription | GeneName_Revised |
|---|---|---:|---:|---|---|---|
| Aegotheles_albertisi | CM078494.1 | 2445494 | 2458490 | AAHN32_012953 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Rissa_tridactyla | NC_071497.1 | 84859140 | 84877649 | LOC128902953 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Columba_livia | NC_088642.1 | 1921455 | 1937181 | LOC135577390 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Patagioenas_fasciata | NC_092560.1 | 83911857 | 83925983 | LOC136114726 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Mergus_octosetaceus | CM072318.1 | 2435566 | 2447305 | V3H86_014483 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |
| Sarcoramphus_papa | CM075626.1 | 93109717 | 93131945 | WM294_016710 | 3-ketoacyl-CoA thiolase%2C mitochondrial-like | ACAA2 |

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
## 8. Use the list of PAR genes to identify the location of the ancestral avian PAR in our short list of species

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
### Add a chromosome label for each blast match. 
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

# Repeat for all 
```
The output was analyzed manually in excel to identify potential PAR misassemblies.

### Case studies: 
Two species (Rhynochetos_jubatus, Zosterops_lateralis) have all PAR genes on Z and W but no surviving regions
Rhynochetos_jubatus: Min pos of a PAR gene is 84475080, max is 85065248
Zosterops_lateralis: Min pos of a PAR gene is 110684, max is 481990

#### Rhynochetos_jubatus
```
PAF=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/Rhynochetos_jubatus_WtoZ.aln.paf
awk -v chr="CM050690.1" -v a=84475080 -v b=85065248 '
BEGIN {
    start = (a < b ? a : b)
    end   = (a > b ? a : b)
}
$1 == chr && $3 >= start && $4 <= end {
    pid = ($10 / $11) * 100
    block = $11

    sum_pid += pid
    sum_block += block
    n++

    if (n == 1 || block > max_block) max_block = block
    if (n == 1 || pid > max_pid) max_pid = pid
}
END {
    if (n > 0) {
        printf "n_alignments\t%d\n", n
        printf "avg_percent_id\t%.2f\n", sum_pid / n
        printf "avg_block_size\t%.2f\n", sum_block / n
        printf "max_block_size\t%d\n", max_block
        printf "max_percent_id\t%.2f\n", max_pid
    } else {
        print "No alignments found in region"
    }
}' "$PAF"
```
Results:
n_alignments    5
avg_percent_id  97.08
avg_block_size  786.60
max_block_size  2522
max_percent_id  100.00

#### Zosterops_lateralis:Chr = OZ246513.1 Min pos of a PAR gene is 110684, max is 481990
Missing the alignment file...
```
PAF=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/Zosterops_lateralis_WtoZ.aln.paf
awk -v chr="OZ246513.1" -v a=110684 -v b=481990 '
BEGIN {
    start = (a < b ? a : b)
    end   = (a > b ? a : b)
}
$1 == chr && $3 >= start && $4 <= end {
    pid = ($10 / $11) * 100
    block = $11

    sum_pid += pid
    sum_block += block
    n++

    if (n == 1 || block > max_block) max_block = block
    if (n == 1 || pid > max_pid) max_pid = pid
}
END {
    if (n > 0) {
        printf "n_alignments\t%d\n", n
        printf "avg_percent_id\t%.2f\n", sum_pid / n
        printf "avg_block_size\t%.2f\n", sum_block / n
        printf "max_block_size\t%d\n", max_block
        printf "max_percent_id\t%.2f\n", max_pid
    } else {
        print "No alignments found in region"
    }
}' "$PAF"
```
Results: No alignments found

