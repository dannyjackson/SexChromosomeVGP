# What two genes are split by the PAR boundary in cetaceans?
cat /data/Wilson_Lab/projects/VGP_Synteny/cetaceans_lifton.txt
Lagenorhynchus_acutus   Bos_taurus
Grampus_griseus Bos_taurus
Inia_geoffrensis        Bos_taurus
Stenella_coeruleoalba   Bos_taurus
Megaptera_novaeangliae  Bos_taurus
Balaenoptera_physalus   Bos_taurus
Mesoplodon_bidens       Bos_taurus


1. Using a curated list of species, identify coordinates of PAR and the boundary genes. Species that meet the curation qualifications will have 3 key essential genomic features:
- Evidence of PAR from alignment
- Evidence of telomere on PAR end of X chromosome
- Evidence from SCINKD that the entire PAR was assembled on the X and Y chromosomes

# Identify coordinates of PAR
Curated species list:
Balaenoptera_physalus
Grampus_griseus
Pseudorca_crassidens
Ovis_canadensis

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/continuous_percentID

grep 'Ovis' surviving_regions.*

Species,PAR_location,PAR_start_qry,PAR_start_ref,PAR_end_qry,PAR_end_ref,PAR_lastgene_qry,PAR_lastgene_ref,NRR_firstgene_qry,NRR_firstgene_ref
qry_end,ref_end
Balaenoptera_physalus,left,0,0,7257636,7193118,
Grampus_griseus,left,0,0,7148355,7106705
Pseudorca_crassidens,rightX_leftY,127639906,0,134794574,7558390
Ovis_canadensis,left,0,0,7072606,7067639
# Identify boundary genes
/data/Wilson_Lab/data/VGP_genomes/
SEXCHROM_CSV="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/sexchrom_accessions.csv"

grep 'Balae' $SEXCHROM_CSV # Balaenoptera_physalus,X,OZ239531.1 -- Balaenoptera_physalus,Y,OZ239532.1
grep 'Gramp' $SEXCHROM_CSV # Grampus_griseus,X,OZ206318.1 -- Grampus_griseus,Y,OZ206335.1
grep 'Pseud' $SEXCHROM_CSV # Pseudorca_crassidens,X,NC_090317.1 -- Pseudorca_crassidens,Y,NC_090318.1
grep 'Ovis' $SEXCHROM_CSV # Ovis_canadensis,Y,NC_091271.1 -- Ovis_canadensis,X,NC_091727.1

From $GFF, output any lines that have "OZ239531.1" in column 1, "gene" in column 3, and col5 is less than 7257636
# Balaenoptera
```
GFF=/data/Wilson_Lab/data/VGP_genomes/Balaenoptera_physalus/Balaenoptera_physalus.gff
# qry
awk -F'\t' '$1=="OZ239531.1" && $3=="gene" && $5<7257636' $GFF
# ref
awk -F'\t' '$1=="OZ239532.1" && $3=="gene" && $5<7193118' $GFF
```
Qry last gene: SHROOM2 (27 lines)
Ref last gene: TBL1X_1 {SHROOM2,GYG2_1,MXRA5,TBL1X}

ID=gene-SHROOM2_1;Dbxref=GeneID:113886919;Name=SHROOM2;gbkey=Gene;gene=SHROOM2;gene_biotype=protein_coding;extra_copy_number=1;source=miniprot

# Grampus
```
GFF=/data/Wilson_Lab/data/VGP_genomes/Grampus_griseus/Grampus_griseus.gff
# qry
awk -F'\t' '$1=="OZ206318.1" && $3=="gene" && $5<7148355' $GFF 
# ref 
awk -F'\t' '$1=="OZ206335.1" && $3=="gene" && $5<7106705' $GFF

```
Qry last gene: SHROOM2 (28 lines)
Ref last gene: ANOS1 {SHROOM2,CRLF2,IL3RA,P2RY8_1,AKAP17A,ZBED1_1,CD99_1,ANOS_1}

# Pseudorca_crassidens
```
# Pseudorca_crassidens,rightX_leftY,127639906,0,134794574,7558390

GFF=/data/Wilson_Lab/data/VGP_genomes/Pseudorca_crassidens/Pseudorca_crassidens.gff
# qry
awk -F'\t' '$1=="NC_090317.1" && $3=="gene" && $5>127639906' $GFF
# ref 
awk -F'\t' '$1=="NC_090318.1" && $3=="gene" && $5<7558390' $GFF

```
Qry last gene: SHROOM2 (53 lines)
Ref last gene: Shroom2-like


Compare GRAMPUS and BALAEN gene lists
# GRAMPUS
GFF=/data/Wilson_Lab/data/VGP_genomes/Grampus_griseus/Grampus_griseus.gff

awk -F'\t' '
$1=="OZ206318.1" && $3=="gene" && $5<7148355 && $9 ~ /gene_biotype=protein_coding/ {
    match($9, /ID=gene-([^;]+)/, a)
    if (a[1] != "") print a[1]
}' "$GFF"

# BALAEN
GFF=/data/Wilson_Lab/data/VGP_genomes/Balaenoptera_physalus/Balaenoptera_physalus.gff

awk -F'\t' '
$1=="OZ239531.1" && $3=="gene" && $5<7257636 && $9 ~ /gene_biotype=protein_coding/ {
    match($9, /ID=gene-([^;]+)/, a)
    if (a[1] != "") print a[1]
}' "$GFF"

# Ovis canadensis
# grep 'Ovis' $SEXCHROM_CSV # Ovis_canadensis,Y,NC_091271.1 -- Ovis_canadensis,X,NC_091727.1
# Ovis_canadensis,left,0,0,7072606,7067639

GFF=/data/Wilson_Lab/data/VGP_genomes/Ovis_canadensis/Ovis_canadensis.gff

awk -F'\t' '
$1=="NC_091727.1" && $3=="gene" && $5<7072606 && $9 ~ /gene_biotype=protein_coding/ {
    match($9, /ID=gene-([^;]+)/, a)
    if (a[1] != "") print a[1]
}' "$GFF"


# Capra hircus
# grep 'Capra' $SEXCHROM_CSV # Capra_hircus,X,CP168640.1 -- Capra_hircus,Y,CP168641.1
# Capra_hircus,left,0,0,7019956,7040093

GFF=/data/Wilson_Lab/data/VGP_genomes/Capra_hircus/Capra_hircus.gff

awk -F'\t' '
$1=="CP168640.1" && $3=="gene" && $5<7040093 && $9 ~ /gene_biotype=protein_coding/ {
    match($9, /ID=gene-([^;]+)/, a)
    if (a[1] != "") print a[1]
}' "$GFF"


# Pseudorca crassidens
# grep 'Pseudorca' $SEXCHROM_CSV # Pseudorca_crassidens,X,NC_090317.1 -- Pseudorca_crassidens,Y,NC_090318.1
# Capra_hircus,left,0,127639906,7558390,end

GFF=/data/Wilson_Lab/data/VGP_genomes/Pseudorca_crassidens/Pseudorca_crassidens.gff

awk -F'\t' '
$1=="NC_090317.1" && $3=="gene" && $5>127639906 && $9 ~ /gene_biotype=protein_coding/ {
    match($9, /ID=gene-([^;]+)/, a)
    if (a[1] != "") print a[1]
}' "$GFF"


awk -F'\t' '
$1=="NC_090317.1" && $3=="gene" && $5>127639906 && $9 ~ /gene_biotype=protein_coding/ {
    match($9, /description=([^;]+)/, a)
    if (a[1] != "") print a[1]
}' "$GFF"