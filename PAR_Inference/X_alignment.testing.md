# modifying alignment parameters until I recover the known PAR boundaries of human genome

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13
~/programs/datasets download genome accession GCF_009914755.1 --include gff3,rna,cds,protein,genome,seq-report

module load samtools
samtools faidx ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna NC_060947.1 > Xchr.fa
samtools faidx ncbi_dataset/data/GCF_009914755.1/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna NC_060948.1 > Ychr.fa


XCHR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13/Xchr.fa
YCHR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13/Ychr.fa


cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment

module load minimap2/2.30

#!/bin/bash

# Parse arguments
while getopts "k:" option; do
  case "${option}" in
    k) K="${OPTARG}" ;;
    *) echo "Usage: $0 [-k KMER]" >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

module load minimap2/2.30

XCHR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13/Xchr.fa
YCHR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/genomes/Homo_sapiens_CHM13/Ychr.fa

minimap2 -x asm5 -c --eqx  -t 3 ${YCHR} ${XCHR} -k ${K} \
 > PAFs/Homo_sapiens_Y_to_X.aln.k${K}.paf


sbatch \
  -c 1 \
  -t 0:30:00 \
  --mem-per-cpu=8G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  minimap_submit.sh -k 19

sbatch \
  -c 1 \
  -t 0:30:00 \
  --mem-per-cpu=8G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  minimap_submit.sh -k 21
  

sbatch \
  -c 1 \
  -t 0:30:00 \
  --mem-per-cpu=8G \
  --mail-user=jacksondan@nih.gov \
  --mail-type=ALL \
  --gres=lscratch:500 \
  minimap_submit.sh -k 23





# Filter and write BED:
# Output BED columns:
# ref_id  ref_start  ref_end  query_id  query_start  query_end  identity_percent  ref_span_len
awk -v thr=0.995 -v minlen=10000 -v OFS='\t' '
  function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
  function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
  {
    # PAF fields:
    # 1 qname 2 qlen 3 qstart 4 qend 5 strand 6 tname 7 tlen 8 tstart 9 tend 10 nmatch 11 alnlen 12 mapq
    tspan = $9 - $8
    de=""; for(i=13;i<=NF;i++) if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }
    id = ident_from_de(de); if (id < 0) id = ident_from_core($10,$11)

    if (id >= thr && tspan >= minlen)
      printf "%s\t%d\t%d\t%s\t%d\t%d\t%.4f\t%d\n", $6, $8, $9, $1, $3, $4, id*100, tspan
  }
' "PAFs/Homo_sapiens_Y_to_X.aln.k23.paf" > BEDs/Homo_sapiens_Y_to_X.aln.k23.id98_5.len10k.refqry.bed

TE_BEDDIR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/TE_annotations
Y_TE_BED=${TE_BEDDIR}/sex_chrs/Homo_sapiens.GCA_009914755.4.TE.Y.bed

sed -i 's/NC_060948.1/chrY/g' BEDs/Homo_sapiens_Y_to_X.aln.k23.id98_5.len10k.refqry.bed
Rscript 2x_par_align_TE_plot.R CHM13v2.0_PAR.chrY.bed ${Y_TE_BED} BEDs/Homo_sapiens_Y_to_X.aln.k23.id98_5.len10k.refqry.bed 62460029 Y_par_align_TE.k23.thr995.png

# thr 0.99
awk -v thr=0.99 -v minlen=10000 -v OFS='\t' '
  function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
  function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
  {
    # PAF fields:
    # 1 qname 2 qlen 3 qstart 4 qend 5 strand 6 tname 7 tlen 8 tstart 9 tend 10 nmatch 11 alnlen 12 mapq
    tspan = $9 - $8
    de=""; for(i=13;i<=NF;i++) if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }
    id = ident_from_de(de); if (id < 0) id = ident_from_core($10,$11)

    if (id >= thr && tspan >= minlen)
      printf "%s\t%d\t%d\t%s\t%d\t%d\t%.4f\t%d\n", $6, $8, $9, $1, $3, $4, id*100, tspan
  }
' "PAFs/Homo_sapiens_Y_to_X.aln.k23.paf" > BEDs/Homo_sapiens_Y_to_X.aln.k23.id98_5.len10k.refqry.bed

TE_BEDDIR=/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/TE_annotations
Y_TE_BED=${TE_BEDDIR}/sex_chrs/Homo_sapiens.GCA_009914755.4.TE.Y.bed

sed -i 's/NC_060948.1/chrY/g' BEDs/Homo_sapiens_Y_to_X.aln.k23.id98_5.len10k.refqry.bed
Rscript 2x_par_align_TE_plot.R CHM13v2.0_PAR.chrY.bed ${Y_TE_BED} BEDs/Homo_sapiens_Y_to_X.aln.k23.id98_5.len10k.refqry.bed 62460029 Y_par_align_TE.k23.thr99.png


# Modifying Simone's threshold cutoff from 0.985 to 0.99 recoups the boundary almost exactly. There are now more gaps, but the end bp of region 1 is identical to the boundary on the CHM13 github (2458320 on Y; 2394410 on X). I can rerun this really easily on all the pafs on Sol really easily so hopefully the results look more clear afterwards
