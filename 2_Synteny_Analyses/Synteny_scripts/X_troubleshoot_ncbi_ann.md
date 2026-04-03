# Prepare genespace
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting
## Have bed / peptide files, broke on genespace
```
Tursiops_truncatus # error
Hoplias_malabaricus # done
Chlamydotis_macqueenii # done
Girardinichthys_multiradiatus # error
Amblyraja_radiata # error
```
## No bed / peptide files
```
Panthera_onca
Pipistrellus_nathusii
Lemur_catta
```
#####################################################################
# Is a specific chromosome causing the error?
#####################################################################
# Processed genome
sinteractive --gres=lscratch:500 --mem=20g
# cd /lscratch/$SLURM_JOBID
export TMPDIR=/lscratch/$SLURM_JOBID

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/chromosomes//

SPECIES=Tursiops_truncatus
source myconda
mamba activate genespace_py3.10

mkdir -p genomes
cp -r /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/${SPECIES}/ ./genomes

## Remove scaffolds from FNA files
f="genomes/${SPECIES}/${SPECIES}.fna"
out="genomes/${SPECIES}.contigs_only.fna"
awk 'BEGIN{IGNORECASE=1}
        /^>/{
        drop = ($0 ~ /(scaffold|unplaced|unlocalized|unloc|mitochon)/)
        }
        !drop{print}
    ' "$f" > "$out"
echo "$f -> $out"


OUT_SUFFIX=".contigs_only.gff"

gff="genomes/${SPECIES}/${SPECIES}.gff"
fasta="genomes/${SPECIES}.contigs_only.fna"
outdir="genomes/${SPECIES}/by_chrom"

mkdir -p "$outdir"

# Split FASTA into one file per chromosome/contig
awk -v outdir="$outdir" '
BEGIN { file="" }
/^>/ {
    hdr = $0
    sub(/^>/, "", hdr)
    sub(/[ \t].*$/, "", hdr)   # keep only first token as seqid
    file = outdir "/" hdr ".fna"
    print $0 > file
    close(file)
    next
}
{
    if (file != "") print $0 >> file
}
' "$fasta"

# For each chromosome FASTA, make matching chromosome GFF
for chrom_fa in "$outdir"/*.fna; do
    chrom=$(basename "$chrom_fa" .fna)
    chrom_gff="$outdir/${chrom}${OUT_SUFFIX}"

    awk -v chrom="$chrom" '
    BEGIN { FS=OFS="\t" }
    /^#/ { print; next }
    $1 == chrom { print }
    ' "$gff" > "$chrom_gff"

    echo "OK: ${SPECIES} -> ${chrom_fa} and ${chrom_gff}"
done

mkdir NC_047034
cd NC_047034


cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/chromosomes/genomes/Tursiops_truncatus/by_chrom/NC_047034.1.contigs_only.gff Tursiops_truncatus.gff
cp /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/chromosomes/genomes/Tursiops_truncatus/by_chrom/NC_047034.1.fna Tursiops_truncatus.fna
cp ../../../raw_genome/genomes/Tursiops_truncatus/Tursiops_truncatus.translated.cds .


library(GENESPACE)

SPECIES <- "Tursiops_truncatus"

parsedPaths <- parse_annotations(
  rawGenomeRepo = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/chromosomes/NC_047034/",
  genomeDirs    = SPECIES,
  genomeIDs     = SPECIES,
  gffString     = "gff",
  faString      = "translated.cds",
  presets       = "ncbi",
  genespaceWd   = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/chromosomes/NC_047034/Tursiops_truncatus/"
)

q()
n

cp ../../../processed_genome/Hoplias_malabaricus/bed/Gallus_gallus.bed ./bed/
cp ../../../processed_genome/Hoplias_malabaricus/peptide/Gallus_gallus.fa ./peptide/








OUT_SUFFIX=".contigs_only.gff"   # output name: ${SPECIES}/${SPECIES}${OUT_SUFFIX}

gff="genomes/${SPECIES}/${SPECIES}.gff"
out="${SPECIES}${OUT_SUFFIX}"
fasta="${SPECIES}.contigs_only.fna"

awk '
BEGIN { FS=OFS="\t" }

## File 1: FASTA -> allowed seqids
FNR==NR {
    if ($0 ~ /^>/) {
    h=$0
    sub(/^>/, "", h)
    sub(/[ \t].*$/, "", h)   # first token only
    keep[h]=1
    }
    next
}

## File 2: GFF -> keep header/comments, and features whose seqid is allowed
/^#/ { print; next }
($1 in keep) { print }
' "$fasta" "$gff" > "$out"

echo "OK: ${SPECIES} -> $out"

# Rename files
mv ${SPECIES}.contigs_only.fna ${SPECIES}.fna
mv ${SPECIES}.contigs_only.gff ${SPECIES}.gff

cds_in=genomes/${SPECIES}/${SPECIES}.cds
cds_out=genomes/${SPECIES}/${SPECIES}.translated.cds
transeq -sequence "${cds_in}" -outseq "${cds_out}"


mkdir -p genomes/${SPECIES}
mv ${SPECIES}* genomes/${SPECIES}/

mamba activate genespace_py3.10

library(GENESPACE)

SPECIES <- "Tursiops_truncatus"

parsedPaths <- parse_annotations(
  rawGenomeRepo = "genomes",
  genomeDirs    = SPECIES,
  genomeIDs     = SPECIES,
  gffString     = "gff",
  faString      = "translated.cds",
  presets       = "ncbi",
  genespaceWd   = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/processed_genome/"
)

# cp ../raw_genome/Hoplias_malabaricus/bed/Gallus_gallus.bed bed/
# cp ../raw_genome/Hoplias_malabaricus/peptide/Gallus_gallus.fa peptide/

mkdir ${SPECIES}
mv bed/ ${SPECIES}
mv peptide/ ${SPECIES}
sbatch genespace.troubleshoot.processed_genome.sh Homo_sapiens ${SPECIES}










#####################################################################
# Prepare each genome and run genespace
#####################################################################

# Processed genome
Output:
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/processed_genome/


sinteractive --gres=lscratch:500
# cd /lscratch/$SLURM_JOBID
export TMPDIR=/lscratch/$SLURM_JOBID

SPECIES=Tursiops_truncatus
source myconda
mamba activate genespace

cp -r /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/${SPECIES}/ ./genomes

## Remove scaffolds from FNA files
f="genomes/${SPECIES}/${SPECIES}.fna"
out="${SPECIES}.contigs_only.fna"
awk 'BEGIN{IGNORECASE=1}
        /^>/{
        drop = ($0 ~ /(scaffold|unplaced|unlocalized|unloc|mitochon)/)
        }
        !drop{print}
    ' "$f" > "$out"
echo "$f -> $out"

## Filter ${SPECIES}/${SPECIES}.gff to just regions kept in file ${SPECIES}/${SPECIES}.contigs_only.fna

OUT_SUFFIX=".contigs_only.gff"   # output name: ${SPECIES}/${SPECIES}${OUT_SUFFIX}

gff="genomes/${SPECIES}/${SPECIES}.gff"
out="${SPECIES}${OUT_SUFFIX}"
fasta="${SPECIES}.contigs_only.fna"

awk '
BEGIN { FS=OFS="\t" }

## File 1: FASTA -> allowed seqids
FNR==NR {
    if ($0 ~ /^>/) {
    h=$0
    sub(/^>/, "", h)
    sub(/[ \t].*$/, "", h)   # first token only
    keep[h]=1
    }
    next
}

## File 2: GFF -> keep header/comments, and features whose seqid is allowed
/^#/ { print; next }
($1 in keep) { print }
' "$fasta" "$gff" > "$out"

echo "OK: ${SPECIES} -> $out"

# Rename files
mv ${SPECIES}.contigs_only.fna ${SPECIES}.fna
mv ${SPECIES}.contigs_only.gff ${SPECIES}.gff

cds_in=genomes/${SPECIES}/${SPECIES}.cds
cds_out=genomes/${SPECIES}/${SPECIES}.translated.cds
transeq -sequence "${cds_in}" -outseq "${cds_out}"


mkdir -p genomes/${SPECIES}
mv ${SPECIES}* genomes/${SPECIES}/

mamba activate genespace_py3.10

library(GENESPACE)

SPECIES <- "Tursiops_truncatus"

parsedPaths <- parse_annotations(
  rawGenomeRepo = "genomes",
  genomeDirs    = SPECIES,
  genomeIDs     = SPECIES,
  gffString     = "gff",
  faString      = "translated.cds",
  presets       = "ncbi",
  genespaceWd   = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/processed_genome/"
)

# cp ../raw_genome/Hoplias_malabaricus/bed/Gallus_gallus.bed bed/
# cp ../raw_genome/Hoplias_malabaricus/peptide/Gallus_gallus.fa peptide/

mkdir ${SPECIES}
mv bed/ ${SPECIES}
mv peptide/ ${SPECIES}
sbatch genespace.troubleshoot.processed_genome.sh Homo_sapiens ${SPECIES}




# try this
## Works for Chlamydotis_macqueenii

parsedPaths2 <- parse_annotations(
  rawGenomeRepo = "genomes",
  genomeDirs = SPECIES,
  genomeIDs  = SPECIES,
  gffString  = "gff",
  faString   = "translated.cds",
  genespaceWd   = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/processed_genome/",
  gffIdColumn = "protein_id",
  headerSep = "protein_id=",
  headerEntryIndex = 2,
  headerStripText = "\\].*",   # remove everything after the protein_id value
)
































Tursiops_truncatus # error
Girardinichthys_multiradiatus # error
Amblyraja_radiata # error

# Raw genome
SPECIES=Tursiops_truncatus
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/raw_genome/

mkdir -p genomes

cp -r /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/${SPECIES}/ ./genomes

cds_in=genomes/${SPECIES}/${SPECIES}.cds
cds_out=genomes/${SPECIES}/${SPECIES}.translated.cds
transeq -sequence "${cds_in}" -outseq "${cds_out}"

cp -r /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Gallus_gallus/ ./genomes

cds_in=genomes/Gallus_gallus/Gallus_gallus.cds
cds_out=genomes/Gallus_gallus/Gallus_gallus.translated.cds
transeq -sequence "${cds_in}" -outseq "${cds_out}"

## Prepare genespace
R

library(GENESPACE)


parsedPaths <- parse_annotations(
  rawGenomeRepo = "genomes",
  genomeDirs    = c("Tursiops_truncatus", "Gallus_gallus"),
  genomeIDs     = c("Tursiops_truncatus", "Gallus_gallus"),
  gffString     = "gff",
  faString      = "translated.cds",
  presets       = "ncbi",
  genespaceWd   = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/Gallus_gallus/troubleshooting/raw_genome/"
)

mkdir ${SPECIES}
mv bed/ ${SPECIES}
mv peptide/ ${SPECIES}

cp ../../sexshared/${SPECIES}/bed/Gal* ${SPECIES}/bed/
cp ../../sexshared/${SPECIES}/peptide/Gal* ${SPECIES}/peptide/

sbatch genespace.troubleshoot.raw_genome.sh Gallus_gallus ${SPECIES}



# Important Findings:
Amblyraja fails in synteny to Homo sapiens if you use the raw genome, but not if you use the filtered one!
All of these errors only occur in some, but not all, pairwise comparisons. It appears that there is a flaw in some genes or something... idk exactly.