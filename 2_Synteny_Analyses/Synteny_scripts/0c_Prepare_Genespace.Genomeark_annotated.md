# Download genome ark annotations
https://genomeark.s3.amazonaws.com/index.html?prefix=downstream_analyses/EGAPx_annotations/

## Identify which species lack an NCBI curated annoation, but have one on genome ark
cd /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists
export GENOMEARK="/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_list_of_annotated_species.txt"
export VGP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD.csv"
```
python3 - <<'PY'
import csv, os, re, sys

genomeark_path = os.environ["GENOMEARK"]
vgp_path       = os.environ["VGP"]
out_path       = os.environ.get("OUT_CSV", "genomeark_species.csv")

def key6(s: str) -> str:
    """
    Normalize IDs like:
      bGulAri2.1        -> GulAri
      mCorTow1_curated  -> CorTow
      aAnoBae1          -> AnoBae
    Rule: drop first letter, take next 6 alnum chars.
    """
    s = s.strip()
    if len(s) < 2:
        return ""
    tail = s[1:]
    tail = re.sub(r'[^A-Za-z0-9]', '', tail)  # remove '.', '_' etc
    return tail[:6] if len(tail) >= 6 else tail

# Read GenomeArk IDs in order (keep full original strings for output)
ids_full = []
seen_full = set()
with open(genomeark_path) as f:
    for line in f:
        x = line.strip()
        if not x:
            continue
        if x not in seen_full:
            ids_full.append(x)
            seen_full.add(x)

# Build VGP lookup by normalized key
key_to_species = {}
with open(vgp_path, newline='') as f:
    r = csv.DictReader(f)
    for row in r:
        vgp_id = (row.get("Assembly.ID") or "").strip()
        if not vgp_id:
            continue

        species = (row.get("Genus_species") or "").strip()
        if not species:
            sci = (row.get("Scientific.Name") or "").strip()
            species = sci.replace(" ", "_") if sci else ""

        if species:
            k = key6(vgp_id)
            if k:
                # If duplicates exist, keep the first encountered (stable)
                key_to_species.setdefault(k, species)

# Write matched to CSV; print unmatched GENOMEARK IDs to stdout
unmatched = []
with open(out_path, "w", newline="") as out_f:
    w = csv.writer(out_f)
    w.writerow(["GenomeArkID", "Species"])
    for full_id in ids_full:
        k = key6(full_id)
        sp = key_to_species.get(k)
        if sp:
            w.writerow([full_id, sp])
        else:
            unmatched.append(full_id)

# Unmatched IDs to stdout (one per line)
for u in unmatched:
    sys.stdout.write(u + "\n")
PY
```
# These aren't matched and require manual curation
```
fHapBur1
fLepSal2
mCanLor1.2
mPhaCin1_curated # may match phaCin_HiC (koala) in VGP
rCarIns
```
# Subset genomeark_species.csv to just those that are missing a gff already
```
export GENOME_ARK_SP="/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_species.csv"
export LACKS_GFF="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/analyses/lacks_gff"

SEX_CHR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv"

python3 - <<'PY'
import csv, os, sys

genome_ark_sp = os.environ["GENOME_ARK_SP"]   # CSV: GenomeArkID,Species
lacks_gff_dir = os.environ["LACKS_GFF"]       # dir with Species-named subdirs
out_path       = os.environ.get("OUT_CSV", "genomeark_species.lacks_NCBI_gff.csv")

# read species present in $LACKS_GFF (directory entries)
have = set()
for name in os.listdir(lacks_gff_dir):
    if name and not name.startswith("."):
        have.add(name)

# filter $GENOME_ARK_SP to rows whose Species is in $LACKS_GFF
with open(genome_ark_sp, newline="") as f:
    r = csv.DictReader(f)
    with open(out_path, "w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(r.fieldnames)
        for row in r:
            sp = row.get("Species", "")
            if sp in have:
                w.writerow([row.get("GenomeArkID",""), sp])
PY
```
## Download gff or gtf files from genome ark
```
ARK="/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_species.lacks_NCBI_gff.csv"
mkdir -p /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations
N=$(($(wc -l < "${ARK}") - 1))
sbatch --array=1-"$N" --export=ALL,ARK="${ARK}" 0d_download_genomeark_gff.sh
```
## If only gtf downloaded, convert to gff
```
ARK="/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_species.lacks_NCBI_gff.csv"
N=$(($(wc -l < "${ARK}") - 1))
sbatch --array=1-"$N" --export=ALL,ARK="${ARK}" 0e_convert_genomeark_gtf_to_gff.sh
```
## Download protein faa files from genome ark
```
ARK="/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_species.lacks_NCBI_gff.csv"
mkdir -p /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations
N=$(($(wc -l < "${ARK}") - 1))
sbatch --array=1-"$N" --export=ALL,ARK="${ARK}" 0f_download_genomeark_faa.sh
```
### Only Inia geoffresnsis lacks one, download it specifically
cd mIniGeo1
source myconda
mamba activate genespace
gffread mIniGeo1.gff -g mIniGeo1.fa -x mIniGeo1.cds
## Make symlink for fna files
ln -sf ${BASE}/genomes/${SPECIES}/ncbi_dataset/data/GCA_051312515.2/GCA_051312515.2_rAnnSte1.2_hap1_genomic.fna ${SP_ID}.fa 
ln -sf ${BASE}/genomes/${SPECIES}/ncbi_dataset/data/GCA_051312515.2/GCA_051312515.2_rAnnSte1.2_hap1_genomic.fna  /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Anniella_stebbinsi/Anniella_stebbinsi.fna
```
# -------- user inputs --------
ARK="${ARK:-/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_species.lacks_NCBI_gff.csv}"
BASE="/data/Wilson_Lab/data/VGP_genomes_phase1"
WORKDIR="/data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations"
# -----------------------------

mkdir -p "${WORKDIR}/logs"
cd "${WORKDIR}"

# Skip header, read CSV lines
tail -n +2 "${ARK}" | while IFS=',' read -r SP_ID SPECIES; do
  # basic validation
  [[ -z "${SP_ID:-}" || -z "${SPECIES:-}" ]] && {
    echo "WARN: skipping malformed line: SP_ID='${SP_ID:-}' SPECIES='${SPECIES:-}'" >&2
    continue
  }

  echo "Processing: SP_ID=${SP_ID} SPECIES=${SPECIES}"

  mkdir -p "${SP_ID}"
  cd "${SP_ID}"

  FNA_DIR="${BASE}/genomes/${SPECIES}/ncbi_dataset/data"

  found=""
  count=0

  # Find matching genome FASTA(s)
  while IFS= read -r f; do
    ((count++))
    if (( count == 1 )); then
      found="$f"
    else
      echo "ERROR: Multiple .fna files found for ${SP_ID} under ${FNA_DIR}:" >&2
      printf '  %s\n' "$found" >&2
      printf '  %s\n' "$f" >&2
      while IFS= read -r rest; do
        printf '  %s\n' "$rest" >&2
      done
      cd "${WORKDIR}"
      continue
    fi
  done < <(find "${FNA_DIR}" -mindepth 1 -maxdepth 3 -type f -name "GC*.fna" 2>/dev/null | sort)

  if (( count == 0 )); then
    echo "ERROR: No .fna files found for ${SP_ID} under ${FNA_DIR}" >&2
    cd "${WORKDIR}"
    continue
  fi

  ln -sf "$found" "${SP_ID}.fa"
  echo "  -> linked ${SP_ID}.fa -> $found"

  cd "${WORKDIR}"
done
```




```

source myconda
mamba activate genespace
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark
 
Rscript 0g_prepare_genespace_genomeark.r /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/ /data/Wilson_Lab/data/VGP_genomes_phase1/genomeark_annotations/
```
## Rename all bed and peptide files to species names, not IDs
```
ARK=/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_species.lacks_NCBI_gff.csv
DIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/bed

## BEDS
# DRY RUN: prints what it would do
awk -F',' 'NR>1 {print $1, $2}' "$ARK" \
| while read -r id sp; do
    src="$DIR/$id.bed"
    dst="$DIR/$sp.bed"
    [[ -e "$src" ]] || continue
    if [[ -e "$dst" ]]; then
      echo "SKIP (dst exists): $src -> $dst"
    else
      echo "mv -v -- '$src' '$dst'"
    fi
  done

# real script
awk -F',' 'NR>1 {print $1, $2}' "$ARK" \
| while read -r id sp; do
    src="$DIR/$id.bed"
    dst="$DIR/$sp.bed"
    [[ -e "$src" ]] || continue
    if [[ -e "$dst" ]]; then
      echo "SKIP (dst exists): $src -> $dst" >&2
    else
      mv -v -- "$src" "$dst"
    fi
  done

## PEPTIDES
DIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/peptide

awk -F',' 'NR>1 {print $1, $2}' "$ARK" \
| while read -r id sp; do
    src="$DIR/$id.fa"
    dst="$DIR/$sp.fa"
    [[ -e "$src" ]] || continue
    if [[ -e "$dst" ]]; then
      echo "SKIP (dst exists): $src -> $dst"
    else
      echo "mv -v -- '$src' '$dst'"
    fi
  done

# real script
awk -F',' 'NR>1 {print $1, $2}' "$ARK" \
| while read -r id sp; do
    src="$DIR/$id.fa"
    dst="$DIR/$sp.fa"
    [[ -e "$src" ]] || continue
    if [[ -e "$dst" ]]; then
      echo "SKIP (dst exists): $src -> $dst" >&2
    else
      mv -v -- "$src" "$dst"
    fi
  done
```
# Reformat bed files to have chromosome numbers, not accession numbers
## Create a list of bed files requiring reformatting
```
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark
ARK=/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_species.lacks_NCBI_gff.csv
mkdir -p ${OUTDIR}/reference_lists
BED_DIR=${OUTDIR}/bed
for f in "$BED_DIR"/*.bed; do
  line=$(grep -m1 -vE '^\s*($|#|track|browser)' "$f" || true)
  [ -z "$line" ] && continue
  chrom=$(awk '{print $1}' <<<"$line")
  if ! [[ "$chrom" =~ ^([0-9]+|X|Y|M|MT)$ ]]; then
    basename "$f" .bed
  fi
done > "${OUTDIR}/reference_lists/species_needing_chrom_reformat.txt"
```
## Create a reformating index for each that needs it
```
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark
LIST="${OUTDIR}/reference_lists/species_needing_chrom_reformat.txt"
FNA_BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"
MAPDIR="${OUTDIR}/reformatting_chr"
mkdir -p "$MAPDIR"

while IFS= read -r SPECIES; do
  [[ -z "${SPECIES// /}" ]] && continue
  [[ "$SPECIES" =~ ^# ]] && continue

  FNA="${FNA_BASE}/${SPECIES}/${SPECIES}.fna"
  MAP="${MAPDIR}/${SPECIES}.chr_reformat.txt"

  [[ ! -s "$FNA" ]] && continue

  awk '
    BEGIN{ OFS="," }
    function trim(s){ gsub(/^[[:space:]]+|[[:space:]]+$/, "", s); return s }
    function up(s){ return toupper(s) }

    /^>/{
      line=$0
      gsub(/\r/,"",line)

      acc=$1; sub(/^>/,"",acc); acc=trim(acc)
      u=up(line)

      # same exclusion as above (optional but recommended)
      if (u ~ /(SCAFFOLD|CONTIG|UNPLACED|UNLOCALIZED|UNLOC|PATCH|ALT_REF_LOCI|ALTERNATE|HAPLOTYP)/) next

      chr=""
      if (match(u, /CHROMOSOME[:[:space:]]+([0-9]+|[WXYZ][0-9]*)/, m)) chr=m[1]
      else if (match(u, /(^|[^A-Z0-9])CHR[[:space:]]*([0-9]+|[WXYZ][0-9]*)/, m)) chr=m[2]
      else if (match(u, /LINKAGE[[:space:]]+GROUP[:[:space:]]+([0-9]+)/, m)) chr=m[1]
      else if (match(u, /(^|[^A-Z0-9])LG[[:space:]]*([0-9]+)/, m)) chr=m[2]

      if (chr!="") print acc, chr
    }
  ' "$FNA" | sort -u > "$MAP"

  [[ -s "$MAP" ]] && echo "Wrote $MAP" >&2
done < "$LIST"
```
# Output list of instances in which chrs in reformating index do not match any chrs in bed
This suggests a mismatch between versions of genome between annotation and assembly
```
LIST="${OUTDIR}/reference_lists/species_needing_chrom_reformat.txt"
INBED_DIR="${OUTDIR}/bed"
OUTBED_DIR="${INBED_DIR}_chrfixed"
MISSING_REPORT="${OUTDIR}/reference_lists/map_chrs_not_in_bed.tsv"

mkdir -p "${OUTBED_DIR}"
: > "${MISSING_REPORT}"  # truncate report

while IFS= read -r SPECIES; do
  [[ -z "${SPECIES// /}" ]] && continue
  [[ "$SPECIES" =~ ^# ]] && continue

  MAP="${MAPDIR}/${SPECIES}.chr_reformat.txt"
  BED="${INBED_DIR}/${SPECIES}.bed"
  OUT="${OUTBED_DIR}/${SPECIES}.bed"

  # ---- 1) Report: chromosomes present in MAP (keys) that do NOT occur in BED col1 ----
  # Output columns: species<TAB>map_chr
  awk -v FS="\t" -v OFS="\t" -v sp="${SPECIES}" '
    # First pass: BED -> collect chromosomes present
    FNR==NR { bed[$1]=1; next }

    # Second pass: MAP -> parse "old,new" and report old not in BED
    {
      split($0,a,",");
      old=a[1];
      if (!(old in bed)) print sp, old;
    }
  ' "${BED}" "${MAP}" >> "${MISSING_REPORT}"

  # (Optional) also echo a per-species summary to stderr
  n_missing=$(awk -v FS="\t" -v sp="${SPECIES}" '$1==sp{c++} END{print c+0}' "${MISSING_REPORT}")
  if (( n_missing > 0 )); then
    echo "[${SPECIES}] MAP keys missing from BED: ${n_missing} (see ${MISSING_REPORT})" >&2
  fi

  # ---- 2) Do the chr reformat on BED using MAP ----
  awk -v FS="\t" -v OFS="\t" '
    FNR==NR { split($0,a,","); map[a[1]]=a[2]; next }
    ($1 in map) { $1=map[$1] }
    { print }
  ' "${MAP}" "${BED}" > "${OUT}"

  [[ -s "${OUT}" ]] && echo "Wrote ${OUT}" >&2
done < "${LIST}"

echo "Missing-chrom report written to: ${MISSING_REPORT}" >&2
```
# Identify all species with mismatched chromosome accession identifiers -- assume that a better annotation will be available before the time of publication
```
awk '{print $1}' ${MISSING_REPORT} | sort -u
SP=Willisornis_vidua
grep ${SP} ${MISSING_REPORT}
awk '{print $1}' bed/${SP}.bed | sort -u 
awk '{print $1}' bed_chrfixed/${SP}.bed | sort -u 

Anniella_stebbinsi  # not chromosome level
Artibeus_lituratus # annotation is missing updated sex chromosomes
Caesio_teres # annotation stops after chr 16... missing 17-23... no sex chrs in annotation
Larus_fuscus # annotation is from a previous version
Microchirus_variegatus # unk sex chr, annotation stops after chr 19 (missing 20-23)
Monodon_monocero # Annotations are on scaffolded assembly, while fasta is on chrs
Myotis_emarginatus # annotation is from a previous version
Myotis_mystacinus # annotation is from a previous version for everything but sex chrs... could actually be okay here? sex chrs are all version OZ*.1
Pipistrellus_nathusii # annotation is from a previous version
Pipistrellus_pygmaeus # annotation is from a previous version
Vespertilio_murinus # annotation is from a previous version
Willisornis_vidua # annotation stops after chr 36... missing 37-39
```
# Remove these species from bed and peptide dirs
```
for SP in \
  Anniella_stebbinsi \
  Artibeus_lituratus \
  Caesio_teres \
  Larus_fuscus \
  Microchirus_variegatus \
  Monodon_monocero \
  Myotis_emarginatus \
  Myotis_mystacinus \
  Pipistrellus_nathusii \
  Pipistrellus_pygmaeus \
  Vespertilio_murinus \
  Willisornis_vidua
do
  rm -v -f "bed/${SP}.bed" "bed_chrfixed/${SP}.bed" "peptide/${SP}.fa"
done
```
## Remove W and Y chromosomes from bed files, and all scaffolds
```
INDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/


INBED="${INDIR}/bed_chrfixed"
INPEP="${INDIR}/peptide"

OUTBED="${OUTDIR}/sexshared/usable_bed"
OUTPEP="${OUTDIR}/sexshared/usable_peptide"

mkdir -p "$OUTBED" "$OUTPEP"

shopt -s nullglob

# 1) Filter BEDs: drop rows where chrom (col 1) is Y*, W*, and drop rows where chrom (col 1) has > 3 total characters
for bed in "$INBED"/*.bed; do
  outbed="${OUTBED}/$(basename "$bed")"

  awk -F'\t' 'BEGIN{OFS="\t"}
    # drop Y, Y1.. and W, W1..
    $1 ~ /^(Y[0-9]*|W[0-9]*)$/ { next }

    # drop if chrom name is longer than 3 chars total
    length($1) > 3 { next }

    { print }
  ' "$bed" > "$outbed"

  echo "Wrote $outbed"
done

# 2) Filter peptide FASTAs to match filtered BEDs (keep only genes still present)
#    - assumes peptide file name matches bed base name: SPECIES.bed <-> SPECIES.fa
for bed in "$OUTBED"/*.bed; do
  base="$(basename "$bed" .bed)"
  pep_in="${INPEP}/${base}.fa"
  pep_out="${OUTPEP}/${base}.fa"

  if [[ ! -s "$pep_in" ]]; then
    echo "WARN: missing/empty peptide FASTA: $pep_in (skipping)" >&2
    continue
  fi

  # IDs to keep = col 4 from filtered BED
  keep_ids="$(mktemp)"
  awk -F'\t' '$4!=""{print $4}' "$bed" | sort -u > "$keep_ids"

  # Subset FASTA by header name (requires seqkit)
  seqkit grep -n -f "$keep_ids" "$pep_in" > "$pep_out"

  rm -f "$keep_ids"
  echo "Wrote $pep_out"
done

```

# Didn't use anything below this line:














```


mkdir -p "$OUTDIR_CHR"

while IFS= read -r SPECIES; do
  # skip blanks/comments
  [[ -z "${SPECIES// /}" ]] && continue
  [[ "$SPECIES" =~ ^# ]] && continue

  FNA="${FNA_BASE}/${SPECIES}/${SPECIES}.fna"
  OUT="${OUTDIR_CHR}/${SPECIES}.chr_reformat.txt"

  if [[ ! -s "$FNA" ]]; then
    echo "WARN: missing/empty $FNA (skipping)" >&2
    continue
  fi

  # Parse headers: accession is token 1; chr captured from either:
  #   - "chromosome <token>"
  #   - "chromosome: <token>"
  # Only emit headers that match one of those patterns.
  awk '
    BEGIN { OFS="," }
    /^>/ {
      acc=$1; sub(/^>/,"",acc)

      chr=""
      if (match($0, /chromosome[[:space:]]+([[:alnum:]_.-]+)/, m)) {
        chr=toupper(m[1])
      } else if (match($0, /chromosome:[[:space:]]+([[:alnum:]_.-]+)/, m)) {
        chr=m[1]
      }

      if (chr != "") print acc, chr
    }
  ' "$FNA" > "$OUT"

  if [[ ! -s "$OUT" ]]; then
    echo "WARN: no chromosome headers found in $FNA (wrote empty $OUT)" >&2
  else
    echo "Wrote $OUT" >&2
  fi

done < "$LIST"
```
## This did not create an index for the following species:
```
reformatting_chr/Sternotherus_odoratus.chr_reformat.txt # not chromosome level, just scaffolded
```
# 
## Recode the chromosome of each bed file
```
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark

LIST="${OUTDIR}/reference_lists/species_needing_chrom_reformat.txt"
MAPDIR="${OUTDIR}/reformatting_chr"

INBED_DIR="${OUTDIR}/bed"
OUTBED_DIR="${INBED_DIR}_chrfixed"

mkdir -p "$OUTBED_DIR"

while IFS= read -r SPECIES; do
  [[ -z "${SPECIES// /}" ]] && continue
  [[ "$SPECIES" =~ ^# ]] && continue

  BED="${INBED_DIR}/${SPECIES}.bed"
  MAP="${MAPDIR}/${SPECIES}.chr_reformat.txt"
  OUT="${OUTBED_DIR}/${SPECIES}.bed"

  if [[ ! -s "$BED" ]]; then
    echo "WARN: missing/empty BED: $BED (skipping)" >&2
    continue
  fi
  if [[ ! -s "$MAP" ]]; then
    echo "WARN: missing/empty MAP: $MAP (skipping)" >&2
    continue
  fi

  awk -v FS="\t" -v OFS="\t" '
    function up(s){ return toupper(s) }

    FNR==NR {
      gsub(/\r$/, "", $0)
      split($0, a, ",")
      acc=a[1]
      chr=up(a[2])

      # Load mapping only for chromosome-like accessions:
      #   CM#######.#  or  NC_########.#
      # Accept chr as:
      #   numeric (1,2,...) OR sex chroms with optional numeric suffix (X, X1.., Y, Y1.., Z, Z1.., W, W1..)
      if ( (acc ~ /^CM[0-9]+\.[0-9]+$/ || acc ~ /^NC_[0-9]+\.[0-9]+$/ || acc ~ /^OZ[0-9]+\.[0-9]+$/ || acc ~ /^OX[0-9]+\.[0-9]+$/ || acc ~ /^OW[0-9]+\.[0-9]+$/ || acc ~ /^OY[0-9]+\.[0-9]+$/ || acc ~ /^LR[0-9]+\.[0-9]+$/ || acc ~ /^PC[0-9]+\.[0-9]+$/ || acc ~ /^OU[0-9]+\.[0-9]+$/ || acc ~ /^OX[0-9]+\.[0-9]+$/ ) &&
           (chr ~ /^[0-9]+$/ || chr ~ /^[WXYZ][0-9]*$/) ) {
        map[acc]=chr
      }
      next
    }

    {
      if ($1 in map) $1 = map[$1]
      print
    }
  ' "$MAP" "$BED" > "$OUT"

  echo "Wrote $OUT" >&2
done < "$LIST"


```
## replace original bed files with reformatted bed files
```
mv ${OUTDIR}/bed_chrfixed/* ${OUTDIR}/bed/
```
## Remove empty or scaffold level beds
```
rm bed/Inia_geoffrensis.bed
rm bed/Sternotherus_odoratus.bed
```

## Remove W and Y chromosomes from bed files, and all scaffolds
```
INDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/


INBED="${INDIR}/bed"
INPEP="${INDIR}/peptide"

OUTBED="${OUTDIR}/sexshared/usable_bed"
OUTPEP="${OUTDIR}/sexshared/usable_peptide"

mkdir -p "$OUTBED" "$OUTPEP"

shopt -s nullglob

# 1) Filter BEDs: drop rows where chrom (col 1) is Y*, W*, and drop rows where chrom (col 1) has > 3 total characters
for bed in "$INBED"/*.bed; do
  outbed="${OUTBED}/$(basename "$bed")"

  awk -F'\t' 'BEGIN{OFS="\t"}
    # drop Y, Y1.. and W, W1..
    $1 ~ /^(Y[0-9]*|W[0-9]*)$/ { next }

    # drop if chrom name is longer than 3 chars total
    length($1) > 3 { next }

    { print }
  ' "$bed" > "$outbed"

  echo "Wrote $outbed"
done

# 2) Filter peptide FASTAs to match filtered BEDs (keep only genes still present)
#    - assumes peptide file name matches bed base name: SPECIES.bed <-> SPECIES.fa
for bed in "$OUTBED"/*.bed; do
  base="$(basename "$bed" .bed)"
  pep_in="${INPEP}/${base}.fa"
  pep_out="${OUTPEP}/${base}.fa"

  if [[ ! -s "$pep_in" ]]; then
    echo "WARN: missing/empty peptide FASTA: $pep_in (skipping)" >&2
    continue
  fi

  # IDs to keep = col 4 from filtered BED
  keep_ids="$(mktemp)"
  awk -F'\t' '$4!=""{print $4}' "$bed" | sort -u > "$keep_ids"

  # Subset FASTA by header name (requires seqkit)
  seqkit grep -n -f "$keep_ids" "$pep_in" > "$pep_out"

  rm -f "$keep_ids"
  echo "Wrote $pep_out"
done

```

# /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/genespace_genomeark/peptide/Anguilla_anguilla.fa
# CP OU OX
# OZ? Larus_fuscus Myotis_emarginatus Myotis_mystacinus Pipistrellus_nathusii Vespertilio_murinus
# scaffolded: Monodon_monocero
SP=Vespertilio_murinus
awk '{print $1}' bed/${SP}.bed | sort -u



## Remove Z and X chromosomes from bed files
```

OUTBED="${OUTDIR}/sexlimited/usable_bed"
OUTPEP="${OUTDIR}/sexlimited/usable_peptide"

mkdir -p "$OUTBED" "$OUTPEP"

shopt -s nullglob

# 1) Filter BEDs: drop rows where chrom (col 1) is Y*, W*
for bed in "$INBED"/*.bed; do
  outbed="${OUTBED}/$(basename "$bed")"

  awk -F'\t' 'BEGIN{OFS="\t"} $1 !~ /^(X[0-9]*|Z[0-9]*)$/ {print}' "$bed" > "$outbed"
  echo "Wrote $outbed"
done

# 2) Filter peptide FASTAs to match filtered BEDs (keep only genes still present)
#    - assumes peptide file name matches bed base name: SPECIES.bed <-> SPECIES.fa
for bed in "$OUTBED"/*.bed; do
  base="$(basename "$bed" .bed)"
  pep_in="${INPEP}/${base}.fa"
  pep_out="${OUTPEP}/${base}.fa"

  if [[ ! -s "$pep_in" ]]; then
    echo "WARN: missing/empty peptide FASTA: $pep_in (skipping)" >&2
    continue
  fi

  # IDs to keep = col 4 from filtered BED
  keep_ids="$(mktemp)"
  awk -F'\t' '$4!=""{print $4}' "$bed" | sort -u > "$keep_ids"

  # Subset FASTA by header name (requires seqkit)
  seqkit grep -n -f "$keep_ids" "$pep_in" > "$pep_out"

  rm -f "$keep_ids"
  echo "Wrote $pep_out"
done
```



# List of all species with annotations
Save as: /data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/genomeark_list_of_annotated_species.txt
```
aAnoBae1
aAscTru1
aDisPic1
aGasCar1
aLepFus1
aManAur1
aMicUni1
aPelIbe1
aPelLes1
aRhiDor1
bAcaChl1
bAegAlb1
bAlcTor1
bAmaOch1
bAmmMai1
bAnsAns1
bAnsBra1
bAraAra1
bArcTri1
bAthNoc1
bAytFer1
bAytMar2
bBraCan1
bBucAby1
bBucCla1
bButBut1
bCalBor7
bCapEur3
bCarCri1
bChlMac1
bCicMag1
bClaHye2
bColMon1
bCotChi1
bCyaCrs1
bCygCol1
bEriRub2
bEudEle1
bFalPun1
bFriCoe1
bGalChl1
bGeoTri1
bGruGru1
bGuaGua1
bGulAri2.1
bGypBar2
bHemCom1
bLarArg3
bLarFus1
bLepDis1
bLonStr1
bLycPyr6
bMerNub1
bMerOct1
bMorBas2
bMorGui1
bNesMay2
bNetRuf1
bNumArq3
bPasSan2
bPelCri1
bPhaAeh10
bPhaSup1
bPhoRub2
bPlaLeu2
bPluApr1
bPodStr1
bPorHoc1
bPsiEch3
bPteGut1
bRhyJub1
bSarPap1
bSphHub1
bSteHir1
bStrDea1
bStrTur12
bTauEry1
bTetUro1.1
bTheCae1
bTroSur1
bVanVan1
bWilVid1
bZonAlb1
bZosLat1
fAbrBra2.1
fAcaSch1
fAmbSpe1
fAmmMar1
fAnaAna1
bAnaPla2
fAntMac1
fAplTae1
fArgSil1
fArrGeo1
fAstCal1
fAulMac1
fAulStu2
fBarBab1
fBarBar1
fBleOce1
fBorAnt1
fCaeTer1
fCenGer3
fCheLab1
fChoSch1
fClaGar1
fColCha1
fCorLav1
fCriAus1
fCypVen1
fDanRer4
fDirArg1
fEchVip8
fEleAnt2
fEleEle2
fEpiLan1
fEpiMut1
fEutGur1
fFunDia1
fGobGob1
fGobNig1
fGymBra2
fGymMic1
fHapBur1
fLatMac1
fLepSal2
fLetNeb1
fLeuLeu2
fLipPho2
fLycPai1
fMelGel1.1
fMenAtl1
fMicKit1
fMicPou1
fMicVar1
fNanAch1
fNanAnt1
fNotCoa1
fNotRos5.1
fOsmMor3
fPhoGun1
fPhoLeu1
fPhoPho1
fPolHol1
fPolLow2
fPolPol2.1
fPorCra3
fPriTyp1
fProBol1
fProPar1
fRhaChi2
fRhiNas1
fSalAlp3
fSalTru1
fScaEry2
fScoMax1
fSilAri3.1
fSpiSpi1
fSprSpr1
fSquCep2.1
fSymMel2
fSymNem1
fSynOce1
fSynPic1
fSynTyp1
fTauAds1
fTauBub2
fTelBon1
fTraTra1
fXenCan1
fZeuFab8
mAcoMin1
mAntPal1
mArtInt1
mArtLit1
mArvNil1
mAseSto1
mBalBor1
mBalPhy1
mCalLat2
mCanLor1.2
mCasCan1
mCorTow1_curated
mCteGun2
mDamDam1
mDasMac1
mDesRot1
mDorCyc1
mDugDug1
mEptNil1.1
mEreDor1
mGloMut1
mGraGri1
mHalGry1
mHetBru1
mHipLar1
mHypAmp2.1
mIniGeo1
mLagAcu1
mLycPic1
mMacLag1
mMarMar1
mMegNov1
mMegSpa1
mMesBid2
mMesMir1
mMicMin1
mMicPen1
mMinSch1
mMolAlv2
mMolNig1
mMonMon1
mMopCon1
mMusAve1
mMusNiv2
mMyoEma1
mMyoMys1
mMyoNat1
mNeoVis2
mPanOnc1
mPerMan1_curated
mPhaCin1_curated
mPipHan1
mPipNat1
mPipPip1
mPipPyg2.1
mPleAur1
mRhiAff1
mRhiHip2
mRhiMic1
mRhiPer1
mRhiTri1
mRhiYon1
mRhyNas2
mRhyPet1
mSaiBol1
mSarHar1
mSciVul1.2
mSmiCra1
mSpeCit3
mSteCoe1.1
mTadBra1
mTalEur1
mTamTet1
mTapInd1
mTriInu1
mTupTan_curated
mVesMur1
mVulVul1
rAnnSte1
rAspTig1
rCarIns
rCycPin1
rDibSmi1
rEreImb1
rFurPar1
rGavGan2
rLepKem1
rLepOli2
rLiaOli1
rMacSuw1
rMyuGeo1
rNatDep1
rNatHel1
rPanTec1
rPodBoc1
rPodCre2.1
rPodErh1
rPodFil1
rPodGai1
rPodLio1
rPodMel1
rPodMur119
rPodPit1
rPodTil1
rPodVau1
rShiCro1
rSteOdo2
rTilSci1
rVipBer3
rVipLat1
rVipUrs1
sRajBra1
```


# follow up on 
bGeoTri1/
fHapBur2_hap1
fHapBur2_hap2
rPodExp1_out/
rPodSic1/
mCorTow1_curated

# not relevant
bRhyJub1_alignments/
rGa/
rLepOli2_alignments/