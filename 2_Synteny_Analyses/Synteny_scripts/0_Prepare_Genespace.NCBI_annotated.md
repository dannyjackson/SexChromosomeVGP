# Genespace installation and environment management
```
mamba env create -f /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/yamls/genespace_py3.10.yml
mamba activate genespace_py3.10

library(devtools)
devtools::install_github("jtlovell/GENESPACE")
```
# Analyses of synteny of sex chromosomes across the VGP Phase 1 genomes
```
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/


mkdir -p ${OUTDIR}/scripts
mkdir -p ${OUTDIR}/analyses
mkdir -p ${OUTDIR}/datafiles
mkdir -p ${OUTDIR}/reference_lists


cd ${OUTDIR}/scripts
```
## Make list of sex chromosomes from genomes with gff
```
python3 0a_filter_sexchrom_to_gff.py ${OUTDIR}/reference_lists
```

## Prepare all bed and peptide files
```

cp -r /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/ ${OUTDIR}/datafiles


HAS_DIR="${OUTDIR}/analyses/has_gff"
LACKS_DIR="${OUTDIR}/analyses/lacks_gff"

mkdir -p "$HAS_DIR" "$LACKS_DIR"

shopt -s nullglob dotglob

for d in ${OUTDIR}/datafiles/symlinks/*/ ; do
  # skip our destination folders if they already exist here
  case "$d" in
    "$HAS_DIR/") continue ;;
    "$LACKS_DIR/") continue ;;
  esac

  # does this subdir contain any .gff files?
  gffs=( "$d"*.gff )

  if ((${#gffs[@]} > 0)); then
    mv -- "$d" "$HAS_DIR/"
    echo "HAS_GFF  -> $d"
  else
    mv -- "$d" "$LACKS_DIR/"
    echo "NO_GFF   -> $d"
  fi
done

# remove all lifton 
find $HAS_DIR -type f -iname "*lifton*" -print
# mv Artibeus_intermedius ../lacks_gff
# mv Artibeus_lituratus ../lacks_gff
mv ${HAS_DIR}/Pseudacris_triseriata ${LACKS_DIR}

```
## Remove scaffolds from FNA files
```
find "${HAS_DIR}" -mindepth 2 -maxdepth 2 -name "*.fna" -print0 |
  while IFS= read -r -d '' f; do
    out="${f%.fna}.contigs_only.fna"
    awk 'BEGIN{IGNORECASE=1}
         /^>/{
           drop = ($0 ~ /(scaffold|unplaced|unlocalized|unloc|mitochon)/)
         }
         !drop{print}
        ' "$f" > "$out"
    echo "$f -> $out"
  done
```
### Fix the ones that are incorrectly filtered
Genomes with irregularities that require manual curation:
```
# scaffolds remain

## remove anything containing '>JA'
Xenentodon_cancila 
Pholidichthys_leucotaenia
Porphyrio_hochstetteri 
Vipera_latastei 
Gastrophryne_carolinensis 
Thomomys_bottae
Fundulus_diaphanus

## remove anything containing '>JB'
Lycaon_pictus
Morphnus_guianensis
Guaruba_guaruba
Trichechus_inunguis
```
# fix these
```
species=(
  Xenentodon_cancila
  Pholidichthys_leucotaenia
  Porphyrio_hochstetteri
  Vipera_latastei
  Gastrophryne_carolinensis
  Thomomys_bottae
  Fundulus_diaphanus
  Lycaon_pictus
  Morphnus_guianensis
  Guaruba_guaruba
  Trichechus_inunguis
)

for sp in "${species[@]}"; do
  f="${HAS_DIR}/${sp}/${sp}.fna"
  out="${HAS_DIR}/${sp}/${sp}.contigs_only.fna"

  if [[ ! -f "$f" ]]; then
    echo "SKIP (missing): $f" >&2
    continue
  fi

  awk 'BEGIN{IGNORECASE=1}
       /^>/{
         drop = ($0 ~ /(scaffold|unplaced|unlocalized|unloc|mitochon|>JB|>JA)/)
       }
       !drop{print}
      ' "$f" > "$out"

  echo "$f -> $out"
done
```

## Filter ${SPECIES}/${SPECIES}.gff to just regions kept in file ${SPECIES}/${SPECIES}.contigs_only.fna
```

OUT_SUFFIX=".contigs_only.gff"   # output name: ${SPECIES}/${SPECIES}${OUT_SUFFIX}

shopt -s nullglob
for fasta in "${HAS_DIR}"/*/*.contigs_only.fna; do
  species="$(basename "$(dirname "$fasta")")"
  gff="${HAS_DIR}/${species}/${species}.gff"
  out="${HAS_DIR}/${species}/${species}${OUT_SUFFIX}"

  if [[ ! -s "$gff" ]]; then
    echo "WARN: missing/empty GFF for ${species}: $gff (skipping)" >&2
    continue
  fi

  awk '
    BEGIN { FS=OFS="\t" }

    # File 1: FASTA -> allowed seqids
    FNR==NR {
      if ($0 ~ /^>/) {
        h=$0
        sub(/^>/, "", h)
        sub(/[ \t].*$/, "", h)   # first token only
        keep[h]=1
      }
      next
    }

    # File 2: GFF -> keep header/comments, and features whose seqid is allowed
    /^#/ { print; next }
    ($1 in keep) { print }
  ' "$fasta" "$gff" > "$out"

  echo "OK: ${species} -> $out"
done

# check that all output files have some length
for gff in ${HAS_DIR}/*/*.contigs_only.gff; do
  species="$(basename "$(dirname "$gff")")"
  lines="$(wc -l < "$gff")"
  printf "%s\t%s\n" "$species" "$lines"
done | sort -k2,2n

# rename gff with contig.only
for gff in ${HAS_DIR}/*/*.contigs_only.gff; do
  species="$(basename "$(dirname "$gff")")"
  echo ${species}
  mv $gff ${HAS_DIR}/${species}/${species}.gff
done 

```
## Translate nucleic acid sequences for all species
```
shopt -s nullglob

for d in "${HAS_DIR}"/*/; do
  SPECIES="$(basename "${d%/}")"
  echo $d
  # skip if any translated.cds already exists
  translated=( "$d"*.translated.cds )
  if (( ${#translated[@]} > 0 )); then
    echo "SKIP (already translated): $SPECIES"
    continue
  fi

  cds_in="${d}${SPECIES}.cds"
  cds_out="${d}${SPECIES}.translated.cds"

  if [[ ! -f "$cds_in" ]]; then
    echo "SKIP (missing input): $SPECIES -> $(basename "$cds_in")"
    continue
  fi

  echo "RUN: $SPECIES"
  (
    transeq -sequence "${cds_in}" -outseq "${cds_out}"
  )
done

```
## Prepare all bed and peptide files
```

Rscript 0b_prepare_genespace.R "$OUTDIR" "$HAS_DIR"

```
# Make a list of all empty beds and peptide files
```
mkdir -p ${OUTDIR}/empty_files/bed ${OUTDIR}/empty_files/peptide

# empty beds/
find ${OUTDIR}/bed -type f -empty -print0 | xargs -0 -r mv -t ${OUTDIR}/empty_files/bed
ls ${OUTDIR}/empty_files/bed | wc -l

# empty peptide/
find ${OUTDIR}/peptide -type f -empty -print0 | xargs -0 -r mv -t ${OUTDIR}/empty_files/peptide
ls ${OUTDIR}/empty_files/peptide | wc -l

mkdir -p ${OUTDIR}/short_files/bed ${OUTDIR}/short_files/peptide

# Make a list of all bed and peptides that are too small (suggests something went wrong)
# short_bed/

find ${OUTDIR}/bed -type f -exec sh -c '
  f="$1"
  [ "$(wc -l < "$f")" -lt 2000 ] && mv -- "$f" ${OUTDIR}/short_files/bed/
' sh {} \;


find "${OUTDIR}/bed" -type f -exec sh -c '
  OUTDIR="$1"
  f="$2"
  [ "$(wc -l < "$f")" -lt 2000 ] && mv -- "$f" "$OUTDIR/short_files/bed/"
' sh "$OUTDIR" {} \;

# short_peptide/
find "${OUTDIR}/peptide" -type f -exec sh -c '
  OUTDIR="$1"
  f="$2"
  [ "$(wc -l < "$f")" -lt 10000 ] && mv -- "$f" "${OUTDIR}/short_files/peptide/"
' sh "$OUTDIR" {} \;

mkdir -p ${OUTDIR}/usable_bed/
mkdir -p ${OUTDIR}/usable_peptide/

mv ${OUTDIR}/bed/* ${OUTDIR}/usable_bed/
mv ${OUTDIR}/peptide/* ${OUTDIR}/usable_peptide/
```


# Create a list of species that did not run properly 
```

ls -1 "${OUTDIR}/empty_files/bed" "${OUTDIR}/short_files/bed" \
  | grep -vE ':$' \
  | sed 's/\.bed$//' | sort -u \
  > "${OUTDIR}/reference_lists/species_error.parse_annotations.txt"
```
# Run parse_annotations on them again, using the protein_id field instead
I troubleshot this on Ammospiza maritima, a random species that did not work with the above parse annotations script. I then applied the same script to all species that did not run properly and, miraculously, this script worked for all of them. Very grateful that I didn't have to individually troubleshoot them all!
${OUTDIR}/reference_lists/species_error.parse_annotations.txt

```
R

library(GENESPACE)

repo <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace//analyses/has_gff"

SPECIES <- readLines("/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/reference_lists/species_error.parse_annotations.txt")
SPECIES <- SPECIES[nzchar(SPECIES)]

parsedPaths2 <- parse_annotations(
  rawGenomeRepo = repo,
  genomeDirs = SPECIES,
  genomeIDs  = SPECIES,
  gffString  = "gff",
  faString   = "translated.cds",
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/",
  gffIdColumn = "protein_id",
  headerSep = "protein_id=",
  headerEntryIndex = 2,
  headerStripText = "\\].*",   # remove everything after the protein_id value
)
```
# Confirm that no empty or short beds remain
```
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/

rm ${OUTDIR}/empty_files/*/*
rm ${OUTDIR}/short_files/*/*


# empty bed/
find ${OUTDIR}/bed -type f -empty -print0 | xargs -0 -r mv -t ${OUTDIR}/empty_files/bed
ls ${OUTDIR}/empty_files/bed

# empty peptide/
find ${OUTDIR}/peptide -type f -empty -print0 | xargs -0 -r mv -t ${OUTDIR}/empty_files/peptide
ls ${OUTDIR}/empty_files/peptide

# short_bed
find ${OUTDIR}/bed -type f -exec sh -c '
  f="$1"
  [ "$(wc -l < "$f")" -lt 2000 ] && mv -- "$f" ${OUTDIR}/short_files/bed/
' sh {} \;
ls ${OUTDIR}/short_files/bed/

# short_peptide/
find ${OUTDIR}/peptide -type f -exec sh -c '
  f="$1"
  [ "$(wc -l < "$f")" -lt 10000 ] && mv -- "$f" ${OUTDIR}/short_files/peptide/
' sh {} \;
ls ${OUTDIR}/short_files/peptide/


mv ${OUTDIR}/bed/* ${OUTDIR}/usable_bed/
mv ${OUTDIR}/peptide/* ${OUTDIR}/usable_peptide
```
# Reformat bed files to have chromosome numbers, not accession numbers
## Create a list of bed files requiring reformatting
```

BED_DIR=${OUTDIR}/usable_bed
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
LIST="${OUTDIR}/reference_lists/species_needing_chrom_reformat.txt"
OUTDIR_CHR="${OUTDIR}/reformatting_chr"
FNA_BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

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

  # Parse headers: accession is token 1; Chr_Num is the integer after 'chromosome '
  # Only emit headers that contain 'chromosome <num>'
  awk '
    BEGIN { OFS="," }
    /^>/ {
      acc=$1; sub(/^>/,"",acc)

      if (match($0, /chromosome[[:space:]]+([[:alnum:]_.-]+)/, m)) {
        chr=toupper(m[1])
        print acc, chr
      }
    }
  ' "$FNA" > "$OUT"

  if [[ ! -s "$OUT" ]]; then
    echo "WARN: no chromosome-number headers found in $FNA (wrote empty $OUT)" >&2
  else
    echo "Wrote $OUT" >&2
  fi

done < "$LIST"
```
Two ran into errors.
WARN: no chromosome-number headers found in /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Lepidochelys_kempii/Lepidochelys_kempii.fna (wrote empty reformatting_chr/Lepidochelys_kempii.chr_reformat.txt)
WARN: no chromosome-number headers found in /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/Lepidochelys_olivacea/Lepidochelys_olivacea.fna (wrote empty reformatting_chr/Lepidochelys_olivacea.chr_reformat.txt)

They have a format of "chromosome: 1" so they need a specific script:
```

OUTDIR_CHR="${OUTDIR}/reformatting_chr"
FNA_BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"
mkdir -p "$OUTDIR_CHR"

for SPECIES in Lepidochelys_olivacea Lepidochelys_kempii; do
  FNA="${FNA_BASE}/${SPECIES}/${SPECIES}.fna"
  OUT="${OUTDIR_CHR}/${SPECIES}.chr_reformat.txt"

  if [[ ! -s "$FNA" ]]; then
    echo "WARN: missing/empty $FNA (skipping)" >&2
    continue
  fi

  awk '
    BEGIN{ OFS="," }
    /^>/{
      acc=$1; sub(/^>/,"",acc)
      if (match($0, /chromosome:[[:space:]]+([[:alnum:]_.-]+)/, m)) {
        print acc, m[1]
      }
    }
  ' "$FNA" > "$OUT"

  echo "Wrote $OUT"
done
```
## Recode the chromosome of each bed file
```
LIST="${OUTDIR}/reference_lists/species_needing_chrom_reformat.txt"
MAPDIR="${OUTDIR}/reformatting_chr"

INBED_DIR="${OUTDIR}/usable_bed"
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
      if ( (acc ~ /^CM[0-9]+\.[0-9]+$/ || acc ~ /^NC_[0-9]+\.[0-9]+$/) &&
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
mv ${OUTDIR}/usable_bed_chrfixed/* ${OUTDIR}/usable_bed/
```


## Remove W and Y chromosomes from bed files
```
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/


INBED="${OUTDIR}/usable_bed"
INPEP="${OUTDIR}/usable_peptide"

OUTBED="${OUTDIR}/sexshared/usable_bed"
OUTPEP="${OUTDIR}/sexshared/usable_peptide"

mkdir -p "$OUTBED" "$OUTPEP"

shopt -s nullglob

# 1) Filter BEDs: drop rows where chrom (col 1) is Y*, W*
for bed in "$INBED"/*.bed; do
  outbed="${OUTBED}/$(basename "$bed")"

  awk -F'\t' 'BEGIN{OFS="\t"} $1 !~ /^(Y[0-9]*|W[0-9]*)$/ {print}' "$bed" > "$outbed"
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

## Remove Z and X chromosomes from bed files
```
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/


INBED="${OUTDIR}/usable_bed"
INPEP="${OUTDIR}/usable_peptide"

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