# VGP Sex Chromosome Synteny Analysis

## Make list of sex chromosomes from genomes with gff
```
python3 filter_sexchrom_to_gff.py
```

## Prepare all bed and peptide files
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/


source myconda
mamba activate genespace


cp -r /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/ .

cd symlinks

HAS_DIR="has_gff"
LACKS_DIR="lacks_gff"

mkdir -p "$HAS_DIR" "$LACKS_DIR"

shopt -s nullglob dotglob

for d in */ ; do
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
cd has_gff
find . -type f -iname "*lifton*" -print
# mv Artibeus_intermedius ../lacks_gff
# mv Artibeus_lituratus ../lacks_gff
mv Pseudacris_triseriata ../lacks_gff

```
## Remove scaffolds from FNA files
```
find . -mindepth 2 -maxdepth 2 -name "*.fna" -print0 |
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
## Check to make sure the chromosome numbers have changed, are reasonable, and that the last chromosome is not a scaffold
```
cd ..

echo "Species;NumChrom;LastChromLine" > chromcounts.contigs_only.txt

for f in has_gff/*/*.contigs_only.fna; do
  [ -f "$f" ] || continue
  species="$(basename "$(dirname "$f")")"

  # Count headers and capture last header
  awk -v sp="$species" '
    BEGIN{n=0; last=""}
    /^>/ {n++; last=$0}
    END{
      if(n==0) {print sp ";0;\"\""}
      else     {print sp ";" n ";\"" last "\""}
    }
  ' "$f" >> chromcounts.contigs_only.txt
done


echo "Species;NumChrom;LastChromLine" > chromcounts.all.txt

for f in has_gff/*/*.fna; do
  [[ "$f" == *.contigs_only.fna ]] && continue
  [ -f "$f" ] || continue

  species="$(basename "$(dirname "$f")")"

  awk -v sp="$species" '
    BEGIN{n=0; last=""}
    /^>/ {n++; last=$0}
    END{
      if(n==0) {print sp ";0;\"\""}
      else     {print sp ";" n ";\"" last "\""}
    }
  ' "$f" >> chromcounts.all.txt
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
cd has_gff
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
  f="./${sp}/${sp}.fna"
  out="./${sp}/${sp}.contigs_only.fna"

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

grep '>' ./Xenentodon_cancila/Xenentodon_cancila.contigs_only.fna
grep '>' ./Pholidichthys_leucotaenia/Pholidichthys_leucotaenia.contigs_only.fna
grep '>' ./Porphyrio_hochstetteri/Porphyrio_hochstetteri.contigs_only.fna
grep '>' ./Vipera_latastei/Vipera_latastei.contigs_only.fna
grep '>' ./Gastrophryne_carolinensis/Gastrophryne_carolinensis.contigs_only.fna
grep '>' ./Thomomys_bottae/Thomomys_bottae.contigs_only.fna
grep '>' ./Fundulus_diaphanus/Fundulus_diaphanus.contigs_only.fna
grep '>' ./Lycaon_pictus/Lycaon_pictus.contigs_only.fna
grep '>' ./Morphnus_guianensis/Morphnus_guianensis.contigs_only.fna
grep '>' ./Guaruba_guaruba/Guaruba_guaruba.contigs_only.fna
grep '>' ./Trichechus_inunguis/Trichechus_inunguis.contigs_only.fna

# check if these made it (likely they don't have annotations)
Note: Only Molosus molossus does
```
## only scaffolds to start with -- a newer GCA is choromosome level, so use that compared to scaffold level GCF
Rousettus_aegyptiacus # replace with GCA_014176215.2 (edited 2025, whereas GCF_014176215.1 was edited 2020)
Molossus_molossus # replace with GCA_014108415.2 (edited 2025, whereas GCF_014108415.1 was edited 2020)
Monodon_monocero # replace with GCA_005190385.4 (edited 2026, whereas GCF_005190385.1 was edited 2019)
Myotis_myotis # replace with GCA_014108235.2 (edited 2026, whereas	GCF_014108235.1 was edited 2019)
Phascolarctos_cinereus # replace with GCA_003287225.3 (edited 2025, whereas GCF_003287225.1 was edited 2022)
Pipistrellus_kuhlii # replace with GCA_014108245.2 (edited 2025, whereas GCF_014108245.1 was edited 2020)
```
## Filter ${SPECIES}/${SPECIES}.gff to just regions kept in file ${SPECIES}/${SPECIES}.contigs_only.fna
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/symlinks/has_gff

FASTA_ROOT="."   # or set to the parent directory that contains ${SPECIES}/${SPECIES}.contigs_only.fna
GFF_ROOT="."     # or set to the parent directory that contains ${SPECIES}/${SPECIES}.gff
OUT_SUFFIX=".contigs_only.gff"   # output name: ${SPECIES}/${SPECIES}${OUT_SUFFIX}

shopt -s nullglob
for fasta in "${FASTA_ROOT}"/*/*.contigs_only.fna; do
  species="$(basename "$(dirname "$fasta")")"
  gff="${GFF_ROOT}/${species}/${species}.gff"
  out="${GFF_ROOT}/${species}/${species}${OUT_SUFFIX}"

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
for gff in */*.contigs_only.gff; do
  species="$(basename "$(dirname "$gff")")"
  lines="$(wc -l < "$gff")"
  printf "%s\t%s\n" "$species" "$lines"
done | sort -k2,2n

# rename gff with contig.only
for gff in */*.contigs_only.gff; do
  species="$(basename "$(dirname "$gff")")"
  echo ${species}
  mv $gff ${species}/${species}.gff
done 

```
## Translate nucleic acid sequences for all species
```
GENOME_DIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/symlinks/has_gff

shopt -s nullglob

for d in "$GENOME_DIR"/*/; do
  SPECIES="$(basename "${d%/}")"

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
    cd "$d"
    transeq -sequence "${SPECIES}.cds" -outseq "${SPECIES}.translated.cds"
  )
done

```
## Prepare all bed and peptide files
```
# mv /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input_backup

mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input

R

library(GENESPACE)

repo <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/symlinks/has_gff"

# get immediate subdir names (species)
SPECIES <- sort(basename(list.dirs(repo, full.names = TRUE, recursive = FALSE)))
parsedPaths <- parse_annotations(
  rawGenomeRepo = repo, 
  genomeDirs = SPECIES,
  genomeIDs = SPECIES,
  gffString = "gff",
  faString = "translated.cds",
  presets = "ncbi", 
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input")

```
# Move empty bed and peptide files to "empty_file" dir for further troubleshooting
```
mkdir -p empty_files/bed empty_files/peptide

# empty beds/
find bed -type f -empty -print0 | xargs -0 -r mv -t empty_files/bed

# empty peptide/
find peptide -type f -empty -print0 | xargs -0 -r mv -t empty_files/peptide


mkdir -p short_files/bed short_files/peptide

# short_bed/
mkdir -p short_files/bed

find bed -type f -exec sh -c '
  f="$1"
  [ "$(wc -l < "$f")" -lt 2000 ] && mv -- "$f" short_files/bed/
' sh {} \;

# short_peptide/
find peptide -type f -exec sh -c '
  f="$1"
  [ "$(wc -l < "$f")" -lt 10000 ] && mv -- "$f" short_files/peptide/
' sh {} \;

mv bed usable_bed
mv peptide usable_peptide
```


# Create a list of species that did not run properly 
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input

ls -1 empty_files/bed short_files/bed \
  | grep -vE '^(empty_files/bed:|short_files/bed:)$' \
  | sed 's/\.bed$//' | sort -u > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/reference_lists/species_error.parse_annotations.txt
```
# Run parse_annotations on them again, using the protein_id field instead
I troubleshot this on Ammospiza maritima, a random species that did not work with the above parse annotations script. I then applied the same script to all species that did not run properly and, miraculously, this script worked for all of them. Very grateful that I didn't have to individually troubleshoot them all!
```
R

library(GENESPACE)

repo <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/symlinks/has_gff"

SPECIES <- readLines("/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/reference_lists/species_error.parse_annotations.txt")
SPECIES <- SPECIES[nzchar(SPECIES)]

parsedPaths2 <- parse_annotations(
  rawGenomeRepo = repo,
  genomeDirs = SPECIES,
  genomeIDs  = SPECIES,
  gffString  = "gff",
  faString   = "translated.cds",
  genespaceWd = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input",
  gffIdColumn = "protein_id",
  headerSep = "protein_id=",
  headerEntryIndex = 2,
  headerStripText = "\\].*",   # remove everything after the protein_id value
)
```
# Confirm that no empty or short beds remain
```
rm empty_files/*/*
rm short_files/*/*


# empty bed/
find bed -type f -empty -print0 | xargs -0 -r mv -t empty_files/bed

# empty peptide/
find peptide -type f -empty -print0 | xargs -0 -r mv -t empty_files/peptide

# short_bed
find bed -type f -exec sh -c '
  f="$1"
  [ "$(wc -l < "$f")" -lt 2000 ] && mv -- "$f" short_files/bed/
' sh {} \;

# short_peptide/
find peptide -type f -exec sh -c '
  f="$1"
  [ "$(wc -l < "$f")" -lt 10000 ] && mv -- "$f" short_files/peptide/
' sh {} \;


mv bed/* usable_bed/
mv peptide/* usable_peptide
```

# Reformat bed files to have chromosome numbers, not accession numbers
## Create a list of bed files requiring reformatting
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input

BED_DIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_bed
for f in "$BED_DIR"/*.bed; do
  line=$(grep -m1 -vE '^\s*($|#|track|browser)' "$f" || true)
  [ -z "$line" ] && continue
  chrom=$(awk '{print $1}' <<<"$line")
  if ! [[ "$chrom" =~ ^([0-9]+|X|Y|M|MT)$ ]]; then
    basename "$f" .bed
  fi
done > "species_needing_chrom_reformat.txt"
```
## Create a reformating index for each that needs it (pseudocode)
```
for each species in "species_needing_chrom_reformat.txt", create a SPECIES.chr_reformat.txt that is:
Accession,Chr_Num

Use files in format /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/${SPECIES}/${SPECIES}.fna to get this information.

For instance, Acanthopagrus_schlegelii:
SPECIES=Acanthopagrus_schlegelii
head /data/Wilson_Lab/data/VGP_genomes_phase1/symlinks/${SPECIES}/${SPECIES}.fna
>CM085472.1 Acanthopagrus schlegelii isolate SG-2024 chromosome 1, whole genome shotgun sequence

Output to file reformatting_chr/Acanthopagrus_schlegelii.chr_reformat.txt:
CM085472.1,1
```
### Code:
```

LIST="species_needing_chrom_reformat.txt"
OUTDIR="reformatting_chr"
FNA_BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

mkdir -p "$OUTDIR"

while IFS= read -r SPECIES; do
  # skip blanks/comments
  [[ -z "${SPECIES// /}" ]] && continue
  [[ "$SPECIES" =~ ^# ]] && continue

  FNA="${FNA_BASE}/${SPECIES}/${SPECIES}.fna"
  OUT="${OUTDIR}/${SPECIES}.chr_reformat.txt"

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

OUTDIR="reformatting_chr"
FNA_BASE="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"
mkdir -p "$OUTDIR"

for SPECIES in Lepidochelys_olivacea Lepidochelys_kempii; do
  FNA="${FNA_BASE}/${SPECIES}/${SPECIES}.fna"
  OUT="${OUTDIR}/${SPECIES}.chr_reformat.txt"

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
LIST="species_needing_chrom_reformat.txt"
MAPDIR="reformatting_chr"

INBED_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_bed"
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
mv usable_bed_chrfixed/* usable_bed/
```
# Filter bed and peptide files to just the X or Z chromosome
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms

SEX_CHR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv"

BED_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_bed"
PEP_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_peptide"

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/bed/Gallus_gallus.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/peptide/Gallus_gallus.fa"

declare -A seen_outbed

tail -n +2 "$SEX_CHR" | while IFS=',' read -r SPECIES CHR_TYPE ACC; do
  ACC="${ACC%$'\r'}"

  out_bed="${CHR_TYPE}/${SPECIES}/bed/${SPECIES}.bed"

  # HARD ERROR if we'd write this same bed file again (violates "one chromosome per bed")
  if [[ -n "${seen_outbed[$out_bed]+x}" ]]; then
    echo "ERROR: ${out_bed} would be written more than once. CSV has duplicate SPECIES+CHR_TYPE or multiple accessions for same chrom." >&2
    echo "       Offending key: SPECIES=${SPECIES} CHR_TYPE=${CHR_TYPE}" >&2
    exit 1
  fi
  seen_outbed["$out_bed"]=1

  mkdir -p "${CHR_TYPE}/${SPECIES}/bed" "${CHR_TYPE}/${SPECIES}/peptide"

  ln -sf "$REF_BED" "${CHR_TYPE}/${SPECIES}/bed/Gallus_gallus_REF.bed"
  ln -sf "$REF_PEP" "${CHR_TYPE}/${SPECIES}/peptide/Gallus_gallus_REF.fa"

  # one chromosome per bed: overwrite (>)
  grep -E "^((${ACC})|(${CHR_TYPE}))([[:space:]]|$)" "${BED_DIR}/${SPECIES}.bed" > "$out_bed"

  # subset full peptide fasta
  # IDs from col 4 (unique, non-empty)
  awk 'BEGIN{OFS="\t"} $4!=""{print $4}' "$out_bed" | sort -u > "${CHR_TYPE}/${SPECIES}/peptide/keep_ids.txt"

  # Subset FASTA by header name
  seqkit grep -n -f "${CHR_TYPE}/${SPECIES}/peptide/keep_ids.txt" \
  "${PEP_DIR}/${SPECIES}.fa" \
  > "${CHR_TYPE}/${SPECIES}/peptide/${SPECIES}.fa"

done
```
## Organize by chromosome type and check that all BED files still contain some lines
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/


mkdir -p sex_shared sex_limited

mv X* sex_shared
mv Z* sex_shared
mv Y* sex_limited
mv W* sex_limited
ROOT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared

find "$ROOT" -type f -name "*.bed" -print0 \
  | xargs -0 -I{} bash -lc 'c=$(sort -u "{}" | wc -l); printf "%10d\t%s\n" "$c" "{}"' \
  | sort -n | head


find "$ROOT" -type f -name "*.fa" -print0 \
  | xargs -0 -I{} bash -lc 'c=$(sort -u "{}" | wc -l); printf "%10d\t%s\n" "$c" "{}"' \
  | sort -n | head
```
### These two do not have sex chromosomes assigned (filter individually)
```
# Narcine_bancroftii (12) and Scyliorhinus_canicula (28)

mkdir -p \
    /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/X/Narcine_bancroftii/ \
    /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/X/Scyliorhinus_canicula/

# filter bed files
grep -E "^((${ACC})|(12))([[:space:]]|$)" /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_bed/Narcine_bancroftii.bed > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/X/Narcine_bancroftii/bed/Narcine_bancroftii.bed

grep -E "^((${ACC})|(28))([[:space:]]|$)" /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_bed/Scyliorhinus_canicula.bed > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/X/Scyliorhinus_canicula/bed/Scyliorhinus_canicula.bed

# filter peptide files
PEP_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_peptide"
## Narcine_bancroftii
out_bed="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/X/Narcine_bancroftii/bed/Narcine_bancroftii.bed"
SPECIES=Narcine_bancroftii

awk 'BEGIN{OFS="\t"} $4!=""{print $4}' "$out_bed" | sort -u > "sex_shared/X/${SPECIES}/peptide/keep_ids.txt"

seqkit grep -n -f "sex_shared/X/${SPECIES}/peptide/keep_ids.txt" \
"${PEP_DIR}/${SPECIES}.fa" \
> "sex_shared/X/${SPECIES}/peptide/${SPECIES}.fa"


## Scyliorhinus_canicula
out_bed="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/X/Scyliorhinus_canicula/bed/Scyliorhinus_canicula.bed"
SPECIES=Scyliorhinus_canicula

awk 'BEGIN{OFS="\t"} $4!=""{print $4}' "$out_bed" | sort -u > "sex_shared/X/${SPECIES}/peptide/keep_ids.txt"

seqkit grep -n -f "sex_shared/X/${SPECIES}/peptide/keep_ids.txt" \
"${PEP_DIR}/${SPECIES}.fa" \
> "sex_shared/X/${SPECIES}/peptide/${SPECIES}.fa"

```
# Run genespace on them all
## Make list of species and chr_type
```
awk -F',' 'NR>1 && $2 ~ /^(X|Z)/ {print $1 "\t" $2}' \
/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv > /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv
```
# Add chicken to each subdir
```

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/X

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/bed/Gallus_gallus.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/peptide/Gallus_gallus.fa"

for SPECIES in */ ; do
  SPECIES="${SPECIES%/}"

  [[ -d "${SPECIES}/bed" && -d "${SPECIES}/peptide" ]] || continue

  ln -sf "$REF_BED" "${SPECIES}/bed/Gallus_gallus_REF.bed"
  ln -sf "$REF_PEP" "${SPECIES}/peptide/Gallus_gallus_REF.fa"
done
```
## Write submission array (save as genespace_to_chicken_array.sh):
```
#!/bin/bash
#SBATCH --job-name=gs_chicken
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

PAIRFILE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv
SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/scripts/genespace.synteny_to_chicken.r
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/sexshared_chickenWG_plots

# Grab the Nth line (array index) and parse "Species Chromosome"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PAIRFILE}")
SPECIES=$(awk '{print $1}' <<< "${line}")
CHR_TYPE=$(awk '{print $2}' <<< "${line}")

WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/${CHR_TYPE}/${SPECIES}

export WORKING_DIR="${WORKDIR}"

cd "${OUTDIR}"

source myconda
mamba activate genespace

Rscript "${SCRIPT}" "${CHR_TYPE}" "${SPECIES}" "${OUTDIR}/"

for r in "${WORKDIR}"/orthofinder/Results*/Orthogroup_Sequences; do
  [ -d "$r" ] || continue
  out="${r%/}.tar.gz"
  tar -C "$(dirname "$r")" -cf - "$(basename "$r")" | bgzip -c > "$out"
done

# Reduce file number generated by genespace
# make globs that match nothing expand to nothing
shopt -s nullglob

# tar+bgzip each Orthogroup_Sequences dir, then delete it
for r in "${WORKDIR}"/orthofinder/Results*/Orthogroup_Sequences; do
  [ -d "$r" ] || continue
  out="${r%/}.tar.gz"
  tar -C "$(dirname "$r")" -cf - "$(basename "$r")" | bgzip -c > "$out" && rm -rf "$r"
done

# 2) One combined tar+bgzip of ALL Single_Copy_Orthologue_Sequences (across Results*), written to ${OUTDIR}/${SPECIES}/
for r in "${WORKDIR}"/orthofinder/Results*/Single_Copy_Orthologue_Sequences; do
  [ -d "$r" ] || continue
  out="${r%/}.tar.gz"
  tar -C "$(dirname "$r")" -cf - "$(basename "$r")" | bgzip -c > "$out" && rm -rf "$r"
done

```
## Submit batch job
```
# create output directory for plots
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/sexshared_chickenWG_plots

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/scripts

N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv)
sbatch --array=1-${N} genespace_to_chicken_array.sh
# sbatch --array=1-1 genespace_to_chicken_array.sh
```
# After running:
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/sex_chroms/sex_shared/sexshared_chickenWG_plots
mkdir pngs pdfs
mv *png pngs
mv *pdf pdfs
```

# Run genespace with the whole genomes 
## Remove W and Y chromosomes from bed files
```

BASE="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input"

INBED="${BASE}/usable_bed"
INPEP="${BASE}/usable_peptide"

OUTBED="${BASE}/usable_beds_noYW"
OUTPEP="${BASE}/usable_peptides_noYW"

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
# Create requisite genespace directory structure
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace

SEX_CHR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv"

BED_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_beds_noYW"
PEP_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_peptides_noYW"

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/bed/Gallus_gallus.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/peptide/Gallus_gallus.fa"

tail -n +2 "$SEX_CHR" | while IFS=',' read -r SPECIES CHR_TYPE ACC; do
  ACC="${ACC%$'\r'}"

  # Make simple per-species structure
  mkdir -p "${SPECIES}/bed" "${SPECIES}/peptide"

  # Symlink species inputs
  ln -sf "${BED_DIR}/${SPECIES}.bed" "${SPECIES}/bed/${SPECIES}.bed"
  ln -sf "${PEP_DIR}/${SPECIES}.fa"  "${SPECIES}/peptide/${SPECIES}.fa"

  # Symlink references (handy to keep alongside each species)
  ln -sf "$REF_BED" "${SPECIES}/bed/Gallus_gallus_REF.bed"
  ln -sf "$REF_PEP" "${SPECIES}/peptide/Gallus_gallus_REF.fa"
done


```
# Confirm contigs but not scaffolds are kept
```
awk '{print $1}' Buteo_buteo/bed/Buteo_buteo.bed | sort -u 
```
# Submission script
```

#!/bin/bash
#SBATCH --job-name=gs_chicken_wg
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

SLURM_ARRAY_TASK_ID=115
PAIRFILE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv
SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/scripts/genespace.synteny_to_chicken.wholegenome.r
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/plots

# Grab the Nth line (array index) and parse "Species Chromosome"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PAIRFILE}")
SPECIES=$(awk '{print $1}' <<< "${line}")
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/${SPECIES}

export WORKING_DIR="${WORKDIR}"

cd "${OUTDIR}"

source myconda
mamba activate genespace

Rscript "${SCRIPT}" "${SPECIES}" "${OUTDIR}/"

# Reduce file number generated by genespace
# make globs that match nothing expand to nothing
shopt -s nullglob

# tar+bgzip each Orthogroup_Sequences dir, then delete it
for r in "${WORKDIR}"/orthofinder/Results*/Orthogroup_Sequences; do
  [ -d "$r" ] || continue
  out="${r%/}.tar.gz"
  tar -C "$(dirname "$r")" -cf - "$(basename "$r")" | bgzip -c > "$out" && rm -rf "$r"
done

# 2) One combined tar+bgzip of ALL Single_Copy_Orthologue_Sequences (across Results*), written to ${OUTDIR}/${SPECIES}/
for r in "${WORKDIR}"/orthofinder/Results*/Single_Copy_Orthologue_Sequences; do
  [ -d "$r" ] || continue
  out="${r%/}.tar.gz"
  tar -C "$(dirname "$r")" -cf - "$(basename "$r")" | bgzip -c > "$out" && rm -rf "$r"
done

```
## Submit batch job
```
# create output directory for plots
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/plots

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/scripts

N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv)
sbatch --array=1-${N} genespace_to_chicken_array.wholegenome.sh
# sbatch --array=1-1 genespace_to_chicken_array.wholegenome.onesp.2.sh
```
## Plot just synteny with sex chrom
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/plots_sexchr

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/scripts

N=$(wc -l < /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv)
sbatch --array=1-${N} genespace_to_chicken_array.wholegenome.sexchr.sh
sbatch --array=1-1 genespace_to_chicken_array.wholegenome.sexchr.oneind.sh
# ...Hipposideros_larvatus; 109

```
### Sbatch script
```
#!/bin/bash
#SBATCH --job-name=gs_chicken_wg
#SBATCH --output=slurm_output/%x_%A_%a.out
#SBATCH --error=slurm_output/%x_%A_%a.err
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --gres=lscratch:500

set -euo pipefail

PAIRFILE=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/species.sexchr.for_genespace_array.csv
SCRIPT=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/scripts/genespace.synteny_to_chicken.wholegenome.sexchr.r
OUTDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/plots_sexchr
SLURM_ARRAY_TASK_ID=20

# Grab the Nth line (array index) and parse "Species Chromosome"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PAIRFILE}")
SPECIES=$(awk '{print $1}' <<< "${line}")
SEXCHR=$(awk '{print $2}' <<< "${line}")
WORKDIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/genespace/${SPECIES}

export WORKING_DIR="${WORKDIR}"

cd "${OUTDIR}"

source myconda
mamba activate genespace

Rscript "${SCRIPT}" "${SPECIES}" "${OUTDIR}/" ${SEXCHR}
```


## Create a table of syntenic regions to chicken
### Sex-shared to chicken WG
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/all/

```
# Combine all *_phasedBlks.csv files into one CSV (one header), adding provenance columns.

library(readr)
library(dplyr)
library(purrr)
library(stringr)

roots <- c("genespace")

files <- unlist(lapply(roots, function(r) {
  if (dir.exists(r)) list.files(r, pattern = "_phasedBlks\\.csv$", recursive = TRUE, full.names = TRUE) else character()
}), use.names = FALSE) %>% sort()

stopifnot(length(files) > 0)

parse_meta <- function(path) {
  # Example: Z/Accipiter_gentilis/riparian/Accipiter_gentilis_phasedBlks.csv
  parts <- strsplit(path, .Platform$file.sep, fixed = TRUE)[[1]]
  tibble(
    topDir = parts[1] %||% NA_character_,
    species = parts[2] %||% NA_character_,
    subdir = parts[3] %||% NA_character_,
    sourceFile = path
  )
}

# Read all, keep union of columns (fills missing with NA), add provenance cols

df <- map_dfr(files, function(f) {
  dat <- read_csv(
    f,
    col_types = cols(.default = col_character()),
    progress = FALSE,
    show_col_types = FALSE
  )
  bind_cols(parse_meta(f), dat)
})

write_csv(df, "all_phasedBlks_combined.csv")
cat("Wrote all_phasedBlks_combined.csv with", nrow(df), "rows from", length(files), "files\n")

# rename instances of "CP162266.1" to X
df$chr1[df$chr1 == "CP162266.1"] <- "X"
df$chr2[df$chr2 == "CP162266.1"] <- "X"
  
# Filter to just rows where a sex chromosome is in either chr 1 or chr 2
df2 <- df %>%
  filter(grepl("[X-Zx-z]", chr1) | grepl("[X-Zx-z]", chr2))

# remove rows that are chicken Z > other species 
df3 <- df2 %>%
  filter(!(genome1 == "Gallus_gallus_REF" & chr1 == "Z" & !grepl("^(Z|X)", chr2))) %>%
  filter(!(genome2 == "Gallus_gallus_REF" & chr2 == "Z" & !grepl("^(Z|X)", chr1)))

VGP_DATAFRAME <- read_csv("/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/VGP_OrdinalList_Phase1Freeze_v1.1_sex_chroms_seqCov_HalfDeep_SCINKD_May8.25_DwnldDec16.25.csv")
# rename Canis lupus baileyi to Canis lupus
VGP_DATAFRAME$Scientific.Name[VGP_DATAFRAME$Scientific.Name == "Canis lupus baileyi"] <- "Canis lupus"

# rename to Guaruba_guaruba to Guaruba guarouba in df
df3$species[df3$species == "Guaruba_guaruba"] <- "Guaruba_guarouba"

# Join df3 and VGP_DATAFRAME on df3$species and VGP_DATAFRAME$Scientific.Name for columns Lineage, Superorder


df4 <- df3 %>%
  mutate(species_join = gsub("_", " ", species)) %>%
  left_join(
    VGP_DATAFRAME %>% 
      transmute(Scientific.Name, Lineage, Superorder,
                Scientific_join = trimws(Scientific.Name)),
    by = c("species_join" = "Scientific_join")
  ) %>%
  select(-species_join)

unmatched <- df4 %>% filter(is.na(Lineage) & is.na(Superorder)) %>% distinct(species)
unmatched


write_csv(df4, "sexchrs.phasedBlks.csv")
unique_chrs <- unique(df4[c("Lineage", "Superorder", "genome1","genome2","chr1","chr2")])
write_csv(unique_chrs, "unique_chrs.csv")

unique_chrs <- read.csv("unique_chrs.csv")

# condense to just one instance of matching, retaining info if a chr is only matched once not twice
gg_name <- "Gallus_gallus_REF"   # change if your column uses a different exact string

out <- unique_chrs %>%
  # keep rows where chicken is in exactly one of the genome columns
  filter(xor(genome1 == gg_name, genome2 == gg_name)) %>%
  mutate(
    Species = if_else(genome1 == gg_name, genome2, genome1),
    SpChr   = if_else(genome1 == gg_name, chr2,   chr1),

    # chicken chromosome when chicken is genome1 -> chicken chr is chr1
    gg_chr1 = if_else(genome1 == gg_name, chr1, NA_character_),

    # chicken chromosome when chicken is genome2 -> chicken chr is chr2
    gg_chr2 = if_else(genome2 == gg_name, chr2, NA_character_)
  ) %>%
  group_by(Species, SpChr) %>%
  summarise(
    gg_chr1 = list(sort(unique(na.omit(gg_chr1)))),
    gg_chr2 = list(sort(unique(na.omit(gg_chr2)))),

    GgChr_both      = list(intersect(gg_chr1[[1]], gg_chr2[[1]])),
    GgChr_chr1_only = list(setdiff(gg_chr1[[1]], gg_chr2[[1]])),
    GgChr_chr2_only = list(setdiff(gg_chr2[[1]], gg_chr1[[1]])),
    .groups = "drop"
  ) %>%
  select(Species, SpChr, GgChr_both, GgChr_chr1_only, GgChr_chr2_only)


out_csv <- out %>%
  mutate(across(starts_with("GgChr"),
                ~ vapply(.x, \(v) paste(v, collapse = ","), character(1))))

readr::write_csv(out_csv, "syntentic_chromosomes.sex_shared.Gallus_gallus_REF.csv")


# Plot it on a phylogeny


# Plot it on the phylogeny
## Blanks for empty?

library(data.table)
library(ggplot2)
library(patchwork)
library(ape)
library(dplyr)
library(tidyr)
library(purrr)


out <- read.csv("syntentic_chromosomes.sex_shared.Gallus_gallus_REF.csv")

tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"

tr_full <- ape::read.tree(tree_file)

tip_key <- tr_full$tip.label

tr_ultra <- chronos(tr_full, lambda = 1)


# out: your tibble with Species and list-cols
# tr_ultra: your ultrametric tree

# union the three list-cols

out_long <- out %>%
  mutate(
    gg_all = map2(
      map2(GgChr_both, GgChr_chr1_only, ~ union(.x, .y)),
      GgChr_chr2_only,
      ~ union(.x, .y)
    )
  ) %>%
  select(Species, gg_all) %>%
  unnest(gg_all) %>%
  filter(!is.na(gg_all), gg_all != "") %>%
  distinct(Species, gg_all) %>%
  rename(GgChr = gg_all) %>%
  # --- NEW: split combos like "22,Z" into separate rows "22" and "Z"
  separate_rows(GgChr, sep = ",") %>%
  mutate(GgChr = str_trim(GgChr)) %>%
  filter(GgChr != "") %>%
  distinct(Species, GgChr)

# order chicken chromosomes: numeric first, then others (Z/W/etc)
chr_levels <- out_long %>%
  distinct(GgChr) %>%
  mutate(
    is_num = suppressWarnings(!is.na(as.integer(GgChr))),
    numval = suppressWarnings(as.integer(GgChr))
  ) %>%
  arrange(desc(is_num), numval, GgChr) %>%
  pull(GgChr)

# full grid for all tips x all chrs
grid_df <- tidyr::expand_grid(
  Species = out_long$Species,
  GgChr   = chr_levels
) %>%
  left_join(out_long %>% mutate(present = 1L), by = c("Species", "GgChr")) %>%
  mutate(present = if_else(is.na(present), 0L, present))





# Plot tree

pdf("tree.pdf", width = 3, height = 10)

# give more right margin for the dot grid + top margin for labels
par(mar = c(2, 2, 6, 12))

# 1) Compute how much horizontal space you need for the dot grid
#    (dx * nchr) plus some padding.
tmp <- plot(tr_ultra, plot = FALSE)
nchr <- length(chr_levels)

# spacing between chromosome columns
dx <- max(tmp$x.lim) * 0.02

# extra horizontal room to the right of the tree (in tree x-units)
grid_width <- dx * (nchr + 6)

# 2) Plot tree but "squish it left" by expanding xlim to the right
plot(tr_ultra,
     cex = 0.1,
     edge.width = 0.25,
     no.margin = TRUE,
     x.lim = c(0, max(tmp$x.lim) + grid_width))

lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tip_y <- lp$yy[1:Ntip(tr_ultra)]
tips  <- tr_ultra$tip.label

# start x position for the grid (a bit to the right of the tree)
x0 <- max(tmp$x.lim) + dx * 3

# full tips/y from the plotted tree
lp   <- get("last_plot.phylo", envir = .PlotPhyloEnv)
tips <- tr_ultra$tip.label
tip_y <- lp$yy[1:Ntip(tr_ultra)]

# subset: tips that appear in grid_df
tips_with_data <- intersect(tips, unique(grid_df$Species))
idx_data <- match(tips_with_data, tips)
y_data <- tip_y[idx_data]

# within each chromosome column:
for (j in seq_along(chr_levels)) {
  chrj <- chr_levels[j]
  xj <- x0 + (j - 1) * dx

  present_species <- grid_df %>%
    dplyr::filter(GgChr == chrj, present == 1L) %>%
    dplyr::pull(Species)

  is_present_data <- tips_with_data %in% present_species

  # light-gray circles ONLY for species present in grid_df
  points(rep(xj, length(y_data)), y_data,
         pch = 21, cex = 0.05, col = "lightgray", bg = "white")

  # black filled circles for present==1 among those species
  points(rep(xj, sum(is_present_data)), y_data[is_present_data],
         pch = 21, cex = 0.05, bg = "black")
}

# 2) Column labels at the top
usr <- par("usr")

# Put labels slightly ABOVE the top of the plotting region
y_lab <- 585
text(
  x = x0 + (seq_along(chr_levels) - 1) * dx,
  y = y_lab,
  labels = chr_levels,
  srt = 90,                 # rotate
  adj = c(0, 0.5),          # anchor at bottom of rotated text
  xpd = NA,                 # allow drawing into margin area
  cex = 0.2                 # tweak as needed
)

dev.off()




































df <- df3
gallus <- "Gallus_gallus_REF"

# Helper to safely compute bp length
bp_len <- function(start, end) {
  s <- suppressWarnings(as.double(start))
  e <- suppressWarnings(as.double(end))
  abs(e - s) + 1
}



bp_len <- function(start, end) {
  s <- suppressWarnings(as.double(start))
  e <- suppressWarnings(as.double(end))
  abs(e - s) + 1
}

# VGP lookup table: "Genus species" -> "Genus_species"
vgp_lut <- VGP_DATAFRAME %>%
  transmute(
    TaxonKey = str_replace_all(str_squish(Scientific.Name), " ", "_"),
    Lineage,
    Superorder
  ) %>%
  distinct(TaxonKey, .keep_all = TRUE)


summary_df <- df %>%
  # only Gallus vs non-Gallus comparisons
  filter((genome1 == gallus & genome2 != gallus) | (genome2 == gallus & genome1 != gallus)) %>%
  mutate(
    Species = if_else(genome1 == gallus, genome2, genome1),
    Chromosome = if_else(genome1 == gallus, chr1, chr2),

    BP_Gallus = if_else(genome1 == gallus,
                        bp_len(startBp1, endBp1),
                        bp_len(startBp2, endBp2)),
    BP_Query = if_else(genome1 == gallus,
                       bp_len(startBp2, endBp2),
                       bp_len(startBp1, endBp1)),

    nHits_Gallus = suppressWarnings(as.double(if_else(genome1 == gallus, nHits1, nHits2))),
    nHits_Query  = suppressWarnings(as.double(if_else(genome1 == gallus, nHits2, nHits1))),

    # join key for VGP (underscored binomial)
    TaxonKey = Species
  ) %>%
  # drop reciprocal duplicates (optional; remove if you want both directions counted)
  distinct(Species, Chromosome, blkID, startBp1, endBp1, startBp2, endBp2, nHits1, nHits2, .keep_all = TRUE) %>%
  group_by(Species, Chromosome, TaxonKey) %>%
  summarise(
    TotalBP_Gallus = sum(BP_Gallus, na.rm = TRUE),
    nHits_Gallus   = sum(nHits_Gallus, na.rm = TRUE),
    TotalBP_Query  = sum(BP_Query, na.rm = TRUE),
    nHits_Query    = sum(nHits_Query, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(vgp_lut, by = "TaxonKey") %>%
  relocate(Lineage, Superorder, .after = Species) %>%
  select(Species, Lineage, Superorder, Chromosome,
         TotalBP_Gallus, nHits_Gallus, TotalBP_Query, nHits_Query) %>%
  arrange(Lineage, Superorder, Species, Chromosome)


# Optional: show which Species didn't match VGP after stripping _REF
unmatched <- summary_df %>%
  filter(is.na(Lineage) | is.na(Superorder)) %>%
  distinct(Species)
if (nrow(unmatched) > 0) {
  cat("\nUnmatched Species (after stripping _REF for lookup):\n")
  print(unmatched)
}
# Umatched:
1 Canis_lupus    
2 Guaruba_guaruba

summary_df <- summary_df %>%
  mutate(
    Lineage = case_when(
      Species == "Canis_lupus" ~ "Mammal",
      Species == "Guaruba_guaruba" ~ "Bird",
      TRUE ~ Lineage
    )
  )

summary_df <- summary_df %>%
  mutate(
    Superorder = case_when(
      Species == "Canis_lupus" ~ "Laurasiatheria",
      Species == "Guaruba_guaruba" ~ "Psittaciformes",
      TRUE ~ Superorder
    )
  )

write_csv(summary_df, "Gallus_vs_species_window_sums_with_lineage.csv")
cat("Wrote Gallus_vs_species_window_sums_with_lineage.csv with", nrow(summary_df), "rows\n")

























```
# Make a plot NOT WORKING YET
```
# Plot synteny blocks on Gallus Z (0..86,044,486 bp) with species on Y

library(readr)
library(dplyr)
library(purrr)
library(stringr)

library(ggplot2)

gallus <- "Gallus_gallus_REF"
chr_target <- "Z"
Z_len <- 86044486

to_num <- function(x) suppressWarnings(as.numeric(x))

synteny_plot_df <- df %>%
  # keep Gallus vs non-Gallus
  filter((genome1 == gallus & genome2 != gallus) | (genome2 == gallus & genome1 != gallus)) %>%
  mutate(
    Species = if_else(genome1 == gallus, genome2, genome1),

    # Gallus chromosome for this row
    GallusChr = if_else(genome1 == gallus, chr1, chr2),

    # Gallus interval for this row
    g_start = if_else(genome1 == gallus, startBp1, startBp2),
    g_end   = if_else(genome1 == gallus, endBp1,   endBp2),

    g_start = to_num(g_start),
    g_end   = to_num(g_end),

    x_start = pmin(g_start, g_end, na.rm = TRUE),
    x_end   = pmax(g_start, g_end, na.rm = TRUE)
  ) %>%
  # only Gallus Z
  filter(as.character(GallusChr) == chr_target) %>%
  # drop reciprocal duplicates (optional but usually desirable)
  distinct(Species, GallusChr, blkID, x_start, x_end, .keep_all = TRUE) %>%
  filter(!is.na(x_start), !is.na(x_end)) %>%
  mutate(Species = factor(Species, levels = sort(unique(as.character(Species)))))

# Ensure A->Z (top to bottom depends on your preference; this is the common "A at top")
synteny_plot_df <- synteny_plot_df %>%
  mutate(Species = factor(as.character(Species), levels = sort(unique(as.character(Species)))))

p <- ggplot(synteny_plot_df, aes(y = Species)) +
  geom_segment(aes(x = x_start, xend = x_end, yend = Species), linewidth = 0.4, alpha = 0.8) +
  coord_cartesian(xlim = c(0, Z_len), clip = "off") +
  scale_y_discrete(limits = levels(synteny_plot_df$Species)) +  # enforce the order
  scale_x_continuous(labels = scales::comma) +
  labs(
    x = "Gallus_gallus Z position (bp)",
    y = "Species",
    title = "Synteny blocks to Gallus_gallus Z"
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 7),
    plot.margin = margin(t = 15, r = 10, b = 15, l = 10)  # <- THIS fixes top/bottom cutoff
  )

p <- p +
  geom_vline(xintercept = 0,  linewidth = 0.2, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = Z_len, linewidth = 0.2, linetype = "dashed", colour = "grey60")

# give the output a touch more vertical room too (optional but helps)
ggsave(
  "synteny_to_gallus_Z.png",
  p,
  width = 10,
  height = max(4, 0.22 * nlevels(synteny_plot_df$Species)), # slightly larger multiplier
  dpi = 300
)


ggsave("synteny_to_gallus_Z.png", p, width = 10, height = max(4, 0.2 * nlevels(synteny_plot_df$Species)), dpi = 300)
```