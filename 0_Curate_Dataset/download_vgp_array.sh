#!/usr/bin/env bash
#SBATCH --job-name=vgp_dl
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_logs/vgp_dl_%A_%a.out
#SBATCH --error=slurm_logs/vgp_dl_%A_%a.err

set -u
set -o pipefail

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  echo "ERROR: This script must be submitted as a SLURM array."
  echo "Example:"
  echo "  sbatch --array=1-N%25 download_vgp_array.sh"
  exit 1
fi

source myconda
conda activate ncbi_datasets

GENOME_DIR="/data/Wilson_Lab/data/VGP_genomes_phase1"
ACCESSION_LIST="/data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/referencelists/VGP_OrdinalList_Phase1Freeze_v1.2_Sept.30.2025_sex_chrs_HalfDeep_SCINKD.csv"
# we still request these via datasets but we'll select files by type after inspecting zips
FILES_TO_DOWNLOAD="gff3,rna,cds,protein,genome,seq-report"

mkdir -p "$GENOME_DIR" slurm_logs

OK_LOG="${GENOME_DIR}/download_ok.${SLURM_ARRAY_JOB_ID}.${SLURM_ARRAY_TASK_ID}.tsv"
FAIL_LOG="${GENOME_DIR}/download_fail.${SLURM_ARRAY_JOB_ID}.${SLURM_ARRAY_TASK_ID}.tsv"

printf "Scientific.Name\tAccessionUsed\tZip\tContains.Genome\tContains.GFF\tContains.RNA\tContains.CDS\tContains.Protein\n" > "$OK_LOG"
printf "Scientific.Name\tAccessionTried\tStatus\tNotes\n" > "$FAIL_LOG"

# --- find column indices from header (1-based) ---
header="$(head -n 1 "$ACCESSION_LIST")"

sci_col=$(
  awk -v FS=',' -v h="$header" 'BEGIN{
    n=split(h,a,FS);
    for(i=1;i<=n;i++) if(a[i]=="Scientific.Name"){print i; exit}
  }'
)
acc_col=$(
  awk -v FS=',' -v h="$header" 'BEGIN{
    n=split(h,a,FS);
    for(i=1;i<=n;i++) if(a[i]=="Accession_num_main_haplotype"){print i; exit}
  }'
)

if [[ -z "${sci_col:-}" || -z "${acc_col:-}" ]]; then
  echo "ERROR: Could not find required columns in header."
  echo "Need: Scientific.Name and Accession_num_main_haplotype"
  exit 1
fi

sanitize_dir() {
  echo "$1" \
    | tr ' ' '_' \
    | sed -E 's/[^A-Za-z0-9_.-]+/_/g; s/_+/_/g; s/^_+|_+$//g'
}

# Download an accession into $zipfile. Return codes:
# 0 = success (zip exists and non-empty)
# 1 = download failed or zip empty
try_download() {
  local accession="$1"
  local zipfile="$2"

  rm -f "$zipfile"  # remove partial from previous failed attempt to keep a clean state
  datasets download genome accession "$accession" \
    --include "$FILES_TO_DOWNLOAD" \
    --filename "$zipfile" \
    --no-progressbar

  local rc=$?
  if [[ $rc -ne 0 || ! -s "$zipfile" ]]; then
    rm -f "$zipfile"
    return 1
  fi
  return 0
}

# Inspect zip contents for presence flags. Outputs a 5-field tab line: genome<TAB>gff<TAB>rna<TAB>cds<TAB>protein
inspect_zip_flags() {
  local zipfile="$1"
  # list contents, tolower to ease matching
  local listing
  listing="$(unzip -l "$zipfile" 2>/dev/null | awk '{print tolower($0)}')"

  local has_genome=0
  local has_gff=0
  local has_rna=0
  local has_cds=0
  local has_protein=0

  # genome fasta patterns: .fna, .fa, _genomic.fna, .fa.gz, .fna.gz
  if echo "$listing" | grep -E '\.fna(\.gz)?$|\.fa(\.gz)?$|_genomic\.fna(\.gz)?$' >/dev/null 2>&1; then
    has_genome=1
  fi

  # gff patterns: .gff, .gff3, gz variants
  if echo "$listing" | grep -E '\.gff3?$|\.gff3?\.gz' >/dev/null 2>&1; then
    has_gff=1
  fi

  # rna/cds/protein: filename contains rna|cds|protein (case-insensitive) and has a fasta-like extension
  if echo "$listing" | grep -E 'rna.*\.(fa|fasta|fna)(\.gz)?|\.rna\.(fa|fasta|fna)(\.gz)?' >/dev/null 2>&1; then
    has_rna=1
  else
    # sometimes sequences are in files named *_rna, or contain "transcript"
    if echo "$listing" | grep -E 'transcript.*\.(fa|fasta|fna)(\.gz)?' >/dev/null 2>&1; then
      has_rna=1
    fi
  fi

  if echo "$listing" | grep -E 'cds.*\.(fa|fasta|fna|faa)(\.gz)?|\.cds\.(fa|fasta|fna|faa)(\.gz)?' >/dev/null 2>&1; then
    has_cds=1
  fi

  if echo "$listing" | grep -E 'protein.*\.(faa|fa|fasta)(\.gz)?|\.protein\.(faa|fa|fasta)(\.gz)?' >/dev/null 2>&1; then
    has_protein=1
  fi

  printf "%d\t%d\t%d\t%d\t%d\n" "$has_genome" "$has_gff" "$has_rna" "$has_cds" "$has_protein"
}

# --- pick the row for this array task ---
line_no=$((SLURM_ARRAY_TASK_ID + 1))
line="$(sed -n "${line_no}p" "$ACCESSION_LIST" || true)"
if [[ -z "${line:-}" ]]; then
  echo "No line at $line_no (array task ${SLURM_ARRAY_TASK_ID}). Exiting."
  exit 0
fi

SCI_NAME="$(echo "$line" | awk -v FPAT='([^,]+)|(\"[^\"]+\")' -v c="$sci_col" '{ gsub(/^"|"$/,"",$c); print $c }')"
ACC="$(echo "$line" | awk -v FPAT='([^,]+)|(\"[^\"]+\")' -v c="$acc_col" '{ gsub(/^"|"$/,"",$c); print $c }')"

SCI_NAME="${SCI_NAME%\"}"; SCI_NAME="${SCI_NAME#\"}"; SCI_NAME="$(echo "$SCI_NAME" | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')"
ACC="${ACC%\"}"; ACC="${ACC#\"}"; ACC="$(echo "$ACC" | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')"

if [[ -z "$SCI_NAME" || -z "$ACC" ]]; then
  echo "Empty Scientific.Name or Accession on line $line_no; exiting."
  exit 0
fi

SPECIES_DIR_NAME="$(sanitize_dir "$SCI_NAME")"
SPECIES_DIR="${GENOME_DIR}/${SPECIES_DIR_NAME}"
mkdir -p "$SPECIES_DIR"

echo "==> $SCI_NAME | original accession: $ACC"

# Parse accession into prefix, numeric id, and version if present.
if [[ "$ACC" =~ ^(GC[AF])_([0-9]+)(\.([0-9]+))?$ ]]; then
  acc_prefix="${BASH_REMATCH[1]}"   # GCA or GCF
  acc_id="${BASH_REMATCH[2]}"       # numeric identifier
  acc_ver="${BASH_REMATCH[4]:-1}"   # version if present; default 1
else
  acc_prefix="$(echo "$ACC" | cut -d'_' -f1)"
  acc_id="$(echo "$ACC" | cut -d'_' -f2 | cut -d'.' -f1)"
  acc_ver="1"
fi

acc_ver=$((acc_ver + 0))
start_ver=$acc_ver

# Variables to track which accession provided which filetype (most recent found)
genome_acc=""
gff_acc=""
rna_acc=""
cds_acc=""
protein_acc=""

genome_zip=""
gff_zip=""
rna_zip=""
cds_zip=""
protein_zip=""

found_genome=0
found_gff=0
found_rna=0
found_cds=0
found_protein=0

tried_list=()

# Loop versions newest -> 1, within each try GCF then GCA
for (( ver=start_ver; ver>=1; ver-- )); do
  for try_pref in GCF GCA; do
    accession_to_try="${try_pref}_${acc_id}.${ver}"
    ZIP="${SPECIES_DIR}/${SPECIES_DIR_NAME}_${accession_to_try}_dataset.zip"

    # avoid re-trying exact same string twice (in case acc_prefix was same as try_pref and user originally provided it later)
    if printf '%s\n' "${tried_list[@]}" | grep -x "$accession_to_try" >/dev/null 2>&1; then
      continue
    fi
    tried_list+=("$accession_to_try")

    echo "Trying $accession_to_try ..."
    if try_download "$accession_to_try" "$ZIP"; then
      # success download; inspect zip contents
      flags="$(inspect_zip_flags "$ZIP")" || flags="$(printf "0\t0\t0\t0\t0\n")"
      IFS=$'\t' read -r has_genome has_gff has_rna has_cds has_protein <<< "$flags"

      notes=""
      [[ $has_genome -eq 1 ]] && notes="${notes}GENOME;"
      [[ $has_gff -eq 1 ]] && notes="${notes}GFF;"
      [[ $has_rna -eq 1 ]] && notes="${notes}RNA;"
      [[ $has_cds -eq 1 ]] && notes="${notes}CDS;"
      [[ $has_protein -eq 1 ]] && notes="${notes}PROTEIN;"

      # Record presence in FAIL/OK logs for traceability per accession tried
      printf "%s\t%s\tDOWNLOAD_OK\t%s\n" "$SCI_NAME" "$accession_to_try" "$notes" >> "$OK_LOG"

      # For each filetype, if we haven't yet found a more recent one, record this accession
      # (we iterate newest -> oldest so first hit is the most recent)
      if [[ $has_genome -eq 1 && $found_genome -eq 0 ]]; then
        found_genome=1
        genome_acc="$accession_to_try"
        genome_zip="$ZIP"
      fi
      if [[ $has_gff -eq 1 && $found_gff -eq 0 ]]; then
        found_gff=1
        gff_acc="$accession_to_try"
        gff_zip="$ZIP"
      fi
      if [[ $has_rna -eq 1 && $found_rna -eq 0 ]]; then
        found_rna=1
        rna_acc="$accession_to_try"
        rna_zip="$ZIP"
      fi
      if [[ $has_cds -eq 1 && $found_cds -eq 0 ]]; then
        found_cds=1
        cds_acc="$accession_to_try"
        cds_zip="$ZIP"
      fi
      if [[ $has_protein -eq 1 && $found_protein -eq 0 ]]; then
        found_protein=1
        protein_acc="$accession_to_try"
        protein_zip="$ZIP"
      fi

      # If we've now found all 5 types, break out early
      if [[ $found_genome -eq 1 && $found_gff -eq 1 && $found_rna -eq 1 && $found_cds -eq 1 && $found_protein -eq 1 ]]; then
        echo "Found all desired filetypes at most recent versions. Stopping search."
        break 2
      fi

      # continue trying older versions until either all types found or versions exhausted
    else
      # download failed for this accession - log it
      printf "%s\t%s\tNO_DOWNLOAD\t-\n" "$SCI_NAME" "$accession_to_try" >> "$FAIL_LOG"
      echo "Could not download $accession_to_try"
    fi
  done
done

# As a last resort: try original ACC string if not already tried
if printf '%s\n' "${tried_list[@]}" | grep -x "$ACC" >/dev/null 2>&1; then
  :
else
  ZIP="${SPECIES_DIR}/${SPECIES_DIR_NAME}_${ACC}_dataset.zip"
  echo "Trying original accession string $ACC as last resort..."
  if try_download "$ACC" "$ZIP"; then
    flags="$(inspect_zip_flags "$ZIP")" || flags="$(printf "0\t0\t0\t0\t0\n")"
    IFS=$'\t' read -r has_genome has_gff has_rna has_cds has_protein <<< "$flags"
    notes=""
    [[ $has_genome -eq 1 ]] && notes="${notes}GENOME;"
    [[ $has_gff -eq 1 ]] && notes="${notes}GFF;"
    [[ $has_rna -eq 1 ]] && notes="${notes}RNA;"
    [[ $has_cds -eq 1 ]] && notes="${notes}CDS;"
    [[ $has_protein -eq 1 ]] && notes="${notes}PROTEIN;"
    printf "%s\t%s\tDOWNLOAD_OK\t%s\n" "$SCI_NAME" "$ACC" "$notes" >> "$OK_LOG"
    if [[ $has_genome -eq 1 && $found_genome -eq 0 ]]; then
      found_genome=1; genome_acc="$ACC"; genome_zip="$ZIP"
    fi
    if [[ $has_gff -eq 1 && $found_gff -eq 0 ]]; then
      found_gff=1; gff_acc="$ACC"; gff_zip="$ZIP"
    fi
    if [[ $has_rna -eq 1 && $found_rna -eq 0 ]]; then
      found_rna=1; rna_acc="$ACC"; rna_zip="$ZIP"
    fi
    if [[ $has_cds -eq 1 && $found_cds -eq 0 ]]; then
      found_cds=1; cds_acc="$ACC"; cds_zip="$ZIP"
    fi
    if [[ $has_protein -eq 1 && $found_protein -eq 0 ]]; then
      found_protein=1; protein_acc="$ACC"; protein_zip="$ZIP"
    fi
  else
    printf "%s\t%s\tNO_DOWNLOAD\t(original accession)\n" "$SCI_NAME" "$ACC" >> "$FAIL_LOG"
  fi
fi

# Report summary and unzip zips that contributed files
echo "Summary for $SCI_NAME:"
echo "  genome:   ${genome_acc:-NOT_FOUND}"
echo "  gff:      ${gff_acc:-NOT_FOUND}"
echo "  rna:      ${rna_acc:-NOT_FOUND}"
echo "  cds:      ${cds_acc:-NOT_FOUND}"
echo "  protein:  ${protein_acc:-NOT_FOUND}"

# For traceability, write a single line to OK_LOG summarizing chosen zips (if any)
# Compose values as accession:zip or NOT_FOUND
fmt_field() {
  local acc="$1"; local zip="$2"
  if [[ -n "$acc" && -n "$zip" ]]; then
    printf "%s:%s" "$acc" "$zip"
  else
    printf "NOT_FOUND"
  fi
}

# if any of the desired types were found, write a consolidated OK line
if [[ $found_genome -eq 1 || $found_gff -eq 1 || $found_rna -eq 1 || $found_cds -eq 1 || $found_protein -eq 1 ]]; then
  printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n" \
    "$SCI_NAME" \
    "$(printf "%s|%s|%s|%s|%s" "$(fmt_field "$genome_acc" "$genome_zip")" "$(fmt_field "$gff_acc" "$gff_zip")" "$(fmt_field "$rna_acc" "$rna_zip")" "$(fmt_field "$cds_acc" "$cds_zip")" "$(fmt_field "$protein_acc" "$protein_zip")")" \
    "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    "$found_genome" "$found_gff" "$found_rna" "$found_cds" "$found_protein" \
    >> "$OK_LOG"
else
  printf "%s\t%s\t%s\n" "$SCI_NAME" "${tried_list[*]}" "NO_DESIRED_FILES_FOUND" >> "$FAIL_LOG"
  echo "No desired files found for any tried accession/version."
  exit 1
fi

# Unzip only the zips that were used to provide at least one desired type
unzipped_any=0
for zip in "$genome_zip" "$gff_zip" "$rna_zip" "$cds_zip" "$protein_zip"; do
  if [[ -n "$zip" && -f "$zip" ]]; then
    echo "Unzipping $zip"
    unzip -o "$zip" -d "$(dirname "$zip")"
    unzipped_any=1
  fi
done

if [[ $unzipped_any -eq 1 ]]; then
  echo "Unzipped zips that supplied desired files."
else
  echo "No zips unzipped (no desired files found in any downloaded zips)."
fi

exit 0
