#!/usr/bin/env bash
set -euo pipefail

: "${TSV:?Need TSV env var set}"
: "${GENOMEDIR:?Need GENOMEDIR env var set}"

OUT="${1:-sexchrom_accessions.csv}"
DEBUG="${DEBUG:-0}"   # set DEBUG=1 for stderr messages

col() {
  awk -F'\t' -v name="$1" '
    NR==1{
      for(i=1;i<=NF;i++) if($i==name){print i; exit}
      exit 1
    }' "$TSV"
}

C_SPECIES=$(col "Species")
C_CHR=$(col "Chromosome name")
C_GB=$(col "GenBank seq accession")
C_RS=$(col "RefSeq seq accession")

echo "Species,Chromosome,Accession" > "$OUT"

TMPDIR=$(mktemp -d)
trap 'rm -rf "$TMPDIR"' EXIT

# Use process substitution (NOT a pipe) so counters work in current shell
n_total=0
n_missing_fai=0
n_no_match=0
n_written=0

while IFS=$'\t' read -r species chr gb rs; do
  ((n_total++)) || true

  fai="${GENOMEDIR}/${species}/${species}.fna.fai"
  cache="${TMPDIR}/${species}.fai.acc"

  # Default to NA unless we prove a match
  acc="NA"

  # If the .fai exists, build/load cache
  if [[ -s "$fai" ]]; then
    if [[ ! -s "$cache" ]]; then
      cut -f1 "$fai" | sort -u > "$cache"
    fi

    # Try exact accession matches against the .fai sequence-name column
    if [[ -n "${gb:-}" && "$gb" != "NA" ]] && grep -qxF "$gb" "$cache"; then
      acc="$gb"
    elif [[ -n "${rs:-}" && "$rs" != "NA" ]] && grep -qxF "$rs" "$cache"; then
      acc="$rs"
    else
      ((n_no_match++)) || true
      if [[ "$DEBUG" == "1" && $n_no_match -le 20 ]]; then
        echo "DEBUG: no match in fai for species='$species' chr='$chr' GB='$gb' RS='$rs' (example fai names: $(head -n 3 "$cache" | paste -sd, -))" >&2
      fi
    fi
  else
    ((n_missing_fai++)) || true
    if [[ "$DEBUG" == "1" && $n_missing_fai -le 20 ]]; then
      echo "DEBUG: missing fai for species='$species' chr='$chr' path='$fai'" >&2
    fi
  fi

  # ALWAYS write one output row per TSV row
  printf '%s,%s,%s\n' "$species" "$chr" "$acc" >> "$OUT"
  ((n_written++)) || true

done < <(
  awk -F'\t' -v OFS='\t' \
      -v cs="$C_SPECIES" -v cc="$C_CHR" -v cgb="$C_GB" -v crs="$C_RS" \
      'NR>1 {print $cs,$cc,$cgb,$crs}' "$TSV"
)

if [[ "$DEBUG" == "1" ]]; then
  echo "DEBUG SUMMARY: total_rows=$n_total written=$n_written missing_fai=$n_missing_fai no_match=$n_no_match" >&2
fi

echo "Wrote: $OUT" >&2

