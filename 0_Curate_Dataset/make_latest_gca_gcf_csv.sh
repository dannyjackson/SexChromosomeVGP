#!/usr/bin/env bash
set -euo pipefail

IN="${1:?usage: $0 input.csv > out.csv}"
echo "Species,Accession_GCF,Accession_GCA"

# exists check via NCBI Datasets CLI
exists() {
  local acc="$1"
  datasets summary genome accession "$acc" --as-json-lines 2>/dev/null \
    | grep -Fq "\"accession\":\"$acc\""
}

# find highest existing version for prefix + id, starting from a seed version
latest_ver() {
  local prefix="$1" id="$2" seed="$3"
  local v="$seed"

  # if seed doesn't exist, walk down to find any
  while (( v >= 1 )) && ! exists "${prefix}_${id}.${v}"; do
    ((v--))
  done
  (( v >= 1 )) || { echo ""; return; }

  # walk up while next exists
  while exists "${prefix}_${id}.$((v+1))"; do
    ((v++))
  done
  echo "${prefix}_${id}.${v}"
}

# read CSV, find needed columns, then process rows
awk -F, '
NR==1{
  for(i=1;i<=NF;i++){
    if($i=="Scientific.Name") sci=i
    if($i=="Accession_num_main_haplotype") acc=i
  }
  next
}
{
  gsub(/^"|"$/,"",$sci); gsub(/^"|"$/,"",$acc)
  print $sci "\t" $acc
}' "$IN" | while IFS=$'\t' read -r sci acc; do
  [[ -n "${sci:-}" && -n "${acc:-}" ]] || continue

  species="$(echo "$sci" | tr ' ' '_' | sed -E 's/[^A-Za-z0-9_]+/_/g; s/_+/_/g; s/^_+|_+$//g')"

  # parse numeric id + seed version from either GCA/GCF input
  id="$(echo "$acc" | sed -E 's/^GC[AF]_([0-9]+)\.?.*$/\1/')"
  seed="$(echo "$acc" | sed -nE 's/^GC[AF]_[0-9]+\.([0-9]+)$/\1/p')"
  [[ -n "${seed:-}" ]] || seed=1

  gcf="$(latest_ver GCF "$id" "$seed")"
  gca="$(latest_ver GCA "$id" "$seed")"

  echo "${species},${gcf},${gca}"
done