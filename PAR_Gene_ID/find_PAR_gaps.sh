#!/usr/bin/env bash

module load bedtools seqkit

PAR_FILE="PAR.species_chr_region.txt"
base="/data/Wilson_Lab/data/VGP_genomes_phase1/symlinks"

outdir="PAR_gap_results"
mkdir -p "$outdir"/per_species "$outdir"/logs

all_gaps="$outdir/PAR.gaps.bed"
summary="$outdir/PAR.gaps.summary.tsv"
regions_bed="$outdir/PAR.regions.bed"

: > "$all_gaps"
: > "$regions_bed"

echo -e "species\tchrom\tPAR_start\tPAR_end\tPAR_len\tn_gaps\ttotal_gap_bp\tmax_gap_bp" > "$summary"

echo "Finding gaps within PAR regions..."

while IFS=, read -r species region_raw; do
    [[ -z "${species:-}" ]] && continue

    # Repair accidental quote typo, if present
    region="$(echo "$region_raw" | sed 's/"/:/g')"

    chrom="${region%%:*}"
    coords="${region#*:}"
    start="${coords%-*}"
    end="${coords#*-}"

    fasta="${base}/${species}/${species}.fna"

    if [[ ! -s "$fasta" ]]; then
        echo "WARNING: missing fasta for $species: $fasta" | tee -a "$outdir/logs/missing_fasta.log"
        continue
    fi

    if [[ "$region" == "$region_raw" ]]; then
        :
    else
        echo "WARNING: repaired malformed region for $species: $region_raw -> $region" | tee -a "$outdir/logs/repaired_regions.log"
    fi

    par_len=$((end - start))

    echo -e "${chrom}\t${start}\t${end}\t${species}" >> "$regions_bed"

    tmp_all_n="$outdir/per_species/${species}.all_N.bed"
    tmp_par="$outdir/per_species/${species}.PAR.gaps.raw.bed"
    tmp_merged="$outdir/per_species/${species}.PAR.gaps.bed"

    echo "  $species  $chrom:$start-$end"

    # Find all N positions/runs in the genome.
    # seqkit locate --bed emits BED-like 0-based intervals.
    seqkit locate --bed -m 0 -i -p "N" "$fasta" \
        > "$tmp_all_n"

    # Keep only N-runs overlapping the PAR interval, trim to PAR boundaries,
    # then merge adjacent/overlapping N intervals.
    awk -v c="$chrom" -v s="$start" -v e="$end" 'BEGIN{OFS="\t"} $1==c && $3>s && $2<e {
        a = ($2 < s ? s : $2)
        b = ($3 > e ? e : $3)
        if (b > a) print $1, a, b
    }' "$tmp_all_n" \
        | sort -k1,1 -k2,2n \
        | bedtools merge -d 1 -i - \
        > "$tmp_merged"

    # Add species name and relative PAR coordinates.
    # Output columns:
    # species chrom gap_start gap_end gap_len PAR_start PAR_end rel_gap_start rel_gap_end
    awk -v sp="$species" -v ps="$start" -v pe="$end" 'BEGIN{OFS="\t"} {
        print sp, $1, $2, $3, $3-$2, ps, pe, $2-ps, $3-ps
    }' "$tmp_merged" >> "$all_gaps"

    n_gaps=$(wc -l < "$tmp_merged" | tr -d ' ')
    total_gap_bp=$(awk '{sum += $3-$2} END{print sum+0}' "$tmp_merged")
    max_gap_bp=$(awk 'BEGIN{max=0} {len=$3-$2; if (len>max) max=len} END{print max+0}' "$tmp_merged")

    echo -e "${species}\t${chrom}\t${start}\t${end}\t${par_len}\t${n_gaps}\t${total_gap_bp}\t${max_gap_bp}" >> "$summary"

done < "$PAR_FILE"

echo "Done."
echo "Gap BED-like table: $all_gaps"
echo "Summary table:      $summary"
echo "PAR regions BED:    $regions_bed"