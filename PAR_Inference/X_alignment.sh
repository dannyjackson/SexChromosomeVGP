# Simone's code for PAR inference via alignment

srun ~/minimap2-2.28_x64-linux/minimap2 -x asm5 -c --eqx -k 23 -t 3 ${species_name}.${small_chr}chrom.fasta ${species_name}.${big_chr}chrom.fasta \
 > /scratch/sgable3/WGA_Benchmarking/minimap2/PAR_analysis/${species_name}_${small_chr}to${big_chr}.aln.paf

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/

base="$in"
infile="${base}.paf"

outbed="${base}.id98_5.len10k.refqry.bed"

# Filter and write BED:
# Output BED columns:
# ref_id  ref_start  ref_end  query_id  query_start  query_end  identity_percent  ref_span_len
awk -v thr=0.985 -v minlen=10000 -v OFS='\t' '
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
' "${base}.paf" > "PAR_annotations/$outbed"

# output continuous percent identity
mkdir -p continuous_percentID

for infile in *.paf; do
  base=$(basename "$infile" .paf)
  outbed="continuous_percentID/${base}.refqry.csv"

  # Write header
  printf "chrom_qry\tlen_qry\tbp_start_qry\tbp_end_qry\tpercent_identity_qry\tchrom_ref\tlen_ref\tbp_start_ref\tbp_end_ref\tpercent_identity_ref\n" \
    > "$outbed"

  # Process file
  awk -v OFS='\t' '
    function ident_from_de(de){ return (de=="" ? -1 : 1.0 - de) }
    function ident_from_core(m,a){ return (a>0 ? m/a : -1) }
    {
      # Extract de:f tag if present
      de=""
      for(i=13;i<=NF;i++)
        if($i ~ /^de:f:/){ split($i,x,":"); de=x[3]; break }

      # Compute identity
      id = ident_from_de(de)
      if (id < 0) id = ident_from_core($10,$11)

      pid = (id >= 0 ? id * 100 : "NA")

      printf "%s\t%d\t%d\t%d\t%s\t%s\t%d\t%d\t%d\t%s\n", \
             $1, $2, $3, $4, pid, \
             $6, $7, $8, $9, pid
    }
  ' "$infile" >> "$outbed"

  echo "Wrote $outbed"
done



# plot continuous percent ID
head /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/datafiles/minimap2/continuous_percentID/H

cd /data/Wilson_Lab/projects/VertebrateSexChr/jacksondan/analyses/PAR_inference/alignment/continuous_percentID

