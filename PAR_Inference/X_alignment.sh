# Simone's code for PAR inference via alignment

srun ~/minimap2-2.28_x64-linux/minimap2 -x asm5 -c --eqx  -t 3 ${species_name}.${small_chr}chrom.fasta ${species_name}.${big_chr}chrom.fasta \
 > /scratch/sgable3/WGA_Benchmarking/minimap2/PAR_analysis/${species_name}_${small_chr}to${big_chr}.aln.paf


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
' "${base}.paf" > "$outbed"