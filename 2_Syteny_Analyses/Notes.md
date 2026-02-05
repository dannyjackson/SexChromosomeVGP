# Slack from Simone on Feb 4 2026
Heyo, Ying Chen generated a chicken to chorus frog alignment for the synteny plot: https://genomeark.s3.amazonaws.com/index.html?prefix=species/Pseudacris_triseriata/temp/

I wasn't sure if you would use an existing alignment or just re-align yourself so I asked for hap2 as well. Chr1 is the sex chrom, the larger chr1 is Y
[10:42 AM]These are the params she used for aligning:
```
REF=GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic.fna

QUERY=GCA_053478255.1_ASM5347825v1_genomic.split.fa 

minimap2 -ax asm20 -L $REF $QUERY -t 10 > chicken_GCF_016700215.2_frog_GCA_053478255.1.sam

samtools sort -@ 10 chicken_GCF_016700215.2_frog_GCA_053478255.1.sam -o chicken_GCF_016700215.2_frog_GCA_053478255.1.sorted.bam

minimap2 -x asm20 -L $REF $QUERY -t 10 > chicken_GCF_016700215.2_frog_GCA_053478255.1.paf
```
