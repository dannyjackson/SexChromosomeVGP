# Inferring PAR from TE
# Pan troglodytes W chromosome
# NC_072422.2 or 	CM054458.2
# 56726463 bp long

cd /scratch/dnjacks4/TE_PAR_Inference/analyses/spot_test/Pan_troglodytes/Y

BEDDIR=/scratch/sgable3/VGP_Projects/RepMask_out/RepMasker_BEDs/
META_FILE=/scratch/dnjacks4/TE_PAR_Inference/reference_lists/VGP_freeze_hap1_combined_sexchroms_seq_reports.tsv

TE_BED=/scratch/dnjacks4/TE_PAR_Inference/datafiles/Pan_troglodytes_GCF_028858775.2.TE.autosomes.bed
Y_TE_BED=/scratch/dnjacks4/TE_PAR_Inference/datafiles/Pan_troglodytes_GCF_028858775.2.TE.Y.bed

grep 'NC_072422.2' ${BEDDIR}/Pan_troglodytes_GCF_028858775.2.repMask.out.bed \
  > $Y_TE_BED

grep -v 'NC_072422.2' ${BEDDIR}/Pan_troglodytes_GCF_028858775.2.repMask.out.bed | grep -v 'NC_072421.2' > $TE_BED

# Autosomes
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2}' "$TE_BED" > autosomes.te.bed

# Y chromosome
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2}' "$Y_TE_BED" > Y.te.bed


# --------------------------
# Y windows (your approach)
# --------------------------
echo -e 'NC_072422.2\t36447294' > ychrom.size
bedtools makewindows -g ychrom.size -w 10000 > Y.10kb.bed

bedtools intersect -a Y.10kb.bed -b Y.te.bed -wa -wb |
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' |
  bedtools groupby -g 1,2,3 -c 4 -o sum > Y.TEbp.10kb.bed

# Add zero-TE W windows
bedtools intersect -a Y.10kb.bed -b Y.TEbp.10kb.bed -v |
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,0}' >> Y.TEbp.10kb.bed


# ------------------------------------------
# Autosomes: infer "chr length" from TE_BED
# then drop the last window per chromosome
# ------------------------------------------

# 1) Infer per-chrom "length" as max end coordinate seen in TE_BED
#    (This is NOT true chromosome length, but enough to make windows;
#     dropping the last window avoids edge/truncation issues.)
awk 'BEGIN{OFS="\t"} {if($3>max[$1]) max[$1]=$3} END{for(c in max) print c, max[c]}' \
  autosomes.te.bed | sort -k1,1 > autosomes.inferred.sizes

# 2) Make 50kb windows for autosomes
bedtools makewindows -g autosomes.inferred.sizes -w 50000 > autosomes.50kb.all.bed

# 3) Drop the last window per chromosome (based on highest start)
#    This is robust even if bedtools made a truncated final window.

sort -k1,1 -k2,2n autosomes.50kb.all.bed |
  awk 'BEGIN{OFS="\t"}
       {chr=$1; start=$2; end=$3
        if(chr==prev_chr){
          print prev_chr, prev_start, prev_end
        }
        prev_chr=chr; prev_start=start; prev_end=end
       }
       END{}' > autosomes.50kb.bed


# 4) Sum TE bp per autosomal window
bedtools intersect -a autosomes.50kb.bed -b autosomes.te.bed -wa -wb |
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' |
  bedtools groupby -g 1,2,3 -c 4 -o sum > autosomes.TEbp.50kb.bed

# 5) Add zero-TE autosomal windows
bedtools intersect -a autosomes.50kb.bed -b autosomes.TEbp.50kb.bed -v |
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,0}' >> autosomes.TEbp.50kb.bed

Rscript ../../Calonectris_borealis/2x_plot_chr_repeats_bar.R $Y_TE_BED 36447294 Y.TE.png



scp "dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/TE_PAR_Inference/analyses/spot_test/Pan_troglodytes/Y/*.png" .


library(readr)
library(dplyr)
library(ggplot2)

# bed: chr  start  end  TEbp
df <- read_tsv(
  "Y.TEbp.10kb.bed",
  col_names = c("chr", "start", "end", "TEbp"),
  col_types = cols(
    chr   = col_character(),
    start = col_double(),
    end   = col_double(),
    TEbp  = col_double()
  )
)

# optional: pick one chromosome if multiple present
# df <- df %>% filter(chr == "NC_072422.2")

df <- df %>%
  arrange(chr, start) %>%
  mutate(
    pos_mb = start / 1e6,
    TE_kb  = TEbp / 1e3
  )

df_step <- df %>%
  group_by(chr) %>%
  arrange(start, .by_group = TRUE) %>%
  reframe(
    pos_mb = c(start/1e6, max(end)/1e6),
    TE_kb  = c(TEbp/1e3, last(TEbp)/1e3)
  )

step_plot <- ggplot(df_step, aes(pos_mb, TE_kb)) +
  geom_step(direction = "hv", linewidth = 0.6) +
  theme_classic() +
  labs(x = "Position (Mb)", y = "TE bp per 10 kb window (kb)") +
  facet_wrap(~ chr, scales = "free_x")

ggsave("Y.10kb.te.png", step_plot, width = 10, height = 5)

step_plot_zoom <- ggplot(df_step, aes(pos_mb, TE_kb)) +
  geom_step(direction = "hv", linewidth = 0.6) +
  theme_classic() +
  labs(
    x = "Position (Mb)",
    y = "TE bp per 50 kb window (kb)"
  ) +
  facet_wrap(~ chr, scales = "free_x") +
  coord_cartesian(ylim = c(0, 50), xlim = c(0,2.5))

ggsave("Y.10kb.te.zoom.png", step_plot_zoom, width = 10, height = 5)
