# Inferring PAR from TE
# Calonectris_borealis W chromosome
# OZ122107.1
# 56726463 bp long

cd /scratch/dnjacks4/TE_PAR_Inference/analyses/spot_test/Calonectris_borealis/W
# Calonecteris borealis

BEDDIR=/scratch/sgable3/VGP_Projects/RepMask_out/RepMasker_BEDs/
META_FILE=/scratch/dnjacks4/TE_PAR_Inference/reference_lists/VGP_freeze_hap1_combined_sexchroms_seq_reports.tsv

grep 'OZ122107.1' ${BEDDIR}/Calonectris_borealis_GCA_964195595.2.repMask.out.bed > /scratch/dnjacks4/TE_PAR_Inference/datafiles/Calonectris_borealis.GCA_964195595.2.TE.W.bed

W_TE_BED=/scratch/dnjacks4/TE_PAR_Inference/datafiles/Calonectris_borealis.GCA_964195595.2.TE.W.bed

# Autosomes
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2}' "$TE_BED" > autosomes.te.bed

# W chromosome
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2}' "$W_TE_BED" > W.te.bed


# --------------------------
# W windows (your approach)
# --------------------------
echo -e 'OZ122107.1\t56726463' > wchrom.size
bedtools makewindows -g wchrom.size -w 50000 > W.50kb.bed

bedtools intersect -a W.50kb.bed -b W.te.bed -wa -wb |
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$7}' |
  bedtools groupby -g 1,2,3 -c 4 -o sum > W.TEbp.50kb.bed

# Add zero-TE W windows
bedtools intersect -a W.50kb.bed -b W.TEbp.50kb.bed -v |
  awk 'BEGIN{OFS="\t"} {print $1,$2,$3,0}' >> W.TEbp.50kb.bed


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

Rscript ../2x_plot_chr_repeats_bar.R $W_TE_BED 56726463 W.TE.png



R

auto <- read.table("autosomes.TEbp.50kb.bed",
                   col.names=c("chr","start","end","TEbp"))
Z <- read.table("W.TEbp.50kb.bed",
                col.names=c("chr","start","end","TEbp"))

auto$prop <- auto$TEbp / 50000
Z$prop <- Z$TEbp / 50000

wilcox.test(Z$prop, auto$prop, alternative="less")

# windowed test
null <- auto$prop

Z$p_empirical <- sapply(Z$prop, function(x) mean(null <= x))
Z$padj <- p.adjust(Z$p_empirical, method="fdr")
sig_Z <- subset(Z, padj < 0.05)

combined <- c(Z$prop, auto$prop)
labels <- c(rep("Z", length(Z$prop)), rep("A", length(auto$prop)))

obs <- mean(Z$prop) - mean(auto$prop)

perm <- replicate(10000, {
  l <- sample(labels)
  mean(combined[l=="Z"]) - mean(combined[l=="A"])
})

mean(perm <= obs)

# Plot it all
library(ggplot2)

# ---------------------------
# Prep (same as before)
# ---------------------------
null <- auto$prop
Z$p_empirical <- sapply(Z$prop, function(x) mean(null <= x))
Z$padj <- p.adjust(Z$p_empirical, method="fdr")
Z$sig <- Z$padj < 0.05

combined <- c(Z$prop, auto$prop)
labels <- c(rep("Z", length(Z$prop)), rep("A", length(auto$prop)))
obs <- mean(Z$prop) - mean(auto$prop)

set.seed(1)
perm <- replicate(10000, {
  l <- sample(labels)
  mean(combined[l=="Z"]) - mean(combined[l=="A"])
})
p_perm <- mean(perm <= obs)

Z$mid <- (Z$start + Z$end) / 2
Z$mb  <- Z$mid / 1e6
Z$mlog10_fdr <- -log10(pmax(Z$padj, 1e-300))

# ---------------------------
# Plot 1: TE proportion along Z
# ---------------------------
p1 <- ggplot(Z, aes(x = mb, y = prop)) +
  geom_line(alpha = 0.25) +
  geom_point(aes(shape = sig), size = 1.6, alpha = 0.8) +
  geom_hline(yintercept = median(auto$prop, na.rm = TRUE), linetype = 2) +
  labs(
    x = "Position on W (Mb)",
    y = "TE proportion per 50 kb window",
    shape = "FDR < 0.05",
    title = "TE density along W (50 kb windows)",
    subtitle = sprintf(
      "Global depletion permutation p = %.3g; dashed line = autosomal median",
      p_perm
    )
  ) +
  theme_classic()

ggsave(
  filename = "W_TE_density_50kb.png",
  plot = p1,
  width = 10,
  height = 4,
  units = "in",
  dpi = 300
)

# ---------------------------------
# Plot 2: significance along Z (FDR)
# ---------------------------------
p2 <- ggplot(Z, aes(x = mb, y = mlog10_fdr)) +
  geom_point(aes(shape = sig), size = 1.6, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(
    x = "Position on W (Mb)",
    y = expression(-log[10]("FDR")),
    shape = "FDR < 0.05",
    title = "Windows with significantly low TE density on W",
    subtitle = "Dashed line = FDR 0.05 threshold"
  ) +
  theme_classic()

ggsave(
  filename = "W_TE_depletion_FDR_50kb.png",
  plot = p2,
  width = 10,
  height = 4,
  units = "in",
  dpi = 300
)


scp "dnjacks4@login.sol.rc.asu.edu:/scratch/dnjacks4/TE_PAR_Inference/analyses/spot_test/Calonectris_borealis/W/*.png" .


