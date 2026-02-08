# Robust ALL-by-AssemblyTech plotting script (phylo-ordered)
library(data.table)
library(ggplot2)
library(patchwork)
library(ape)


tree_file <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/roadies_v1.1.16b.numbers.scientific.nwk"

SPECIES <- c("Homo_sapiens", "Gallus_gallus", "Anolis_sagrei", "Podarcis_raffonei", "Hyla_sarda", "Hoplias_malabaricus", "Narcine_bancroftii")

tr_full <- ape::read.tree(tree_file)

tip_key <- tr_full$tip.label


# Keep only tips in SPECIES (using the key)
keep_idx  <- tip_key %in% SPECIES
tr_pruned <- keep.tip(tr_full, tr_full$tip.label[keep_idx])

# Recompute key after pruning (important)
tip_key <- tip_key[keep_idx]

# Desired order for the CURRENT tip labels
desired <- tr_pruned$tip.label[order(match(tip_key, SPECIES))]

# Rotate nodes to try to achieve that tip order (topology-preserving)
tr_ordered <- rotateConstr(tr_pruned, desired)

# Quick check: tip order you’ll get when plotting
tr_ordered$tip.label

# --- plot (root/“past” on left, tips/“present” on right) and save to PDF ---


# --- plot: branch lengths + tips aligned to the right ---

# Force tips to align by making the tree ultrametric (tips all same distance from root)
# This preserves topology but may slightly adjust terminal/internal lengths to fit.
tr_ultra <- chronos(tr_ordered, lambda = 1)

pdf("species_ordered_tree_ultrametric_aligned.pdf", width = 3, height = 10)
par(mar = c(1, 1, 1, 1))
plot(
  tr_ultra,
  type = "phylogram",
  use.edge.length = TRUE,
  direction = "rightwards",
  show.tip.label = TRUE,
  cex = 0.8,
  no.margin = TRUE,
  edge.width = 2,
  lwd = 2
)
dev.off()
pdf("tree_branch_lengths.pdf", width=3, height=10)
plot(tr_ordered,
     type="phylogram",
     direction="rightwards",
     use.edge.length=TRUE,
     show.tip.label=TRUE,
     cex=0.8,
     no.margin = TRUE,
     edge.width = 2,
     lwd = 2)
dev.off()