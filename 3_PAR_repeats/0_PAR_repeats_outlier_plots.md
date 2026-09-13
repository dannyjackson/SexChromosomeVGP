# Analyze the relationship between PAR size and gene count to identify trends
## Set up environment
```
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/PAR_repeats/mammals


export REPEAT_DIR=/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/repeats
export GORILLA_REP=$REPEAT_DIR/Gorilla_gorilla_GCA_029281585.3_hap1.repeatMasker.XY.bed
export SHREW_REP=$REPEAT_DIR/Rhynchocyon_petersi_GCA_043290085.1_hap1.repeatMasker.XY.bed
export SIAMANG_REP=$REPEAT_DIR/Symphalangus_syndactylus_GCF_028878055.3_hap1.repeatMasker.XY.bed

export TOGA_DIR="/data/Wilson_Lab/data/TOGA2_Hiller/Homo_sapiens_hg38"
export GORILLA_GENES=$TOGA_DIR/Gorilla_gorilla__western_gorilla__HLgorGor7__GCA_029281585.3/query_annotation.gtf.gz
export SHREW_GENES=$TOGA_DIR/Rhynchocyon_petersi__Black_and_rufous_elephant_shrew__HLrhyPet1A__GCA_043290085.1/query_annotation.gtf.gz
export SIAMANG_GENES=$TOGA_DIR/Symphalangus_syndactylus__siamang__HLsymSyn4__GCA_028878055.3/query_annotation.gtf.gz

export GORILLA_PAR_GENES="Gorilla_gorilla,chr1_pat_hsa1:0-12039748"
export GORILLA_PAR_REP="Gorilla_gorilla,CM055469.2:0-12039748"
export GORILLA_PAR_PAF="Gorilla_gorilla,NC_073247.2:0-12039748"

export SHREW_PAR_GENES="Rhynchocyon_petersi,CM091802:0-20182068"
export SHREW_PAR_REP="Rhynchocyon_petersi,CM091802.1:0-20182068"
export SHREW_PAR_PAF="Rhynchocyon_petersi,CM091802.1:0-20182068"

export SIAMANG_PAR_GENES="Symphalangus_syndactylus,chrX_hap1:0-16994066"
export SIAMANG_PAR_REP="Symphalangus_syndactylus,NC_072447.2:0-16994066"
export SIAMANG_PAR_PAF="Symphalangus_syndactylus,NC_072447.2:0-16994066"

export PAF_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/continuous_percentID"
export GORILLA_PAF=$PAF_DIR/Gorilla_gorilla_YtoX.aln.refqry.csv
export SIAMANG_PAF=$PAF_DIR/Symphalangus_syndactylus_YtoX.aln.refqry.csv
export SHREW_PAF=$PAF_DIR/Rhynchocyon_petersi_YtoX.aln.refqry.csv

export PAF_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/datafiles/minimap2/"
export GORILLA_PAF=$PAF_DIR/Gorilla_gorilla_YtoX.aln.paf
export SIAMANG_PAF=$PAF_DIR/Symphalangus_syndactylus_YtoX.aln.paf
export SHREW_PAF=$PAF_DIR/Rhynchocyon_petersi_YtoX.aln.paf
```
## Analyze and plot using R
```
Rscript 1a_plot_repeats_PAR.R GORILLA
Rscript 1a_plot_repeats_PAR.R SIAMANG
Rscript 1a_plot_repeats_PAR.R SHREW

Rscript 1b_plot_repeats_PAR.ZNF.R GORILLA
Rscript 1b_plot_repeats_PAR.ZNF.R SIAMANG
Rscript 1b_plot_repeats_PAR.ZNF.R SHREW

Rscript 1c_plot_repeats_PAR.each_class.R GORILLA
Rscript 1c_plot_repeats_PAR.each_class.R SIAMANG
Rscript 1c_plot_repeats_PAR.each_class.R SHREW

Rscript 1d_anova.R
```