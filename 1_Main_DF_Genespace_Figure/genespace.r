library(GENESPACE)

wd <- "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/formatted_annotations"

gpar <- init_genespace(
  wd = wd,
  path2mcscanx = "/vf/users/jacksondan/conda/envs/genespace/bin/")

out <- run_genespace(gpar, overwrite = F)