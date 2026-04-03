# Create requisite genespace directory structure
```
SEX_CHR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv"

BED_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_beds_noYW"
PEP_DIR="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/XZ_Synteny_Gallus_gallus/genespace_input/usable_peptides_noYW"

REF_BED="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/bed/Gallus_gallus.bed"
REF_PEP="/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/VGP_DF_Figure/Option2/peptide/Gallus_gallus.fa"

tail -n +2 "$SEX_CHR" | while IFS=',' read -r SPECIES CHR_TYPE ACC; do
  ACC="${ACC%$'\r'}"

  # Make simple per-species structure
  mkdir -p "${SPECIES}/bed" "${SPECIES}/peptide"

  # Symlink species inputs
  ln -sf "${BED_DIR}/${SPECIES}.bed" "${SPECIES}/bed/${SPECIES}.bed"
  ln -sf "${PEP_DIR}/${SPECIES}.fa"  "${SPECIES}/peptide/${SPECIES}.fa"

  # Symlink references (handy to keep alongside each species)
  ln -sf "$REF_BED" "${SPECIES}/bed/Gallus_gallus_REF.bed"
  ln -sf "$REF_PEP" "${SPECIES}/peptide/Gallus_gallus_REF.fa"
done
```