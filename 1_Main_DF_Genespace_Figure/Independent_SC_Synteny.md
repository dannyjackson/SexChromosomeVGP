# Plot Sex Chr Synteny to Chicken for subset of VGP genomes

source myconda
mamba activate genespace_py3.10

cat > species.txt <<'EOF'
Carcharodon carcharias
Coilia mystus
Hoplias malabaricus
Argentina silus
Aulostomus maculatus
Girardinichthys multiradiatus
Echiichthys vipera
Gasterosteus aculeatus
Hyla sarda
Ornithorhynchus anatinus
Homo sapiens
Dibamus smithi
Tiliqua scincoides
Shinisaurus crocodilurus
Anniella stebinsi
Bradypodion pumilum
Furcifer pardalis
Cyclura pinguis
Vipera latastei
Podarcis bocagei
Gallus gallus
EOF

Rscript plot_genespace.r Gallus_gallus
