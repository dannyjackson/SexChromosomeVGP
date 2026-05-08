# Goal:
Create a plot that has chicken chromosomes on the X axis and unique SC systems on the Y. For each unique SC on the Y, there will be a horizontal line from start-stop of each syntenic block for a given clade. This will be min and max of all species that have that SC.
# Dataset:
## Sharks and Rays
```
Carcharodon_carcharias
Narcine_bancroftii
Hypanus_sabinus
Mobula_birostris
```
## Ray finned fish
```
Colia_mystus
Hoplias_malabaricus
Argentina_silus
Aulostomus_maculatus
Girardinichthys_multiradiatus
Echiichthys_vipera
Gasterosteus_aculeatus
```
## Mammals
```
# Monotremes
Tachyglossus_aculeatus
Ornithorhynchus_anatinus
# Therian
Sarcophilus_harrisii
```
## Lepidosaurs
```
Cyclura_pinguis
Vipera_latastei
Podarcis_raffonei
Shinisaurus_crocodilurus
Anniella_stebbinsi
Furcifer_pardalis
Anolis_sagrei
```
## Turtles
```
```
## Amphibians
```
# Pseudacris_triseriata
Hyla_sarda
# Ambystoma_mexicanum_x_Amblystoma_tigrinum
```
## Birds
```
Gallus_gallus
Taeniopygia_guttata
```
# Estimate chromosome lengths of chicken genome
```
# Set working directory

mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Vertebrate_Sex_Chromosomes
cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Vertebrate_Sex_Chromosomes

module load samtools
samtools faidx /data/Wilson_Lab/data/VGP_genomes_phase1/genomes/Gallus_gallus/ncbi_dataset/data/GCF_016700215.2/GCF_016700215.2_bGalGal1.pat.whiteleghornlayer.GRCg7w_genomic.fna
../Gallus_gallus/combined_phaseblks/
```
# Create gg_chr_len.csv
```
cat > gg_chr_len.csv <<'EOF'
1,196449156
2,149539284
3,110642502
4,90861225
5,59506338
6,36220557
7,36382834
8,29578256
9,23733309
10,20453248
11,19638187
12,20119077
13,17905061
14,15331188
15,12703657
16,2706039
17,11092391
18,11623896
19,10455293
20,14265659
21,6970754
22,4686657
23,6253421
24,6478339
25,3067737
26,5349051
27,5228753
28,5437364
29,726478
30,755666
31,2457334
32,125424
33,3839931
34,3469343
35,554126
36,358375
37,157853
38,667312
39,177356
Z,86044486
EOF
```
# set up environment
```
mkdir -p /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Vertebrate_Sex_Chromosomes

cd /data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/analyses/Genespace/All_Vertebrate_Sex_Chromosomes
```
# Create Ancestral_SC_list.txt
```
cat > Ancestral_SC_list.txt <<'EOF'
Carcharodon_carcharias,Chondricthyes
Scyliorhinus_canicula,
Heterodontus_francisci,Chondricthyes
Colia_mystus,Actinopterygii
Hoplias_malabaricus,Actinopterygii
Argentina_silus,Actinopterygii
Aulostomus_maculatus,Actinopterygii
Girardinichthys_multiradiatus,Actinopterygii
Echiichthys_vipera,Actinopterygii
Gasterosteus_aculeatus,Actinopterygii
Tachyglossus_aculeatus,Mammalia
Ornithorhynchus_anatinus,Mammalia
Sarcophilus_harrisii,Mammalia
Monodelphis_domestica,Mammalia
Dromiciops_gliroides,Mammalia
Cyclura_pinguis,Lepidosauria
Vipera_latastei,Lepidosauria
Podarcis_siculus,Lepidosauria
Podarcis_muralis,Lepidosauria
Podarcis_liolepis,Lepidosauria
Podarcis_bocagei,Lepidosauria
Podarcis_vaucheri,Lepidosauria
Podarcis_pityusensis,Lepidosauria
Podarcis_tiliguerta,Lepidosauria
Podarcis_raffoei,Lepidosauria
Podarcis_filfolensis,Lepidosauria
Podarcis_melisellensis,Lepidosauria
Podarcis_gaigeae,Lepidosauria
Podarcis_cretensis,Lepidosauria
Podarcis_erhardii,Lepidosauria
Shinisaurus_crocodilurus,Lepidosauria
Anniella_stebbinsi,Lepidosauria
Furcifer_pardalis,Lepidosauria
Anolis_sagrei,Lepidosauria
Hyla_sarda,Amphibia
Taeniopygia_guttata,Aves
Anser_anser,Aves
Larus_fuscus,Aves
EOF
```
# Save this as palette.Ancestral.txt
```
cat > palette.Ancestral.txt <<'EOF'
Monotremes,#E69F00
Therian_mammals,#E69F00
Aves,#00796B
Cyclura_pinguis,#80CBC4
Vipera_latastei,#80CBC4
Lacertids,#80CBC4
Shinisaurus_crocodilurus,#80CBC4
Anniella_stebbinsi,#80CBC4
Furcifer_pardalis,#80CBC4
Anolis_sagrei,#80CBC4
Hyla_sarda,#984EA3
Colia_mystus,#56B4E9
Hoplias_malabaricus,#56B4E9
Argentina_silus,#56B4E9
Aulostomus_maculatus,#56B4E9
Girardinichthys_multiradiatus,#56B4E9
Echiichthys_vipera,#56B4E9
Gasterosteus_aculeatus,#56B4E9
Chondricthyes,#0072B2
EOF
```
# Subset combined_phased_blks
```
awk -F',' '
BEGIN {
    OFS=","
}

# species list: Species,Clade
FNR==NR {
    gsub(/\r/, "", $1)
    gsub(/\r/, "", $2)

    if ($1 != "") {
        want[$1] = 1
        clade[$1] = $2
        order[++n] = $1
    }
    next
}

# header row of CSV
FNR==1 {
    for (i=1; i<=NF; i++) {
        gsub(/\r/, "", $i)
        col[$i] = i
    }
    next
}

{
    sp = $(col["partner"])

    if (!(sp in want)) next
    if ($(col["genome1"]) == "Gallus_gallus") next
    if ($(col["genome2"]) != "Gallus_gallus") next
    if ($(col["chr1"]) != "X" && $(col["chr1"]) != "X1" && $(col["chr1"]) != "Z" && $(col["chr1"]) != "Z1") next

    region = $(col["chr2"]) ":" $(col["startBp2"]) "-" $(col["endBp2"])

    key = sp SUBSEP region
    if (!(key in seen)) {
        seen[key] = 1
        if (sp in regions) {
            regions[sp] = regions[sp] ";" region
        } else {
            regions[sp] = region
        }
    }
}

END {
    print "Species","Clade","Gg_regions"
    for (i=1; i<=n; i++) {
        sp = order[i]
        print sp, clade[sp], regions[sp]
    }
}
' Ancestral_SC_list.txt ../Gallus_gallus/combined_phaseblks/Gallus_gallus_sexchrs.phasedBlks.csv > gg_regions.Ancestral.csv
```
# Plot the data
```
Rscript plot_gg_regions_by_clade.R gg_regions.Ancestral.csv gg_chr_len.csv palette.Ancestral.txt Ancestral.clade.png

Rscript plot_gg_regions_by_taxa.R gg_regions.Ancestral.csv gg_chr_len.csv palette.Ancestral.txt Ancestral.taxa.png
```
# Repeat with neo sex chromosomes
# Save this as palette.Neo.txt
```
Chondricthyes,#4797ba
Monotremes,#ffcc66
Therian_mammals,#ffcc66
Aves,#5fb2a6
```
# Subset combined_phased_blks
```
awk -F',' '
BEGIN {
    OFS=","
}

# species list: Species,Clade
FNR==NR {
    gsub(/\r/, "", $1)
    gsub(/\r/, "", $2)

    if ($1 != "") {
        want[$1] = 1
        clade[$1] = $2
        order[++n] = $1
    }
    next
}

# header row of CSV
FNR==1 {
    for (i=1; i<=NF; i++) {
        gsub(/\r/, "", $i)
        col[$i] = i
    }
    next
}

{
    sp = $(col["partner"])

    if (!(sp in want)) next
    if ($(col["genome1"]) == "Gallus_gallus") next
    if ($(col["genome2"]) != "Gallus_gallus") next
    if ($(col["chr1"]) != "X" && $(col["chr1"]) != "X1" && $(col["chr1"]) != "X2" && $(col["chr1"]) != "X3" && $(col["chr1"]) != "X4" && $(col["chr1"]) != "X5" && $(col["chr1"]) != "X6" && $(col["chr1"]) != "Z" && $(col["chr1"]) != "Z1" && $(col["chr1"]) != "Z2" && $(col["chr1"]) != "Z3" && $(col["chr1"]) != "Z4" && $(col["chr1"]) != "Z5" && $(col["chr1"]) != "Z6") next

    region = $(col["chr2"]) ":" $(col["startBp2"]) "-" $(col["endBp2"])

    key = sp SUBSEP region
    if (!(key in seen)) {
        seen[key] = 1
        if (sp in regions) {
            regions[sp] = regions[sp] ";" region
        } else {
            regions[sp] = region
        }
    }
}

END {
    print "Species","Clade","Gg_regions"
    for (i=1; i<=n; i++) {
        sp = order[i]
        print sp, clade[sp], regions[sp]
    }
}
' Neo_SC_list.txt ../Gallus_gallus/combined_phaseblks/Gallus_gallus_sexchrs.phasedBlks.csv > gg_regions.Neo.csv
```
# Plot the data
```
Rscript plot_gg_regions_by_clade.R gg_regions.Neo.csv gg_chr_len.csv palette.Neo.txt Neo.clade.png

Rscript plot_gg_regions_by_taxa.R gg_regions.Neo.csv gg_chr_len.csv palette.Neo.txt Neo.taxa.png

Rscript plot_gg_regions_by_clade.Neo_and_Ancestral.R \
  gg_regions.Neo.csv \
  gg_regions.Ancestral.csv \
  gg_chr_len.csv \
  palette.Neo.txt \
  palette.Ancestral.txt \
  Neo_and_Ancestral.clade.png
```