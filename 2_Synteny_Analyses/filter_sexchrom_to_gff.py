#!/usr/bin/env python3
import csv

sex_path = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.csv"
gff_path = "/data/Wilson_Lab/data/VGP_genomes_phase1/reference_lists/gff_file_list.tsv"
out_path = "/data/Wilson_Lab/projects/VGP_Phase_1_Sex_Chr_Project/jacksondan/referencelists/sexchrom_accessions.with_gff.csv"

# 1) Read species present in the GFF list (TSV)
keep_species = set()
with open(gff_path, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        sp = (row.get("Species") or "").strip()
        if sp:
            keep_species.add(sp)

# 2) Filter sexchrom_accessions.csv (CSV) by Species
with open(sex_path, newline="") as fin, open(out_path, "w", newline="") as fout:
    r = csv.DictReader(fin)
    w = csv.DictWriter(fout, fieldnames=r.fieldnames)
    w.writeheader()

    kept = 0
    total = 0
    for row in r:
        total += 1
        sp = (row.get("Species") or "").strip()
        if sp in keep_species:
            w.writerow(row)
            kept += 1

print(f"Done. Kept {kept}/{total} rows -> {out_path}")

