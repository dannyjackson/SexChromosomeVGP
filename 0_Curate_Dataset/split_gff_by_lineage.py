#!/usr/bin/env python3
import csv, argparse, os, re, sys

def safe_name(s):
    # produce a safe filename component
    return re.sub(r'[^A-Za-z0-9._-]+', '_', s).strip('_')

def main():
    p = argparse.ArgumentParser(description="Subset GFF list into files per Lineage using REFERENCE CSV's Genus_species -> Lineage mapping.")
    p.add_argument("--ref", required=True, help="REFERENCE_FILE CSV")
    p.add_argument("--gff", required=True, help="GFF_LIST (one species per line, e.g. Genus_species)")
    p.add_argument("--outdir", default="gff_by_lineage", help="output directory")
    p.add_argument("--gencol", default="Genus_species", help="CSV column to match (default: Genus_species)")
    p.add_argument("--lineagecol", default="Lineage", help="CSV column to use for grouping (default: Lineage)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # read reference CSV, build mapping Genus_species -> Lineage
    mapping = {}
    with open(args.ref, newline='') as fh:
        reader = csv.DictReader(fh)
        cols = reader.fieldnames
        if args.gencol not in cols or args.lineagecol not in cols:
            print("ERROR: CSV does not contain required columns. Found columns: {}".format(cols), file=sys.stderr)
            sys.exit(2)
        for r in reader:
            key = r.get(args.gencol) or ""
            lineage = r.get(args.lineagecol) or ""
            key = key.strip()
            lineage = lineage.strip() if lineage is not None else ""
            if key:
                mapping[key] = lineage

    # prepare filehandles dictionary
    handles = {}
    unmapped = open(os.path.join(args.outdir, "UNMAPPED.txt"), "w")

    def get_handle(lineage):
        if lineage is None or lineage == "":
            return unmapped
        name = safe_name(lineage)
        fname = os.path.join(args.outdir, f"Lineage_{name}.txt")
        if fname not in handles:
            handles[fname] = open(fname, "w")
        return handles[fname]

    # read gff list and dispatch
    with open(args.gff) as fh:
        for line in fh:
            sp = line.strip()
            if not sp:
                continue
            # exact match first; allow case-insensitive fallback
            lineage = mapping.get(sp)
            if lineage is None:
                # try case-insensitive match
                match = next((v for k,v in mapping.items() if k.lower()==sp.lower()), None)
                lineage = match
            if lineage is None:
                # also try replacing underscores/spaces
                variant = sp.replace('_',' ')
                lineage = mapping.get(variant)
            if lineage is None:
                # last resort: try matching genus only (not recommended, but helpful)
                genus = sp.split('_')[0] if '_' in sp else None
                if genus:
                    match = next((v for k,v in mapping.items() if k.startswith(genus + "_") or k.startswith(genus + " ")), None)
                    lineage = match
            fh_out = get_handle(lineage)
            fh_out.write(sp + "\n")

    # close all handles
    for h in handles.values():
        h.close()
    unmapped.close()
    print("Done. Output files in:", os.path.abspath(args.outdir))

if __name__ == "__main__":
    main()
