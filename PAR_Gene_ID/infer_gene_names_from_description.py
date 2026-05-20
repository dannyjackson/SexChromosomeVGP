#!/usr/bin/env python3

import argparse
import re
from collections import defaultdict


def is_loc_gene(gene_name):
    return gene_name.startswith("LOC")


def normalize_description(desc):
    """
    Normalize descriptions so entries like:
        anosmin-1-like
        anosmin 1
    can match.

    This intentionally removes common suffixes like '-like',
    lowercases, converts punctuation to spaces, and collapses whitespace.
    """
    desc = desc.lower().strip()

    # Decode common URL-style comma from GFF-derived descriptions
    desc = desc.replace("%2c", ",")

    # Remove trailing "-like" or " like"
    desc = re.sub(r"[-\s]+like$", "", desc)

    # Convert punctuation/separators to spaces
    desc = re.sub(r"[-_/(),;:]+", " ", desc)

    # Collapse whitespace
    desc = re.sub(r"\s+", " ", desc).strip()

    return desc


def read_table(path):
    rows = []

    with open(path, "r") as handle:
        first = handle.readline().rstrip("\n")
        fields = first.split("\t")

        has_header = (
            len(fields) >= 6
            and fields[0] == "Species"
            and fields[1] == "Chromosome"
        )

        if has_header:
            header = fields
        else:
            header = [
                "Species",
                "Chromosome",
                "StartPos",
                "StopPos",
                "GeneName",
                "GeneDescription",
            ]
            rows.append(fields)

        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))

    return header, rows


def main():
    parser = argparse.ArgumentParser(
        description="Infer LOC gene names from matching descriptions in x_par_genes.tsv."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="x_par_genes.tsv",
        help="Input TSV file. Default: x_par_genes.tsv",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="x_par_genes.with_inferred_names.tsv",
        help="Output TSV with InferredGeneName column.",
    )
    parser.add_argument(
        "-f",
        "--failed",
        default="x_par_genes.failed_to_rename.tsv",
        help="Output TSV of LOC genes that could not be renamed.",
    )
    args = parser.parse_args()

    header, rows = read_table(args.input)

    try:
        gene_col = header.index("GeneName")
        desc_col = header.index("GeneDescription")
    except ValueError:
        raise SystemExit(
            "ERROR: Input must contain GeneName and GeneDescription columns, "
            "or use the default 6-column format."
        )

    # Dictionary:
    # normalized_description -> set of known non-LOC gene names
    desc_to_gene_names = defaultdict(set)

    for row in rows:
        if len(row) <= max(gene_col, desc_col):
            continue

        gene = row[gene_col].strip()
        desc = row[desc_col].strip()

        if gene and desc and not is_loc_gene(gene):
            norm_desc = normalize_description(desc)
            if norm_desc:
                desc_to_gene_names[norm_desc].add(gene)

    inferred_rows = []
    failed_rows = []

    for row in rows:
        # Pad short rows so output stays rectangular
        if len(row) < len(header):
            row = row + [""] * (len(header) - len(row))

        gene = row[gene_col].strip()
        desc = row[desc_col].strip()
        norm_desc = normalize_description(desc)

        inferred = gene
        status = "original"

        if is_loc_gene(gene):
            candidates = sorted(desc_to_gene_names.get(norm_desc, []))

            if len(candidates) == 1:
                inferred = candidates[0]
                status = "renamed"
            elif len(candidates) > 1:
                inferred = ";".join(candidates)
                status = "ambiguous"
                failed_rows.append(row + [inferred, status, norm_desc])
            else:
                inferred = ""
                status = "failed"
                failed_rows.append(row + [inferred, status, norm_desc])

        inferred_rows.append(row + [inferred, status])

    with open(args.output, "w") as out:
        out.write("\t".join(header + ["InferredGeneName", "InferenceStatus"]) + "\n")
        for row in inferred_rows:
            out.write("\t".join(row) + "\n")

    with open(args.failed, "w") as out:
        out.write(
            "\t".join(
                header
                + ["InferredGeneName", "InferenceStatus", "NormalizedGeneDescription"]
            )
            + "\n"
        )
        for row in failed_rows:
            out.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()