#!/usr/bin/env python3

import sys
import re
import pandas as pd

# -----------------------------
# Regex patterns (case-insensitive)
# -----------------------------
PAT_CLR = re.compile(r'\b(CLR|RSII|Mecat|FALCON)\b', re.IGNORECASE)
PAT_HIFI = re.compile(r'\bHi[- ]?Fi\b', re.IGNORECASE)   # HiFi, Hi-Fi, Hi Fi
PAT_HIC = re.compile(r'\bHi[- ]?C\sphasing\b', re.IGNORECASE)
PAT_TRIO = re.compile(r'\btrio\b', re.IGNORECASE)       # trio, Trio-binning, TrioCanu
PAT_SOLO = re.compile(r'\bsolo\b', re.IGNORECASE)
PAT_CANU = re.compile(r'\bCanu\b', re.IGNORECASE)
PAT_HIFIASM = re.compile(r'\bhifiasm\b', re.I)
PAT_PACBIO = re.compile(r'\bpacbio\b', re.I)
# -----------------------------
# Classification logic
# Priority:
# 1) CLR
# 2) HiFi_trio_HiC
# 3) HiFi_Solo
# -----------------------------
def classify_row(row):
    # Combine all columns into one searchable string
    txt = " ".join(str(v) for v in row.values if pd.notna(v))

    # 1) CLR
    if PAT_CLR.search(txt):
        return "CLR"

    # 2 & 3) HiFi-based

    is_hifi_like = PAT_HIFI.search(txt) or (PAT_PACBIO.search(txt) and PAT_HIFIASM.search(txt))

    if is_hifi_like:
        if PAT_SOLO.search(txt) or PAT_CANU.search(txt):
            return "HiFi_Solo"
        if PAT_TRIO.search(txt) or PAT_HIC.search(txt):
            return "HiFi_trio_HiC"

    return "Unclassified"

# -----------------------------
# Main
# -----------------------------
def main():
    if len(sys.argv) != 3:
        sys.stderr.write(
            "Usage: python3 classify_seq_method.py <input.tsv> <output.tsv>\n"
        )
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    # Read TSV as strings
    df = pd.read_csv(infile, sep="\t", dtype=str)

    # Apply classification
    df["Classification"] = df.apply(classify_row, axis=1)

    # Write TSV
    df.to_csv(outfile, sep="\t", index=False)

    # -----------------------------
    # Print counts per group
    # -----------------------------
    print("\nClassification counts:")
    counts = df["Classification"].value_counts().sort_index()
    for cls, n in counts.items():
        print(f"{cls}\t{n}")

    print(f"\nTotal rows:\t{len(df)}")

if __name__ == "__main__":
    main()
