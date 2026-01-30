#!/usr/bin/env python3
import argparse
import json
import os
import sys
from glob import glob

def safe_get(d, path, default=""):
    """Get nested dict value by dotted path, return default if missing."""
    cur = d
    for key in path.split("."):
        if isinstance(cur, dict) and key in cur:
            cur = cur[key]
        else:
            return default
    return cur if cur is not None else default

def normalize(s):
    if s is None:
        return ""
    s = str(s).strip()
    # collapse whitespace
    return " ".join(s.split())

def find_reports(root):
    pattern = os.path.join(root, "*", "ncbi_dataset", "data", "assembly_data_report.jsonl")
    return sorted(glob(pattern))

def species_from_path(report_path):
    # .../genomes/<SPECIES>/ncbi_dataset/data/assembly_data_report.jsonl
    return os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(report_path))))

def parse_jsonl(report_path):
    records = []
    with open(report_path, "r", encoding="utf-8") as f:
        for line_no, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"WARNING: JSON parse error in {report_path}:{line_no}: {e}", file=sys.stderr)
    return records

def extract_methods(records):
    """
    Return (sequencing_methods, assembly_methods, assembly_comments)
    as sorted unique lists.
    """
    seq_set = set()
    asm_set = set()
    cmt_set = set()

    for rec in records:
        # Sequencing
        seq = safe_get(rec, "assemblyInfo.sequencingTech", "")
        if not seq:
            seq = safe_get(rec, "sequencingTech", "")
        seq = normalize(seq)
        if seq:
            seq_set.add(seq)

        # Assembly method
        asm = safe_get(rec, "assemblyInfo.assemblyMethod", "")
        if not asm:
            asm = safe_get(rec, "assemblyMethod", "")
        asm = normalize(asm)
        if asm:
            asm_set.add(asm)

        # Assembly comments
        cmt = safe_get(rec, "assemblyInfo.comments", "")
        if not cmt:
            cmt = safe_get(rec, "comments", "")
        cmt = normalize(cmt)
        if cmt:
            cmt_set.add(cmt)

    return sorted(seq_set), sorted(asm_set), sorted(cmt_set)

def main():
    ap = argparse.ArgumentParser(
        description="Scan VGP genomes directory for NCBI assembly_data_report.jsonl and output a TSV summary."
    )
    ap.add_argument(
        "--root",
        default="/data/Wilson_Lab/data/VGP_genomes_phase1/genomes",
        help="Root genomes directory (default: /data/Wilson_Lab/data/VGP_genomes_phase1/genomes)",
    )
    ap.add_argument(
        "--out",
        default="species_sequencing_assembly_methods.tsv",
        help="Output TSV path (default: species_sequencing_assembly_methods.tsv)",
    )
    args = ap.parse_args()

    reports = find_reports(args.root)
    if not reports:
        print(f"ERROR: No assembly_data_report.jsonl found under {args.root}", file=sys.stderr)
        sys.exit(1)

    rows = []
    for rp in reports:
        sp = species_from_path(rp)
        records = parse_jsonl(rp)
        seq_list, asm_list, cmt_list = extract_methods(records)

        seq = "; ".join(seq_list) if seq_list else "NA"
        asm = "; ".join(asm_list) if asm_list else "NA"
        cmt = "; ".join(cmt_list) if cmt_list else "NA"

        rows.append((sp, seq, asm, cmt))

    rows.sort(key=lambda x: x[0].lower())

    with open(args.out, "w", encoding="utf-8") as out:
        out.write("Species\tSequencing_Method\tAssembly_Method\tAssembly_Comments\n")
        for sp, seq, asm, cmt in rows:
            out.write(f"{sp}\t{seq}\t{asm}\t{cmt}\n")

    print(f"Wrote {len(rows)} rows to {args.out}", file=sys.stderr)

if __name__ == "__main__":
    main()
