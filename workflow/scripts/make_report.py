#!/usr/bin/env python3
"""Summarise the run: how many peptides mapped, and a per-peptide feature map
listing which methods produced sequence for each peptide."""
import argparse
import csv
import os


def count_fasta(path):
    if not os.path.exists(path):
        return 0
    with open(path) as fh:
        return sum(1 for l in fh if l.startswith(">"))


def peptides_in_fasta(path):
    ids = set()
    if not os.path.exists(path):
        return ids
    with open(path) as fh:
        for l in fh:
            if l.startswith(">"):
                ids.add(l[1:].split("|")[0].split()[0])
    return ids


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--methods-dir", required=True)
    ap.add_argument("--out-summary", required=True)
    ap.add_argument("--out-featmap", required=True)
    args = ap.parse_args()

    md = args.methods_dir
    files = {
        "method1_amp": f"{md}/method1_amp.fa",
        "method2_propeptide": f"{md}/method2_propeptide.fa",
        "method3_pro_plus_amp": f"{md}/method3_pro_plus_amp.fa",
        "method5_locus": f"{md}/method5_locus.fa",
    }
    present = {k: peptides_in_fasta(v) for k, v in files.items()}

    # Summary counts
    with open(args.out_summary, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["stage", "count"])
        w.writerow(["peptides_mapped_to_precursor", sum(1 for _ in open(args.map)) - 1])
        w.writerow(["regions_defined", sum(1 for _ in open(args.regions)) - 1])
        for k, v in files.items():
            w.writerow([f"{k}_records", count_fasta(v)])
        w.writerow(["method4_enzyme_records", count_fasta(f"{md}/method4_enzyme.fa")])

    # Per-peptide feature map
    with open(args.regions) as fh, open(args.out_featmap, "w", newline="") as out:
        r = csv.DictReader(fh, delimiter="\t")
        w = csv.writer(out, delimiter="\t")
        w.writerow(["peptide_id", "precursor_id", "gene_symbol", "precursor_len",
                    "signal_end", "mature_start", "mature_end", "pro_start",
                    "pro_end", "has_m1", "has_m2", "has_m3", "has_m5"])
        for row in r:
            pep = row["peptide_id"]
            w.writerow([row["peptide_id"], row["precursor_id"], row["gene_symbol"],
                        row["precursor_len"], row["signal_end"],
                        row["mature_start"], row["mature_end"],
                        row["pro_start"], row["pro_end"],
                        int(pep in present["method1_amp"]),
                        int(pep in present["method2_propeptide"]),
                        int(pep in present["method3_pro_plus_amp"]),
                        int(pep in present["method5_locus"])])


if __name__ == "__main__":
    main()
