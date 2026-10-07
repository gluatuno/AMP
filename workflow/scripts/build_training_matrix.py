#!/usr/bin/env python3
"""Assemble the labelled feature matrix for the frog AMP genome model.

Positives  = frog AMP loci (mature sequence + genomic features from the mapping
             / annotation / cleavage steps).
Negatives  = non-AMP genomic windows (sample_negatives.py).

Both are passed through frog_features.assemble_row so columns are identical and
comparable, then written as one matrix with a `label` column (1=AMP, 0=not).
"""
import argparse
import csv
import sys

from frog_features import ALL_FEATURES, assemble_row


def load(path, label):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            feats = assemble_row(r)
            rows.append((r.get("id", r.get("peptide_id", "?")),
                         r.get("species", "."), label, feats))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--positives", required=True)
    ap.add_argument("--negatives", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = load(args.positives, 1) + load(args.negatives, 0)
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["id", "species", "label"] + ALL_FEATURES)
        for rid, sp, label, feats in rows:
            w.writerow([rid, sp, label] + [f"{feats[k]:.5g}" for k in ALL_FEATURES])
    npos = sum(1 for _ in rows if _[2] == 1)
    sys.stderr.write(f"build_training_matrix: {len(rows)} rows "
                     f"({npos} positive, {len(rows) - npos} negative).\n")


if __name__ == "__main__":
    main()
