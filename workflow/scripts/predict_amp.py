#!/usr/bin/env python3
"""Score peptides for AMP likelihood with the deployable (M1) model.

Sequence-only: takes a peptide FASTA and outputs an AMP probability per record,
directly comparable to peptide-only tools. (The full-context M1-5 model is for
the genome scanner, where the gene context is available.)
"""
import argparse
import csv
import pickle
import sys

from method_features import compute


def read_fasta(path):
    recs, sid, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    recs.append((sid, "".join(buf)))
                sid, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if sid is not None:
        recs.append((sid, "".join(buf)))
    return recs


def score_rows(model, cols, X):
    if model["kind"] == "sklearn":
        import numpy as np
        return list(model["clf"].predict_proba(np.array(X))[:, 1])
    return model["clf"].predict_proba(X)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threshold", type=float, default=0.5)
    args = ap.parse_args()

    with open(args.model, "rb") as fh:
        model = pickle.load(fh)
    cols = model["columns"]

    recs = read_fasta(args.fasta)
    X = [[compute({"mature_seq": seq}).get(c, 0.0) for c in cols]
         for _sid, seq in recs]
    probs = score_rows(model, cols, X) if X else []

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["id", "length", "amp_probability", "call"])
        for (sid, seq), p in zip(recs, probs):
            w.writerow([sid, len(seq), f"{p:.4f}",
                        "AMP" if p >= args.threshold else "non-AMP"])
    sys.stderr.write(f"predict_amp: scored {len(recs)} sequences.\n")


if __name__ == "__main__":
    main()
