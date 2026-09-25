#!/usr/bin/env python3
"""Re-score candidate cleavage sites with the trained model / PWM.

Replaces scan_cleavage.py's placeholder score with a calibrated one:
  * logistic-regression probability when model.pkl is present (sklearn), else
  * PWM log-odds squashed to (0,1) via a logistic.
Sites whose full P4-P4' window runs off a terminus get score 0.
"""
import argparse
import csv
import os
import pickle
import sys

from cleavage_features import raw_window, featurize, score_pwm, logistic


def read_fasta(path):
    seqs, sid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf)
                sid, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if sid is not None:
        seqs[sid] = "".join(buf)
    return seqs


def load_pwm(path):
    if not os.path.exists(path):
        return None
    pwm = []
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        aas = header[1:]
        for row in r:
            pwm.append({a: float(v) for a, v in zip(aas, row[1:])})
    return pwm


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--precursors", required=True)
    ap.add_argument("--model", help="model.pkl (optional)")
    ap.add_argument("--pwm", help="pwm.tsv (fallback)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = read_fasta(args.precursors)
    pwm = load_pwm(args.pwm) if args.pwm else None

    clf = None
    if args.model and os.path.exists(args.model):
        try:
            import numpy as np  # noqa: F401
            with open(args.model, "rb") as fh:
                clf = pickle.load(fh).get("clf")
        except Exception as e:  # sklearn/numpy missing or unpickle fail
            sys.stderr.write(f"score_sites: model unusable ({e}); using PWM.\n")

    def score(window):
        if window is None:
            return 0.0
        if clf is not None:
            try:
                import numpy as np
                return float(clf.predict_proba(np.array([featurize(window)]))[0, 1])
            except Exception:
                pass
        if pwm is not None:
            s = score_pwm(window, pwm)
            return logistic(s) if s is not None else 0.0
        return 0.0

    with open(args.candidates) as fh, open(args.out, "w", newline="") as out:
        r = csv.DictReader(fh, delimiter="\t")
        cols = r.fieldnames
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        n = 0
        for row in r:
            seq = seqs.get(row["precursor_id"], "")
            win = raw_window(seq, int(row["cut_after_aa"])) if seq else None
            row["score"] = f"{score(win):.4f}"
            w.writerow(row)
            n += 1
    sys.stderr.write(f"score_sites: rescored {n} sites "
                     f"(model={'logreg' if clf is not None else 'pwm'}).\n")


if __name__ == "__main__":
    main()
