#!/usr/bin/env python3
"""Train a cleavage-site scorer from the known junctions.

Positives  = P4-P4' windows at true signal->pro and pro->mature bonds.
Negatives  = windows at random bonds away from any true boundary.

Always builds a log-odds PWM (pure Python). If scikit-learn is available it
also trains a logistic-regression classifier on the featurised windows and
reports cross-validated AUC; otherwise it reports the PWM's own AUC. Artefacts
(model.pkl, pwm.tsv, metrics.tsv) are consumed by score_sites.py.
"""
import argparse
import csv
import os
import pickle
import random
import sys

from cleavage_features import (AA, WINDOW_LEN, FLANK, raw_window, featurize,
                               compute_bg, build_pwm, score_pwm)


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


def auc(scores_pos, scores_neg):
    """Mann-Whitney U -> AUC (pure Python)."""
    if not scores_pos or not scores_neg:
        return float("nan")
    wins = 0.0
    for p in scores_pos:
        for n in scores_neg:
            wins += 1.0 if p > n else (0.5 if p == n else 0.0)
    return wins / (len(scores_pos) * len(scores_neg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--precursors", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--neg-per-pos", type=int, default=3)
    ap.add_argument("--tol", type=int, default=2)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    seqs = read_fasta(args.precursors)
    truth = {}   # precursor -> list of true bonds
    with open(args.regions) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = row["precursor_id"]
            bonds = []
            if int(row["mature_start"]) > 1:
                bonds.append(int(row["mature_start"]) - 1)
            if int(row["signal_end"]) > 0:
                bonds.append(int(row["signal_end"]))
            truth.setdefault(pid, []).extend(bonds)

    pos_windows, neg_windows = [], []
    for pid, bonds in truth.items():
        seq = seqs.get(pid)
        if not seq:
            continue
        for b in bonds:
            w = raw_window(seq, b)
            if w:
                pos_windows.append(w)
        # negatives: random bonds not near any true bond
        avoid = set()
        for b in bonds:
            avoid.update(range(b - args.tol, b + args.tol + 1))
        choices = [k for k in range(FLANK, len(seq) - FLANK) if k not in avoid]
        rng.shuffle(choices)
        for k in choices[:args.neg_per_pos * max(1, len(bonds))]:
            w = raw_window(seq, k)
            if w:
                neg_windows.append(w)

    if not pos_windows:
        sys.stderr.write("train_cleavage_model: no positive windows; writing empty artefacts.\n")

    # PWM (always)
    bg = compute_bg(seqs.values())
    pwm = build_pwm(pos_windows, bg)
    with open(os.path.join(args.outdir, "pwm.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["position"] + list(AA))
        for i, col in enumerate(pwm):
            w.writerow([i + 1] + [f"{col[a]:.3f}" for a in AA])

    pwm_pos = [s for s in (score_pwm(w, pwm) for w in pos_windows) if s is not None]
    pwm_neg = [s for s in (score_pwm(w, pwm) for w in neg_windows) if s is not None]
    pwm_auc = auc(pwm_pos, pwm_neg)

    # Optional ML classifier
    model_kind, ml_auc = "none", float("nan")
    try:
        import numpy as np
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import cross_val_score
        X = [featurize(w) for w in pos_windows + neg_windows]
        y = [1] * len(pos_windows) + [0] * len(neg_windows)
        if len(set(y)) == 2 and len(y) >= 10:
            Xn, yn = np.array(X), np.array(y)
            clf = LogisticRegression(max_iter=1000, class_weight="balanced")
            k = min(5, sum(yn == 1), sum(yn == 0))
            if k >= 2:
                ml_auc = float(cross_val_score(clf, Xn, yn, cv=k,
                                                scoring="roc_auc").mean())
            clf.fit(Xn, yn)
            with open(os.path.join(args.outdir, "model.pkl"), "wb") as fh:
                pickle.dump({"type": "logreg", "clf": clf}, fh)
            model_kind = "logreg"
    except ImportError:
        sys.stderr.write("scikit-learn/numpy not available; PWM-only scorer.\n")

    with open(os.path.join(args.outdir, "metrics.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["n_positives", len(pos_windows)])
        w.writerow(["n_negatives", len(neg_windows)])
        w.writerow(["pwm_auc", f"{pwm_auc:.3f}"])
        w.writerow(["model_kind", model_kind])
        w.writerow(["ml_cv_auc", f"{ml_auc:.3f}" if ml_auc == ml_auc else "NA"])
    sys.stderr.write(
        f"train_cleavage_model: pos={len(pos_windows)} neg={len(neg_windows)} "
        f"pwm_auc={pwm_auc:.3f} model={model_kind} ml_auc={ml_auc}\n")


if __name__ == "__main__":
    main()
