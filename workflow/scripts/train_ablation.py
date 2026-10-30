#!/usr/bin/env python3
"""Method 1-5 ablation + ensemble for the context-aware AMP predictor.

Trains a classifier on each CUMULATIVE feature set (M1, M1-2, ... M1-5),
cross-validates ROC-AUC, and reports how much each added gene-context Method
improves over peptide-only (M1). Then trains:
  * a deployable sequence-only model (M1) for the CLI,
  * a full-context model (M1-5) for the genome scanner,
  * a soft-voting ensemble of the cumulative models.

Uses scikit-learn when available; otherwise the pure-Python logistic model in
ml_lite (so the ablation always runs).
"""
import argparse
import csv
import os
import pickle
import sys

from method_features import compute, cumulative_columns, LEVELS


def load(dataset):
    rows, y = [], []
    with open(dataset) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(compute(r))
            y.append(int(r["label"]))
    return rows, y


def matrix(rows, cols):
    return [[row.get(c, 0.0) for c in cols] for row in rows]


def _sklearn():
    try:
        import numpy  # noqa
        import sklearn  # noqa
        return True
    except ImportError:
        return False


def cv_auc(X, y, cv, seed, use_sk):
    if use_sk:
        import numpy as np
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.model_selection import cross_val_score
        k = min(cv, sum(y), len(y) - sum(y))
        if k < 2:
            return float("nan")
        clf = RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                     random_state=seed)
        return float(cross_val_score(clf, np.array(X), np.array(y),
                                     cv=k, scoring="roc_auc").mean())
    from ml_lite import kfold_auc
    return kfold_auc(X, y, k=cv, seed=seed)


def fit_model(X, y, seed, use_sk):
    if use_sk:
        import numpy as np
        from sklearn.ensemble import RandomForestClassifier
        clf = RandomForestClassifier(n_estimators=300, class_weight="balanced",
                                     random_state=seed).fit(np.array(X), np.array(y))
        return {"kind": "sklearn", "clf": clf}
    from ml_lite import LogReg
    return {"kind": "ml_lite", "clf": LogReg(seed=seed).fit(X, y)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--cv", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    rows, y = load(args.dataset)
    if len(set(y)) < 2:
        sys.stderr.write("train_ablation: need both classes.\n")
        sys.exit(2)
    use_sk = _sklearn()
    example = rows[0]

    # --- ablation over cumulative Method sets ---
    results, prev = [], None
    for lv in LEVELS:
        cols = cumulative_columns(lv, example)
        X = matrix(rows, cols)
        a = cv_auc(X, y, args.cv, args.seed, use_sk)
        delta = (a - prev) if (prev is not None and a == a and prev == prev) else float("nan")
        results.append((lv, len(cols), a, delta))
        prev = a

    with open(os.path.join(args.outdir, "ablation.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["method_set", "n_features", "cv_roc_auc", "delta_vs_prev"])
        for lv, nf, a, d in results:
            w.writerow([f"M1..{lv}" if lv != "M1" else "M1", nf,
                        f"{a:.3f}" if a == a else "NA",
                        f"{d:+.3f}" if d == d else "NA"])

    # --- deployable (M1 only) + full (M1-6) models ---
    m1_cols = cumulative_columns("M1", example)
    full_cols = cumulative_columns(LEVELS[-1], example)
    m1_model = fit_model(matrix(rows, m1_cols), y, args.seed, use_sk)
    full_model = fit_model(matrix(rows, full_cols), y, args.seed, use_sk)
    for name, mdl, cols in [("model_M1", m1_model, m1_cols),
                            ("model_full", full_model, full_cols)]:
        with open(os.path.join(args.outdir, name + ".pkl"), "wb") as fh:
            pickle.dump({**mdl, "columns": cols}, fh)

    with open(os.path.join(args.outdir, "metrics.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["metric", "value"])
        w.writerow(["backend", "sklearn_rf" if use_sk else "ml_lite_logreg"])
        w.writerow(["n_positive", sum(y)])
        w.writerow(["n_negative", len(y) - sum(y)])
        w.writerow(["auc_M1_peptide_only", f"{results[0][2]:.3f}"])
        w.writerow([f"auc_full_{results[-1][0]}", f"{results[-1][2]:.3f}"])
        gain = results[-1][2] - results[0][2]
        w.writerow(["total_gain_auc", f"{gain:+.3f}" if gain == gain else "NA"])
    sys.stderr.write(f"train_ablation: M1={results[0][2]:.3f} "
                     f"M1-5={results[-1][2]:.3f} backend="
                     f"{'sklearn' if use_sk else 'ml_lite'}\n")


if __name__ == "__main__":
    main()
