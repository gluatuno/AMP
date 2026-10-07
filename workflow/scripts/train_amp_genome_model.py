#!/usr/bin/env python3
"""Train the frog AMP genome model from the labelled feature matrix.

With scikit-learn: RandomForest or LogisticRegression, cross-validated AUC and
feature importances, model saved to disk. Without it: a pure-Python baseline
(single-feature AUC ranking) so the step always produces an interpretable
report. The model later scores candidate genomic loci genome-wide.
"""
import argparse
import csv
import pickle
import sys

from frog_features import ALL_FEATURES


def load_matrix(path):
    ids, species, X, y = [], [], [], []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ids.append(r["id"])
            species.append(r.get("species", "."))
            X.append([float(r[k]) for k in ALL_FEATURES])
            y.append(int(r["label"]))
    return ids, species, X, y


def auc(pos, neg):
    if not pos or not neg:
        return float("nan")
    wins = sum((1.0 if p > n else 0.5 if p == n else 0.0)
               for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--classifier", default="random_forest")
    ap.add_argument("--cv", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    import os
    os.makedirs(args.outdir, exist_ok=True)
    ids, species, X, y = load_matrix(args.matrix)
    npos, nneg = sum(y), len(y) - sum(y)

    metrics = [("n_total", len(y)), ("n_positive", npos), ("n_negative", nneg)]

    trained = False
    try:
        import numpy as np
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.model_selection import cross_val_score
        from sklearn.preprocessing import StandardScaler
        from sklearn.pipeline import make_pipeline
        if npos >= 5 and nneg >= 5:
            Xn, yn = np.array(X), np.array(y)
            if args.classifier == "logreg":
                clf = make_pipeline(StandardScaler(),
                                    LogisticRegression(max_iter=1000,
                                                       class_weight="balanced"))
            else:
                clf = RandomForestClassifier(n_estimators=300,
                                             class_weight="balanced",
                                             random_state=args.seed)
            k = min(args.cv, npos, nneg)
            if k >= 2:
                cvauc = cross_val_score(clf, Xn, yn, cv=k, scoring="roc_auc").mean()
                metrics.append(("cv_auc", round(float(cvauc), 3)))
            clf.fit(Xn, yn)
            with open(os.path.join(args.outdir, "model.pkl"), "wb") as fh:
                pickle.dump({"type": args.classifier, "clf": clf,
                             "features": ALL_FEATURES}, fh)
            # feature importances (RF) or |coef| (logreg)
            imp = getattr(clf, "feature_importances_", None)
            if imp is None and hasattr(clf, "named_steps"):
                lr = clf.named_steps.get("logisticregression")
                imp = abs(lr.coef_[0]) if lr is not None else None
            if imp is not None:
                with open(os.path.join(args.outdir, "feature_importance.tsv"),
                          "w", newline="") as fh:
                    w = csv.writer(fh, delimiter="\t")
                    w.writerow(["feature", "importance"])
                    for f, v in sorted(zip(ALL_FEATURES, imp),
                                       key=lambda t: -t[1]):
                        w.writerow([f, f"{v:.4f}"])
            trained = True
            metrics.append(("model_kind", args.classifier))
    except ImportError:
        sys.stderr.write("scikit-learn/numpy absent; pure-Python baseline.\n")

    if not trained:
        # univariate AUC per feature (which single signal separates AMP loci)
        with open(os.path.join(args.outdir, "univariate_auc.tsv"), "w",
                  newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["feature", "auc"])
            ranked = []
            for j, f in enumerate(ALL_FEATURES):
                pos = [X[i][j] for i in range(len(y)) if y[i] == 1]
                neg = [X[i][j] for i in range(len(y)) if y[i] == 0]
                ranked.append((f, auc(pos, neg)))
            for f, a in sorted(ranked, key=lambda t: -(t[1] if t[1] == t[1] else 0)):
                w.writerow([f, f"{a:.3f}" if a == a else "NA"])
        metrics.append(("model_kind", "baseline_univariate_auc"))

    with open(os.path.join(args.outdir, "metrics.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in metrics:
            w.writerow([k, v])
    sys.stderr.write(f"train_amp_genome_model: {dict(metrics)}\n")


if __name__ == "__main__":
    main()
