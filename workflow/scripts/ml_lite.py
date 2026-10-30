#!/usr/bin/env python3
"""Minimal pure-Python ML so the predictor trains/evaluates without sklearn.

Standardised logistic regression (L2, gradient descent) + k-fold ROC-AUC. Used
as a fallback; train_ablation.py prefers scikit-learn when it is installed.
"""
import math
import random


def auc(pos, neg):
    """Mann-Whitney U -> ROC-AUC."""
    if not pos or not neg:
        return float("nan")
    wins = sum((1.0 if p > n else 0.5 if p == n else 0.0)
               for p in pos for n in neg)
    return wins / (len(pos) * len(neg))


def _standardizer(X):
    m = len(X[0])
    means = [sum(row[j] for row in X) / len(X) for j in range(m)]
    stds = []
    for j in range(m):
        var = sum((row[j] - means[j]) ** 2 for row in X) / max(1, len(X) - 1)
        stds.append(math.sqrt(var) or 1.0)
    return means, stds


def _apply(X, means, stds):
    return [[(row[j] - means[j]) / stds[j] for j in range(len(means))] for row in X]


class LogReg:
    def __init__(self, l2=1.0, lr=0.1, epochs=300, seed=0):
        self.l2, self.lr, self.epochs, self.seed = l2, lr, epochs, seed
        self.w, self.b = None, 0.0
        self.means = self.stds = None

    def fit(self, X, y):
        self.means, self.stds = _standardizer(X)
        Xs = _apply(X, self.means, self.stds)
        n, m = len(Xs), len(Xs[0])
        self.w, self.b = [0.0] * m, 0.0
        # class weighting for imbalance
        npos = sum(y) or 1
        nneg = (n - sum(y)) or 1
        wpos, wneg = n / (2 * npos), n / (2 * nneg)
        for _ in range(self.epochs):
            gw, gb = [0.0] * m, 0.0
            for i in range(n):
                z = self.b + sum(self.w[j] * Xs[i][j] for j in range(m))
                p = 1.0 / (1.0 + math.exp(-max(-30, min(30, z))))
                wgt = wpos if y[i] == 1 else wneg
                err = wgt * (p - y[i])
                for j in range(m):
                    gw[j] += err * Xs[i][j]
                gb += err
            for j in range(m):
                gw[j] = gw[j] / n + self.l2 * self.w[j] / n
                self.w[j] -= self.lr * gw[j]
            self.b -= self.lr * gb / n
        return self

    def predict_proba(self, X):
        Xs = _apply(X, self.means, self.stds)
        out = []
        for row in Xs:
            z = self.b + sum(self.w[j] * row[j] for j in range(len(row)))
            out.append(1.0 / (1.0 + math.exp(-max(-30, min(30, z)))))
        return out


def kfold_auc(X, y, k=5, seed=0, **kw):
    idx = list(range(len(X)))
    random.Random(seed).shuffle(idx)
    folds = [idx[i::k] for i in range(k)]
    aucs = []
    for f in folds:
        test = set(f)
        Xtr = [X[i] for i in idx if i not in test]
        ytr = [y[i] for i in idx if i not in test]
        Xte = [X[i] for i in f]
        yte = [y[i] for i in f]
        if len(set(ytr)) < 2 or len(set(yte)) < 2:
            continue
        clf = LogReg(**kw).fit(Xtr, ytr)
        p = clf.predict_proba(Xte)
        pos = [p[i] for i in range(len(yte)) if yte[i] == 1]
        neg = [p[i] for i in range(len(yte)) if yte[i] == 0]
        a = auc(pos, neg)
        if a == a:
            aucs.append(a)
    return sum(aucs) / len(aucs) if aucs else float("nan")
