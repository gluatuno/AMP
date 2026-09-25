#!/usr/bin/env python3
"""Shared featurisation + PWM scoring for the cleavage-site model.

Pure Python (no hard numpy/sklearn dependency) so the same window encoding is
used by train_cleavage_model.py and score_sites.py and can run anywhere.
"""
import math

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}
FLANK = 4
WINDOW_LEN = 2 * FLANK

# Kyte-Doolittle hydropathy (for a coarse physchem feature)
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}
SMALL = set("AGSCT")
AROMATIC = set("FWY")
BASIC = set("KR")


def raw_window(seq, bond, flank=FLANK):
    """P(flank)..P(flank)' window (length 2*flank) around a bond that is
    C-terminal to 1-based residue `bond`. Returns None if it runs off an end."""
    if bond - flank < 0 or bond + flank > len(seq):
        return None
    return seq[bond - flank: bond + flank]


def featurize(window):
    """One-hot (20 x WINDOW_LEN) + a few physicochemical features -> list[float]."""
    feats = []
    for ch in window:
        row = [0.0] * len(AA)
        if ch in AA_INDEX:
            row[AA_INDEX[ch]] = 1.0
        feats.extend(row)
    p1 = window[FLANK - 1] if len(window) >= FLANK else "X"       # P1 residue
    p1p = window[FLANK] if len(window) > FLANK else "X"           # P1' residue
    cation = sum(1 for c in window if c in BASIC) / len(window)
    hydro = sum(KD.get(c, 0.0) for c in window) / len(window)
    feats.extend([
        cation,
        hydro,
        1.0 if p1 in BASIC else 0.0,
        1.0 if p1 in SMALL else 0.0,
        1.0 if p1p in AROMATIC else 0.0,
    ])
    return feats


def compute_bg(seqs):
    counts = {a: 1 for a in AA}   # Laplace pseudocount
    for s in seqs:
        for c in s:
            if c in AA_INDEX:
                counts[c] += 1
    tot = sum(counts.values())
    return {a: counts[a] / tot for a in AA}


def build_pwm(pos_windows, bg):
    """Log-odds PWM (list over positions of {aa: log2(p/bg)})."""
    L = WINDOW_LEN
    pwm = []
    for i in range(L):
        counts = {a: 0.5 for a in AA}   # pseudocount
        n = 0
        for w in pos_windows:
            if len(w) == L and w[i] in AA_INDEX:
                counts[w[i]] += 1
                n += 1
        denom = (n + 0.5 * len(AA))
        pwm.append({a: math.log2((counts[a] / denom) / bg[a]) for a in AA})
    return pwm


def score_pwm(window, pwm):
    if len(window) != len(pwm):
        return None
    return sum(pwm[i].get(window[i], 0.0) for i in range(len(window)))


def logistic(x):
    try:
        return 1.0 / (1.0 + math.exp(-x))
    except OverflowError:
        return 0.0 if x < 0 else 1.0
