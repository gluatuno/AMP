#!/usr/bin/env python3
"""Feature computation for the frog AMP genome model.

Sequence features are computed identically for positives (known mature AMPs /
their precursors) and negatives (translated ORFs from random genomic windows),
so the two are comparable. Genomic features (exon count, locus length, GC, cut
sites) are passed through from the mapping/annotation steps when available.

Pure Python; reuses the hydropathy scale from cleavage_features.
"""
from cleavage_features import KD

BASIC = set("KR")
ACIDIC = set("DE")
AROMATIC = set("FWY")

# Canonical feature column order (missing values filled with 0.0).
SEQ_FEATURES = [
    "len", "cationicity", "acidic_frac", "net_charge_per_res",
    "hydropathy", "aromatic_frac", "cys_frac", "gly_frac", "boman_like",
]
GENOMIC_FEATURES = [
    "precursor_len", "n_exons", "locus_len_bp", "cds_len_bp",
    "intron_frac", "gc_locus", "signal_len", "pro_len",
    "has_dibasic", "n_cutsites",
]
ALL_FEATURES = SEQ_FEATURES + GENOMIC_FEATURES


def seq_features(seq):
    seq = (seq or "").upper()
    n = len(seq)
    if n == 0:
        return {k: 0.0 for k in SEQ_FEATURES}
    basic = sum(seq.count(x) for x in BASIC)
    acidic = sum(seq.count(x) for x in ACIDIC)
    hyd = sum(KD.get(c, 0.0) for c in seq) / n
    return {
        "len": float(n),
        "cationicity": basic / n,
        "acidic_frac": acidic / n,
        "net_charge_per_res": (basic - acidic) / n,
        "hydropathy": hyd,
        "aromatic_frac": sum(1 for c in seq if c in AROMATIC) / n,
        "cys_frac": seq.count("C") / n,
        "gly_frac": seq.count("G") / n,
        # crude Boman-like index proxy: mean negative hydropathy (binding potential)
        "boman_like": -hyd,
    }


def assemble_row(fields):
    """fields: dict possibly containing 'mature_seq' and any GENOMIC_FEATURES.
    Returns an ordered feature dict over ALL_FEATURES."""
    row = {k: 0.0 for k in ALL_FEATURES}
    sf = seq_features(fields.get("mature_seq", ""))
    row.update(sf)
    for k in GENOMIC_FEATURES:
        if k in fields and fields[k] not in (None, ""):
            try:
                row[k] = float(fields[k])
            except (TypeError, ValueError):
                row[k] = 0.0
    return row
