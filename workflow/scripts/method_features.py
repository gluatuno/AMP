#!/usr/bin/env python3
"""Per-Method feature blocks for the context-aware AMP predictor.

The reference figure defines five nested views of an AMP's gene context:
  M1 mature peptide            (what peptide-only tools use)
  M2 + pro-peptide             (processing / anionic pro-region)
  M3 + full precursor          (signal peptide + pro + mature)
  M4 + processing-enzyme       (is a convertase/protease nearby?)
  M5 + whole locus             (transporter...AMP neighbourhood)

M1-M3 are computable from sequence alone (usable at predict time with no
genome); M4-M5 are genomic-context fields passed through from the pipeline
(used in the ablation and by the genome scanner).

Pure Python; reuses the Kyte-Doolittle scale from cleavage_features.
"""
from cleavage_features import KD

AA = "ACDEFGHIKLMNPQRSTVWY"
BASIC, ACIDIC, AROMATIC = set("KR"), set("DE"), set("FWY")


def _frac(seq, chars):
    return sum(seq.count(c) for c in chars) / len(seq) if seq else 0.0


def _seqstats(seq, prefix):
    seq = (seq or "").upper()
    n = len(seq)
    if n == 0:
        base = {f"{prefix}_len": 0.0, f"{prefix}_charge": 0.0,
                f"{prefix}_cationic": 0.0, f"{prefix}_acidic": 0.0,
                f"{prefix}_hydropathy": 0.0, f"{prefix}_aromatic": 0.0,
                f"{prefix}_cys": 0.0, f"{prefix}_gly": 0.0, f"{prefix}_boman": 0.0}
        return base
    hyd = sum(KD.get(c, 0.0) for c in seq) / n
    return {
        f"{prefix}_len": float(n),
        f"{prefix}_charge": (_frac(seq, BASIC) - _frac(seq, ACIDIC)),
        f"{prefix}_cationic": _frac(seq, BASIC),
        f"{prefix}_acidic": _frac(seq, ACIDIC),
        f"{prefix}_hydropathy": hyd,
        f"{prefix}_aromatic": _frac(seq, AROMATIC),
        f"{prefix}_cys": _frac(seq, "C"),
        f"{prefix}_gly": _frac(seq, "G"),
        f"{prefix}_boman": -hyd,
    }


def _aa_composition(seq, prefix):
    seq = (seq or "").upper()
    n = len(seq) or 1
    return {f"{prefix}_aa_{a}": seq.count(a) / n for a in AA}


def _passthrough(fields, keys):
    out = {}
    for k in keys:
        v = fields.get(k, "")
        try:
            out[k] = float(v) if v not in ("", None) else 0.0
        except (TypeError, ValueError):
            out[k] = 0.0
    return out


# Genomic-context fields supplied by the pipeline (absent -> 0).
M4_KEYS = ["enzyme_present", "enzyme_dist_kb", "cut_score", "n_strong_cuts"]
M5_KEYS = ["locus_len_bp", "intron_frac", "gc_locus", "n_exons",
           "transporter_present", "gene_density", "neighbor_amp_count"]


def compute(fields):
    """fields: dict with mature_seq, pro_seq, precursor_seq, signal_len,
    has_signal, and any M4/M5 genomic keys. Returns the full feature dict."""
    feats = {}
    # M1 mature
    feats.update(_seqstats(fields.get("mature_seq", ""), "m1"))
    feats.update(_aa_composition(fields.get("mature_seq", ""), "m1"))
    # M2 pro-peptide
    feats.update(_seqstats(fields.get("pro_seq", ""), "m2pro"))
    feats["m2_has_dibasic"] = float(fields.get("has_dibasic", 0) or 0)
    # M3 precursor (signal + pro + mature)
    feats.update(_seqstats(fields.get("precursor_seq", ""), "m3prec"))
    feats["m3_has_signal"] = float(fields.get("has_signal", 0) or 0)
    feats["m3_signal_len"] = float(fields.get("signal_len", 0) or 0)
    # M4 enzyme context (genomic passthrough)
    feats.update(_passthrough(fields, M4_KEYS))
    # M5 locus context (genomic passthrough)
    feats.update(_passthrough(fields, M5_KEYS))
    return feats


def _cols_for(prefixes, example):
    return [k for k in example if any(k.startswith(p) for p in prefixes)]


# Which feature columns belong to each Method block (prefix-based).
BLOCK_PREFIXES = {
    "M1": ["m1_", "m1_aa_"],
    "M2": ["m2pro_", "m2_"],
    "M3": ["m3prec_", "m3_"],
    "M4": M4_KEYS,
    "M5": M5_KEYS,
}
LEVELS = ["M1", "M2", "M3", "M4", "M5"]


def cumulative_columns(level, example):
    """Columns for the cumulative set M1..level (e.g. level='M3' -> M1+M2+M3)."""
    idx = LEVELS.index(level)
    prefixes = []
    for lv in LEVELS[:idx + 1]:
        prefixes += BLOCK_PREFIXES[lv]
    cols, seen = [], set()
    for k in example:
        if k in seen:
            continue
        if k in prefixes or any(k.startswith(p) for p in prefixes if p.endswith("_")):
            cols.append(k)
            seen.add(k)
    return cols
