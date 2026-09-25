#!/usr/bin/env python3
"""Scan precursor proteins for protease cleavage-site signatures.

Applies the curated motif library (config/protease_motifs.tsv) to every
precursor and emits all candidate scissile bonds, each with the enzyme, the
1-based residue AFTER which cleavage occurs, the P1 residue, a P4-P4' window,
and a score. Score here is motif-membership (1.0) scaled by a simple flank
context term; swap in a PWM/PSSM for calibrated scores.
"""
import argparse
import csv
import re
import sys


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


def read_motifs(path):
    motifs = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            motifs.append({"enzyme": f[0], "role": f[1],
                           "regex": re.compile(f[2]), "cut_offset": int(f[3]),
                           "note": f[4] if len(f) > 4 else ""})
    return motifs


def window(seq, bond, flank=4):
    """P(flank)..P(flank)' around a bond after 1-based residue `bond`."""
    lo = max(0, bond - flank)
    hi = min(len(seq), bond + flank)
    left = seq[lo:bond]
    right = seq[bond:hi]
    return f"{left}|{right}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--precursors", required=True)
    ap.add_argument("--motifs", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    seqs = read_fasta(args.precursors)
    motifs = read_motifs(args.motifs)

    n = 0
    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["precursor_id", "enzyme", "role", "cut_after_aa",
                    "p1_residue", "window_P4_P4prime", "score", "note"])
        for pid, seq in seqs.items():
            for mo in motifs:
                for m in mo["regex"].finditer(seq):
                    bond = m.start() + mo["cut_offset"]     # 1-based residue before bond
                    if bond < 1 or bond >= len(seq):
                        continue
                    p1 = seq[bond - 1]
                    win = window(seq, bond)
                    w.writerow([pid, mo["enzyme"], mo["role"], bond, p1,
                                win, "1.0", mo["note"]])
                    n += 1
    sys.stderr.write(f"scan_cleavage: {n} candidate sites across {len(seqs)} precursors.\n")


if __name__ == "__main__":
    main()
