#!/usr/bin/env python3
"""Derive signal / pro-peptide / mature amino-acid coordinates per precursor.

For each peptide->precursor pair:
  * locate the mature peptide inside its precursor (exact, else best window)
  * read the signal-peptide cleavage site from deepsig
  * define the pro-region as everything between the signal peptide and mature

Outputs one row per peptide with 1-based inclusive aa coordinates.
"""
import argparse
import csv
import sys


def read_fasta(path):
    seqs, sid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf)
                sid = line[1:].split()[0]
                buf = []
            elif line:
                buf.append(line)
    if sid is not None:
        seqs[sid] = "".join(buf)
    return seqs


def find_mature(precursor: str, peptide: str):
    """Return 1-based (start, end) of peptide within precursor.

    Exact substring first; otherwise the min-mismatch window (peptides may be
    slightly variant isoforms of the precursor)."""
    if not peptide or not precursor or len(peptide) > len(precursor):
        return None
    idx = precursor.find(peptide)
    if idx >= 0:
        return idx + 1, idx + len(peptide)
    idxu = precursor.upper().find(peptide.upper())
    if idxu >= 0:
        return idxu + 1, idxu + len(peptide)
    # fall back to best sliding window (Hamming distance)
    best_i, best_mm = -1, len(peptide) + 1
    L = len(peptide)
    for i in range(0, len(precursor) - L + 1):
        mm = sum(1 for a, b in zip(precursor[i:i + L], peptide) if a != b)
        if mm < best_mm:
            best_mm, best_i = mm, i
            if mm == 0:
                break
    # accept only if reasonably similar (<=20% mismatch)
    if best_i >= 0 and best_mm <= max(1, int(0.2 * L)):
        return best_i + 1, best_i + L
    return None


def read_signalp(path):
    """deepsig GFF-like output -> {accession: cleavage_position}."""
    sig = {}
    try:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 5:
                    continue
                acc, ftype = f[0], f[2].lower()
                if "signal" in ftype:
                    try:
                        sig[acc] = int(float(f[4]))  # end = cleavage site
                    except ValueError:
                        pass
    except FileNotFoundError:
        sys.stderr.write(f"NOTE: signalp file {path} missing; signal_end=0.\n")
    return sig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--precursors", required=True)
    ap.add_argument("--peptides", required=True)
    ap.add_argument("--map", required=True)
    ap.add_argument("--signalp", required=True)
    ap.add_argument("--dibasic", default="KR,RR,RK,KK")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    precursors = read_fasta(args.precursors)
    peptides = read_fasta(args.peptides)
    sig = read_signalp(args.signalp)
    motifs = [m.strip() for m in args.dibasic.split(",") if m.strip()]

    rows = []
    with open(args.map) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for m in r:
            pid, prec_id = m["peptide_id"], m["precursor_id"]
            prec = precursors.get(prec_id)
            pep = peptides.get(pid)
            if prec is None or pep is None:
                continue
            loc = find_mature(prec, pep)
            if loc is None:
                sys.stderr.write(f"skip {pid}: mature not locatable in {prec_id}\n")
                continue
            mstart, mend = loc
            signal_end = sig.get(prec_id, 0)
            pro_start = signal_end + 1
            pro_end = mstart - 1
            if pro_end < pro_start:
                pro_start = pro_end = 0  # no pro-region

            # nearest upstream dibasic motif (informational)
            upstream = prec[max(0, mstart - 3):mstart - 1]
            dibasic = upstream if upstream in motifs else "."

            rows.append([pid, prec_id, m.get("gene_symbol", "."), len(prec),
                         signal_end, mstart, mend, pro_start, pro_end, dibasic])

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["peptide_id", "precursor_id", "gene_symbol", "precursor_len",
                    "signal_end", "mature_start", "mature_end",
                    "pro_start", "pro_end", "dibasic_hint"])
        w.writerows(rows)
    sys.stderr.write(f"Defined regions for {len(rows)} peptides.\n")


if __name__ == "__main__":
    main()
