#!/usr/bin/env python3
"""De-novo cleavage-motif discovery from the pipeline's own precursors.

Using the known boundaries in regions_aa.tsv, extract P4-P4' windows around
each signal->pro and pro->mature junction, then:
  * always compute a position-frequency matrix (PFM) + consensus (pure Python),
  * additionally run STREME/MEME (if installed) for a statistical motif.

Because the true boundaries are known, the discovered motif can be compared to
the curated specificity in protease_motifs.tsv (validation / novelty).
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import Counter

FLANK = 4
POS_LABELS = [f"P{i}" for i in range(FLANK, 0, -1)] + [f"P{i}'" for i in range(1, FLANK + 1)]


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


def window_at(seq, bond):
    """Full-length P4..P4' window around a bond after 1-based residue `bond`,
    or None if it runs off either end."""
    if bond - FLANK < 0 or bond + FLANK > len(seq):
        return None
    return seq[bond - FLANK: bond + FLANK]


def collect_windows(precursors, regions):
    pro, sig = [], []
    with open(regions) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            seq = precursors.get(row["precursor_id"])
            if not seq:
                continue
            ms = int(row["mature_start"])
            se = int(row["signal_end"])
            w = window_at(seq, ms - 1)          # pro -> mature bond
            if w:
                pro.append(w)
            if se > 0:
                w2 = window_at(seq, se)         # signal -> pro bond
                if w2:
                    sig.append(w2)
    return pro, sig


def pfm(windows):
    if not windows:
        return []
    L = len(windows[0])
    cols = []
    for i in range(L):
        c = Counter(w[i] for w in windows if len(w) == L)
        tot = sum(c.values()) or 1
        top = c.most_common(3)
        cols.append([(aa, f"{n / tot:.2f}") for aa, n in top])
    return cols


def write_fasta(path, windows, tag):
    with open(path, "w") as out:
        for i, w in enumerate(windows):
            out.write(f">{tag}_{i}\n{w}\n")


def run_streme(fa, outdir):
    tool = shutil.which("streme") or shutil.which("meme")
    if not tool:
        return "not_installed"
    try:
        if tool.endswith("streme"):
            subprocess.run([tool, "--p", fa, "--protein", "--nmotifs", "3",
                            "--oc", outdir], check=True,
                           capture_output=True, text=True)
        else:  # meme
            subprocess.run([tool, fa, "-protein", "-nmotifs", "3", "-mod", "zoops",
                            "-oc", outdir], check=True,
                           capture_output=True, text=True)
        return "ok"
    except subprocess.CalledProcessError as e:
        sys.stderr.write(e.stderr or str(e))
        return "failed"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--precursors", required=True)
    ap.add_argument("--regions", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    precursors = read_fasta(args.precursors)
    pro, sig = collect_windows(precursors, args.regions)

    pro_fa = os.path.join(args.outdir, "promature_windows.fa")
    sig_fa = os.path.join(args.outdir, "signalpro_windows.fa")
    write_fasta(pro_fa, pro, "promature")
    write_fasta(sig_fa, sig, "signalpro")

    # PFM report (pure Python; always produced)
    with open(os.path.join(args.outdir, "motif_report.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["junction", "position", "top_residues(freq)"])
        for junction, wins in (("pro_mature", pro), ("signal_pro", sig)):
            for lab, col in zip(POS_LABELS, pfm(wins)):
                w.writerow([junction, lab, ";".join(f"{aa}:{fr}" for aa, fr in col)])

    # Statistical motif (optional; STREME/MEME if available)
    status = run_streme(pro_fa, os.path.join(args.outdir, "streme_promature"))
    sys.stderr.write(
        f"discover_cutsites: pro->mature n={len(pro)}, signal->pro n={len(sig)}; "
        f"STREME/MEME={status}\n")


if __name__ == "__main__":
    main()
