#!/usr/bin/env python3
"""Assemble per-AMP positive feature rows for the frog model.

Combines the outputs of the per-species mapping/annotation/cleavage steps into
one row per mapped AMP: sequence fields + genomic features (exon count, locus
span, intron fraction, GC, signal/pro lengths, cut-site count). Reuses the
miniprot CDS model via gmap.
"""
import argparse
import csv
import sys

from gmap import parse_gff


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


def gc_fraction(seq):
    seq = seq.upper()
    n = len(seq) or 1
    return (seq.count("G") + seq.count("C")) / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", required=True)     # peptides.tsv (id, seq, ...)
    ap.add_argument("--regions", required=True)      # regions_aa.tsv
    ap.add_argument("--gff", required=True)          # miniprot gff
    ap.add_argument("--cleavage", default="")        # candidate_sites_scored.tsv
    ap.add_argument("--genome", default="")          # genome FASTA (for GC)
    ap.add_argument("--species", default=".")
    ap.add_argument("--cut-score", type=float, default=0.5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    pep = {}
    with open(args.peptides) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            pep[r["peptide_id"]] = r["seq"]

    models = parse_gff(args.gff)
    genome = read_fasta(args.genome) if args.genome else {}

    # cut-site counts per precursor (score >= threshold)
    cuts = {}
    if args.cleavage:
        with open(args.cleavage) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                try:
                    sc = float(r.get("score", 0) or 0)
                except ValueError:
                    sc = 0.0
                if sc >= args.cut_score:
                    cuts[r["precursor_id"]] = cuts.get(r["precursor_id"], 0) + 1

    cols = ["id", "species", "mature_seq", "precursor_len", "n_exons",
            "locus_len_bp", "cds_len_bp", "intron_frac", "gc_locus",
            "signal_len", "pro_len", "has_dibasic", "n_cutsites"]
    n = 0
    with open(args.regions) as fh, open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(cols)
        for r in csv.DictReader(fh, delimiter="\t"):
            pid, prec = r["peptide_id"], r["precursor_id"]
            rec = models.get(prec)
            if rec is None:
                continue
            cds = rec["cds"]
            gstarts = [s for s, _e, _p in cds]
            gends = [e for _s, e, _p in cds]
            locus_len = max(gends) - min(gstarts) + 1
            cds_len = sum(e - s + 1 for s, e, _p in cds)
            intron_frac = 1 - cds_len / locus_len if locus_len else 0.0
            gc = 0.0
            if rec["chrom"] in genome:
                gc = gc_fraction(genome[rec["chrom"]][min(gstarts) - 1:max(gends)])
            signal_len = int(r.get("signal_end", 0) or 0)
            ps, pe = int(r.get("pro_start", 0) or 0), int(r.get("pro_end", 0) or 0)
            pro_len = (pe - ps + 1) if (ps > 0 and pe >= ps) else 0
            has_dibasic = 1 if r.get("dibasic_hint", ".") not in (".", "") else 0
            w.writerow([pid, args.species, pep.get(pid, ""),
                        r.get("precursor_len", 0), len(cds), locus_len, cds_len,
                        f"{intron_frac:.3f}", f"{gc:.3f}", signal_len, pro_len,
                        has_dibasic, cuts.get(prec, 0)])
            n += 1
    sys.stderr.write(f"locus_features: wrote {n} positive feature rows.\n")


if __name__ == "__main__":
    main()
