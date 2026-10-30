#!/usr/bin/env python3
"""Assemble the labelled dataset for the context-aware AMP predictor.

Positives: AMP records with sequence fields (mature/pro/precursor, signal) and
genomic-context fields (M4/M5), as produced by the genome pipeline.
Negatives: non-AMP proteins (--neg-fasta) and/or composition-matched decoys
made by shuffling the positives' mature sequences (--shuffle-per-pos).

Negatives get empty pro/precursor-processing and zeroed genomic context, which
is biologically the point: a real AMP sits in an AMP-like gene context, decoys
do not. (For a harder negative set, pass curated non-secreted proteins that
have been through the same pipeline.)
"""
import argparse
import csv
import random
import sys

FIELDS = ["id", "label", "mature_seq", "pro_seq", "precursor_seq",
          "has_signal", "signal_len", "has_dibasic",
          "enzyme_present", "enzyme_dist_kb", "cut_score", "n_strong_cuts",
          "locus_len_bp", "intron_frac", "gc_locus", "n_exons",
          "transporter_present", "gene_density", "neighbor_amp_count"]


def read_fasta(path):
    recs, sid, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    recs.append((sid, "".join(buf)))
                sid, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if sid is not None:
        recs.append((sid, "".join(buf)))
    return recs


def blank(fields):
    return {k: fields.get(k, "") for k in FIELDS}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--records", required=True, help="positive records TSV")
    ap.add_argument("--neg-fasta", default="", help="non-AMP protein FASTA")
    ap.add_argument("--shuffle-per-pos", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    rows = []
    with open(args.records) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rec = blank(r)
            rec["id"] = r.get("id", r.get("peptide_id", f"pos{len(rows)}"))
            rec["label"] = 1
            rows.append(rec)
    npos = len(rows)

    # composition-matched shuffled decoys
    for i in range(args.shuffle_per_pos):
        for r in rows[:npos]:
            seq = list(r["mature_seq"])
            rng.shuffle(seq)
            s = "".join(seq)
            neg = {k: "" for k in FIELDS}
            neg.update({"id": f"{r['id']}_shuf{i}", "label": 0,
                        "mature_seq": s, "precursor_seq": s})
            rows.append(neg)

    # external non-AMP negatives
    if args.neg_fasta:
        for sid, seq in read_fasta(args.neg_fasta):
            neg = {k: "" for k in FIELDS}
            neg.update({"id": f"neg_{sid}", "label": 0,
                        "mature_seq": seq, "precursor_seq": seq})
            rows.append(neg)

    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in FIELDS})
    nneg = len(rows) - npos
    sys.stderr.write(f"make_dataset: {npos} positives, {nneg} negatives.\n")


if __name__ == "__main__":
    main()
