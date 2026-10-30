#!/usr/bin/env python3
"""Bridge: build the predictor's positive records from genome-pipeline outputs.

Pulls M1-M3 sequence fields from the precursor/regions tables and M4/M5 genomic
context from the cleavage and locus-feature outputs (optional; absent -> 0).
Produces one positive record per mapped AMP in the dataset schema.
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
                sid, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if sid is not None:
        seqs[sid] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", required=True)      # id<TAB>seq (or peptide_id,seq)
    ap.add_argument("--regions", required=True)       # regions_aa.tsv
    ap.add_argument("--precursors", required=True)    # precursors.faa
    ap.add_argument("--cleavage", default="")         # candidate_sites_scored.tsv (M4)
    ap.add_argument("--locus", default="")            # locus features tsv (M5)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # peptides
    pep = {}
    with open(args.peptides) as fh:
        r = csv.DictReader(fh, delimiter="\t")
        key = "peptide_id" if "peptide_id" in (r.fieldnames or []) else "id"
        sal = "seq" if "seq" in (r.fieldnames or []) else None
        for row in r:
            pep[row[key]] = row.get(sal, "") if sal else ""
    precs = read_fasta(args.precursors)

    # M4 cleavage: best convertase score + strong-cut count per precursor
    cut_score, n_cuts = {}, {}
    if args.cleavage:
        with open(args.cleavage) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                p = row["precursor_id"]
                try:
                    sc = float(row.get("score", 0) or 0)
                except ValueError:
                    sc = 0.0
                cut_score[p] = max(cut_score.get(p, 0.0), sc)
                if sc >= 0.5:
                    n_cuts[p] = n_cuts.get(p, 0) + 1

    # M5 locus features keyed by peptide id
    locus = {}
    if args.locus:
        with open(args.locus) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                locus[row.get("id", row.get("peptide_id", ""))] = row

    cols = ["id", "mature_seq", "pro_seq", "precursor_seq", "has_signal",
            "signal_len", "has_dibasic", "enzyme_present", "enzyme_dist_kb",
            "cut_score", "n_strong_cuts", "locus_len_bp", "intron_frac",
            "gc_locus", "n_exons", "transporter_present", "gene_density",
            "neighbor_amp_count"]
    n = 0
    with open(args.regions) as fh, open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in csv.DictReader(fh, delimiter="\t"):
            pid, prec_id = r["peptide_id"], r["precursor_id"]
            prec = precs.get(prec_id, "")
            se = int(r.get("signal_end", 0) or 0)
            ms = int(r.get("mature_start", 0) or 0)
            pro = prec[se:ms - 1] if (prec and ms > se + 1) else ""
            lo = locus.get(pid, {})
            w.writerow({
                "id": pid,
                "mature_seq": pep.get(pid, ""),
                "pro_seq": pro,
                "precursor_seq": prec,
                "has_signal": 1 if se > 0 else 0,
                "signal_len": se,
                "has_dibasic": 1 if r.get("dibasic_hint", ".") not in (".", "") else 0,
                "enzyme_present": 0,            # fill from enzyme-locus table if available
                "enzyme_dist_kb": 0,
                "cut_score": f"{cut_score.get(prec_id, 0.0):.4f}",
                "n_strong_cuts": n_cuts.get(prec_id, 0),
                "locus_len_bp": lo.get("locus_len_bp", 0),
                "intron_frac": lo.get("intron_frac", 0),
                "gc_locus": lo.get("gc_locus", 0),
                "n_exons": lo.get("n_exons", 0),
                "transporter_present": 0,
                "gene_density": 0,
                "neighbor_amp_count": 0,
            })
            n += 1
    sys.stderr.write(f"assemble_positives: wrote {n} positive records.\n")


if __name__ == "__main__":
    main()
