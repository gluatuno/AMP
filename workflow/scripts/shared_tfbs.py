#!/usr/bin/env python3
"""Scan AMP and enzyme promoters for TF binding motifs, then report the motifs
SHARED between each co-regulated AMP<->enzyme pair (a candidate common
regulatory signature), plus a genome-mapped TFBS BED track.

Pure-Python IUPAC-consensus scan of both strands. Swap in JASPAR PWMs + FIMO
for calibrated scoring; this gives a deterministic, dependency-free baseline.
"""
import argparse
import csv
import sys
from collections import defaultdict

IUPAC = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"), "B": set("CGT"), "D": set("AGT"),
    "H": set("ACT"), "V": set("ACG"), "N": set("ACGT"),
}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def motif_sets(consensus):
    return [IUPAC.get(c, set("ACGT")) for c in consensus.upper()]


def revcomp_sets(sets):
    comp = [{COMP.get(b, "N") for b in s} for s in sets]
    return list(reversed(comp))


def read_motifs(path):
    motifs = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                motifs.append((f[0], motif_sets(f[1]),
                               f[2] if len(f) > 2 else ""))
    return motifs


def read_promoters(path):
    """Return list of (gene, role, chrom, ps, pe, strand, seq)."""
    out, hdr, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if hdr:
                    out.append(hdr + ["".join(buf)])
                parts = line[1:].split("|")
                hdr = [parts[0], parts[1], parts[2], int(parts[3]),
                       int(parts[4]), parts[5]]
                buf = []
            elif line:
                buf.append(line)
    if hdr:
        out.append(hdr + ["".join(buf)])
    return out


def scan(seq, sets, max_mm):
    """Yield (start, matched) for windows matching within max_mm mismatches."""
    w = len(sets)
    for i in range(0, len(seq) - w + 1):
        mm = 0
        ok = True
        for j in range(w):
            if seq[i + j] not in sets[j]:
                mm += 1
                if mm > max_mm:
                    ok = False
                    break
        if ok:
            yield i, seq[i:i + w]


def to_genome(ps, pe, strand, local, w):
    """Map a local promoter offset -> genomic BED (start0, end0)."""
    if strand == "+":
        g1 = ps + local                # 1-based genomic start
        return g1 - 1, g1 - 1 + w
    else:                              # stored seq is revcomp of genome
        g_end = pe - local            # 1-based
        return g_end - w, g_end


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--promoters", required=True)
    ap.add_argument("--motifs", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--max-mismatch", type=int, default=1)
    ap.add_argument("--out-hits", required=True)
    ap.add_argument("--out-shared", required=True)
    ap.add_argument("--out-bed", required=True)
    args = ap.parse_args()

    motifs = read_motifs(args.motifs)
    proms = read_promoters(args.promoters)

    tf_by_gene = defaultdict(set)     # gene -> set(tf)
    hits_rows, bed_rows = [], []
    for gene, role, chrom, ps, pe, strand, seq in proms:
        seq = seq.upper()
        for tf, sets, _note in motifs:
            for fwd_sets, orient in ((sets, "+"), (revcomp_sets(sets), "-")):
                for local, matched in scan(seq, fwd_sets, args.max_mismatch):
                    tf_by_gene[gene].add(tf)
                    s0, e0 = to_genome(ps, pe, strand, local, len(fwd_sets))
                    hits_rows.append([gene, role, tf, orient, local, chrom, s0, e0, matched])
                    bed_rows.append([chrom, s0, e0, f"{tf}|{gene}|{role}", ".", "."])

    with open(args.out_hits, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["gene", "role", "tf", "motif_orient", "promoter_offset",
                    "chrom", "start0", "end0", "matched_seq"])
        w.writerows(hits_rows)

    with open(args.out_bed, "w") as out:
        for r in bed_rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    # shared signatures per co-regulated pair
    with open(args.pairs) as fh, open(args.out_shared, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["amp_gene", "partner_gene", "shared_TFs", "n_shared",
                    "amp_only_TFs", "partner_only_TFs"])
        for row in csv.DictReader(fh, delimiter="\t"):
            a, p = row["amp_gene"], row["partner_gene"]
            ta, tp = tf_by_gene.get(a, set()), tf_by_gene.get(p, set())
            shared = sorted(ta & tp)
            w.writerow([a, p, ",".join(shared) or ".", len(shared),
                        ",".join(sorted(ta - tp)) or ".",
                        ",".join(sorted(tp - ta)) or "."])
    sys.stderr.write(f"shared_tfbs: {len(hits_rows)} hits across {len(proms)} promoters.\n")


if __name__ == "__main__":
    main()
