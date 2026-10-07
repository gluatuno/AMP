#!/usr/bin/env python3
"""Sample non-AMP genomic windows as negatives for the frog AMP model.

Baseline strategy: draw random windows from the genome (length matched to the
positive AMP loci), excluding windows that overlap known AMP loci, translate
the window, and emit the same feature schema as the positives. Swap in
'non_amp_gene' sampling (random annotated genes that are not AMPs) for a harder
negative set.

Reads the genome as plain FASTA (no samtools dependency) so it is testable on a
small genome; for GRCh-scale frog genomes prefer streaming / faidx.
"""
import argparse
import csv
import random
import sys

STOP = {"TAA", "TAG", "TGA"}
CODON = {
    # minimal standard table (enough for ORF translation of negatives)
}
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
for i, aa in enumerate(_AAS):
    CODON["".join(_BASES[(i >> s) & 3] for s in (4, 2, 0))] = aa


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


def read_bed(path):
    byc = {}
    if not path:
        return byc
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            f = line.split("\t")
            if len(f) >= 3:
                byc.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return byc


def overlaps(chrom, s, e, avoid):
    for (a, b) in avoid.get(chrom, []):
        if s < b and a < e:
            return True
    return False


def translate(dna):
    dna = dna.upper()
    best = ""
    for frame in range(3):
        prot, cur = [], []
        for i in range(frame, len(dna) - 2, 3):
            aa = CODON.get(dna[i:i + 3], "X")
            if aa == "*":
                if len("".join(cur)) > len(best):
                    best = "".join(cur)
                cur = []
            else:
                cur.append(aa)
        if len("".join(cur)) > len(best):
            best = "".join(cur)
    return best


def gc(dna):
    dna = dna.upper()
    n = len(dna) or 1
    return (dna.count("G") + dna.count("C")) / n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("--avoid-bed", default="")
    ap.add_argument("--n", type=int, default=100)
    ap.add_argument("--win-min", type=int, default=300)
    ap.add_argument("--win-max", type=int, default=3000)
    ap.add_argument("--species", default=".")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    genome = read_fasta(args.genome)
    avoid = read_bed(args.avoid_bed)
    chroms = [c for c in genome if len(genome[c]) > args.win_max]
    if not chroms:
        chroms = list(genome)
    rng = random.Random(args.seed)

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["id", "species", "mature_seq", "precursor_len", "n_exons",
                    "locus_len_bp", "cds_len_bp", "intron_frac", "gc_locus",
                    "signal_len", "pro_len", "has_dibasic", "n_cutsites"])
        made, tries = 0, 0
        while made < args.n and tries < args.n * 50:
            tries += 1
            chrom = rng.choice(chroms)
            L = len(genome[chrom])
            wlen = rng.randint(args.win_min, min(args.win_max, L))
            s = rng.randint(0, max(0, L - wlen))
            e = s + wlen
            if overlaps(chrom, s, e, avoid):
                continue
            dna = genome[chrom][s:e]
            prot = translate(dna)
            has_dibasic = 1 if any(prot[i:i + 2] in ("KR", "RR", "KK", "RK")
                                   for i in range(len(prot) - 1)) else 0
            w.writerow([f"neg_{chrom}_{s}", args.species, prot, len(prot), 1,
                        wlen, wlen, 0.0, f"{gc(dna):.3f}", 0, 0, has_dibasic, 0])
            made += 1
    sys.stderr.write(f"sample_negatives: wrote {made} negatives.\n")


if __name__ == "__main__":
    main()
