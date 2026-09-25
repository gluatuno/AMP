#!/usr/bin/env python3
"""Extract promoter sequences for each AMP gene and its processing-enzyme gene,
and record the co-regulated pairs.

Promoter = [TSS - upstream, TSS + downstream] on the gene's strand (minus-strand
promoters are reverse-complemented). Genomic origin is kept in the FASTA header
so downstream TFBS hits can be mapped back to the genome.
"""
import argparse
import csv
import re
import subprocess
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def index_genes(gff):
    genes = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"(?:^|;)Name=([^;]+)", f[8])
            if m:
                genes[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
    return genes


def faidx(genome, chrom, start, end):
    start = max(1, start)
    region = f"{chrom}:{start}-{end}"
    out = subprocess.run(["samtools", "faidx", genome, region],
                         capture_output=True, text=True, check=True).stdout
    return "".join(l for l in out.splitlines() if not l.startswith(">"))


def promoter_window(coord, up, down):
    chrom, gstart, gend, strand = coord
    if strand == "+":
        tss = gstart
        return chrom, tss - up, tss + down, strand
    else:
        tss = gend
        return chrom, tss - down, tss + up, strand


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True)
    ap.add_argument("--enzyme-tsv", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--upstream", type=int, default=1000)
    ap.add_argument("--downstream", type=int, default=200)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-pairs", required=True)
    args = ap.parse_args()

    genes = index_genes(args.gff)

    amp_genes = set()
    with open(args.map) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g = row.get("gene_symbol", ".")
            if g and g != ".":
                amp_genes.add(g)

    pairs = []          # (amp_gene, enzyme_gene, role)
    enzyme_genes = {}   # gene -> role
    with open(args.enzyme_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene, role = row["gene"], row["role"]
            if gene == "." or row.get("status") not in ("ok", None):
                # still record the pairing intent even if extraction failed
                pass
            enzyme_genes[gene] = role
            for fam in (row.get("requested_by_families", "") or "").split(","):
                fam = fam.strip()
                if fam and fam in amp_genes:
                    pairs.append((fam, gene, role))

    # write pairs
    with open(args.out_pairs, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["amp_gene", "partner_gene", "partner_role"])
        for p in sorted(set(pairs)):
            w.writerow(list(p))

    # extract promoters for the union of genes
    wanted = {g: "amp" for g in amp_genes}
    for g, role in enzyme_genes.items():
        if g != ".":
            wanted.setdefault(g, role)

    n = 0
    with open(args.out_fasta, "w") as out:
        for gene, role in sorted(wanted.items()):
            if gene not in genes:
                sys.stderr.write(f"promoter: {gene} not in GFF; skipped\n")
                continue
            chrom, ps, pe, strand = promoter_window(genes[gene],
                                                    args.upstream, args.downstream)
            try:
                seq = faidx(args.genome, chrom, ps, pe)
            except subprocess.CalledProcessError:
                sys.stderr.write(f"promoter: faidx failed for {gene}\n")
                continue
            if strand == "-":
                seq = revcomp(seq)
            out.write(f">{gene}|{role}|{chrom}|{max(1, ps)}|{pe}|{strand}\n{seq}\n")
            n += 1
    sys.stderr.write(f"extract_promoters: {n} promoters, {len(set(pairs))} pairs.\n")


if __name__ == "__main__":
    main()
