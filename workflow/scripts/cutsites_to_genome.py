#!/usr/bin/env python3
"""Project cleavage sites onto the genome as a 1-bp BED track.

Each candidate scissile bond (from scan_cleavage.py) is mapped through the
miniprot CDS model to the genomic base at the 3' end of its P1 codon, yielding
an annotatable "cut-site" track (the signature 'presented in the genome').
"""
import argparse
import csv
import sys

from gmap import parse_gff, coding_units, map_bond_point


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    models = parse_gff(args.gff)
    cache = {}
    n, miss = 0, 0
    with open(args.candidates) as fh, open(args.out, "w") as out:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = row["precursor_id"]
            rec = models.get(pid)
            if rec is None:
                miss += 1
                continue
            if pid not in cache:
                cache[pid] = coding_units(rec)
            units, total = cache[pid]
            bond = int(row["cut_after_aa"])
            pt = map_bond_point(units, total, bond)
            if pt is None:
                miss += 1
                continue
            s0, e0, strand = pt
            name = f"{row['enzyme']}|{pid}|cutAfter{bond}|{row.get('p1_residue','.')}"
            out.write(f"{rec['chrom']}\t{s0}\t{e0}\t{name}\t{row.get('score','.')}\t{strand}\n")
            n += 1
    sys.stderr.write(f"cutsites_to_genome: wrote {n} sites, {miss} unmapped.\n")


if __name__ == "__main__":
    main()
