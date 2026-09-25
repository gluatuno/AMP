#!/usr/bin/env python3
"""Convert per-precursor amino-acid regions into genomic BED intervals.

Reads the miniprot GFF (spliced protein->genome alignment) and the aa region
table, then writes BED files for:
  Method 1  mature (AMP)              -- exonic
  Method 2  pro-peptide               -- exonic
  Method 3  pro-peptide + AMP         -- exonic
  Method 5  locus span (whole mRNA)   -- single interval, introns included

Amino-acid position p (1-based) occupies coding nucleotides [(p-1)*3, p*3).
We walk the CDS segments in translation order (strand-aware, first-segment
phase applied) to translate those coding offsets back to genomic coordinates.
"""
import argparse
import csv
import os
import sys

from gmap import parse_gff, coding_units, map_aa_interval


def write_bed(path, records):
    with open(path, "w") as out:
        for (chrom, s0, e0, name, strand) in records:
            out.write(f"{chrom}\t{s0}\t{e0}\t{name}\t.\t{strand}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    models = parse_gff(args.gff)
    os.makedirs(args.outdir, exist_ok=True)
    m1, m2, m3, m5 = [], [], [], []

    with open(args.regions) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            prec = row["precursor_id"]
            rec = models.get(prec)
            if rec is None:
                sys.stderr.write(f"no genome model for {prec}; skipping {row['peptide_id']}\n")
                continue
            units, total = coding_units(rec)
            chrom, strand = rec["chrom"], rec["strand"]
            pep = row["peptide_id"]
            ms, me = int(row["mature_start"]), int(row["mature_end"])
            ps, pe = int(row["pro_start"]), int(row["pro_end"])

            def emit(dst, a, b, tag):
                for i, (s0, e0) in enumerate(map_aa_interval(units, total, a, b)):
                    dst.append((chrom, s0, e0,
                                f"{pep}|{tag}|{prec}|ex{i+1}", strand))

            emit(m1, ms, me, "AMP")                       # Method 1
            if ps > 0 and pe >= ps:
                emit(m2, ps, pe, "PRO")                   # Method 2
                emit(m3, ps, me, "PRO_AMP")               # Method 3 (pro+mature)
            else:
                emit(m3, ms, me, "PRO_AMP")               # no pro-region -> mature

            # Method 5: whole-locus genomic envelope (introns included)
            gstarts = [s for (s, _e, _p) in rec["cds"]]
            gends = [e for (_s, e, _p) in rec["cds"]]
            m5.append((chrom, min(gstarts) - 1, max(gends),
                       f"{pep}|LOCUS|{prec}", strand))

    write_bed(os.path.join(args.outdir, "method1_amp.bed"), m1)
    write_bed(os.path.join(args.outdir, "method2_propeptide.bed"), m2)
    write_bed(os.path.join(args.outdir, "method3_pro_plus_amp.bed"), m3)
    write_bed(os.path.join(args.outdir, "method5_locus.bed"), m5)
    sys.stderr.write(
        f"BED written: m1={len(m1)} m2={len(m2)} m3={len(m3)} m5={len(m5)}\n")


if __name__ == "__main__":
    main()
