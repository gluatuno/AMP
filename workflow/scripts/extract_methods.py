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
import re
import sys


def parse_gff(path):
    """Return {protein_id: best_mrna} where best_mrna =
    {chrom, strand, score, cds:[(start,end,phase)...]} (1-based inclusive)."""
    mrna = {}   # mrna_id -> record
    cds = {}    # mrna_id -> list
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _, ftype, start, end, score, strand, phase, attrs = f[:9]
            if ftype == "mRNA":
                mid = _attr(attrs, "ID")
                target = _attr(attrs, "Target")
                prot = target.split()[0] if target else None
                try:
                    sc = float(score)
                except ValueError:
                    sc = 0.0
                mrna[mid] = {"protein": prot, "chrom": chrom, "strand": strand,
                             "score": sc}
            elif ftype == "CDS":
                pid = _attr(attrs, "Parent")
                ph = 0 if phase in (".", "") else int(phase)
                cds.setdefault(pid, []).append((int(start), int(end), ph))

    best = {}
    for mid, rec in mrna.items():
        prot = rec["protein"]
        if prot is None or mid not in cds:
            continue
        rec = dict(rec)
        rec["cds"] = cds[mid]
        cur = best.get(prot)
        if cur is None or rec["score"] > cur["score"]:
            best[prot] = rec
    return best


def _attr(attrs, key):
    m = re.search(rf"(?:^|;){re.escape(key)}=([^;]+)", attrs)
    return m.group(1) if m else None


def coding_units(rec):
    """Order CDS in translation order and return coding units:
    list of (cds_off_start, cds_off_end, seg) where seg carries a mapper
    from coding offset -> genomic 1-based coordinate."""
    strand = rec["strand"]
    segs = sorted(rec["cds"], key=lambda s: s[0], reverse=(strand == "-"))
    units = []
    off = 0
    for i, (gs, ge, ph) in enumerate(segs):
        seg_len = ge - gs + 1
        skip = ph if i == 0 else 0     # apply phase only to first coding segment
        usable = seg_len - skip
        if usable <= 0:
            continue
        if strand == "+":
            base0 = gs + skip           # genomic coord of this segment's first coding base
            mapper = ("+", base0)
        else:
            base0 = ge - skip
            mapper = ("-", base0)
        units.append((off, off + usable - 1, mapper))
        off += usable
    return units, off  # off == total coding length


def map_aa_interval(units, total, a, b):
    """aa [a,b] 1-based inclusive -> list of genomic (start0, end0) BED intervals."""
    nt_start = (a - 1) * 3
    nt_end = b * 3 - 1
    nt_end = min(nt_end, total - 1)
    if nt_start > nt_end:
        return []
    out = []
    for (us, ue, (strand, base0)) in units:
        if ue < nt_start or us > nt_end:
            continue
        lo = max(nt_start, us) - us      # offset within this segment
        hi = min(nt_end, ue) - us
        if strand == "+":
            g_lo = base0 + lo            # 1-based
            g_hi = base0 + hi
            out.append((g_lo - 1, g_hi))         # BED 0-based half-open
        else:
            g_hi = base0 - lo            # decreasing genomic with offset
            g_lo = base0 - hi
            out.append((g_lo - 1, g_hi))
    return sorted(out)


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
