#!/usr/bin/env python3
"""Shared amino-acid -> genome coordinate mapping via a miniprot CDS model.

Used by extract_methods.py (feature blocks) and cutsites_to_genome.py
(cleavage-site track). All genomic coordinates are 1-based inclusive on input
(GFF) and returned as 0-based half-open (BED) intervals.
"""
import re


def parse_gff(path):
    """Return {protein_id: best_mrna} where best_mrna =
    {protein, chrom, strand, score, cds:[(start,end,phase)...]}."""
    mrna, cds = {}, {}
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
                mrna[mid] = {"protein": prot, "chrom": chrom,
                             "strand": strand, "score": sc}
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
    """CDS in translation order -> (units, total_coding_len).
    units = list of (cds_off_start, cds_off_end, (strand, base0))."""
    strand = rec["strand"]
    segs = sorted(rec["cds"], key=lambda s: s[0], reverse=(strand == "-"))
    units, off = [], 0
    for i, (gs, ge, ph) in enumerate(segs):
        seg_len = ge - gs + 1
        skip = ph if i == 0 else 0
        usable = seg_len - skip
        if usable <= 0:
            continue
        base0 = (gs + skip) if strand == "+" else (ge - skip)
        units.append((off, off + usable - 1, (strand, base0)))
        off += usable
    return units, off


def map_aa_interval(units, total, a, b):
    """aa [a,b] (1-based inclusive) -> list of genomic (start0, end0) BED intervals."""
    nt_start = (a - 1) * 3
    nt_end = min(b * 3 - 1, total - 1)
    if nt_start > nt_end:
        return []
    out = []
    for (us, ue, (strand, base0)) in units:
        if ue < nt_start or us > nt_end:
            continue
        lo = max(nt_start, us) - us
        hi = min(nt_end, ue) - us
        if strand == "+":
            out.append((base0 + lo - 1, base0 + hi))
        else:
            out.append((base0 - hi - 1, base0 - lo))
    return sorted(out)


def map_bond_point(units, total, k):
    """Genomic 1-bp point for the peptide bond C-terminal to residue k
    (1-based). Returns (start0, end0, strand) or None."""
    ivs = map_aa_interval(units, total, k, k)
    if not ivs:
        return None
    # strand from the first unit
    strand = units[0][2][0]
    if strand == "+":
        e0 = max(e for _s, e in ivs)     # 3' end of codon k
        return (e0 - 1, e0, "+")
    else:
        s0 = min(s for s, _e in ivs)     # 3' end of codon k on minus strand
        return (s0, s0 + 1, "-")
