#!/usr/bin/env python3
"""Route each frog AMP to the genome it should be mapped against.

Rule (per the 'few genomed species' strategy):
  1. If the APD organism names a species/genus that has a genome -> map there.
  2. else if the family/name implies a genus (frog_families.tsv) whose genus
     has a genome -> map there.
  3. else -> status 'no_genome' (set aside; NOT force-mapped to a distant ref).

Emits an assignment table and one mapped-peptide FASTA per target species.
"""
import argparse
import csv
import os
import re
import sys
from collections import Counter


def read_species(path):
    """Return (genus->species_row, species_name->row). has_genome=yes only."""
    by_genus, by_species = {}, {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            row = {"species": f[0], "taxid": f[1], "genus": f[2],
                   "accession": f[3], "assembly": f[4],
                   "has_genome": f[5], "notes": f[6]}
            if row["has_genome"].lower() in ("yes", "partial"):
                by_species[row["species"].lower()] = row
                by_genus.setdefault(row["genus"].lower(), row)
    return by_genus, by_species


def read_families(path):
    fams = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                fams.append((re.compile(f[0], re.IGNORECASE), f[1]))
    return fams


def resolve(organism, family, name, by_genus, by_species, fams):
    org = (organism or "").strip()
    if org:
        if org.lower() in by_species:
            return by_species[org.lower()], "organism_species"
        genus = org.split()[0].lower()
        if genus in by_genus:
            return by_genus[genus], "organism_genus"
    hay = f"{family} {name}"
    for rx, genus in fams:
        if rx.search(hay):
            g = genus.lower()
            if g in by_genus:
                return by_genus[g], "family_genus"
            return None, f"family->{genus}:no_genome"
    return None, "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--peptides", required=True)   # peptides.tsv from fetch
    ap.add_argument("--species-table", required=True)
    ap.add_argument("--family-table", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    by_genus, by_species = read_species(args.species_table)
    fams = read_families(args.family_table)
    os.makedirs(os.path.join(args.outdir, "by_species"), exist_ok=True)

    per_species = {}   # species -> list of (pid, seq)
    status_counts = Counter()
    rows = []
    with open(args.peptides) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sp, how = resolve(r.get("organism", ""), r.get("family", ""),
                              r.get("name", ""), by_genus, by_species, fams)
            if sp is not None:
                status = "mapped"
                per_species.setdefault(sp["species"], []).append(
                    (r["peptide_id"], r["seq"]))
                rows.append([r["peptide_id"], r.get("organism", ""),
                             r.get("family", ""), sp["species"], sp["taxid"],
                             sp["accession"], status, how])
            else:
                status = "no_genome" if how.endswith("no_genome") else how
                rows.append([r["peptide_id"], r.get("organism", ""),
                             r.get("family", ""), ".", ".", ".", "unmapped", how])
            status_counts[status] += 1

    with open(os.path.join(args.outdir, "assignments.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["peptide_id", "organism", "family", "assigned_species",
                    "taxid", "accession", "status", "how"])
        w.writerows(rows)

    for sp, peps in per_species.items():
        safe = sp.replace(" ", "_")
        with open(os.path.join(args.outdir, "by_species", safe + ".faa"), "w") as out:
            for pid, seq in peps:
                out.write(f">{pid}\n{seq}\n")

    with open(os.path.join(args.outdir, "assignment_summary.tsv"), "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["category", "count"])
        for sp, peps in sorted(per_species.items()):
            w.writerow([f"species:{sp}", len(peps)])
        for st, c in sorted(status_counts.items()):
            w.writerow([f"status:{st}", c])
    sys.stderr.write(
        f"assign_species: mapped {status_counts['mapped']} peptides to "
        f"{len(per_species)} species; {sum(status_counts.values())} total.\n")


if __name__ == "__main__":
    main()
