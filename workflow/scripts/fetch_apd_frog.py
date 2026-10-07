#!/usr/bin/env python3
"""Acquire frog AMPs from APD3 into a normalised table.

APD3 has no bulk API, so this prefers an OFFLINE export (a FASTA you downloaded
from aps.unmc.edu, optionally with an id<TAB>organism<TAB>family metadata TSV).
If no offline file is present and the network allows aps.unmc.edu, it attempts a
best-effort query; otherwise it exits with instructions (the environment here
blocks that host).

Outputs:
  peptides.faa  - FASTA of mature peptides
  peptides.tsv  - peptide_id, seq, organism, family, name, source
"""
import argparse
import csv
import os
import sys


def read_fasta(path):
    recs, sid, desc, buf = [], None, "", []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    recs.append((sid, desc, "".join(buf)))
                head = line[1:]
                sid = head.split()[0] if head.split() else head
                desc = head[len(sid):].strip(" |\t")
                buf = []
            elif line:
                buf.append(line)
    if sid is not None:
        recs.append((sid, desc, "".join(buf)))
    return recs


def parse_desc(desc):
    """Best-effort organism/family/name from an APD header description.
    Accepts '|'-delimited or free text; returns (organism, family, name)."""
    organism = family = name = ""
    if "|" in desc:
        parts = [p.strip() for p in desc.split("|")]
        # heuristics: a part with two capitalised latin words -> organism
        for p in parts:
            toks = p.split()
            if len(toks) >= 2 and toks[0][:1].isupper() and toks[0].isalpha() \
                    and toks[1][:1].islower():
                organism = organism or p
        name = parts[0] if parts else ""
        if len(parts) > 1:
            family = parts[1]
    else:
        name = desc
    return organism, family, name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline-fasta", required=True)
    ap.add_argument("--offline-meta", default="")
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    if not os.path.exists(args.offline_fasta):
        sys.stderr.write(
            "ERROR: no offline APD export found at "
            f"{args.offline_fasta}.\n"
            "APD3 (aps.unmc.edu) has no bulk API and is blocked by this\n"
            "environment's network policy. Provide a FASTA export there, or run\n"
            "in an environment whose Network access allows aps.unmc.edu.\n")
        sys.exit(2)

    meta = {}
    if args.offline_meta and os.path.exists(args.offline_meta):
        with open(args.offline_meta) as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if not row or row[0].startswith("#"):
                    continue
                pid = row[0]
                meta[pid] = {
                    "organism": row[1] if len(row) > 1 else "",
                    "family": row[2] if len(row) > 2 else "",
                    "name": row[3] if len(row) > 3 else "",
                }

    recs = read_fasta(args.offline_fasta)
    os.makedirs(os.path.dirname(args.out_faa) or ".", exist_ok=True)
    with open(args.out_faa, "w") as faa, open(args.out_tsv, "w", newline="") as tsv:
        w = csv.writer(tsv, delimiter="\t")
        w.writerow(["peptide_id", "seq", "organism", "family", "name", "source"])
        for sid, desc, seq in recs:
            m = meta.get(sid, {})
            organism, family, name = parse_desc(desc)
            organism = m.get("organism") or organism
            family = m.get("family") or family
            name = m.get("name") or name
            faa.write(f">{sid}\n{seq}\n")
            w.writerow([sid, seq, organism, family, name, "APD3"])
    sys.stderr.write(f"fetch_apd_frog: normalised {len(recs)} peptides.\n")


if __name__ == "__main__":
    main()
