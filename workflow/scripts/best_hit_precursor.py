#!/usr/bin/env python3
"""Pick the best human precursor for each mature AMP peptide.

Reads DIAMOND blastp outfmt-6 with the trailing `stitle` and `full_sseq`
columns, applies identity + coverage filters, and writes:
  * a peptide -> precursor mapping table
  * a FASTA of the (deduplicated) precursor proteins

Coverage is measured on the *query* (the mature peptide): a real precursor
should contain almost the whole mature peptide.
"""
import argparse
import csv
import sys

COLS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        "stitle", "full_sseq"]


def parse_gene_symbol(stitle: str) -> str:
    """Extract GN= gene symbol from a UniProt FASTA title if present."""
    for tok in stitle.split():
        if tok.startswith("GN="):
            return tok[3:]
    return "."


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blast", required=True)
    ap.add_argument("--min-pident", type=float, default=90.0)
    ap.add_argument("--min-cov", type=float, default=0.80)
    ap.add_argument("--out-map", required=True)
    ap.add_argument("--out-fasta", required=True)
    args = ap.parse_args()

    best = {}  # qseqid -> row dict (highest bitscore passing filters)
    with open(args.blast) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(COLS):
                # tolerate stitle containing tabs by re-joining the tail
                parts = parts[:len(COLS) - 1] + ["\t".join(parts[len(COLS) - 1:])]
            row = dict(zip(COLS, parts))
            pident = float(row["pident"])
            qcov = int(row["length"]) / max(int(row["qlen"]), 1)
            if pident < args.min_pident or qcov < args.min_cov:
                continue
            bs = float(row["bitscore"])
            cur = best.get(row["qseqid"])
            if cur is None or bs > float(cur["bitscore"]):
                row["_qcov"] = f"{qcov:.3f}"
                best[row["qseqid"]] = row

    if not best:
        sys.stderr.write("WARNING: no peptide passed the precursor filters. "
                         "Loosen --min-pident / --min-cov in config.\n")

    # Write mapping table
    with open(args.out_map, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["peptide_id", "precursor_id", "gene_symbol", "pident",
                    "query_cov", "precursor_len", "evalue", "stitle"])
        for q, r in sorted(best.items()):
            w.writerow([q, r["sseqid"], parse_gene_symbol(r["stitle"]),
                        r["pident"], r["_qcov"], r["slen"], r["evalue"],
                        r["stitle"]])

    # Write unique precursor FASTA (one record per precursor_id).
    seen = {}
    for r in best.values():
        seen.setdefault(r["sseqid"], r["full_sseq"].replace("-", ""))
    with open(args.out_fasta, "w") as out:
        for sid, seq in seen.items():
            out.write(f">{sid}\n")
            for i in range(0, len(seq), 60):
                out.write(seq[i:i + 60] + "\n")

    sys.stderr.write(f"Mapped {len(best)} peptides to {len(seen)} precursors.\n")


if __name__ == "__main__":
    main()
