#!/usr/bin/env python3
"""Method 4: resolve trans-acting processing enzymes (and transporters) for
each AMP precursor from a curated table, then extract their genomic loci.

Human AMPs are single-gene precursors, so the figure's "Enzyme" and
"Transporter" blocks are SEPARATE loci elsewhere in the genome. This script
maps each precursor's gene symbol to enzyme/transporter genes via regex rules,
looks those genes up in the Ensembl GFF3, and extracts their DNA with samtools.
"""
import argparse
import csv
import re
import subprocess
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_table(path):
    rules = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            rules.append({"regex": f[0], "role": f[1], "gene": f[2],
                          "note": f[3] if len(f) > 3 else ""})
    return rules


def index_genes(gff):
    """Ensembl GFF3 'gene' features -> {SYMBOL: (chrom, start, end, strand)}."""
    genes = {}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"(?:^|;)Name=([^;]+)", f[8])
            if not m:
                continue
            genes[m.group(1)] = (f[0], int(f[3]), int(f[4]), f[6])
    return genes


def faidx(genome, chrom, start, end):
    region = f"{chrom}:{start}-{end}"
    out = subprocess.run(["samtools", "faidx", genome, region],
                         capture_output=True, text=True, check=True).stdout
    return "".join(l for l in out.splitlines() if not l.startswith(">"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--map", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--table", required=True)
    ap.add_argument("--out-fasta", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    rules = read_table(args.table)

    # Which enzyme/transporter genes does each precursor family require?
    requested = {}  # (role, gene) -> set(families)
    with open(args.map) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            fam = row.get("gene_symbol", ".") or "."
            title = row.get("stitle", "")
            hay = f"{fam} {title}"
            matched = False
            for rule in rules:
                if rule["regex"] == ".*":
                    continue
                if re.search(rule["regex"], hay, re.IGNORECASE):
                    requested.setdefault((rule["role"], rule["gene"]), set()).add(fam)
                    matched = True
            if not matched:  # fall back to default rules
                for rule in rules:
                    if rule["regex"] == ".*":
                        requested.setdefault((rule["role"], rule["gene"]), set()).add(fam)

    genes = index_genes(args.gff)
    fa_out = open(args.out_fasta, "w")
    tsv_out = csv.writer(open(args.out_tsv, "w", newline=""), delimiter="\t")
    tsv_out.writerow(["role", "gene", "chrom", "start", "end", "strand",
                      "requested_by_families", "status"])

    for (role, gene), fams in sorted(requested.items()):
        fam_str = ",".join(sorted(f for f in fams if f and f != "."))
        if gene == "." or gene not in genes:
            tsv_out.writerow([role, gene, ".", ".", ".", ".", fam_str,
                              "not_found_in_gff"])
            continue
        chrom, start, end, strand = genes[gene]
        try:
            seq = faidx(args.genome, chrom, start, end)
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"samtools faidx failed for {gene}: {e}\n")
            tsv_out.writerow([role, gene, chrom, start, end, strand, fam_str,
                              "extract_failed"])
            continue
        if strand == "-":
            seq = revcomp(seq)
        fa_out.write(f">{gene}|{role}|{chrom}:{start}-{end}|{strand} families={fam_str}\n")
        for i in range(0, len(seq), 60):
            fa_out.write(seq[i:i + 60] + "\n")
        tsv_out.writerow([role, gene, chrom, start, end, strand, fam_str, "ok"])

    fa_out.close()
    sys.stderr.write(f"Method 4: processed {len(requested)} trans-acting loci.\n")


if __name__ == "__main__":
    main()
