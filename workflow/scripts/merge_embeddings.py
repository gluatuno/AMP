#!/usr/bin/env python3
"""Join M6 embedding TSV(s) onto the dataset by id, producing the feature table
the ablation trains on. Rows without an embedding (e.g. negatives lacking a
locus for m6g) get zeros for those columns.
"""
import argparse
import csv
import sys


def read_emb(path):
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t")
        header = next(r)
        cols = header[1:]
        table = {row[0]: dict(zip(cols, row[1:])) for row in r if row}
    return cols, table


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--emb", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    emb_cols, tables = [], []
    for p in args.emb:
        cols, tbl = read_emb(p)
        emb_cols += cols
        tables.append((cols, tbl))

    with open(args.dataset) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
        base_cols = rows[0].keys() if rows else []

    out_cols = list(base_cols) + emb_cols
    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            rid = r.get("id", "")
            for cols, tbl in tables:
                vals = tbl.get(rid, {})
                for c in cols:
                    r[c] = vals.get(c, "0")
            w.writerow(r)
    sys.stderr.write(f"merge_embeddings: added {len(emb_cols)} embedding columns "
                     f"to {len(rows)} rows.\n")


if __name__ == "__main__":
    main()
