#!/usr/bin/env python3
"""Predict signal / pro / mature boundaries from cleavage signatures ALONE
(no known mature peptide), then benchmark against the known boundaries.

Model:
  signal_end     = deepsig signal-peptide cleavage
  pro->mature    = best convertase/protease cut C-terminal to the signal,
                   ranked by enzyme priority + model/PWM score + cationicity of
                   the resulting mature candidate + a plausible-length term
  mature_end     = C-terminus, UNLESS a downstream cut leaves a low-cationicity
                   C-terminal pro-segment (then trim to that cut)

Emits a predicted regions table that supersets define_regions.py's schema
(extra trailing columns), so it can feed extract_methods.py for precursors with
no known mature. When a ground-truth regions table is given, writes a
precision/recall report.
"""
import argparse
import csv
import sys

ROLE_PRIORITY = {"convertase": 3.0, "serine": 1.5, "metallo": 1.0, "aspartic": 1.0}


def read_fasta(path):
    seqs, sid, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = "".join(buf)
                sid, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if sid is not None:
        seqs[sid] = "".join(buf)
    return seqs


def read_signalp(path):
    sig = {}
    try:
        with open(path) as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) >= 5 and "signal" in f[2].lower():
                    try:
                        sig[f[0]] = int(float(f[4]))
                    except ValueError:
                        pass
    except FileNotFoundError:
        pass
    return sig


def cationicity(seg):
    if not seg:
        return 0.0
    return sum(seg.count(x) for x in "KR") / len(seg)


def net_anionic(seg):
    """True if the segment is net negatively charged (D/E > K/R fraction)."""
    if not seg:
        return False
    acidic = sum(seg.count(x) for x in "DE") / len(seg)
    return acidic > cationicity(seg)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--precursors", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--signalp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--truth", help="regions_aa.tsv ground truth (optional)")
    ap.add_argument("--accuracy", help="output precision/recall report")
    ap.add_argument("--min-mature", type=int, default=5)
    ap.add_argument("--max-mature", type=int, default=120)
    ap.add_argument("--tol", type=int, default=2)
    ap.add_argument("--no-cterm-pro", action="store_true",
                    help="disable C-terminal pro-segment trimming")
    ap.add_argument("--min-cterm-tail", type=int, default=3,
                    help="shortest C-terminal pro-segment worth trimming")
    args = ap.parse_args()

    seqs = read_fasta(args.precursors)
    sig = read_signalp(args.signalp)

    # group candidate bonds per precursor: (bond, role, model_score)
    cand = {}
    with open(args.candidates) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                sc = float(row.get("score", "0") or 0)
            except ValueError:
                sc = 0.0
            cand.setdefault(row["precursor_id"], []).append(
                (int(row["cut_after_aa"]), row["role"], sc))

    predictions = {}   # precursor -> dict
    for pid, seq in seqs.items():
        se = sig.get(pid, 0)
        L = len(seq)
        # --- (b) score-aware pro->mature selection ---
        best, best_score = None, -1e9
        for bond, role, msc in cand.get(pid, []):
            if bond <= se:
                continue
            mature = seq[bond:]
            mlen = len(mature)
            if mlen < args.min_mature or mlen > args.max_mature:
                continue
            score = (ROLE_PRIORITY.get(role, 1.0)
                     + 3.0 * msc
                     + 2.0 * cationicity(mature)
                     + (0.5 if args.min_mature <= mlen <= 60 else 0.0))
            if score > best_score:
                best_score, best = score, bond
        mstart = (best + 1) if best is not None else (se + 1)

        # --- (c) optional C-terminal pro-segment trimming ---
        # AMP C-terminal pro-segments are characteristically ANIONIC; trim to
        # the most C-terminal cut that leaves a net-anionic tail (removing the
        # acidic pro while keeping the mature as long as possible).
        mend, mend_src = L, "terminus"
        if not args.no_cterm_pro and mstart <= L:
            chosen = None
            for bond, role, msc in sorted(cand.get(pid, []), key=lambda t: t[0]):
                if bond <= mstart or (L - bond) < args.min_cterm_tail:
                    continue
                mature = seq[mstart - 1:bond]
                if not (args.min_mature <= len(mature) <= args.max_mature):
                    continue
                tail = seq[bond:]
                if net_anionic(tail) and cationicity(mature) > cationicity(tail):
                    chosen = bond          # keep most C-terminal qualifying cut
            if chosen is not None:
                mend, mend_src = chosen, "cterm_cut"

        predictions[pid] = {
            "signal_end": se,
            "pro_start": se + 1 if mstart - 1 >= se + 1 else 0,
            "pro_end": mstart - 1 if mstart - 1 >= se + 1 else 0,
            "mature_start": mstart,
            "mature_end": mend,
            "precursor_len": L,
            "sel_score": best_score if best is not None else 0.0,
            "mature_end_source": mend_src,
        }

    with open(args.out, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["peptide_id", "precursor_id", "gene_symbol", "precursor_len",
                    "signal_end", "mature_start", "mature_end",
                    "pro_start", "pro_end", "dibasic_hint",
                    "sel_score", "mature_end_source"])
        for pid, p in sorted(predictions.items()):
            w.writerow([pid + "_pred", pid, ".", p["precursor_len"],
                        p["signal_end"], p["mature_start"], p["mature_end"],
                        p["pro_start"], p["pro_end"], ".",
                        f"{p['sel_score']:.3f}", p["mature_end_source"]])

    # Benchmark against ground truth
    if args.truth and args.accuracy:
        truth = {}
        with open(args.truth) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                truth[row["precursor_id"]] = int(row["mature_start"])
        rows, exact, covered, n = [], 0, 0, 0
        for pid, true_ms in truth.items():
            if pid not in predictions:
                continue
            n += 1
            pred_ms = predictions[pid]["mature_start"]
            hit = abs(pred_ms - true_ms) <= args.tol
            any_cand = any(abs(b - true_ms) <= args.tol
                           for b, _r, _s in cand.get(pid, []))
            exact += hit
            covered += any_cand
            rows.append([pid, true_ms, pred_ms, int(hit), int(any_cand)])
        with open(args.accuracy, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["precursor_id", "true_mature_start",
                        "pred_mature_start", "boundary_hit", "scanner_covered"])
            w.writerows(rows)
            w.writerow([])
            w.writerow(["# precursors_evaluated", n])
            w.writerow(["# boundary_accuracy(+/-%dres)" % args.tol,
                        f"{exact / n:.3f}" if n else "NA"])
            w.writerow(["# scanner_recall(+/-%dres)" % args.tol,
                        f"{covered / n:.3f}" if n else "NA"])
        sys.stderr.write(
            f"predict_boundaries: n={n} accuracy={exact}/{n} "
            f"scanner_recall={covered}/{n}\n")


if __name__ == "__main__":
    main()
