#!/usr/bin/env python3
"""Method 6: embed sequences with a foundation model into fixed-length vectors.

Emits a TSV: id  <prefix>_0 ... <prefix>_{D-1}  (prefix m6p for protein, m6g for
DNA), which merge_embeddings.py joins onto the dataset so the ablation picks up
M6 automatically.

Backends:
  stub  - deterministic k-mer-hashed embedding (no deps, no download);
          sequence-dependent, so the plumbing + ablation are testable offline.
  esm2  - protein language model via `fair-esm` (needs torch; downloads weights).
  nt    - DNA language model (Nucleotide Transformer / DNABERT) via
          `transformers` (needs torch; downloads weights).
Real backends need network + ideally a GPU; use them on HCC, not in a sandbox.
"""
import argparse
import hashlib
import math
import sys


def read_fasta(path):
    recs, sid, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    recs.append((sid, "".join(buf)))
                sid, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if sid is not None:
        recs.append((sid, "".join(buf)))
    return recs


def _l2(v):
    n = math.sqrt(sum(x * x for x in v)) or 1.0
    return [x / n for x in v]


def stub_embed(seq, dim, k):
    """Deterministic k-mer-hashed, L2-normalised embedding (no dependencies)."""
    seq = (seq or "").upper()
    vec = [0.0] * dim
    if len(seq) >= k:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i + k]
            h = int(hashlib.md5(kmer.encode()).hexdigest(), 16)
            vec[h % dim] += 1.0
    return _l2(vec)


def esm2_embed(seqs, model_name):
    import torch, esm  # noqa
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    bc = alphabet.get_batch_converter()
    model.eval()
    layer = model.num_layers
    out = {}
    with torch.no_grad():
        for sid, seq in seqs:
            _, _, toks = bc([(sid, seq[:1022])])
            rep = model(toks, repr_layers=[layer])["representations"][layer]
            out[sid] = rep[0, 1:len(seq) + 1].mean(0).tolist()  # mean-pool residues
    return out


def nt_embed(seqs, model_name):
    import torch
    from transformers import AutoTokenizer, AutoModel
    tok = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
    model = AutoModel.from_pretrained(model_name, trust_remote_code=True).eval()
    out = {}
    with torch.no_grad():
        for sid, seq in seqs:
            enc = tok(seq, return_tensors="pt", truncation=True, max_length=1000)
            h = model(**enc).last_hidden_state
            out[sid] = h[0].mean(0).tolist()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--kind", choices=["protein", "dna"], required=True)
    ap.add_argument("--backend", default="stub",
                    choices=["stub", "esm2", "nt", "auto"])
    ap.add_argument("--prefix", default="")
    ap.add_argument("--dim", type=int, default=32, help="stub/output dim")
    ap.add_argument("--model-name", default="")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    prefix = args.prefix or ("m6p" if args.kind == "protein" else "m6g")
    k = 2 if args.kind == "protein" else 3
    seqs = read_fasta(args.fasta)

    backend = args.backend
    emb = None
    if backend in ("esm2", "nt", "auto"):
        try:
            if backend == "esm2" or (backend == "auto" and args.kind == "protein"):
                mn = args.model_name or "esm2_t12_35M_UR50D"
                emb = esm2_embed(seqs, mn)
            elif backend == "nt" or (backend == "auto" and args.kind == "dna"):
                mn = args.model_name or "InstaDeepAI/nucleotide-transformer-500m-human-ref"
                emb = nt_embed(seqs, mn)
        except Exception as e:  # missing torch/transformers or no network
            if args.backend == "auto":
                sys.stderr.write(f"embed_sequences: {backend} unavailable ({e}); "
                                 "falling back to stub.\n")
            else:
                sys.stderr.write(f"embed_sequences: backend '{backend}' failed: {e}\n"
                                 "Install torch + fair-esm/transformers and provide "
                                 "network/GPU, or use --backend stub.\n")
                sys.exit(2)
    if emb is None:
        emb = {sid: stub_embed(seq, args.dim, k) for sid, seq in seqs}

    dim = len(next(iter(emb.values()))) if emb else args.dim
    with open(args.out, "w") as out:
        out.write("id\t" + "\t".join(f"{prefix}_{i}" for i in range(dim)) + "\n")
        for sid, _seq in seqs:
            v = emb.get(sid, [0.0] * dim)
            out.write(sid + "\t" + "\t".join(f"{x:.5g}" for x in v) + "\n")
    sys.stderr.write(f"embed_sequences: {len(seqs)} sequences -> {dim}-d "
                     f"{prefix} ({backend}).\n")


if __name__ == "__main__":
    main()
