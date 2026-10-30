# Context-aware AMP predictor

A better AMP prediction tool by using **gene context**, not just the peptide.

Most AMP predictors (AMPScanner, amPEP, CAMP…) score the **mature peptide
sequence** alone — that is **Method 1** in the reference figure. But a real AMP
sits in a characteristic gene context: a precursor with a signal peptide and an
(often anionic) **pro-peptide**, frequently beside a processing **enzyme** and a
**transporter**, as one locus. This workflow encodes that context as features
and measures how much it improves prediction.

```
Method 1  mature peptide            baseline (what current tools use)
Method 2  + pro-peptide             processing / anionic pro-region signal
Method 3  + full precursor          signal peptide + pro + mature
Method 4  + processing enzyme       convertase/protease cut score & proximity
Method 5  + whole locus             exon structure, GC, transporter/enzyme nearby
```

Run: `snakemake -s workflow/Snakefile_predict --use-conda --cores 4`

## What it delivers

1. **Method 1-5 ablation** (`results_predict/model/ablation.tsv`): cross-validated
   ROC-AUC for each cumulative feature set (M1, M1-2, … M1-5) with the delta each
   added Method contributes — i.e. *how much gene context beats peptide-only*.
2. **Deployable classifier** (`model/model_M1.pkl`): sequence-only, scores any
   peptide FASTA (`predict_amp.py`) — directly comparable to existing tools, no
   genome needed at predict time.
3. **Full-context model** (`model/model_full.pkl`): M1-5, for the **stage-2
   genome scanner** (below).

## Two prediction stages (per your choice: classifier first, scanner later)

- **Stage 1 — classifier (this workflow).** Input a peptide/precursor, output an
  AMP probability. M1-M3 features are sequence-derivable, so the deployable model
  needs no genome at predict time.
- **Stage 2 — genome scanner (next).** Slide over an annotated genome, build
  M1-M5 features per candidate locus (reusing the genome pipeline + this
  `model_full.pkl`), and call AMP loci. The full-context model is trained here so
  it is ready to plug in.

## Inputs

| Config key | What |
|------------|------|
| `dataset.positive_records` | AMP records (schema below). Build from the genome pipeline with `assemble_positives.py`. |
| `dataset.neg_fasta` | **Real non-AMP proteins** — the primary negatives. Required for a meaningful peptide-only model. |
| `dataset.shuffle_per_pos` | Composition-matched shuffled decoys — a **hard control**, not training data (see below). |
| `predict.query_fasta` | Peptides to score with the deployable model. |

`assemble_positives.py` builds `positives.tsv` from the genome pipeline's
`regions_aa.tsv` + `precursors.faa` + peptides (M1-M3), plus the cleavage scored
sites (M4) and locus features (M5) when available.

## Negative design (important)

- **Train** the peptide-only (M1) model against **real non-AMP proteins**
  (different amino-acid composition). This is what lets composition features
  discriminate.
- **Composition-matched shuffled decoys** preserve length and composition, so by
  construction M1 **cannot** separate them (expect AUC ≈ 0.5). That is not a
  failure — it is the rigorous test: if the model beats shuffled decoys, the
  signal must come from the **context Methods (2-5)**, not just composition.
- Evaluate generalisation with grouped CV (e.g. leave-one-family-out) so near-
  duplicate AMP families don't leak between train and test.

## Dataset schema (`positives.tsv` / dataset)

`id, mature_seq, pro_seq, precursor_seq, has_signal, signal_len, has_dibasic,`
`enzyme_present, enzyme_dist_kb, cut_score, n_strong_cuts,`
`locus_len_bp, intron_frac, gc_locus, n_exons, transporter_present,`
`gene_density, neighbor_amp_count`

M1-M3 come from sequence; M4-M5 are the genomic-context fields (0 when unknown).

## Implementation notes

- **Backend:** scikit-learn RandomForest when available; a pure-Python logistic
  model (`ml_lite.py`) otherwise, so the ablation always runs.
- **Tested offline:** feature blocks, ablation ordering/deltas, dataset assembly,
  the pipeline bridge, and the CLI all pass on synthetic data (M1 AUC 1.0 vs real
  non-AMP negatives; cationic AMP scored 0.99 vs non-AMP controls ~0.15).
- **Scope:** this is a research/benchmark pipeline. For a published tool, use a
  large curated positive set (APD/DBAASP) and a matched non-AMP negative set,
  and report grouped-CV ROC/PR against peptide-only baselines.
