# Frog AMP → genome training-set pipeline (initial project)

Build a labelled training set for a **genome-wide frog AMP predictor**: collect
natural frog antimicrobial peptides from APD3, map each to the genome of its
**source species**, extract the gene's genomic sequence + features, and assemble
a matrix where AMP loci are positives.

Run: `snakemake -s workflow/Snakefile_frog --use-conda --cores 8`

## The central design choice: species vs. genome

APD "frog" AMPs come from **dozens of species**, most without a sequenced
genome. You can only obtain genomic coordinates for an AMP whose species has a
genome — and a genome predictor must train on positives that genuinely live in
the genome. So this project uses the **few-genomed-species** strategy: keep only
AMPs whose source species has a genome; AMPs from genome-less genera
(*Phyllomedusa*, *Pelophylax*, *Hyla*, *Odorrana*…) are **reported and set
aside**, never force-mapped onto a distant reference.

Target species (`config/frog_species.tsv`, accessions editable/verifiable):

| Species | Assembly | AMP families |
|---------|----------|--------------|
| *Xenopus tropicalis* | UCB_Xtro_10.0 | magainin, PGLa, CPF, XPF |
| *Xenopus laevis* | Xenopus_laevis_v10.1 | magainins (original) |
| *Rana temporaria* | aRanTem1.1 | brevinin, temporin, esculentin, ranatuerin |
| *Nanorana parkeri* | ASM93562v1 | — |
| *Bufo bufo* | aBufBuf1.1 | buforin |
| *Bombina bombina* | aBomBom1 (verify) | bombinin, maximin |

## Flow

```
APD3 frog AMPs (fetch_apd_frog.py; offline export in data/frog/)
   │
 route to genomed species            assign_species.py
   │   (organism → genus → family hint; genome-less → set aside)
   ▼  per species:
 peptide → precursor protein         diamond blastp vs species proteome
 precursor → genome (spliced)        miniprot → CDS/exon model
 annotate                            deepsig signal peptide + pro/mature
 cleavage sites                      scan + PWM/ML score (reused module)
 locus features                      locus_features.py (exons, span, GC, …)
   │
 positives (AMP loci) ──┐
 negatives (random non-AMP windows)  sample_negatives.py
   │                    │
   ▼                    ▼
 feature matrix         build_training_matrix.py
   ▼
 AMP genome model       train_amp_genome_model.py  (RandomForest / baseline)
```

Each step reuses the existing per-step scripts (`best_hit_precursor`,
`define_regions`, the cleavage `scan`/`train`/`score` module, `gmap`) against
each species' own genome + proteome.

## Features used for model construction

Per locus (`frog_features.py`), computed identically for positives and
negatives so they are comparable:

- **Sequence:** length, cationicity (K+R), acidic fraction, net charge,
  hydropathy, aromatic/Cys/Gly fractions, Boman-like index.
- **Genomic:** precursor length, exon count, locus span (bp), CDS length,
  intron fraction, GC of the locus, signal-peptide & pro-region lengths,
  dibasic-site presence, cleavage-site count.

Positives = mapped AMP loci; negatives = random genomic windows not overlapping
AMP loci (translated; `random_orf`/`random_window`/`non_amp_gene` strategies).
The trained model then scores candidate loci genome-wide (next step).

## Outputs (`results_frog/`)

```
amp/peptides.tsv, amp/assignments.tsv        # what mapped where, what was set aside
map/<species>/…                              # per-species precursor, genome map, regions, cleavage
features/<species>.positives.tsv, .loci.bed  # positive feature rows + AMP loci
model/feature_matrix.tsv                     # labelled matrix (1=AMP locus, 0=not)
model/metrics.tsv, feature_importance.tsv    # CV AUC + what drives the model
```

## Status & caveats

- **Tested offline:** species routing, locus-feature assembly, negative
  sampling, matrix build, and the baseline trainer all pass on synthetic data.
- **Needs network:** APD export and NCBI `datasets` genome downloads — this
  environment blocks `aps.unmc.edu` and `api.ncbi.nlm.nih.gov`. Open the
  environment's Network access (or run elsewhere), then the DAG runs end to end.
- **Mature-peptide mapping:** frog AMPs are short; mapping goes via the
  **precursor** (preprodermaseptin/prepromagainin style), not the bare mature
  peptide, which miniprot could not place reliably.
- **Negative design matters:** random-window negatives make signal-peptide /
  exon-structure features look perfectly discriminative; use `non_amp_gene`
  negatives (other annotated genes) for an honest model, and evaluate with
  **leave-one-species-out** CV to measure genome-transfer.
- **Accessions** in `frog_species.tsv` should be confirmed before a run
  (assembly versions change).
