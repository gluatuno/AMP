# AMP genomics pipeline

Recover **genomic DNA and genomic features** for human antimicrobial peptides
(AMPs) starting from their **mature amino-acid sequences**.

Input: `data/human_amp.fasta` — 171 mature AMP peptides (APD-style IDs such as
`AP00307` = cathelicidin/LL-37, `AP00451` = a defensin, `AP00509` = histatin).

## The idea

Mature AMPs are the *processed* products of larger precursors. To get back to
DNA and to the feature blocks in the reference figure, the pipeline walks
backward:

```
mature peptide (aa)
   → full-length precursor protein        (DIAMOND blastp vs human proteome)
   → genome, splice-aware                 (miniprot → CDS/exon model on GRCh38)
   → genomic DNA per feature block         (aa coords → genomic BED → bedtools)
   → precursor feature annotation          (deepsig signal peptide + pro/mature)
   → trans-acting enzyme / transporter     (curated table → Ensembl GFF3 loci)
```

### Why two tracks (important biology)

The figure's layout — **Transporter · Enzyme · Pro-peptide · AMP** in a row — is
the classic *microbial* AMP/bacteriocin **biosynthetic gene cluster (operon)**.
**Human** host-defense peptides are organised differently: one gene encodes a
single precursor (**signal peptide → pro-region → mature AMP**), while the
processing enzyme and any transporter are **separate genes elsewhere**. So:

* **Methods 1–3 & 5 (cis):** derived from the *same* precursor gene.
* **Method 4 (trans):** the enzyme (and transporter) are looked up as separate
  loci via `config/processing_enzymes.tsv` (edit this table freely).

> If your real input is *microbial*, the right tools are **antiSMASH** and
> **BAGEL4**, which call the operon directly. This repo implements the human
> track you selected.

## Figure "Methods" → pipeline outputs

| Method | Figure target        | Output file                                   | How |
|:------:|----------------------|-----------------------------------------------|-----|
| 1 | AMP (mature)              | `results/methods/method1_amp.fa`              | exonic DNA of the mature-peptide codons |
| 2 | Pro-peptide              | `results/methods/method2_propeptide.fa`       | exonic DNA between signal peptide and mature |
| 3 | Pro-peptide + AMP         | `results/methods/method3_pro_plus_amp.fa`     | exonic DNA of pro-region through mature |
| 4 | Enzyme (trans-acting)    | `results/methods/method4_enzyme.fa` (+`.tsv`) | processing-enzyme/transporter gene loci |
| 5 | Transporter…AMP locus     | `results/methods/method5_locus.fa`            | whole precursor gene envelope (introns incl.) |

Plus reports: `results/report/pipeline_summary.tsv` and
`results/report/precursor_feature_map.tsv` (per-peptide aa coordinates + which
methods yielded sequence).

## Cleavage-site signatures (enzyme "cutting sites" in the genome)

Proteases act on the precursor **protein**, but their cleavage motifs are
encoded in the genomic CDS — so a "cutting site" can be found as a sequence
signature and projected back to a genomic coordinate. This module does that in
three parts:

1. **Find** — `scan_cleavage.py` scans every precursor with a curated protease
   motif library (`config/protease_motifs.tsv`: furin/PCSK dibasic, KLK5/7/14,
   PRTN3, MMP7, cathepsin-D), and `discover_cutsites.py` learns the motifs
   *de novo* from your own signal→pro and pro→mature junctions (position-
   frequency matrix always; STREME/MEME when the `cleavage` env is used).
2. **Predict** — `predict_boundaries.py` calls signal/pro/mature boundaries
   from the signatures **alone** (no known mature peptide needed), so the
   pipeline generalises to novel precursors. Because the true boundaries are
   known here, it also writes a **precision/recall** report
   (`boundary_accuracy.tsv`: boundary accuracy + scanner recall).
3. **Annotate** — `cutsites_to_genome.py` maps each cleavage bond through the
   miniprot CDS model to a genomic base, producing a 1-bp **cut-site BED track**
   (`cutsites.bed`) you can load in a browser next to Methods 1–5.

The candidate sites are then **rescored** by a trained model
(`train_cleavage_model.py` → `score_sites.py`): a log-odds **PWM** learned from
the true junctions, plus a **logistic-regression** classifier over one-hot +
physicochemical window features when scikit-learn is available (reports
cross-validated AUC in `model/metrics.tsv`; PWM-only otherwise). The boundary
predictor also trims a **C-terminal pro-segment** when a downstream cut leaves a
net-anionic tail — so mature ends are called correctly for precursors whose
pro-region is C-terminal, not just N-terminal.

Outputs land in `results/cleavage/`:
`candidate_sites.tsv`, `candidate_sites_scored.tsv`, `denovo/motif_report.tsv`,
`model/{pwm.tsv,metrics.tsv,model.pkl}`, `predicted_regions.tsv` (feeds
`extract_methods.py` for novel inputs; carries `mature_end_source` =
`terminus`/`cterm_cut`), `boundary_accuracy.tsv`, `cutsites.bed`.

## Cis-regulatory layer — shared promoter signatures

Beyond the protein-encoded cut sites, the AMP gene and its processing enzyme
may be **co-regulated at the DNA level**. This module (`rules/regulatory.smk`)
extracts the promoter (`[TSS-upstream, TSS+downstream]`) of every AMP gene and
its paired enzyme/transporter gene (`extract_promoters.py`), scans both strands
for TF binding motifs (`config/tf_motifs.tsv`: VDR/RXR — the cathelicidin
vitamin-D element — NF-κB, C/EBP, AP-1, STAT, ISRE, HNF4), and reports the
motifs **shared between each co-regulated pair** (`shared_tfbs.py`), plus a
genome-mapped TFBS BED track. Outputs in `results/regulatory/`:
`shared_signatures.tsv` (per AMP↔enzyme pair: shared / AMP-only / enzyme-only
TFs), `tfbs_hits.tsv`, `tfbs.bed`.

The scan uses a deterministic IUPAC-consensus match (both strands, ≤N
mismatches) — a dependency-free baseline; drop in JASPAR PWMs + FIMO for
calibrated scoring. Keep consensus motifs ≥7 bp (shorter ones match too often).

## Install

```bash
conda env create -f environment.yml
conda activate amp-genomics
# (Snakemake also creates per-rule envs from workflow/envs/ when run with --use-conda)
```

## Run

```bash
# Dry run — see the DAG without executing
snakemake -n --use-conda

# Full run (downloads GRCh38 + Ensembl GFF3 + UniProt human on first use)
snakemake --use-conda --cores 8
```

Pin releases and tune thresholds in `config/config.yaml`
(`diamond_min_pident`, `min_query_cov`, Ensembl release, etc.).

## Layout

```
config/
  config.yaml                 # paths, reference URLs, thresholds
  processing_enzymes.tsv      # curated AMP-family → enzyme/transporter genes (Method 4)
  protease_motifs.tsv         # curated protease cleavage-site signatures
  tf_motifs.tsv               # curated TF binding motifs (cis-regulatory layer)
workflow/
  Snakefile
  rules/       references, prep, homology, genome_map, annotate, methods,
               cleavage, regulatory
  scripts/     best_hit_precursor, define_regions, gmap, extract_methods,
               enzyme_loci, scan_cleavage, discover_cutsites,
               cleavage_features, train_cleavage_model, score_sites,
               predict_boundaries, cutsites_to_genome,
               extract_promoters, shared_tfbs, make_report
  envs/        core.yaml, annotate.yaml, cleavage.yaml
data/human_amp.fasta          # input peptides
resources/                    # downloaded references (gitignored)
results/                      # outputs (gitignored)
```

## Notes & limitations

* **Multi-exon features** (Methods 1–3) emit one FASTA record per exon
  (`...|ex1`, `...|ex2`); concatenate by shared name for the full coding
  sequence. Most short mature AMPs fall within a single exon.
* **Signal-peptide caller:** `deepsig` (eukaryotic model). Swap in SignalP 6.0
  or Phobius in `workflow/envs/annotate.yaml` if you have a licence.
* **Coordinate model** assumes miniprot's best-scoring, in-frame alignment per
  precursor; frameshifted/partial alignments are approximated and flagged.
* **Method 4** enzyme assignments in `processing_enzymes.tsv` are a curated
  starting point (e.g. PRTN3/KLK5 for cathelicidin, MMP7 for enteric
  α-defensins, FURIN as a generic convertase) — refine per your biology.
* Optional domain annotation (InterProScan/Pfam) is stubbed in
  `environment.yml`; enable if you want family-level calls.
* **Cleavage motifs** in `protease_motifs.tsv` are specificity-based starting
  points; the de-novo step and the trained PWM/ML scorer refine them from your
  data. The boundary predictor now trims C-terminal pro-segments when the tail
  is net-anionic (`mature_end_source=cterm_cut`).
* **Cut-site model:** logistic-regression when scikit-learn is present, else a
  log-odds PWM; both trained on the known junctions with random-bond negatives.
  Small precursor sets give optimistic AUC — validate on held-out families.
* **TFBS scan** is a consensus baseline; co-regulation calls are hypotheses to
  confirm with JASPAR PWMs + FIMO and expression/ChIP evidence.
