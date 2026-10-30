# Predictor inputs

- `positives.tsv` — AMP positive records (dataset schema). Build from the genome
  pipeline with `workflow/scripts/assemble_positives.py`, or supply directly.
- `nonamp.faa` — real non-AMP proteins (primary negatives; different composition).
  e.g. a random UniProt sample excluding known AMPs.
- `query.faa` — peptides to score with the deployable model.

See `README_predict.md` for the schema and negative-design notes.
