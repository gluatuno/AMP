# Frog AMP input (APD3)

Place the APD3 export here so the pipeline can read it offline:

- `apd_frog.fasta` — mature peptide FASTA exported from APD
  (https://aps.unmc.edu/ → search by amphibian source organisms).
- `apd_frog_meta.tsv` *(optional)* — `peptide_id<TAB>organism<TAB>family<TAB>name`
  to give each peptide a source species (used to pick the genome). Without it,
  species is inferred from the peptide name/family via `config/frog_families.tsv`.

APD3 has no bulk API and is blocked by this environment's network policy, so the
download step cannot run here; export once and drop the file in this folder, or
run the pipeline where `aps.unmc.edu` is reachable.
