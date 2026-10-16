# Running on HCC (Holland Computing Center)

The pipeline runs on HCC's SLURM clusters (e.g. **Swan**). HCC compute nodes
have outbound internet, so the Ensembl/UniProt/NCBI downloads work there.

## 0. One-time: log in and get the code

```bash
ssh <user>@swan.unl.edu          # your HCC cluster login host
git clone <this-repo-url> amp && cd amp
git checkout claude/amp-genomics-pipeline-7xvfpl
```

## 1. Put big data on $WORK, keep code on $HOME

`$HOME` has a small, backed-up quota; `$WORK` is large Lustre scratch (not
backed up, periodically purged). Genomes + results belong on `$WORK`:

```bash
mkdir -p "$WORK/amp/resources" "$WORK/amp/results"
rm -rf resources results                       # remove the in-repo placeholders
ln -s "$WORK/amp/resources" resources
ln -s "$WORK/amp/results"   results
```

(Alternatively, edit `resources_dir:` and `outdir:` in `config/config.yaml`
to absolute `$WORK` paths instead of symlinking.)

## 2. Create the conda environment

```bash
module purge
module load anaconda                 # provides conda/mamba on HCC
mamba env create -f environment.yml  # one-time; ~10 min
conda activate amp-genomics
```

## 3. Edit the two cluster-specific values

In **`hcc/profiles/slurm/config.yaml`** (and/or the `#SBATCH` lines in
`hcc/run.slurm` / `hcc/controller.slurm`):

- `slurm_partition` — your partition (Swan default: `batch`).
- `slurm_account` — only if your HCC group requires `-A <group>` (check with
  `sacctmgr show assoc user=$USER format=account,partition`).

Resource defaults are tuned for a **human-scale genome**; lower `mem_mb` /
`runtime` for smaller genomes.

## 4. Dry run (always do this first)

```bash
snakemake -n -p                      # shows the DAG without running anything
```

## 5. Submit — pick one style

**A. Single allocation (simplest):** everything in one job.
```bash
sbatch hcc/run.slurm
```

**B. SLURM-executor profile (scales best):** a small controller job submits
each rule as its own sized sbatch job.
```bash
sbatch hcc/controller.slurm
# or, interactively on the login node (lightweight controller):
# snakemake --profile hcc/profiles/slurm
```

## 6. Monitor

```bash
squeue -u $USER                      # your running/queued jobs
tail -f logs/amp-*.out               # controller / run log
ls results/                          # outputs appear here (-> $WORK)
```

## Notes & tips

- **First run downloads references** (GRCh38 ~3 GB + indexes, UniProt, GFF3)
  into `resources/` on `$WORK`. Budget disk and the `download_*` walltimes.
- **Resume after a killed job:** just resubmit — `--rerun-incomplete` picks up
  where it stopped; finished outputs are not recomputed.
- **Per-rule conda envs** (`workflow/envs/*.yaml`) are built on first use with
  `--use-conda`; this is automatic in both submit scripts.
- **Don't run heavy steps on the login node.** The profile controller is light
  (it only submits/waits), so login-node or a tiny job is fine for it.
- **Partition/cluster names change** over time at HCC — confirm against the
  current HCC docs if `batch`/`swan.unl.edu` no longer match.
