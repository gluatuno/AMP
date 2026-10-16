#!/bin/bash
# Single-allocation run: the whole pipeline runs inside ONE job with N cores.
# Simplest way to run on HCC; good for this pipeline's size. For large fan-out
# use the SLURM-executor profile instead (sbatch hcc/controller.slurm).
#
# Submit from the repo root:   sbatch hcc/run.slurm
#
#SBATCH --job-name=amp-pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G                    # miniprot whole-genome index is the peak
#SBATCH --time=12:00:00
#SBATCH --partition=batch            # HCC Swan default; change if needed
# #SBATCH --account=YOUR_GROUP        # uncomment + set if your group requires it
#SBATCH --output=logs/amp-%j.out
#SBATCH --error=logs/amp-%j.err

set -euo pipefail
mkdir -p logs

module purge
module load anaconda                 # HCC: provides conda/mamba
conda activate amp-genomics

cd "$SLURM_SUBMIT_DIR"

# -p prints shell commands; --rerun-incomplete recovers from a killed job.
snakemake --use-conda --conda-frontend mamba \
    --cores "${SLURM_CPUS_PER_TASK:-16}" \
    --rerun-incomplete --keep-going -p "$@"
