#!/bin/bash
# Controller job for the SLURM-executor profile: this small job runs the
# snakemake process, which SUBMITS each rule as its own sbatch job. Keep it
# lightweight but give it a long walltime (it waits for the children).
#
# Submit from the repo root:   sbatch hcc/controller.slurm
#
#SBATCH --job-name=amp-ctl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=24:00:00
#SBATCH --partition=batch
# #SBATCH --account=YOUR_GROUP        # uncomment + set if your group requires it
#SBATCH --output=logs/amp-ctl-%j.out
#SBATCH --error=logs/amp-ctl-%j.err

set -euo pipefail
mkdir -p logs

module purge
module load anaconda            # HCC: provides conda/mamba
conda activate amp-genomics

cd "$SLURM_SUBMIT_DIR"

# Snakemake submits child jobs per the profile in hcc/profiles/slurm/.
snakemake --profile hcc/profiles/slurm "$@"
