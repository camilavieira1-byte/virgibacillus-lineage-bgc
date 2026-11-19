#!/usr/bin/env bash
#SBATCH --job-name=virgi_gtdb
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/gtdbtk_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/gtdbtk_%j.err
set -euo pipefail

# Activate environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc

# Point GTDB-Tk to your home directory database
export GTDBTK_DATA_PATH="$HOME/db/gtdbtk_db"

BASE="/PATH/TO/virgibacillus_v2"

# Run classification on the QC-passed isolates
gtdbtk classify_wf \
  --genome_dir $BASE/02_isolates/fna \
  -x fna \
  --out_dir $BASE/04_gtdbtk \
  --cpus 16 \
  --pplacer_cpus 1
