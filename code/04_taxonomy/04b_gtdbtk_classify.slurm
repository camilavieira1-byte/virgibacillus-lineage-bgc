#!/usr/bin/env bash
#SBATCH --job-name=virgi_gtdb2
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/gtdbtk2_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/gtdbtk2_%j.err
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc
export GTDBTK_DATA_PATH="$HOME/db/gtdbtk_db"

BASE=/PATH/TO/virgibacillus_v2

gtdbtk classify \
  --genome_dir $BASE/02_isolates/fna \
  --align_dir $BASE/04_gtdbtk \
  --out_dir $BASE/04_gtdbtk \
  -x fna \
  --cpus 16 \
  --pplacer_cpus 1 \
  --scratch_dir $BASE/04_gtdbtk/scratch

echo "GTDB-Tk classify finished."
