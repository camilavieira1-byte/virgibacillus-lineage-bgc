#!/usr/bin/env bash
#SBATCH --job-name=virgi_gtdb3
#SBATCH --cpus-per-task=8
#SBATCH --mem=300G
#SBATCH --time=12:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/gtdbtk3_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/gtdbtk3_%j.err
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc
export GTDBTK_DATA_PATH="$HOME/db/gtdbtk_db"
BASE=/PATH/TO/virgibacillus_v2

gtdbtk classify \
  --genome_dir $BASE/02_isolates/fna \
  --align_dir $BASE/04_gtdbtk_b \
  --out_dir $BASE/04_gtdbtk_b \
  -x fna \
  --cpus 8 \
  --pplacer_cpus 1
echo "GTDB-Tk classify (no scratch) finished."
