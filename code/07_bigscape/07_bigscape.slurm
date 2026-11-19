#!/usr/bin/env bash
#SBATCH --job-name=virgi_bigscape
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/bigscape_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/bigscape_%j.err
set -euo pipefail

module purge
module load big-scape/2.0.3

BASE=/PATH/TO/virgibacillus_v2
PFAM=$HOME/db/pfam/Pfam-A.hmm
IN=$BASE/07_bigscape/input_fixed
echo "Input regions: $(ls $IN | wc -l) | Pfam: $PFAM"

for MODE in glocal auto; do
  echo "===== alignment mode: $MODE ====="
  rm -rf $BASE/07_bigscape/out_$MODE
  bigscape cluster \
    --input-dir $IN \
    --input-mode flat \
    --output-dir $BASE/07_bigscape/out_$MODE \
    --pfam-path $PFAM \
    --gcf-cutoffs 0.3,0.5,0.7 \
    --alignment-mode $MODE \
    --cores ${SLURM_CPUS_PER_TASK} \
    --no-trees \
    --label virgi_$MODE \
    --verbose
done
echo "BiG-SCAPE finished."
