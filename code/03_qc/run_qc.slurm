#!/usr/bin/env bash
#SBATCH --job-name=virgi_qc
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/qc_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/qc_%j.err
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc

THREADS=${SLURM_CPUS_PER_TASK:-16}
BASE=/PATH/TO/virgibacillus_v2
FNA=$BASE/02_isolates/fna
OUT=$BASE/03_qc
DB=$(find $HOME/db -name "uniref100.KO.1.dmnd" | head -n1)
echo "Threads: $THREADS | CheckM2 DB: $DB"

seqkit stats -a -T -j "$THREADS" $FNA/*.fna > $OUT/seqkit_stats.tsv

quast.py -t "$THREADS" --min-contig 300 --no-plots --no-html --no-icarus \
  -o $OUT/quast $FNA/*.fna

checkm2 predict --threads "$THREADS" --input $FNA -x fna \
  --output-directory $OUT/checkm2 --database_path "$DB" --force

echo "QC finished."
