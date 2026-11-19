#!/usr/bin/env bash
#SBATCH --job-name=virgi_asm_slow
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/asmslow_%A_%a.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/asmslow_%A_%a.err
set -euo pipefail

module purge
module load antismash/8.0.4

BASE=/PATH/TO/virgibacillus_v2
LIST=$BASE/05_ani/genome_list.txt
FA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")
ACC=$(basename "$FA" .fna | sed 's/^usr_//')
OUT=$BASE/06_antismash/$ACC
mkdir -p $BASE/06_antismash

if [ -f "$OUT/${ACC}.json" ]; then echo "$ACC already done, skipping"; exit 0; fi
rm -rf "$OUT"

LONG=""
if antismash --help 2>&1 | grep -q -- "--allow-long-headers"; then LONG="--allow-long-headers"; fi

echo "Task ${SLURM_ARRAY_TASK_ID}: $ACC ($FA)"
antismash "$FA" \
  --taxon bacteria \
  --genefinding-tool prodigal \
  --cpus ${SLURM_CPUS_PER_TASK} \
  --output-dir "$OUT" \
  --output-basename "$ACC" \
  --cb-general --cb-knownclusters --cb-subclusters \
  --pfam2go --verbose --logfile $BASE/06_antismash/${ACC}.log \
  $LONG

echo "Done $ACC"
