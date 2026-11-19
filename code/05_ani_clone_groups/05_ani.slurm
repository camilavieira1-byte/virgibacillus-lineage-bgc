#!/usr/bin/env bash
#SBATCH --job-name=virgi_ani
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/ani_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/ani_%j.err
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc

BASE=/PATH/TO/virgibacillus_v2
mkdir -p $BASE/05_ani

# genomes that passed QC (CheckM2 names have no usr_ prefix)
awk -F'\t' 'NR>1 && $2>=90 && $3<5 {print "'$BASE'/02_isolates/fna/usr_"$1".fna"}' \
  $BASE/03_qc/checkm2/quality_report.tsv > $BASE/05_ani/genome_list.txt
echo "Genomes in list: $(wc -l < $BASE/05_ani/genome_list.txt)"
head -n 2 $BASE/05_ani/genome_list.txt

skani triangle -l $BASE/05_ani/genome_list.txt -t 16 -E -s 80 \
  -o $BASE/05_ani/skani_edges.tsv
echo "skani finished."
