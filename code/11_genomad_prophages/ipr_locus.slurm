#!/usr/bin/env bash
#SBATCH --job-name=ipr_locus
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/ipr_locus_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/ipr_locus_%j.err

set -euo pipefail
cd ~/virgibacillus_v2/11_genomad

module purge
module load interproscan/5.61-93.0

interproscan.sh -i locus_AUQA_12-22kb.faa -f tsv -o locus_AUQA.iprscan.tsv -dp -cpu 4
