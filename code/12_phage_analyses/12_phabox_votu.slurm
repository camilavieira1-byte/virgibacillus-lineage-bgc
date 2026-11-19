#!/usr/bin/env bash
#SBATCH --job-name=virgi_phabox_votu
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/phabox_votu_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/phabox_votu_%j.err
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate phabox_env
BASE=/PATH/TO/virgibacillus_v2/12_phages
DB=/PATH/TO/project/phabox_db_v2_2
rm -rf $BASE/phabox_votu_ani $BASE/phabox_votu_aai
phabox2 --task votu --dbdir $DB --outpth $BASE/phabox_votu_ani --contigs $BASE/prophages_short.fna --threads 8 --mode ANI --ani 95 --tcov 85
phabox2 --task votu --dbdir $DB --outpth $BASE/phabox_votu_aai --contigs $BASE/prophages_short.fna --threads 8 --mode AAI
for d in ani aai; do echo "== $d"; find $BASE/phabox_votu_$d -maxdepth 3 -name "*.tsv" | head; done
