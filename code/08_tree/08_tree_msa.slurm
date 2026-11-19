#!/usr/bin/env bash
#SBATCH --job-name=virgi_msa
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/msa_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/msa_%j.err
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc
export GTDBTK_DATA_PATH="$HOME/db/gtdbtk_db"
BASE=/PATH/TO/virgibacillus_v2
OUT=$BASE/08_tree
mkdir -p $OUT

gtdbtk identify --genome_dir $BASE/02_isolates/fna -x fna --out_dir $OUT --cpus 16
gtdbtk align --identify_dir $OUT --out_dir $OUT --skip_gtdb_refs --cpus 16

echo "Genomes in MSA: $(zcat $OUT/align/gtdbtk.bac120.user_msa.fasta.gz | grep -c '>')"
zcat $OUT/align/gtdbtk.bac120.user_msa.fasta.gz > $OUT/bac120_user_msa.faa
fasttree -wag -gamma $OUT/bac120_user_msa.faa > $OUT/virgi_bac120.fasttree.nwk
echo "Tree done"
