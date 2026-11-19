#!/usr/bin/env bash
#SBATCH --job-name=virgi_genomad
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/genomad_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/genomad_%j.err

set -euo pipefail

module purge
module load genomad/1.7.0

BASE="/PATH/TO/virgibacillus_v2"
IN_DIR="$BASE/02_isolates/fna"
OUT_DIR="$BASE/11_genomad"
DB="$HOME/db/genomad_db"

mkdir -p $OUT_DIR

# Adjust *.fna to *.fasta or *.fa if your assembly extensions are different
echo "Starting geNomad for $(ls $IN_DIR/*.fna 2>/dev/null | wc -l) genomes..."

for FASTA in $IN_DIR/*.fna; do
    # Skip if no files found
    [ -e "$FASTA" ] || continue 
    
    GENOME_ID=$(basename "$FASTA" .fna)
    echo "Processing $GENOME_ID..."
    
    # Run geNomad end-to-end pipeline
    genomad end-to-end \
        --cleanup \
        --threads ${SLURM_CPUS_PER_TASK} \
        "$FASTA" \
        "$OUT_DIR/$GENOME_ID" \
        "$DB"
done

echo "geNomad finished."
