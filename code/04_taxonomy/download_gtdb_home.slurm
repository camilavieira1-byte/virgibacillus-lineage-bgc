#!/usr/bin/env bash
#SBATCH --job-name=dl_gtdb
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/dl_gtdb_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/dl_gtdb_%j.err
set -euo pipefail

# Initialize conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc

# Force GTDB-Tk to use your home directory
mkdir -p $HOME/db/gtdbtk_db
export GTDBTK_DATA_PATH="$HOME/db/gtdbtk_db"

echo "Temporarily disabling wget SSL checks..."
echo "check_certificate = off" > ~/.wgetrc

echo "Starting download to $GTDBTK_DATA_PATH..."
download-db.sh

echo "Restoring wget security..."
rm ~/.wgetrc

echo "Verifying installation..."
gtdbtk check_install

echo "Done."
