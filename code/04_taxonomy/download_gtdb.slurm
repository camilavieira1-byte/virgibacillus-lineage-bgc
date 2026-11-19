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

echo "Temporarily disabling wget SSL checks..."
echo "check_certificate = off" > ~/.wgetrc

echo "Starting download..."
download-db.sh

echo "Restoring wget security..."
rm ~/.wgetrc

echo "Verifying installation..."
gtdbtk check_install

echo "Done."
