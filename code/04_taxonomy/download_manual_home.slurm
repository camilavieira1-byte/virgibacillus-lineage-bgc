#!/usr/bin/env bash
#SBATCH --job-name=dl_gtdb
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/dl_gtdb_manual_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/dl_gtdb_manual_%j.err
set -euo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc

# 1. Setup exact home paths
mkdir -p $HOME/db/gtdbtk_db
export GTDBTK_DATA_PATH="$HOME/db/gtdbtk_db"
cd $HOME/db/gtdbtk_db

# 2. Download manually (bypassing download-db.sh)
echo "check_certificate = off" > ~/.wgetrc

echo "Downloading archive to $GTDBTK_DATA_PATH..."
wget -c https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz

echo "Extracting archive..."
tar -xvzf gtdbtk_r232_data.tar.gz -C $HOME/db/gtdbtk_db --strip-components=1

echo "Cleaning up..."
rm gtdbtk_r232_data.tar.gz
rm ~/.wgetrc

echo "Verifying installation..."
gtdbtk check_install

echo "Done."
