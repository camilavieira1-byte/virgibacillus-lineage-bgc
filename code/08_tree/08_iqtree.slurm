#!/usr/bin/env bash
#SBATCH --job-name=virgi_iqtree
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/iqtree_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/iqtree_%j.err
set -euo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc
BASE=/PATH/TO/virgibacillus_v2
cd $BASE/08_tree
IQ=$(command -v iqtree2 || command -v iqtree)
OG=$(python - << 'PY'
import pandas as pd
s = pd.read_csv("/PATH/TO/virgibacillus_v2/04_gtdbtk_b/gtdb_classification.tsv", sep="\t")
print("usr_" + s[s.genus.fillna("").str.startswith("Paenibacillus")].accession.iloc[0])
PY
)
echo "Outgroup tip: $OG | IQ-TREE: $IQ"
$IQ -s bac120_user_msa.faa -m MFP -bb 1000 -alrt 1000 -T AUTO --threads-max 16 -o "$OG" --prefix virgi_iqtree
echo "IQ-TREE finished."
