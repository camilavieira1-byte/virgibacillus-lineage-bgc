#!/usr/bin/env bash
#SBATCH --job-name=dl_pfam
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/dl_pfam_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/dl_pfam_%j.err
set -euo pipefail

module purge
module load hmmer/3.3.2

D=$HOME/db/pfam
URL=https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
mkdir -p $D && cd $D

echo "== download"
wget -c -t 5 --timeout=60 $URL/Pfam-A.hmm.gz
wget -t 3 --timeout=60 $URL/Pfam.version.gz  || echo "no version file"
wget -t 3 --timeout=60 $URL/md5_checksums     || echo "no checksum file"

echo "== integrity"
gzip -t Pfam-A.hmm.gz && echo "gzip OK"
if [ -s md5_checksums ]; then
  grep " Pfam-A.hmm.gz$" md5_checksums | md5sum -c - || { echo "MD5 MISMATCH"; exit 1; }
fi

echo "== unpack and press"
gunzip -f Pfam-A.hmm.gz
[ "$(stat -c %s Pfam-A.hmm)" -gt 1000000000 ] || { echo "Pfam-A.hmm looks too small"; exit 1; }
hmmpress -f Pfam-A.hmm
zcat Pfam.version.gz 2>/dev/null | head -n 3 || true
ls -lh $D
echo "Pfam ready."
