#!/usr/bin/env bash
#SBATCH --job-name=virgi_pan
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/PATH/TO/virgibacillus_v2/logs/pan_%j.out
#SBATCH --error=/PATH/TO/virgibacillus_v2/logs/pan_%j.err
set -euo pipefail

BASE=/PATH/TO/virgibacillus_v2
PAN=$BASE/10_pangenome
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate qc
rm -rf $PAN/genes $PAN/proteins
mkdir -p $PAN/genes $PAN/proteins
echo "Prodigal: $(prodigal -v 2>&1 | head -n 1)"

run_prodigal() {
  a=$1
  prodigal -q -g 11 -p single -i "$PAN/fna/$a.fna" \
    -a "$PAN/genes/$a.faa" -d "$PAN/genes/$a.ffn" -f gff -o "$PAN/genes/$a.gff"
  sed '/^>/! s/\*//g' "$PAN/genes/$a.faa" > "$PAN/proteins/$a.faa"
}
export -f run_prodigal; export PAN
xargs -a $PAN/pan_set.txt -P ${SLURM_CPUS_PER_TASK:-16} -I{} bash -c 'run_prodigal {}'

for a in $(cat $PAN/pan_set.txt); do printf "%s\t%s\n" $a $(grep -c '>' $PAN/proteins/$a.faa); done | sort -k2,2n > $PAN/cds_counts.tsv
awk '{print $2}' $PAN/cds_counts.tsv | sort -n | awk '{a[NR]=$1} END{print "CDS per genome - genomes:",NR,"| min:",a[1],"| median:",a[int((NR+1)/2)],"| max:",a[NR]}'

conda deactivate
module purge
module load orthofinder/2.5.5
orthofinder -f $PAN/proteins -t ${SLURM_CPUS_PER_TASK:-16} -a 4 -M msa -A mafft -T fasttree
echo "OrthoFinder done"
ls $PAN/proteins/OrthoFinder/*/
