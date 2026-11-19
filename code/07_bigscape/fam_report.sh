#!/usr/bin/env bash
# usage: fam_report.sh FAM_00030 [cutoff, default 0.5]
set -euo pipefail
FAM=${1:?usage: fam_report.sh FAM_xxxxx [cutoff]}
C=${2:-0.5}
OUT=/PATH/TO/virgibacillus_v2/07_bigscape/out_glocal/output_files
GBK=/PATH/TO/virgibacillus_v2/07_bigscape/input_fixed

tmp=$(mktemp); trap 'rm -f $tmp' EXIT

# members of the family: bin <tab> region name
for f in $OUT/*_c${C}/*/*_clustering_c${C}.tsv; do
  bin=$(basename "$(dirname "$f")")
  awk -F'\t' -v fam="$FAM" -v bin="$bin" '$6==fam {print bin"\t"$2}' "$f"
done > $tmp

if [ ! -s $tmp ]; then echo "$FAM not found at cutoff $C"; exit 1; fi

echo "##### $FAM at cutoff $C"
echo "Bin(s):";        cut -f1 $tmp | sort | uniq -c
echo "Regions:        $(wc -l < $tmp)"
echo "Distinct genomes: $(cut -f2 $tmp | cut -d_ -f1,2 | sort -u | wc -l)"

echo; echo "antiSMASH product(s) of the members:"
for m in $(cut -f2 $tmp); do
  awk '/^     region/{f=1;next} f&&/^     [^ ]/{exit}
       f&&/\/product=/{gsub(/.*\/product="|"$/,""); print}' "$GBK/$m.gbk"
done | sort | uniq -c | sort -rn

echo; echo "Contig edge (True = possibly fragmented cluster):"
for m in $(cut -f2 $tmp); do
  awk '/^     region/{f=1;next} f&&/^     [^ ]/{exit}
       f&&/\/contig_edge=/{gsub(/.*\/contig_edge="|"$/,""); print}' "$GBK/$m.gbk"
done | sort | uniq -c

echo; echo "Region length (bp): min / median / max"
for m in $(cut -f2 $tmp); do head -1 "$GBK/$m.gbk" | awk '{print $3}'; done \
  | sort -n | awk '{a[NR]=$1} END{print a[1]" / "a[int((NR+1)/2)]" / "a[NR]}'

echo; echo "Where these same regions fall at each cutoff (count, family):"
for c in 0.3 0.5 0.7; do
  echo "-- cutoff $c"
  cat $OUT/*_c${c}/*/*_clustering_c${c}.tsv \
    | awk -F'\t' 'NR==FNR{m[$1]=1;next} ($2 in m){print $6}' <(cut -f2 $tmp) - \
    | sort | uniq -c | sort -rn
done
