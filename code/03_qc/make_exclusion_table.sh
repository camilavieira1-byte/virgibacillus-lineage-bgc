#!/usr/bin/env bash
set -euo pipefail
cd ~/virgibacillus_v2/03_qc
awk -F'\t' 'BEGIN{OFS="\t"; print "accession","completeness","contamination","reason"}
NR>1 && !($2>=90 && $3<5){
  r="";
  if($2<90) r=r"low_completeness;";
  if($3>=5) r=r"high_contamination;";
  print $1,$2,$3,r
}' checkm2/quality_report.tsv > excluded_qc.tsv
awk -F'\t' 'NR>1 && ($2>=90 && $3<5){print $1}' checkm2/quality_report.tsv > ../00_metadata/pass_qc.txt
echo "Excluded: $(tail -n +2 excluded_qc.tsv | wc -l) | Passed: $(wc -l < ../00_metadata/pass_qc.txt)"
cat excluded_qc.tsv
