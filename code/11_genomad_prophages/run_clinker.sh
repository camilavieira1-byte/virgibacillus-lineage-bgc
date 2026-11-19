#!/bin/bash
cd ~/virgibacillus_v2/11_genomad/visualizations_gbk
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate qc

# one plot per family (families with >= 2 members)
for fam in RiPP_FAM_00262 RiPP_FAM_00264 terpene_FAM_00030; do
    n=$(ls ${fam}__*.gbk 2>/dev/null | wc -l)
    [ "$n" -ge 2 ] || continue
    clinker ${fam}__*.gbk -i 0.3 -p ../clinker_${fam}.html \
        -o ../clinker_${fam}_identity.csv -dl "," 2>&1 | grep -v -e "parent gene" -e "conda entry point"
done

# all hits together, to see if any phage region is shared across families
clinker *.gbk -i 0.3 -p ../clinker_all_near.html 2>&1 | grep -v -e "parent gene" -e "conda entry point"
