import pandas as pd
import glob
import os

base = "/PATH/TO/virgibacillus_v2/11_genomad"

# 1. Aggregate Viruses
virus_files = glob.glob(f"{base}/usr_*/usr_*_summary/*_virus_summary.tsv")
v_df = pd.concat([pd.read_csv(f, sep="\t").assign(genome=os.path.basename(f).split("_virus_")[0]) for f in virus_files], ignore_index=True)
v_df.to_csv(f"{base}/all_viruses.tsv", sep="\t", index=False)

# 2. Aggregate Plasmids
plasmid_files = glob.glob(f"{base}/usr_*/usr_*_summary/*_plasmid_summary.tsv")
p_df = pd.concat([pd.read_csv(f, sep="\t").assign(genome=os.path.basename(f).split("_plasmid_")[0]) for f in plasmid_files], ignore_index=True)
p_df.to_csv(f"{base}/all_plasmids.tsv", sep="\t", index=False)

# 3. Summarize per genome
summary = pd.DataFrame({
    "viruses": v_df.groupby("genome").size(),
    "proviruses": v_df[v_df["topology"] == "Provirus"].groupby("genome").size(),
    "plasmids": p_df.groupby("genome").size()
}).fillna(0).astype(int)

summary.to_csv(f"{base}/mobility_summary_per_genome.tsv", sep="\t")

print(f"Total viruses: {len(v_df)} ({v_df[v_df.topology == 'Provirus'].shape[0]} proviruses)")
print(f"Total plasmids: {len(p_df)}")
print("\nTop 10 genomes by total mobile elements:")
summary["total"] = summary["viruses"] + summary["plasmids"]
print(summary.sort_values("total", ascending=False).head(10).to_string())
