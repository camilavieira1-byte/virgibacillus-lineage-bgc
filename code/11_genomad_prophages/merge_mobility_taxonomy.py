import pandas as pd

base = "/PATH/TO/virgibacillus_v2"

# 1. Load geNomad mobility summary
mob = pd.read_csv(f"{base}/11_genomad/mobility_summary_per_genome.tsv", sep="\t")

# Recreate 'total' if it is missing from the saved TSV
if "total" not in mob.columns:
    mob["total"] = mob.get("viruses", 0) + mob.get("plasmids", 0)

# 2. Load GTDB-Tk taxonomy summary
gtdb = pd.read_csv(f"{base}/04_gtdbtk_b/gtdbtk.bac120.summary.tsv", sep="\t")

# 3. Extract Species name
def extract_species(classification):
    if pd.isna(classification):
        return "Unknown"
    parts = classification.split(";")
    for p in parts:
        if p.startswith("s__"):
            species = p.replace("s__", "")
            return species if species else "Unclassified Species"
    return "Unclassified Species"

gtdb["species"] = gtdb["classification"].apply(extract_species)

# 4. Merge datasets
merged = mob.merge(gtdb[["user_genome", "species"]], left_on="genome", right_on="user_genome", how="left")
merged = merged.drop(columns=["user_genome"])

# Reorder columns
cols = ["genome", "species", "viruses", "proviruses", "plasmids", "total"]
merged = merged[cols].sort_values("total", ascending=False)

# Save merged table
merged.to_csv(f"{base}/11_genomad/mobility_with_taxonomy.tsv", sep="\t", index=False)

print("Top 10 genomes with GTDB species taxonomy:")
print(merged.head(10).to_string(index=False))

print("\nSummary of Mobile Elements per Species (Mean per genome):")
species_summary = merged.groupby("species").agg(
    n_genomes=("genome", "count"),
    mean_viruses=("viruses", "mean"),
    mean_proviruses=("proviruses", "mean"),
    mean_plasmids=("plasmids", "mean"),
    mean_total=("total", "mean")
).round(2).sort_values("n_genomes", ascending=False)

print(species_summary.to_string())
