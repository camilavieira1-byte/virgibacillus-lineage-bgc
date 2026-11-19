import pandas as pd, numpy as np
base = "/PATH/TO/virgibacillus_v2"
m = pd.read_csv(f"{base}/06_antismash/mibig_best_per_region.tsv", sep="\t")
m["gbk"] = m.accession + "_" + m.file.str.replace(".gbk", "", regex=False)
rec = pd.read_csv(f"{base}/09_matrix/gcf_records_c0.5.tsv", sep="\t")
x = rec.merge(m[["gbk", "best_similarity", "mibig", "mibig_desc"]], on="gbk", how="left")
print("Regions in the BiG-SCAPE table without a MIBiG row:", int(x.best_similarity.isna().sum()))
x["best_similarity"] = x.best_similarity.fillna(0)
def top(s):
    s = s.dropna(); return s.value_counts().index[0] if len(s) else ""
g = x.groupby("family").agg(regions=("gbk", "size"), genomes=("genome", "nunique"), max_sim=("best_similarity", "max"),
                            product=("product", top), mibig_desc=("mibig_desc", top))
g["status"] = pd.cut(g.max_sim, [-1, 0, 49, 79, 101], labels=["no hit", "low (1-49%)", "medium (50-79%)", "high (>=80%)"])
g.to_csv(f"{base}/06_antismash/gcf_novelty_c0.5.tsv", sep="\t")
print("\nGene cluster families (cutoff 0.5): %d" % len(g))
print("All families:"); print(g.status.value_counts().to_string())
print("\nFamilies present in >= 2 genomes (%d):" % int((g.genomes >= 2).sum())); print(g[g.genomes >= 2].status.value_counts().to_string())
print("\nFamilies in one genome only (%d):" % int((g.genomes == 1).sum())); print(g[g.genomes == 1].status.value_counts().to_string())
top12 = g["product"].value_counts().head(12).index
print("\nFamilies by product x status (top 12 products):")
print(pd.crosstab(g[g["product"].isin(top12)]["product"], g[g["product"].isin(top12)].status).to_string())
print("\nFamilies with a hit at >= 50%% similarity:")
print(g[g.max_sim >= 50].sort_values("genomes", ascending=False).head(25)[["product", "regions", "genomes", "max_sim", "mibig_desc"]].to_string())
h = x[x.best_similarity >= 50].groupby(["product", "mibig_desc"]).size().sort_values(ascending=False).head(15)
print("\nRegions at >= 50%: product x known cluster"); print(h.to_string())
