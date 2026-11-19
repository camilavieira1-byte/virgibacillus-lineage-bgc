import pandas as pd, numpy as np
base = "/PATH/TO/virgibacillus_v2"
m = pd.read_csv(f"{base}/06_antismash/mibig_best_per_region.tsv", sep="\t")
m["gbk"] = m.accession + "_" + m.file.str.replace(".gbk", "", regex=False)
rec = pd.read_csv(f"{base}/09_matrix/gcf_records_c0.5.tsv", sep="\t")
x = rec.merge(m[["gbk", "best_similarity", "mibig_desc"]], on="gbk", how="left")
x["best_similarity"] = x.best_similarity.fillna(0); x["ge50"] = x.best_similarity >= 50
rows = []
for fam, d in x.groupby("family"):
    i = d.best_similarity.idxmax()
    rows.append({"family": fam, "regions": len(d), "genomes": d.genome.nunique(), "product": d["product"].value_counts().index[0],
                 "median_sim": d.best_similarity.median(), "max_sim": d.best_similarity.max(),
                 "pct_members_ge50": round(100 * d.ge50.mean()), "best_member_desc": d.loc[i, "mibig_desc"] if d.loc[i, "best_similarity"] > 0 else ""})
g = pd.DataFrame(rows).set_index("family")
g["status"] = np.select([g.pct_members_ge50 >= 50, g.max_sim >= 50, g.max_sim > 0], ["known", "partial", "weak only"], "no hit")
g.to_csv(f"{base}/06_antismash/gcf_novelty_c0.5_v2.tsv", sep="\t")
order = ["known", "partial", "weak only", "no hit"]
print("Families:", len(g))
print("All families:", g.status.value_counts().reindex(order).to_dict())
print("Families in >= 2 genomes (%d):" % int((g.genomes >= 2).sum()), g[g.genomes >= 2].status.value_counts().reindex(order).fillna(0).astype(int).to_dict())
print("Families in 1 genome (%d):" % int((g.genomes == 1).sum()), g[g.genomes == 1].status.value_counts().reindex(order).fillna(0).astype(int).to_dict())
print("\nProduct x status (families, top 12 products):")
top = g["product"].value_counts().head(12).index
print(pd.crosstab(g[g["product"].isin(top)]["product"], g[g["product"].isin(top)].status).reindex(columns=order).fillna(0).astype(int).to_string())
print("\nFamilies in >= 5 genomes (name from the best-matching member):")
pd.set_option("display.width", 220)
print(g[g.genomes >= 5].sort_values("genomes", ascending=False)[["product", "genomes", "median_sim", "max_sim", "pct_members_ge50", "status", "best_member_desc"]].to_string())
