import pandas as pd, numpy as np
base = "/PATH/TO/virgibacillus_v2"
d = pd.read_csv(f"{base}/06_antismash/mibig_hit_detail.tsv", sep="\t")
m = pd.read_csv(f"{base}/06_antismash/mibig_best_per_region.tsv", sep="\t")
x = m.merge(d[["accession", "contig", "region", "n_hits", "known_size"]], on=["accession", "contig", "region"], how="left")
x["n_hits"] = x.n_hits.fillna(0)
labs = ["no hit", "weak (<50%)", "short match (>=50%, <3 genes)", "known (>=50%, >=3 genes)"]
x["strict"] = np.select([x.best_similarity == 0, x.best_similarity < 50, (x.best_similarity >= 50) & (x.n_hits < 3)], labs[:3], labs[3])
print("Regions:", len(x)); print(x.strict.value_counts().reindex(labs).to_string())
print("\nWithout ectoine-containing products (%d regions):" % int((~x["product"].str.contains("ectoine")).sum()))
print(x[~x["product"].str.contains("ectoine")].strict.value_counts().reindex(labs).to_string())
print("\nComplete regions only:"); print(x[~x.contig_edge].strict.value_counts().reindex(labs).to_string())
x["gbk"] = x.accession + "_" + x.file.str.replace(".gbk", "", regex=False)
rec = pd.read_csv(f"{base}/09_matrix/gcf_records_c0.5.tsv", sep="\t").merge(x[["gbk", "strict"]], on="gbk", how="left")
def fam_status(s):
    k = (s == labs[3]).mean()
    if k >= .5: return "known"
    if k > 0: return "partly known"
    if (s == labs[2]).any(): return "short match only"
    if (s == labs[1]).any(): return "weak only"
    return "no hit"
f = rec.groupby("family").strict.agg(fam_status)
order = ["known", "partly known", "short match only", "weak only", "no hit"]
print("\nFamilies (cutoff 0.5): %d" % len(f)); print(f.value_counts().reindex(order).fillna(0).astype(int).to_string())
fp = rec.groupby("family")["product"].agg(lambda s: s.value_counts().index[0]); top = fp.value_counts().head(12).index
print("\nProduct x family status (top 12 products):")
print(pd.crosstab(fp[fp.isin(top)], f[fp.isin(top).index[fp.isin(top)]]).reindex(columns=order).fillna(0).astype(int).to_string())
x.to_csv(f"{base}/06_antismash/mibig_strict_per_region.tsv", sep="\t", index=False)
