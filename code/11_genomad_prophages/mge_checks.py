import pandas as pd, numpy as np, glob
from scipy.stats import mannwhitneyu
base = "/PATH/TO/virgibacillus_v2"; gd = f"{base}/11_genomad"
vals = []
for f in glob.glob(f"{gd}/usr_*/usr_*_find_proviruses/*_provirus.tsv"):
    try: d = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
    except Exception: continue
    if "integrases" in d.columns: vals += d.integrases.tolist()
s = pd.Series(vals)
print("== integrase column (raw): provirus rows %d | not NA/empty: %d" % (len(s), int((~s.isin(["NA", ""])).sum())))
print(s[~s.isin(["NA", ""])].head(5).tolist())

g = pd.read_csv(f"{gd}/mge_per_genome.tsv", sep="\t", dtype={"in_A": str, "in_B": str})
cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str)
cg = cg[cg.in_A == "True"][["accession", "species"]]
g = g.merge(cg, on="accession", how="inner")
g["group"] = np.where(g.in_B == "True", "Virgibacillus (B)", np.where(g.genus.fillna("").str.startswith("Virgibacillus"), "other Virgibacillus_X", "other genera"))
g["contig"] = np.where(g.Contig_N50 >= 500e3, "N50 >= 500 kb", "N50 < 500 kb")
g["has"] = g.n_prophage > 0
print("\n== representatives of the clade: %d" % len(g))
print(g.groupby(["group", "contig"]).agg(reps=("accession", "size"), pct_with_prophage=("has", lambda x: round(100 * x.mean())),
      median_prophages=("n_prophage", "median"), prophages=("n_prophage", "sum"), at_contig_end=("n_prophage_edge", "sum"),
      virus_contigs_per_genome=("n_virus_contigs", "mean")).round(2).to_string())
h = g[g.Contig_N50 >= 500e3]; b = h[h.in_B == "True"].n_prophage; o = h[h.in_B != "True"].n_prophage
print("\nN50 >= 500 kb: B (n=%d) vs others (n=%d): Mann-Whitney p = %.2g" % (len(b), len(o), mannwhitneyu(b, o)[1]))
print("\n== by GTDB genus (genera with >= 3 representatives)")
t = g.groupby("genus").agg(reps=("accession", "size"), pct_with_prophage=("has", lambda x: round(100 * x.mean())), median_prophages=("n_prophage", "median"), median_N50_kb=("Contig_N50", lambda x: round(x.median() / 1e3)))
print(t[t.reps >= 3].sort_values("median_prophages", ascending=False).to_string())
print("\n== Virgibacillus (B) by species")
bb = g[g.in_B == "True"]
print(bb.groupby("species").agg(reps=("accession", "size"), pct_with_prophage=("has", lambda x: round(100 * x.mean())), median_prophages=("n_prophage", "median"), median_N50_kb=("Contig_N50", lambda x: round(x.median() / 1e3))).to_string())
