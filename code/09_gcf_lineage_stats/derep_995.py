import pandas as pd, re
base = "/PATH/TO/virgibacillus_v2"
cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str)
g = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str)
ed = pd.read_csv(f"{base}/05_ani/skani_edges.tsv", sep="\t")
acc = lambda s: re.search(r"(GCA_\d+\.\d)", s).group(1)
ed["a"], ed["b"] = ed.Ref_file.map(acc), ed.Query_file.map(acc)
ed["af"] = ed[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
par = {a: a for a in g.accession}
def f(x):
    while par[x] != x: par[x] = par[par[x]]; x = par[x]
    return x
for r in ed[(ed.ANI >= 99.5) & (ed.af >= 50)].itertuples():
    if r.a in par and r.b in par: par[f(r.a)] = f(r.b)
cg["grp995"] = cg.accession.map({a: f(a) for a in g.accession})
cg["comp"] = pd.to_numeric(cg.Completeness); cg["n50"] = pd.to_numeric(cg.Contig_N50)
best = cg.sort_values(["comp", "n50"], ascending=False).groupby("grp995").head(1)
out = cg[cg.accession.isin(best.accession)].drop(columns=["comp", "n50", "grp995"])
out.to_csv(f"{base}/09_matrix/analysis_table_clonegroups_995.tsv", sep="\t", index=False)
b = pd.read_csv(f"{base}/10_pangenome/b_reps.tsv", sep="\t", dtype=str); b = b[b.accession.isin(out.accession)]
b.to_csv(f"{base}/10_pangenome/b_reps_995.tsv", sep="\t", index=False)
b.accession.to_csv(f"{base}/10_pangenome/b_reps_995.txt", index=False, header=False)
print("Representatives at 99.5%%: clade A' %d | B %d (expected 74 and 17)" % ((out.in_A == "True").sum(), (out.in_B == "True").sum()))
print(out[out.in_B == "True"].species.value_counts().to_dict())
