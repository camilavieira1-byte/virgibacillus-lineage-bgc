import pandas as pd, numpy as np
from Bio import Phylo
rng = np.random.default_rng(7)
base = "/PATH/TO/virgibacillus_v2"
cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str)
cg = cg[cg.in_B == "True"].copy()
t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick")
tips = {x.name: x for x in t.get_terminals()}
cg = cg[cg.accession.map(lambda a: "usr_" + a in tips)].reset_index(drop=True)
reps = cg.accession.tolist(); n = len(reps)
P = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        P[i, j] = P[j, i] = t.distance(tips["usr_" + reps[i]], tips["usr_" + reps[j]])

M, prod, sets = {}, {}, {}
for c in ("0.3", "0.5", "0.7"):
    M[c] = pd.read_csv(f"{base}/09_matrix/gcf_matrix_c{c}_all.tsv", sep="\t", index_col=0).loc[reps]
    rec = pd.read_csv(f"{base}/09_matrix/gcf_records_c{c}.tsv", sep="\t")
    prod[c] = rec.groupby("family")["product"].agg(lambda s: s.value_counts().index[0])
    sets[c] = {f: set(M[c].index[M[c][f] == 1]) for f in M[c].columns}

X = M["0.5"]
cand = [f for f in X.columns if 5 <= X[f].sum() <= 0.4 * n and prod["0.5"].get(f, "") != "terpene-precursor"]
print("Representatives:", n, "| candidate GCFs (5 to %d carriers, no terpene-precursor): %d" % (int(0.4 * n), len(cand)))

NP = 4999
def mpd(idx):
    k = len(idx); return P[np.ix_(idx, idx)][np.triu_indices(k, 1)].mean()
null = {k: np.array([mpd(rng.choice(n, k, replace=False)) for _ in range(NP)]) for k in sorted({int(X[f].sum()) for f in cand})}

def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); r = np.empty(len(p))
    q = p[o] * len(p) / np.arange(1, len(p) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
    r[o] = np.minimum(q, 1); return r

def match(carr, c):
    best = (None, 0.0)
    for f, s in sets[c].items():
        j = len(s & carr) / len(s | carr)
        if j > best[1]: best = (f, j)
    return "%s(J=%.2f)" % best

rows = []
for f in cand:
    idx = np.where(X[f].values == 1)[0]; k = len(idx); o = mpd(idx); nl = null[k]
    carr = set(X.index[X[f] == 1])
    sub = cg.iloc[idx]
    rows.append({"family": f, "product": prod["0.5"].get(f, ""), "carriers": k,
                 "n_genera": sub.genus.nunique(), "SES": (o - nl.mean()) / nl.std(),
                 "p_clustered": (np.sum(nl <= o) + 1) / (NP + 1), "p_dispersed": (np.sum(nl >= o) + 1) / (NP + 1),
                 "c0.3": match(carr, "0.3"), "c0.7": match(carr, "0.7")})
d = pd.DataFrame(rows)
d["fdr_clustered"] = bh(d.p_clustered); d["fdr_dispersed"] = bh(d.p_dispersed)
d.to_csv(f"{base}/09_matrix/gcf_scatter_test_B.tsv", sep="\t", index=False)

print("\nGCFs significantly clustered on the tree (FDR<0.05): %d" % int((d.fdr_clustered < 0.05).sum()))
print("GCFs significantly MORE dispersed than random (FDR<0.05): %d" % int((d.fdr_dispersed < 0.05).sum()))
print("GCFs indistinguishable from random placement: %d" % int(((d.fdr_clustered >= 0.05) & (d.fdr_dispersed >= 0.05)).sum()))
pd.set_option("display.width", 250, "display.max_colwidth", 60)
print("\nHighest SES (most spread across the tree):")
print(d.sort_values("SES", ascending=False).head(12)[["family", "product", "carriers", "n_genera", "SES", "p_clustered", "fdr_clustered", "c0.3", "c0.7"]].round(3).to_string(index=False))
print("\nBeta-lactone GCFs:"); print(d[d["product"].str.contains("betalactone")][["family", "product", "carriers", "n_genera", "SES", "fdr_clustered", "c0.3", "c0.7"]].round(3).to_string(index=False))
