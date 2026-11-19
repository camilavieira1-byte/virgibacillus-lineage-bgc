import pandas as pd, numpy as np
from Bio import Phylo
rng = np.random.default_rng(1)
base = "/PATH/TO/virgibacillus_v2"
cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str)
cg = cg[cg.in_B == "True"]
t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick")
tips = {x.name: x for x in t.get_terminals()}
cg = cg[cg.accession.map(lambda a: "usr_" + a in tips)].reset_index(drop=True)
reps = cg.accession.tolist(); n = len(reps)
P = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        P[i, j] = P[j, i] = t.distance(tips["usr_" + reps[i]], tips["usr_" + reps[j]])
rec = pd.read_csv(f"{base}/09_matrix/gcf_records_c0.5.tsv", sep="\t")
fam_prod = rec.groupby("family")["product"].agg(lambda s: s.value_counts().index[0])
X = pd.read_csv(f"{base}/09_matrix/gcf_matrix_c0.5_all.tsv", sep="\t", index_col=0).loc[reps]
X = X.loc[:, X.sum() > 0]
print("Clone-group representatives in Virgibacillus (B):", n, "| GCFs present:", X.shape[1])
print("GCFs by number of representatives carrying them:", pd.cut(X.sum(), [0, 1, 2, 4, 9, 100], labels=["1", "2", "3-4", "5-9", ">=10"]).value_counts().sort_index().to_dict())
print("GCFs per representative: mean %.1f, range %d-%d" % (X.sum(1).mean(), X.sum(1).min(), X.sum(1).max()))
keep = [f for f in X.columns if X[f].sum() >= 2 and fam_prod.get(f, "") != "terpene-precursor"]
print("Kept for distances (>=2 carriers, no terpene-precursor):", len(keep))
Bm = X[keep].values.astype(int)
inter = Bm @ Bm.T; sz = Bm.sum(1); union = sz[:, None] + sz[None, :] - inter
D = 1 - np.where(union > 0, inter / np.maximum(union, 1), 0); np.fill_diagonal(D, 0)
iu = np.triu_indices(n, 1)
def rk(A):
    r = pd.Series(A[iu]).rank().values; R = np.zeros_like(A); R[iu] = r; return R + R.T
RD, RP = rk(D), rk(P)
obs = np.corrcoef(RD[iu], RP[iu])[0, 1]; c = 0
for _ in range(4999):
    p = rng.permutation(n)
    if np.corrcoef(RD[iu], RP[np.ix_(p, p)][iu])[0, 1] >= obs: c += 1
print("Mantel (Spearman), GCF Jaccard distance vs tree distance: r = %.3f, p = %.4f" % (obs, (c + 1) / 5000))
