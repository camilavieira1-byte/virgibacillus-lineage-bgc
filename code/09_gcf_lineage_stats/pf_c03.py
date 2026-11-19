import pandas as pd, numpy as np
from Bio import Phylo
rng = np.random.default_rng(1)
base = "/PATH/TO/virgibacillus_v2"

cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str)
cg = cg[cg.in_A == "True"].copy()
rec = pd.read_csv(f"{base}/09_matrix/gcf_records_c0.3.tsv", sep="\t")
M = pd.read_csv(f"{base}/09_matrix/gcf_matrix_c0.3_all.tsv", sep="\t", index_col=0)

t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick")
tips = {x.name: x for x in t.get_terminals()}
cg = cg[cg.accession.map(lambda a: "usr_" + a in tips)].reset_index(drop=True)
reps = cg.accession.tolist(); n = len(reps)
print("Clone-group representatives in clade A' and in the tree:", n)

P = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n):
        P[i, j] = P[j, i] = t.distance(tips["usr_" + reps[i]], tips["usr_" + reps[j]])

X = M.loc[reps]
fam_prod = rec.groupby("family")["product"].agg(lambda s: s.value_counts().index[0])
prev = X.mean().sort_values(ascending=False)
print("\nMost widespread GCFs (share of representatives) and their antiSMASH product:")
print(pd.DataFrame({"prevalence": prev.head(6).round(2), "product": [fam_prod.get(f, "") for f in prev.head(6).index]}).to_string())
tp = [f for f in X.columns if fam_prod.get(f, "") == "terpene-precursor"]
keep = [f for f in X.columns if X[f].sum() >= 2 and f not in tp]
print("\nGCFs total: %d | terpene-precursor families removed: %d | private (1 genome) removed: %d | kept: %d" % (
    X.shape[1], len(tp), int(((X.sum() < 2) & ~X.columns.isin(tp)).sum()), len(keep)))
B = X[keep].values.astype(int)
print("Representatives with no kept GCF:", int((B.sum(1) == 0).sum()))
inter = B @ B.T; sz = B.sum(1); union = sz[:, None] + sz[None, :] - inter
D = 1 - np.where(union > 0, inter / np.maximum(union, 1), 0); np.fill_diagonal(D, 0)

iu = np.triu_indices(n, 1)
def rankmat(A):
    r = pd.Series(A[iu]).rank().values; R = np.zeros_like(A); R[iu] = r; return R + R.T
def mantel(A, Bm, nperm=4999):
    RA, RB = rankmat(A), rankmat(Bm)
    obs = np.corrcoef(RA[iu], RB[iu])[0, 1]; c = 0
    for _ in range(nperm):
        p = rng.permutation(n)
        if np.corrcoef(RA[iu], RB[np.ix_(p, p)][iu])[0, 1] >= obs: c += 1
    return obs, (c + 1) / (nperm + 1)
r, p = mantel(D, P)
print("\n== Mantel (Spearman): GCF Jaccard distance vs patristic distance: r = %.3f, p = %.4f" % (r, p))

def pseudoF(Dm, lab):
    m = len(lab); D2 = Dm ** 2
    sst = D2[np.triu_indices(m, 1)].sum() / m; ssw = 0; labs = pd.unique(lab)
    for g in labs:
        idx = np.where(lab == g)[0]; k = len(idx)
        if k > 1: ssw += D2[np.ix_(idx, idx)][np.triu_indices(k, 1)].sum() / k
    kk = len(labs)
    if kk < 2 or m <= kk or ssw == 0: return np.nan, np.nan
    return ((sst - ssw) / (kk - 1)) / (ssw / (m - kk)), (sst - ssw) / sst
def permanova(Dm, lab, strata=None, nperm=999):
    lab = np.asarray(lab); F0, R2 = pseudoF(Dm, lab); c = 0
    for _ in range(nperm):
        if strata is None: pl = rng.permutation(lab)
        else:
            pl = lab.copy()
            for s in pd.unique(strata):
                idx = np.where(strata == s)[0]; pl[idx] = rng.permutation(lab[idx])
        if pseudoF(Dm, pl)[0] >= F0: c += 1
    return F0, R2, (c + 1) / (nperm + 1)

gen = cg.genus.values
F, R2, pv = permanova(D, gen)
print("== PERMANOVA by GTDB genus (%d genera, %d reps): R2 = %.3f, F = %.2f, p = %.4f" % (len(set(gen)), n, R2, F, pv))
for var in ("habitat", "region"):
    lab = cg[var].values; cnt = pd.Series(lab).value_counts()
    ok = [k for k in cnt.index if k != "unresolved" and cnt[k] >= 5]
    sel = np.isin(lab, ok); sub = D[np.ix_(sel, sel)]
    print("\n== %s: groups used (>=5 clone groups, resolved only): %s | n = %d" % (var, {k: int(cnt[k]) for k in ok}, int(sel.sum())))
    F1, R1, p1 = permanova(sub, lab[sel])
    F2, R2b, p2 = permanova(sub, lab[sel], strata=gen[sel])
    print("   plain PERMANOVA:             R2 = %.3f, p = %.4f" % (R1, p1))
    print("   permuted within genus:       R2 = %.3f, p = %.4f   (lineage held fixed)" % (R2b, p2))
    print("   genera with >1 %s label: %d" % (var, int((pd.DataFrame({"g": gen[sel], "l": lab[sel]}).groupby("g").l.nunique() > 1).sum())))

mean_all = P[iu].mean(); rows = []
for f in keep:
    idx = np.where(X[f].values == 1)[0]; k = len(idx)
    if k >= 3:
        rows.append((f, k, P[np.ix_(idx, idx)][np.triu_indices(k, 1)].mean() / mean_all, fam_prod.get(f, "")))
d = pd.DataFrame(rows, columns=["family", "carriers", "dist_ratio", "product"])
d["pattern"] = pd.cut(d.dist_ratio, [0, 0.5, 0.8, 10], labels=["clade-restricted (<0.5)", "intermediate", "scattered (>0.8)"])
print("\n== GCFs with >=3 carriers, by how far apart the carriers sit on the tree (carrier distance / mean distance)")
print(d.pattern.value_counts().to_string())
print("\nMost scattered GCFs (>=5 carriers):")
print(d[d.carriers >= 5].sort_values("dist_ratio", ascending=False).head(10).round(2).to_string(index=False))
d.to_csv(f"{base}/09_matrix/gcf_phylo_pattern.tsv", sep="\t", index=False)
