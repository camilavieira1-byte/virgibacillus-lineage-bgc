import pandas as pd, numpy as np, os
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
rng = np.random.default_rng(5)
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/13_figures"; os.makedirs(out, exist_ok=True)
cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str).fillna("")
cg = cg[cg.in_A == "True"].copy()
t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick"); tips = {x.name: x for x in t.get_terminals()}
cg = cg[cg.accession.map(lambda a: "usr_" + a in tips)].reset_index(drop=True)
reps = cg.accession.tolist(); n = len(reps)
P = np.zeros((n, n))
for i in range(n):
    for j in range(i + 1, n): P[i, j] = P[j, i] = t.distance(tips["usr_" + reps[i]], tips["usr_" + reps[j]])
rec = pd.read_csv(f"{base}/09_matrix/gcf_records_c0.5.tsv", sep="\t")
fam_prod = rec.groupby("family")["product"].agg(lambda s: s.value_counts().index[0])
X = pd.read_csv(f"{base}/09_matrix/gcf_matrix_c0.5_all.tsv", sep="\t", index_col=0).loc[reps]
keep = [f for f in X.columns if X[f].sum() >= 2 and fam_prod.get(f, "") != "terpene-precursor"]
Bm = X[keep].values.astype(int); inter = Bm @ Bm.T; sz = Bm.sum(axis=1); union = sz[:, None] + sz[None, :] - inter
D = 1 - np.where(union > 0, inter / np.maximum(union, 1), 0); np.fill_diagonal(D, 0)
iu = np.triu_indices(n, 1)
def rk(M):
    r = pd.Series(M[iu]).rank().values; R = np.zeros_like(M); R[iu] = r; return R + R.T
RD, RP = rk(D), rk(P); rho = np.corrcoef(RD[iu], RP[iu])[0, 1]; print("pairs of genomes sharing no kept GCF (Jaccard distance = 1): %.0f%% of %d" % (100 * np.mean(D[iu] >= 0.999), len(iu[0])))
hit = sum(np.corrcoef(RD[iu], RP[np.ix_(p, p)][iu])[0, 1] >= rho for p in (rng.permutation(n) for _ in range(4999))); pval = (hit + 1) / 5000
gen = cg.genus.values
def r2(Dm, lab):
    m = len(lab); D2 = Dm ** 2; sst = D2[np.triu_indices(m, 1)].sum() / m; ssw = 0.0
    for g in pd.unique(lab):
        idx = np.where(lab == g)[0]; k = len(idx)
        if k > 1: ssw += D2[np.ix_(idx, idx)][np.triu_indices(k, 1)].sum() / k
    return (sst - ssw) / sst
R2 = r2(D, gen); ng = len(set(gen)); R2a = 1 - (1 - R2) * (n - 1) / (n - ng)
J = np.eye(n) - 1.0 / n; Bc = -0.5 * J @ (D ** 2) @ J; w, v = np.linalg.eigh(Bc); o = np.argsort(w)[::-1]; w, v = w[o], v[:, o]
co = v[:, :2] * np.sqrt(np.maximum(w[:2], 0)); var = 100 * w[:2] / w[w > 0].sum()
sp = f"{base}/09_matrix/gcf_scatter_test.tsv"
if not os.path.exists(sp): sp = f"{base}/09_matrix/gcf_scatter_test_Aprime.tsv"
S = pd.read_csv(sp, sep="\t"); sig = S.fdr_clustered < 0.05; disp = S.fdr_dispersed < 0.05; rnd = ~sig & ~disp
print("representatives: %d | GCFs kept for distances: %d | Mantel r = %.3f (p = %.4f) | genus PERMANOVA R2 = %.3f (adjusted %.3f, %d genera)" % (n, len(keep), rho, pval, R2, R2a, ng))
print("scatter test file %s: %d GCFs | clustered %d | more dispersed %d | indistinguishable from random %d" % (os.path.basename(sp), len(S), sig.sum(), disp.sum(), rnd.sum()))
cmap = LinearSegmentedColormap.from_list("pl", ["#F6E8EF", "#C07BA0", "#7A2F5E"])
fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.4), gridspec_kw={"width_ratios": [1.15, 1.1, 1.1]})
hb = ax[0].hexbin(np.maximum(P[iu], 1e-4), D[iu], gridsize=(40, 22), xscale="log", bins="log", cmap=cmap, mincnt=1, linewidths=0.1); plt.colorbar(hb, ax=ax[0], pad=0.01, shrink=0.7, label="pairs of genomes (log scale)")
xq = np.logspace(-3.6, np.log10(P[iu].max()), 16); xq[0] = 0.0; xm, ym = [], []
for lo, hi in zip(xq[:-1], xq[1:]):
    m = (P[iu] >= lo) & (P[iu] < hi)
    if m.sum() >= 15: xm.append(np.median(np.maximum(P[iu][m], 1e-4))); ym.append(np.median(D[iu][m]))
ax[0].plot(xm, ym, color="#3B4252", lw=2, label="median per distance bin"); ax[0].legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.16))
sub_a = "Mantel r = %.2f, p %s | PERMANOVA by genus: R2 = %.2f (adjusted %.2f)" % (rho, "< 0.001" if pval <= 0.001 else "= %.3f" % pval, R2, R2a)
ax[0].set_xlabel("Distance between genomes on the tree"); ax[0].set_ylabel("Difference in GCF content (Jaccard)"); ax[0].set_title("A. GCF content follows the phylogeny\n" + sub_a, fontsize=8.5)
SP = {"Virgibacillus salexigens": "#E8999B", "Virgibacillus salarius": "#F4D6D0", "Virgibacillus pantothenticus": "#F7C99A", "Virgibacillus proomii": "#A0467A", "Virgibacillus chiguensis": "#C07BA0", "Virgibacillus dokdonensis": "#EBD0DE"}
groups = {}
for i in range(n):
    if cg.in_B[i] == "True": key = (SP.get(cg.species[i], "#CCCCCC"), cg.species[i].replace("Virgibacillus ", "V. "))
    elif cg.genus[i].startswith("Virgibacillus"): key = ("#8FA1B3", "other Virgibacillus_X")
    else: key = ("#3B4252", "other genera")
    groups.setdefault(key, []).append(i)
for (col, lab), idx in groups.items(): ax[1].scatter(co[idx, 0], co[idx, 1], s=46, c=col, edgecolors="#444444", linewidths=0.5, label="%s (%d)" % (lab, len(idx)), zorder=3)
ax[1].set_xlabel("PCoA 1 (%.0f%%)" % var[0]); ax[1].set_ylabel("PCoA 2 (%.0f%%)" % var[1]); ax[1].legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2); ax[1].set_title("B. Ordination of GCF content", fontsize=10)
ax[2].axhline(0, color="#999999", lw=0.8, ls="--", zorder=1)
ax[2].scatter(S.carriers[sig], S.SES[sig], s=40, c="#7A2F5E", edgecolors="white", linewidths=0.4, label="carriers clustered on the tree (FDR < 0.05): %d" % sig.sum(), zorder=3)
ax[2].scatter(S.carriers[~sig], S.SES[~sig], s=56, c="#F4B183", edgecolors="#7A2F5E", linewidths=1.0, label="not different from random placement: %d" % (~sig).sum(), zorder=4)
OFF = {"FAM_00100": (8, 6), "FAM_00032": (8, -14), "FAM_00163": (8, 8), "FAM_00107": (8, -14)}
for _, r in S[~sig].iterrows(): ax[2].annotate("%s (%s)" % (r.family, r["product"]), (r.carriers, r.SES), xytext=OFF.get(r.family, (8, 6)), textcoords="offset points", fontsize=6.5, arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5))
ax[2].set_ylim(S.SES.min() - 0.8, S.SES.max() + 1.3)
ax[2].set_xlabel("Genomes carrying the GCF"); ax[2].set_ylabel("Standardised spread of carriers (SES)\nnegative = carriers are close relatives"); ax[2].legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14)); ax[2].set_title("C. Do a GCF's carriers sit together on the tree? (%d GCFs)" % len(S), fontsize=10)
for a in ax: a.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig(f"{out}/fig_lineage_structure.png", dpi=250, bbox_inches="tight"); plt.savefig(f"{out}/fig_lineage_structure.pdf", bbox_inches="tight"); print("wrote", f"{out}/fig_lineage_structure.png")
