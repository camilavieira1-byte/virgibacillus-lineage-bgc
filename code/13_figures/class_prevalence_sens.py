import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact, mannwhitneyu
rng = np.random.default_rng(2)
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/13_figures"
cg = pd.read_csv(f"{base}/09_matrix/" + os.environ.get("TABLE", "analysis_table_clonegroups.tsv"), sep="\t", dtype=str).fillna("")
cg = cg[cg.in_A == "True"].copy(); cg["n50"] = pd.to_numeric(cg.Contig_N50, errors="coerce")
cg = cg[cg.n50 >= float(os.environ.get("MIN_N50", "0"))]
G = ["Virgibacillus (GTDB)", "other Virgibacillus_X", "other genera"]; COL = dict(zip(G, ["#7A2F5E", "#E8999B", "#8FA1B3"]))
SHORT = {G[0]: "Virgibacillus\n(GTDB)", G[1]: "other\nVirgibacillus_X", G[2]: "other\ngenera"}
cg["group"] = np.where(cg.in_B == "True", G[0], np.where(cg.genus.str.startswith("Virgibacillus"), G[1], G[2]))
reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t"); reg = reg[reg.accession.isin(cg.accession)]
tok = reg.assign(tok=reg["product"].str.split("+")).explode("tok")
pres = pd.crosstab(tok.accession, tok.tok).clip(upper=1).reindex(cg.accession).fillna(0)
classes = [c for c in pres.columns if c != "terpene-precursor" and pres[c].sum() >= 8]
isB = (cg.group == G[0]).values; rows = []
for c in classes:
    a = pres[c].values.astype(bool); tab = [[a[isB].sum(), (~a[isB]).sum()], [a[~isB].sum(), (~a[~isB]).sum()]]
    rows.append({"class": c, **{g: 100 * a[(cg.group == g).values].mean() for g in G}, "others_all": 100 * a[~isB].mean(), "p": fisher_exact(tab)[1]})
R = pd.DataFrame(rows); o = np.argsort(R.p.values); q = R.p.values[o] * len(R) / np.arange(1, len(R) + 1); q = np.minimum.accumulate(q[::-1])[::-1]
R["fdr"] = np.nan; R.loc[R.index[o], "fdr"] = np.minimum(q, 1)
R["diff"] = R[G[0]] - R.others_all; R = R.sort_values("diff", ascending=True).reset_index(drop=True); print(R.round(4).to_string())
nop = reg[reg["product"] != "terpene-precursor"]
cg["n_nopre"] = cg.accession.map(nop.groupby("accession").size()).fillna(0); cg["n_nopre_complete"] = cg.accession.map(nop[~nop.contig_edge].groupby("accession").size()).fillna(0)
def box(ax, d, col, title):
    data = [d[d.group == g][col].values for g in G]
    bp = ax.boxplot(data, widths=0.55, patch_artist=True, showfliers=False, medianprops=dict(color="#000000"))
    for patch, g in zip(bp["boxes"], G): patch.set_facecolor(COL[g]); patch.set_alpha(0.55)
    for i, v in enumerate(data, 1): ax.scatter(i + rng.uniform(-0.16, 0.16, len(v)), v + rng.uniform(-0.12, 0.12, len(v)), s=9, c="#444444", alpha=0.5, zorder=3)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["%s\nn = %d" % (SHORT[g], len(v)) for g, v in zip(G, data)], fontsize=7.5)
    pv = mannwhitneyu(data[0], np.concatenate(data[1:]))[1]; ax.set_title("%s\nVirgibacillus vs the rest: Mann-Whitney p = %.1e" % (title, pv), fontsize=9); return pv
fig, ax = plt.subplots(1, 3, figsize=(16, 5.8), gridspec_kw={"width_ratios": [1.5, 1, 1]})
for i, r in R.iterrows():
    lo, hi = min(r[g] for g in G), max(r[g] for g in G); ax[0].plot([lo, hi], [i, i], color="#CCCCCC", lw=1.5, zorder=1)
for g, mk in zip(G, ["o", "o", "o"]): ax[0].scatter(R[g], range(len(R)), s=48, c=COL[g], edgecolors="#444444", linewidths=0.4, label=g, zorder=3)
ax[0].set_yticks(range(len(R))); ax[0].set_yticklabels([c + (" *" if f < 0.05 else "") for c, f in zip(R["class"], R.fdr)], fontsize=8)
ax[0].set_xlabel("% of genomes carrying the class"); ax[0].set_title("A. BGC classes in Virgibacillus and its relatives\n(* Fisher test, Virgibacillus vs the rest, FDR < 0.05)", fontsize=9); ax[0].legend(fontsize=7, frameon=False, loc="lower right")
p1 = box(ax[1], cg, "n_nopre", "B. BGC regions per genome\n(terpene-precursor excluded; all representatives)"); ax[1].set_ylabel("Regions per genome")
sub = cg[cg.n50 >= 500e3]; p2 = box(ax[2], sub, "n_nopre_complete", "C. Complete regions per genome\n(contiguous assemblies, N50 >= 500 kb)")
for a in ax: a.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); T = os.environ.get("TAG", "main"); plt.savefig(f"{out}/fig_class_prevalence_{T}.png", dpi=250, bbox_inches="tight"); plt.savefig(f"{out}/fig_class_prevalence_{T}.pdf", bbox_inches="tight")
print("representatives: %d (%s)" % (len(cg), cg.group.value_counts().to_dict())); print("classes in >= 8 genomes: %d | significant (FDR < 0.05): %s" % (len(R), list(R[R.fdr < 0.05]["class"])))
print("p (regions, all): %.2e [earlier 4.8e-07] | p (complete, N50 >= 500 kb, n = %d): %.2e [earlier 5.4e-03]" % (p1, len(sub), p2)); print("wrote", f"{out}/fig_class_prevalence.png")

print(cg.groupby("group")[["n_nopre", "n_nopre_complete"]].agg(["count", "median", "mean"]).round(2).to_string())
print(sub.groupby("group")["n_nopre_complete"].agg(["count", "median", "mean"]).round(2).to_string())
