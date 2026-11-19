import pandas as pd, numpy as np, itertools, os
from scipy.stats import spearmanr, mannwhitneyu
from Bio import Phylo
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = "/PATH/TO/virgibacillus_v2"; ph = f"{base}/12_phages"; out = f"{base}/13_figures"
P = pd.read_csv(f"{ph}/prophages.tsv", sep="\t").merge(pd.read_csv(f"{ph}/prophage_idmap.tsv", sep="\t"), on="pid")
C = pd.read_csv(f"{ph}/clusters_phabox_aai.tsv", sep="\t", header=None, names=["pid", "c"])
reps = set(pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str).accession)
A = pd.read_csv(f"{ph}/phabox_votu_aai/midfolder/self_results_aai.tsv", sep="\t"); A = A[A["query"] != A["target"]]
aai = {}
for r in A.itertuples(index=False):
    k = frozenset((r.query, r.target)); aai[k] = max(aai.get(k, 0), r.aai)
t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick"); tips = {x.name: x for x in t.get_terminals()}
spp = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna("").set_index("accession").species.to_dict()
cache = {}
def hd(a, b):
    k = frozenset((a, b))
    if k not in cache: cache[k] = t.distance(tips["usr_" + a], tips["usr_" + b])
    return cache[k]
D = P[P.accession.isin(reps)].merge(C, on="pid"); rows = []
for c, d in D.groupby("c"):
    d = d.drop_duplicates("short")
    if len(d) < 3 or d.accession.nunique() < 3: continue
    hosts, ids, med = d.accession.tolist(), d.short.tolist(), d.length.median() / 1e3
    for i, j in itertools.combinations(range(len(ids)), 2):
        k = frozenset((ids[i], ids[j]))
        if k in aai and hosts[i] != hosts[j]:
            si, sj = spp.get(hosts[i], ""), spp.get(hosts[j], ""); rows.append((c, aai[k], hd(hosts[i], hosts[j]), bool(si) and si == sj, med))
df = pd.DataFrame(rows, columns=["group", "aai", "host_dist", "same_species", "median_kb"])
df = df[df.groupby("group")["group"].transform("size") >= 3]
sizes = df.groupby("group").size().sort_values(ascending=False); big = sizes.index[0]; top = [g for g in sizes.index if sizes[g] >= 6][:5]
rho_all = spearmanr(df.host_dist, df.aai)[0]; nb = df[df.group != big]; rho_nb = spearmanr(nb.host_dist, nb.aai)[0]
print("groups: %d | pairs: %d | pooled rho = %.3f (earlier run: -0.582) | without %s: rho = %.3f (earlier: -0.408)" % (df.group.nunique(), len(df), rho_all, big, rho_nb))
PAL = ["#7A2F5E", "#D9697A", "#F4B183", "#8C6BB1", "#4C78A8"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.8), gridspec_kw={"width_ratios": [2.2, 1]})
rest = df[~df.group.isin(top)]
a1.scatter(np.maximum(rest.host_dist, 1e-5), rest.aai, s=22, c="#BBBBBB", alpha=0.8, label="other groups (%d pairs)" % len(rest), zorder=2)
for g, col in zip(top, PAL):
    d = df[df.group == g]; rho = spearmanr(d.host_dist, d.aai)[0]
    a1.scatter(np.maximum(d.host_dist, 1e-5), d.aai, s=26, c=col, edgecolors="white", linewidths=0.3, alpha=0.9, zorder=3, label="%s: %d pairs, median %.0f kb, rho = %.2f" % (g, len(d), d.median_kb.iloc[0], rho))
a1.set_xscale("log"); a1.set_xlabel("Host divergence (patristic distance, log scale)"); a1.set_ylabel("Prophage AAI (%)")
a1.text(0.02, 0.03, "All pairs: rho = %.2f (n = %d)\nWithout %s: rho = %.2f" % (rho_all, len(df), big, rho_nb), transform=a1.transAxes, fontsize=8)
a1.legend(fontsize=7, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2); a1.set_title("A. Pairs of prophages within PhaBOX2 AAI groups", fontsize=10)
same, diff = df[df.same_species].aai.values, df[~df.same_species].aai.values
bp = a2.boxplot([same, diff], widths=0.5, patch_artist=True, showfliers=False, medianprops=dict(color="#000000"))
for patch, col in zip(bp["boxes"], ["#EBD0DE", "#C07BA0"]): patch.set_facecolor(col)
rng = np.random.default_rng(0)
for i, v in enumerate([same, diff], 1): a2.scatter(i + rng.uniform(-0.15, 0.15, len(v)), v, s=8, c="#444444", alpha=0.5, zorder=3)
a2.set_xticks([1, 2]); a2.set_xticklabels(["same host species\n(n = %d)" % len(same), "different host species\n(n = %d)" % len(diff)], fontsize=8); a2.set_ylabel("Prophage AAI (%)")
a2.set_title("B. Median %.1f vs %.1f (Mann-Whitney p = %.1e)" % (np.median(same), np.median(diff), mannwhitneyu(same, diff)[1]), fontsize=9)
for a in (a1, a2): a.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig(f"{out}/fig_codivergence.png", dpi=250, bbox_inches="tight"); plt.savefig(f"{out}/fig_codivergence.pdf", bbox_inches="tight")
print("wrote", f"{out}/fig_codivergence.png")
