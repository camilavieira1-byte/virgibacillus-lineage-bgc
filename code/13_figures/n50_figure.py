import pandas as pd, numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/13_figures"
g = pd.read_csv(f"{base}/06_antismash/bgc_per_genome.tsv", sep="\t"); reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t")
g["n_complete"] = g.n_bgc - g.n_edge; labs = ["< 50 kb", "50-100 kb", "100-500 kb", "> 500 kb"]
g["bin"] = pd.cut(g.Contig_N50, [0, 50e3, 100e3, 500e3, 1e9], labels=labs)
agg = g.groupby("bin", observed=True).agg(genomes=("accession", "size"), regions=("n_bgc", "sum"), edge=("n_edge", "sum")); agg["edge_pct"] = 100 * agg.edge / agg.regions
mm = g.groupby("bin", observed=True)[["n_bgc", "n_complete"]].agg(["mean", "sem"]); overall = 100 * g.n_edge.sum() / g.n_bgc.sum()
tk = reg.assign(tok=reg["product"].str.split("+")).explode("tok"); cl = tk.groupby("tok").agg(regions=("region", "size"), edge=("contig_edge", "mean")); cl = cl[cl.regions >= 10].sort_values("edge")
fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.9), gridspec_kw={"width_ratios": [1, 1.1, 1.1]})
ramp = ["#7A2F5E", "#A0467A", "#E8999B", "#F4D6D0"]; x = np.arange(len(agg))
ax[0].bar(x, agg.edge_pct, color=ramp, edgecolor="#777777", linewidth=0.4)
for i, (b, r) in enumerate(agg.iterrows()): ax[0].text(i, max(r.edge_pct, overall) + 2.5, "%.0f%%\n%d genomes\n%d regions" % (r.edge_pct, r.genomes, r.regions), ha="center", va="bottom", fontsize=7.5)
ax[0].axhline(overall, color="#3B4252", ls="--", lw=1, zorder=0, label="all genomes: %.0f%%" % overall); ax[0].legend(fontsize=7.5, frameon=False, loc="upper right")
ax[0].set_xticks(x); ax[0].set_xticklabels(agg.index.astype(str), fontsize=8); ax[0].set_ylim(0, 112); ax[0].set_xlabel("Assembly N50"); ax[0].set_ylabel("% of BGC regions on a contig edge"); ax[0].set_title("A. Fragmented assemblies truncate BGCs", fontsize=10)
w = 0.38
ax[1].bar(x - w / 2, mm[("n_bgc", "mean")], w, yerr=mm[("n_bgc", "sem")], color="#D9697A", edgecolor="#777777", linewidth=0.4, capsize=2, label="all regions")
ax[1].bar(x + w / 2, mm[("n_complete", "mean")], w, yerr=mm[("n_complete", "sem")], color="#7A2F5E", edgecolor="#777777", linewidth=0.4, capsize=2, label="complete regions only")
ax[1].set_xticks(x); ax[1].set_xticklabels(agg.index.astype(str), fontsize=8); ax[1].set_xlabel("Assembly N50"); ax[1].set_ylabel("Mean BGC regions per genome (+/- SEM)"); ax[1].set_ylim(0, 12); ax[1].legend(fontsize=7.5, frameon=False, loc="upper right"); ax[1].set_title("B. Neither count is reliable in fragmented assemblies", fontsize=10)
ax[2].barh(range(len(cl)), 100 * cl.edge, color="#7A2F5E", edgecolor="#777777", linewidth=0.4)
for i, (c, r) in enumerate(cl.iterrows()): ax[2].text(100 * r.edge + 1, i, "%.0f%% (n = %d)" % (100 * r.edge, r.regions), va="center", fontsize=7)
ax[2].set_yticks(range(len(cl))); ax[2].set_yticklabels(cl.index, fontsize=8); ax[2].set_xlim(0, 115); ax[2].set_xlabel("% of regions on a contig edge"); ax[2].set_title("C. Edge share by class (classes with >= 10 regions)", fontsize=10)
for a in ax: a.spines[["top", "right"]].set_visible(False)
plt.tight_layout(); plt.savefig(f"{out}/fig_supp_n50_edges.png", dpi=250, bbox_inches="tight"); plt.savefig(f"{out}/fig_supp_n50_edges.pdf", bbox_inches="tight")
print(agg.round(1).to_string()); print("overall edge share: %.1f%% [earlier 27.7%%]" % overall); print(cl.assign(edge_pct=(100 * cl.edge).round(0))[["regions", "edge_pct"]].tail(5).to_string()); print("wrote", f"{out}/fig_supp_n50_edges.png")
