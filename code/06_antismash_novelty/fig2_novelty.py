import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = "/PATH/TO/virgibacillus_v2/06_antismash"
m = pd.read_csv(f"{base}/mibig_best_per_region.tsv", sep="\t")
n_g = pd.read_csv(f"{base}/bgc_per_genome.tsv", sep="\t").shape[0]

top = m["product"].value_counts().head(14).index
n_other = m[~m["product"].isin(top)]["product"].nunique()
m["cls"] = m["product"].where(m["product"].isin(top), "other (%d products)" % n_other)
order = ["no hit", "low (1-49%)", "medium (50-79%)", "high (>=80%)"]

ct = pd.crosstab(m.cls, m.bin).reindex(columns=order).fillna(0)
ct = ct.loc[ct.sum(axis=1).sort_values().index]
pct = (m.groupby("cls").accession.nunique() / n_g * 100).reindex(ct.index)

# Updated iTOL Palette for the 4 similarity bins
cols = ["#e6e6e6", "#f4b183", "#d9697a", "#7a2f5e"]

fig, ax = plt.subplots(1, 2, figsize=(13, 6.5), sharey=True, gridspec_kw={"width_ratios": [2, 1]})
left = 0

for c, col in zip(order, cols):
    ax[0].barh(ct.index, ct[c], left=left, color=col, label=c); left = left + ct[c].values
ax[0].set_xlabel("BGC regions"); ax[0].set_title("A. Regions per product and similarity to MIBiG"); ax[0].legend(frameon=False, title="Best MIBiG similarity")

# Updated iTOL Palette for the percentage bars (Plum/Purple)
ax[1].barh(ct.index, pct, color="#975283"); ax[1].set_xlabel("%% of genomes (n = %d) with the product" % n_g); ax[1].set_title("B. Distribution across genomes")

fig.text(0.01, 0.01, "Similarity = % of genes of the closest MIBiG cluster with a hit (KnownClusterBlast, antiSMASH 8.0.4); hybrids shown under their product label.", fontsize=8)
ax[0].set_xlim(0, ct.sum(axis=1).max() * 1.05)
plt.tight_layout(); plt.savefig(f"{base}/fig2_novelty.png", dpi=300)
print("Wrote fig2_novelty.png")
