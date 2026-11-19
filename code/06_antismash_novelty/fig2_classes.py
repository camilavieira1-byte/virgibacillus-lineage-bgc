import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
base = "/PATH/TO/virgibacillus_v2/06_antismash"
b = pd.read_csv(f"{base}/bgc_regions.tsv", sep="\t")
g = pd.read_csv(f"{base}/bgc_per_genome.tsv", sep="\t")
n_genomes = len(g)
e = b.assign(cls=b["product"].str.split("+")).explode("cls")
t = e.groupby("cls").agg(n_regions=("region", "size"), n_genomes=("accession", "nunique"), n_edge=("contig_edge", "sum"))
t["n_complete"] = t.n_regions - t.n_edge
t["pct_genomes"] = (100 * t.n_genomes / n_genomes).round(1)
t = t.sort_values("n_regions", ascending=False)
t.to_csv(f"{base}/class_summary.tsv", sep="\t")
print("Genomes in analysis:", n_genomes, "| with >=1 BGC:", int((g.n_bgc > 0).sum())); print(t.head(20).to_string())

top = t.head(20).iloc[::-1]
fig, ax = plt.subplots(1, 2, figsize=(13, 6.5), sharey=True)
ax[0].barh(top.index, top.n_complete, label="complete region")
ax[0].barh(top.index, top.n_edge, left=top.n_complete, label="on contig edge", color="#c9c9c9")
ax[0].set_xlabel("BGC regions"); ax[0].legend(frameon=False); ax[0].set_title("A. Regions per class")
ax[1].barh(top.index, top.pct_genomes, color="#4c78a8")
ax[1].set_xlabel("% of genomes carrying the class"); ax[1].set_title("B. Distribution across genomes")
fig.text(0.01, 0.01, "Hybrid regions are counted once in each of their classes.", fontsize=8)
plt.tight_layout(); plt.savefig(f"{base}/fig2_class_frequency.png", dpi=300)
print("Wrote fig2_class_frequency.png")
