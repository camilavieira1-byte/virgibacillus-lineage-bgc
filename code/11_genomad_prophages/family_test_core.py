import re, sys
import numpy as np, pandas as pd
from scipy.stats import binomtest

base = "/PATH/TO/virgibacillus_v2"
f = sys.argv[1] if len(sys.argv) > 1 else "bgc_virus_distance_core.tsv"
W = 25000
norm = lambda s: re.sub(r"\.\d+$", "", str(s))

d = pd.read_csv(f, sep="\t")
d["near"] = d["dist"] <= W
p0 = d["near"].mean()
print(f"{f}: baseline near-virus rate (<= {W//1000} kb): {d['near'].sum()}/{len(d)} = {p0:.4f}")

g = (d.dropna(subset=["fam"]).groupby("fam")
       .agg(n=("near", "size"), k=("near", "sum"), prev=("prev", "first"), genomes=("genome", "nunique")))
g = g[g["n"] >= 3].copy()
g["p"] = [binomtest(int(r.k), int(r.n), p0, alternative="greater").pvalue for r in g.itertuples()]
g = g.sort_values("p")
m = len(g)
g["q_BH"] = (g["p"] * m / np.arange(1, m + 1))[::-1].cummin()[::-1].clip(upper=1)
print(f"\nfamilies tested: {m}")
print(g[g["k"] >= 1].head(15).to_string())

vall = pd.read_csv(f"{base}/11_genomad/all_viruses.tsv", sep="\t")
vall["contig"] = vall["seq_name"].str.split("|").str[0].map(norm)
x = vall["seq_name"].str.extract(r"\|provirus_(\d+)_(\d+)")
vall["pstart"] = pd.to_numeric(x[0]).fillna(0); vall["pend"] = pd.to_numeric(x[1]).fillna(np.inf)
rows = []
for r in d[d["near"]].itertuples():
    s = vall[vall["contig"] == r.contig].copy()
    s["gap"] = np.maximum(0, np.maximum(s["pstart"] - r.end, r.start - s["pend"]))
    v = s.sort_values("gap").iloc[0]
    rows.append(dict(genome=r.genome, product=r.product, fam=r.fam, dist=int(r.dist),
                     core_len=int(r.end - r.start), prov_len=int(v["length"]),
                     score=v["virus_score"], hallmarks=v["n_hallmarks"]))
print("\n=== near hits (core coordinates) ===")
print(pd.DataFrame(rows).sort_values(["fam", "dist"]).to_string(index=False))
