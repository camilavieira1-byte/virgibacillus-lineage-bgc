import pandas as pd
base = "/PATH/TO/virgibacillus_v2"
b = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t")
g = pd.read_csv(f"{base}/06_antismash/bgc_per_genome.tsv", sep="\t")
g["n_complete"] = g.n_bgc - g.n_edge

print("== 1) Edge fraction by BGC class (classes with >=10 regions)")
e = b.assign(cls=b["product"].str.split("+")).explode("cls")
t = e.groupby("cls").agg(regions=("region", "size"), on_edge=("contig_edge", "sum"))
t["edge_pct"] = (100 * t.on_edge / t.regions).round(1)
print(t[t.regions >= 10].sort_values("edge_pct", ascending=False).to_string())
print("\nMedian region length (bp), complete vs edge:")
print(b.groupby("contig_edge").length_bp.median().to_string())

print("\n== 2) Genomes and BGCs retained at each N50 threshold")
rows = []
for thr in (0, 50e3, 100e3, 200e3, 500e3):
    x = g[g.Contig_N50 >= thr]
    rows.append({"min_N50_kb": int(thr / 1e3), "genomes": len(x), "bgc_all": int(x.n_bgc.sum()),
                 "bgc_complete": int(x.n_complete.sum()),
                 "mean_complete_per_genome": round(x.n_complete.mean(), 2)})
print(pd.DataFrame(rows).to_string(index=False))

print("\n== 2b) 'Edge BGCs removed, all genomes kept': mean complete BGCs per genome by N50")
g["N50_bin"] = pd.cut(g.Contig_N50, [0, 50e3, 100e3, 500e3, 1e9], labels=["<50k", "50-100k", "100-500k", ">500k"])
print(g.groupby("N50_bin", observed=True).agg(genomes=("accession", "size"),
      mean_all=("n_bgc", "mean"), mean_complete=("n_complete", "mean")).round(2).to_string())

print("\n== 3) Is assembly quality related to habitat? (habitat_proposed is a draft, not curated)")
m = pd.read_csv(f"{base}/00_metadata/metadata_curation_template.tsv", sep="\t", usecols=["accession", "habitat_proposed"])
h = g.merge(m, on="accession", how="left")
h["habitat_proposed"] = h.habitat_proposed.fillna("").replace("", "unresolved")
h["low_N50"] = h.Contig_N50 < 100e3
print(h.groupby("habitat_proposed").agg(genomes=("accession", "size"),
      median_N50_kb=("Contig_N50", lambda s: round(s.median() / 1e3)),
      pct_below_100kb=("low_N50", lambda s: round(100 * s.mean()))).to_string())
