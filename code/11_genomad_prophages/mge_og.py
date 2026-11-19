import pandas as pd, numpy as np, glob
base = "/PATH/TO/virgibacillus_v2"; gd = f"{base}/11_genomad"
mge = pd.read_csv(f"{gd}/mge_regions.tsv", sep="\t"); pro = mge[mge.kind == "prophage"]
pro_by = {k: d[["start", "end"]].values for k, d in pro.groupby(["accession", "contig"])}
b = pd.read_csv(f"{base}/10_pangenome/b_reps.tsv", sep="\t", dtype=str)
res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]
og = pd.read_csv(f"{res}/Orthogroups/Orthogroups.tsv", sep="\t", index_col=0, dtype=str)
cat = pd.read_csv(f"{base}/10_pangenome/pangenome_categories.tsv", sep="\t", index_col=0).category
inpro = {}
for acc in b.accession:
    for line in open(f"{base}/10_pangenome/proteins/{acc}.faa"):
        if line.startswith(">"):
            p = line[1:].strip().split(" # "); gid = p[0].split()[0]; ctg = gid.rsplit("_", 1)[0]
            mid = (int(p[1]) + int(p[2])) / 2; arr = pro_by.get((acc, ctg))
            inpro[gid] = bool(arr is not None and ((arr[:, 0] <= mid) & (arr[:, 1] >= mid)).any())
rows = []
for ogname, r in og.iterrows():
    g = [x.strip() for acc in b.accession if acc in og.columns and isinstance(r[acc], str) for x in r[acc].split(",")]
    g = [x for x in g if x in inpro]
    if g: rows.append((ogname, cat.get(ogname), len(g), sum(inpro[x] for x in g)))
d = pd.DataFrame(rows, columns=["og", "category", "copies", "in_prophage"]).dropna(subset=["category"])
d["frac"] = d.in_prophage / d.copies
t = d.groupby("category").agg(orthogroups=("og", "size"), with_any_copy_in_prophage=("in_prophage", lambda s: int((s > 0).sum())),
    mostly_in_prophage_ge50pct=("frac", lambda s: int((s >= .5).sum())))
t["pct_any"] = (100 * t.with_any_copy_in_prophage / t.orthogroups).round(1); t["pct_mostly"] = (100 * t.mostly_in_prophage_ge50pct / t.orthogroups).round(1)
print(t.reindex(["core", "soft core", "shell", "cloud"]).dropna().to_string())
acc = d[d.category.isin(["shell", "cloud"])]
print("\nShell + cloud orthogroups: %d | with any copy in a prophage: %d (%.1f%%) | mostly in prophages: %d (%.1f%%)" % (
    len(acc), int((acc.in_prophage > 0).sum()), 100 * (acc.in_prophage > 0).mean(), int((acc.frac >= .5).sum()), 100 * (acc.frac >= .5).mean()))
d.to_csv(f"{gd}/orthogroup_prophage.tsv", sep="\t", index=False)
