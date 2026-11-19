import pandas as pd, re
base = "/PATH/TO/virgibacillus_v2"
e = pd.read_csv(f"{base}/05_ani/skani_edges.tsv", sep="\t")
print("Columns:", e.columns.tolist())
acc = lambda s: re.search(r"(GCA_\d+\.\d)", s).group(1)
e["a"], e["b"] = e["Ref_file"].map(acc), e["Query_file"].map(acc)
e["af"] = e[["Align_fraction_ref", "Align_fraction_query"]].min(axis=1)
genomes = [acc(l) for l in open(f"{base}/05_ani/genome_list.txt")]
meta = pd.read_csv(f"{base}/00_metadata/metadata_curation_template.tsv", sep="\t").set_index("accession")

def clusters(thr, min_af=50):
    p = {g: g for g in genomes}
    def f(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    for r in e[(e.ANI >= thr) & (e.af >= min_af)].itertuples():
        p[f(r.a)] = f(r.b)
    out = {}
    for g in genomes: out.setdefault(f(g), []).append(g)
    return list(out.values())

for thr in (95, 99, 99.9):
    cl = clusters(thr)
    sizes = sorted((len(c) for c in cl), reverse=True)
    print(f"\nANI >= {thr}%: {len(cl)} clusters from {len(genomes)} genomes | singletons: {sum(s==1 for s in sizes)} | largest sizes: {sizes[:10]}")
    if thr == 99.9:
        for c in sorted(cl, key=len, reverse=True)[:6]:
            if len(c) < 3: continue
            sub = meta.reindex(c)
            print(f"  size {len(c)} | countries: {sub.country.fillna('NA').value_counts().to_dict()} | habitat: {sub.habitat_proposed.fillna('NA').value_counts().to_dict()}")
