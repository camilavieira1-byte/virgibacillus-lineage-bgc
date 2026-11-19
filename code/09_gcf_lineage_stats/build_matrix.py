import pandas as pd, glob, os
from math import comb
base = "/PATH/TO/virgibacillus_v2"
out = f"{base}/09_matrix"
run = glob.glob(f"{base}/07_bigscape/out_glocal/output_files/virgi_glocal_*_full.network")[0].replace("_full.network", "")

staged = sorted(os.path.basename(p)[:-4] for p in glob.glob(f"{base}/07_bigscape/input_fixed/*.gbk"))
genomes = pd.read_csv(f"{base}/06_antismash/bgc_per_genome.tsv", sep="\t").accession.tolist()
reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t")
reg["gbk"] = reg.accession + "_" + reg.file.str.replace(".gbk", "", regex=False)
edge = dict(zip(reg.gbk, reg.contig_edge)); prod = dict(zip(reg.gbk, reg["product"]))
genome_of = lambda g: "_".join(g.split("_")[:2])
print("Staged regions:", len(staged), "| genomes:", len(genomes), "| staged regions without antiSMASH match:", sum(1 for g in staged if g not in edge))

def ari(a, b):
    ct = pd.crosstab(a, b).values
    s = sum(comb(int(x), 2) for x in ct.ravel())
    sa = sum(comb(int(x), 2) for x in ct.sum(1)); sb = sum(comb(int(x), 2) for x in ct.sum(0))
    n = comb(int(ct.sum()), 2); e = sa * sb / n; m = (sa + sb) / 2
    return (s - e) / (m - e)

maps = {}
for c in ("0.3", "0.5", "0.7"):
    files = glob.glob(f"{run}_c{c}/*/*_clustering_c{c}.tsv")
    df = pd.concat([pd.read_csv(f, sep="\t").assign(bin=os.path.basename(os.path.dirname(f))) for f in files])
    dup = df.GBK.duplicated().sum()
    df = df.drop_duplicates("GBK")
    rec = pd.DataFrame({"gbk": staged})
    rec["genome"] = rec.gbk.map(genome_of)
    rec["family"] = rec.gbk.map(dict(zip(df.GBK, df.Family)))
    rec["bin"] = rec.gbk.map(dict(zip(df.GBK, df.bin)))
    rec["singleton"] = rec.family.isna()
    rec.loc[rec.singleton, "family"] = ["SING_%04d" % i for i in range(1, int(rec.singleton.sum()) + 1)]
    rec.loc[rec.singleton, "bin"] = "singleton"
    rec["edge"] = rec.gbk.map(edge).fillna(False).astype(bool)
    rec["product"] = rec.gbk.map(prod)
    rec.to_csv(f"{out}/gcf_records_c{c}.tsv", sep="\t", index=False)
    maps[c] = dict(zip(df.GBK, df.Family))

    print(f"\n===== cutoff {c} =====")
    print("regions in families: %d | singletons added back: %d | duplicated in several bins: %d" % ((~rec.singleton).sum(), rec.singleton.sum(), dup))
    print("Singleton regions on contig edge: %.1f%% | non-singleton regions on contig edge: %.1f%%" % (
        100 * rec[rec.singleton].edge.mean(), 100 * rec[~rec.singleton].edge.mean()))
    for tag, sub in (("all", rec), ("complete", rec[~rec.edge])):
        M = pd.crosstab(sub.genome, sub.family).clip(upper=1).reindex(genomes).fillna(0).astype(int)
        M.to_csv(f"{out}/gcf_matrix_c{c}_{tag}.tsv", sep="\t")
        k = M.sum(1); g = M.sum(0)
        multi = (~rec.singleton & ~rec.edge if tag == "complete" else ~rec.singleton)
        print("[%s] GCFs: %d | per genome: mean %.1f, median %.0f, range %d-%d | genomes with 0 GCFs: %d" % (
            tag, M.shape[1], k.mean(), k.median(), k.min(), k.max(), int((k == 0).sum())))
        if tag == "all":
            bins = pd.cut(g, [0, 1, 2, 4, 9, 19, 1000], labels=["1", "2", "3-4", "5-9", "10-19", ">=20"])
            print("   GCFs by number of genomes carrying them:", bins.value_counts().reindex(["1", "2", "3-4", "5-9", "10-19", ">=20"]).to_dict())
    print("Regions by bin:", rec.bin.value_counts().to_dict())

print("\n===== stability of families between cutoffs (adjusted Rand index) =====")
common = sorted(set(maps["0.3"]) & set(maps["0.5"]) & set(maps["0.7"]))
for a, b in (("0.3", "0.5"), ("0.5", "0.7"), ("0.3", "0.7")):
    print("c%s vs c%s on %d shared regions: ARI = %.3f" % (a, b, len(common), ari(pd.Series([maps[a][x] for x in common]), pd.Series([maps[b][x] for x in common]))))

ann = pd.read_csv(f"{run}_c0.5/record_annotations.tsv", sep="\t")
print("\n===== BiG-SCAPE Category x Class (c0.5), top 25 =====")
print(ann.groupby("Category").Class.value_counts().head(25).to_string())
