import json, glob, pandas as pd
base = "/PATH/TO/virgibacillus_v2"
rows = []
for f in sorted(glob.glob(f"{base}/06_antismash/GCA_*/GCA_*.json")):
    acc = f.split("/")[-2]; d = json.load(open(f))
    for r in d["records"]:
        kc = r.get("modules", {}).get("antismash.modules.clusterblast", {}).get("knowncluster")
        if not kc: continue
        for res in kc.get("results", []):
            best = (0, "", "", "")
            for item in res.get("ranking", []):
                if isinstance(item, (list, tuple)) and len(item) >= 2:
                    info, sc = item[0], item[1]; s = sc.get("similarity", 0) or 0
                    if s > best[0]: best = (s, info.get("accession", ""), info.get("description", ""), info.get("cluster_type", ""))
            rows.append({"accession": acc, "contig": r.get("id"), "region": int(res["region_number"]),
                         "best_similarity": best[0], "mibig": best[1], "mibig_desc": best[2], "mibig_type": best[3]})
m = pd.DataFrame(rows)
reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t"); reg["region"] = reg.region.astype(int)
m = reg.merge(m, on=["accession", "contig", "region"], how="left")
m["best_similarity"] = m.best_similarity.fillna(0)
m["bin"] = pd.cut(m.best_similarity, [-1, 0, 49, 79, 101], labels=["no hit", "low (1-49%)", "medium (50-79%)", "high (>=80%)"])
m.to_csv(f"{base}/06_antismash/mibig_best_per_region.tsv", sep="\t", index=False)
print("Regions:", len(m)); print(m.bin.value_counts().to_string())
print("\nBy product (top 12), % of regions in each bin:")
top = m["product"].value_counts().head(12).index
print((pd.crosstab(m[m["product"].isin(top)]["product"], m[m["product"].isin(top)].bin, normalize="index") * 100).round(0).to_string())
print("\nComplete regions only:"); print(m[~m.contig_edge].bin.value_counts().to_string())
h = m[m.best_similarity >= 50].groupby(["mibig", "mibig_desc"]).agg(regions=("region", "size"), genomes=("accession", "nunique"), median_sim=("best_similarity", "median")).sort_values("regions", ascending=False)
print("\nMost frequent known clusters at >= 50% similarity:"); print(h.head(15).to_string())
