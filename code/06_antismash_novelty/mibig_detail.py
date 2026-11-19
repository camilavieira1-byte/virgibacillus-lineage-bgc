import json, glob, pandas as pd, numpy as np
base = "/PATH/TO/virgibacillus_v2"
rows = []
for f in sorted(glob.glob(f"{base}/06_antismash/GCA_*/GCA_*.json")):
    acc = f.split("/")[-2]; d = json.load(open(f))
    for r in d["records"]:
        kc = r.get("modules", {}).get("antismash.modules.clusterblast", {}).get("knowncluster")
        if not kc: continue
        for res in kc.get("results", []):
            best = None
            for it in res.get("ranking", []):
                if isinstance(it, (list, tuple)) and len(it) >= 2 and (best is None or (it[1].get("similarity") or 0) > (best[1].get("similarity") or 0)): best = it
            if best:
                info, sc = best
                rows.append({"accession": acc, "contig": r.get("id"), "region": int(res["region_number"]), "mibig": info.get("accession"),
                             "mibig_desc": info.get("description"), "similarity": sc.get("similarity"), "n_hits": sc.get("hits"),
                             "core_hits": sc.get("core_gene_hits"), "known_size": len(info.get("proteins", []))})
m = pd.DataFrame(rows); m.to_csv(f"{base}/06_antismash/mibig_hit_detail.tsv", sep="\t", index=False)
m["calc"] = 100 * m.n_hits / m.known_size
print("Regions with a ranked hit:", len(m))
print("Similarity vs 100*hits/known-cluster size: median abs difference = %.2f, share within 1.5 points = %.0f%%" % (
    (m.similarity - m.calc).abs().median(), 100 * ((m.similarity - m.calc).abs() <= 1.5).mean()))
h = m[m.similarity >= 50]
print("\nRegions with similarity >= 50: %d | with >= 3 matched genes: %d | with >= 5: %d" % (len(h), int((h.n_hits >= 3).sum()), int((h.n_hits >= 5).sum())))
g = h.groupby(["mibig", "mibig_desc"]).agg(regions=("region", "size"), median_sim=("similarity", "median"), matched_genes=("n_hits", "median"), known_cluster_genes=("known_size", "median")).sort_values("regions", ascending=False)
pd.set_option("display.width", 200); print("\nKnown clusters behind the >= 50% hits:"); print(g.head(15).to_string())
