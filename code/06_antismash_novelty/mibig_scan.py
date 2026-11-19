import json, glob, pandas as pd
base = "/PATH/TO/virgibacillus_v2"
rows, example = [], None
for f in sorted(glob.glob(f"{base}/06_antismash/GCA_*/GCA_*.json")):
    acc = f.split("/")[-2]
    d = json.load(open(f))
    for r in d["records"]:
        kc = r.get("modules", {}).get("antismash.modules.clusterblast", {}).get("knowncluster")
        if not kc: continue
        for res in kc.get("results", []):
            th = res.get("total_hits", 0)
            rows.append({"accession": acc, "contig": r.get("id"), "region": int(res.get("region_number")), "total_hits": th})
            if th > 0 and example is None: example = (acc, r.get("id"), res)
df = pd.DataFrame(rows)
reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t"); reg["region"] = reg.region.astype(int)
m = reg.merge(df, on=["accession", "contig", "region"], how="left")
m.to_csv(f"{base}/06_antismash/mibig_hits_per_region.tsv", sep="\t", index=False)
print("Regions in antiSMASH table: %d | with a knowncluster result: %d | unmatched: %d" % (len(reg), int(m.total_hits.notna().sum()), int(m.total_hits.isna().sum())))
m = m[m.total_hits.notna()].copy(); m["has_hit"] = m.total_hits > 0
print("Regions with >= 1 significant MIBiG hit: %d of %d (%.1f%%)" % (m.has_hit.sum(), len(m), 100 * m.has_hit.mean()))
print("Hits per region among those with hits: median %.0f, range %d-%d" % (m[m.has_hit].total_hits.median(), m[m.has_hit].total_hits.min(), m[m.has_hit].total_hits.max()))
print("\nShare with a hit, contig edge vs complete:"); print(m.groupby("contig_edge").has_hit.agg(["size", "mean"]).round(3).to_string())
t = m.groupby("product").has_hit.agg(regions="size", pct_with_hit=lambda s: round(100 * s.mean(), 1)).sort_values("regions", ascending=False).head(15)
print("\nBy antiSMASH product (top 15 by regions):"); print(t.to_string())
if example:
    acc, rid, res = example
    print("\n=== first region with hits: %s %s region %s | total_hits %s" % (acc, rid, res.get("region_number"), res.get("total_hits")))
    rk = res.get("ranking", [])
    print("ranking: %s of length %d" % (type(rk).__name__, len(rk)))
    if rk: print(json.dumps(rk[0], default=str)[:1800])
