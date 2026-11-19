import pandas as pd, numpy as np, glob, os
base = "/PATH/TO/virgibacillus_v2"; gd = f"{base}/11_genomad"
def contig_lengths(fa):
    L, cur = {}, None
    for line in open(fa):
        if line.startswith(">"): cur = line[1:].split()[0]; L[cur] = 0
        elif cur: L[cur] += len(line.strip())
    return L
rows, gen, pl_cols = [], [], None
for d in sorted(glob.glob(f"{gd}/usr_*")):
    name = os.path.basename(d); acc = name[4:]
    L = contig_lengths(f"{base}/02_isolates/fna/{name}.fna")
    v = pd.read_csv(f"{d}/{name}_summary/{name}_virus_summary.tsv", sep="\t")
    pvf = f"{d}/{name}_find_proviruses/{name}_provirus.tsv"
    if os.path.exists(pvf) and os.path.getsize(pvf) > 0:
        p = pd.read_csv(pvf, sep="\t", dtype=str).set_index("seq_name")
        edge = p.in_seq_edge.eq("True").to_dict(); integ = p.integrases.notna().to_dict()
    else: edge, integ = {}, {}
    for r in v.itertuples(index=False):
        sn = r.seq_name
        if "|provirus_" in sn:
            contig = sn.split("|provirus_")[0]; s, e = (int(x) for x in str(r.coordinates).split("-")); kind = "prophage"
        else:
            contig = sn; s, e = 1, int(r.length); kind = "virus_contig"
        cl = L.get(contig, np.nan)
        rows.append({"accession": acc, "kind": kind, "contig": contig, "start": s, "end": e, "length": e - s + 1, "contig_len": cl,
                     "dist_to_contig_end": min(s - 1, cl - e) if cl == cl else np.nan,
                     "in_seq_edge": edge.get(sn, np.nan), "has_integrase": integ.get(sn, np.nan),
                     "virus_score": r.virus_score, "n_hallmarks": r.n_hallmarks, "n_genes": r.n_genes, "taxonomy": r.taxonomy})
    pf = f"{d}/{name}_summary/{name}_plasmid_summary.tsv"; npl = 0
    if os.path.exists(pf) and os.path.getsize(pf) > 0:
        ps = pd.read_csv(pf, sep="\t"); npl = len(ps)
        if pl_cols is None and npl: pl_cols = ps.columns.tolist()
    gen.append({"accession": acc, "genome_bp": sum(L.values()), "n_contigs": len(L), "plasmid_calls": npl})
m = pd.DataFrame(rows)
for c in ("in_seq_edge", "has_integrase"): m[c] = m[c].astype(float)
m.to_csv(f"{gd}/mge_regions.tsv", sep="\t", index=False)
pp = m[m.kind == "prophage"].groupby("accession").agg(n_prophage=("length", "size"), prophage_bp=("length", "sum"),
        n_prophage_edge=("in_seq_edge", "sum"), n_prophage_integrase=("has_integrase", "sum"))
vc = m[m.kind == "virus_contig"].groupby("accession").size().rename("n_virus_contigs")
g = pd.DataFrame(gen).merge(pp, on="accession", how="left").merge(vc, on="accession", how="left").fillna(0)
g["prophage_pct"] = 100 * g.prophage_bp / g.genome_bp
t = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str)[["accession", "genus", "in_A", "in_B", "clone_group", "Contig_N50"]]
t["Contig_N50"] = pd.to_numeric(t.Contig_N50)
g = g.merge(t, on="accession", how="left"); g.to_csv(f"{gd}/mge_per_genome.tsv", sep="\t", index=False)

P = m[m.kind == "prophage"]; q = g[g.Contig_N50.notna()].copy()
print("Genomes: %d | prophage regions: %d | whole-contig virus calls: %d" % (len(g), len(P), int((m.kind == "virus_contig").sum())))
print("Genomes with >= 1 prophage: %d | per genome: median %.0f (IQR %.0f-%.0f) | median %% of genome in prophage: %.2f" % (
    int((g.n_prophage > 0).sum()), g.n_prophage.median(), g.n_prophage.quantile(.25), g.n_prophage.quantile(.75), g.prophage_pct.median()))
print("Prophage length, kb: median %.0f (IQR %.0f-%.0f) | with integrase: %.0f%% | touching a contig end: %.0f%% | Caudoviricetes: %.0f%%" % (
    P.length.median() / 1e3, P.length.quantile(.25) / 1e3, P.length.quantile(.75) / 1e3, 100 * P.has_integrase.mean(), 100 * P.in_seq_edge.mean(),
    100 * P.taxonomy.fillna("").str.contains("Caudoviricetes").mean()))
print("Plasmid summary columns:", pl_cols, "| plasmid calls: %d in %d genomes" % (g.plasmid_calls.sum(), int((g.plasmid_calls > 0).sum())))
q["N50_bin"] = pd.cut(q.Contig_N50, [0, 100e3, 500e3, 1e9], labels=["<100k", "100-500k", ">500k"])
print("\nProphage detection by assembly N50 (134 QC-passing genomes):")
print(q.groupby("N50_bin", observed=True).agg(genomes=("accession", "size"), mean_prophages=("n_prophage", "mean"),
      median_prophage_pct=("prophage_pct", "median"), edge_prophages=("n_prophage_edge", "sum"), all_prophages=("n_prophage", "sum")).round(2).to_string())
reps = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str).accession
r = g[g.accession.isin(reps) & (g.Contig_N50 >= 500e3)].copy()
r["group"] = np.where(r.in_B == "True", "Virgibacillus (B)", np.where(r.genus.fillna("").str.startswith("Virgibacillus"), "other Virgibacillus_X", "other genera"))
print("\nClone-group representatives with N50 >= 500 kb:")
print(r.groupby("group").agg(reps=("accession", "size"), median_prophages=("n_prophage", "median"), median_prophage_pct=("prophage_pct", "median"),
      median_with_integrase=("n_prophage_integrase", "median")).round(2).to_string())
