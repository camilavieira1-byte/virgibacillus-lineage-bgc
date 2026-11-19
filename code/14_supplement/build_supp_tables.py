import pandas as pd, numpy as np, glob, os, re, shutil, sys
from Bio import Phylo
base = "/PATH/TO/virgibacillus_v2"; out = f"{base}/14_supplement"; M = f"{base}/09_matrix"
os.makedirs(f"{out}/tables_csv", exist_ok=True); print("python:", sys.executable); sheets = {}
def rd(p, **k): return pd.read_csv(p, sep="\t", **k)
def cgt(): return rd(f"{M}/analysis_table_clonegroups.tsv", dtype=str).fillna("")
def add(name, fn):
    try: r = fn(); sheets[name] = r; print("ok    %-27s %6d rows x %3d cols" % (name, len(r), r.shape[1]))
    except Exception as e: print("FAIL  %-27s %s: %s" % (name, type(e).__name__, str(e)[:220]))
def s1():
    d = None
    for fn in ("metadata_raw_flagged.tsv", "metadata_raw.tsv"):
        p = f"{base}/00_metadata/{fn}"
        if os.path.exists(p): d = rd(p, dtype=str); print("      metadata file:", fn, d.shape); break
    if d is None: raise Exception("no metadata file; 00_metadata has: " + ", ".join(sorted(os.listdir(f"{base}/00_metadata"))[:12]))
    ac = next((c for c in d.columns if c.lower() in ("accession", "assembly_accession", "assembly accession")), None)
    if ac is None: raise Exception("no accession column; columns: " + ", ".join(d.columns[:15]))
    d = d.rename(columns={ac: "accession"})
    q = rd(f"{base}/03_qc/checkm2/quality_report.tsv"); q["accession"] = q.Name.str.replace(r"^usr_", "", regex=True).str.replace(r"\.fna$", "", regex=True)
    q = q.drop(columns=[c for c in ("Name", "Additional_Notes") if c in q.columns]); q.columns = [c if c == "accession" else "checkm2_" + c for c in q.columns]
    gt = rd(f"{base}/04_gtdbtk_b/gtdb_classification.tsv"); gt["accession"] = gt.user_genome.str.replace(r"^usr_", "", regex=True)
    gt = gt[["accession"] + [c for c in gt.columns if c not in ("accession", "user_genome", "Completeness", "Contamination", "pass_qc", "is_virgi")]]
    gt.columns = [c if c == "accession" else "gtdb_" + c for c in gt.columns]
    d = d.merge(q, on="accession", how="left").merge(gt, on="accession", how="left")
    comp = pd.to_numeric(d.checkm2_Completeness, errors="coerce"); con = pd.to_numeric(d.checkm2_Contamination, errors="coerce")
    d["status"] = np.where(comp.isna(), "metagenome-assembled genome (excluded)", np.where((comp >= 90) & (con < 5), "passed QC", "failed QC"))
    print("      status:", d.status.value_counts().to_dict(), "(expected 134 passed, 4 failed, 12 MAGs)")
    G = rd(f"{M}/analysis_table_genomes.tsv", dtype=str); lab = {r: "CG%03d" % (i + 1) for i, r in enumerate(G.clone_group.value_counts().index)}
    G["clone_group_ANI99.9"] = G.clone_group.map(lab); G["clone_group_size"] = G.groupby("clone_group").accession.transform("size")
    keep = ["accession", "clone_group_ANI99.9", "clone_group_size"] + [c for c in ("in_A", "in_B") if c in G.columns]
    d = d.merge(G[keep].rename(columns={"in_A": "in_Virgibacillus_clade", "in_B": "in_GTDB_Virgibacillus"}), on="accession", how="left")
    d["representative_ANI99.9"] = d.accession.isin(set(cgt().accession)); d["representative_ANI99.5"] = d.accession.isin(set(rd(f"{M}/analysis_table_clonegroups_995.tsv", dtype=str).accession))
    ok = d[d.status == "passed QC"]
    print("      GTDB Virgibacillus species among QC-passing genomes:", ok[ok.gtdb_genus == "Virgibacillus"].gtdb_species.value_counts().to_dict(), "(Results 1 says 31, 8, 8, 6, 4, 1)")
    print("      genomes in genera named Virgibacillus*:", int(d.gtdb_genus.fillna("").str.startswith("Virgibacillus").sum()), "of", int(d.gtdb_genus.notna().sum()), "| clone groups:", G.clone_group.nunique(), "| largest group:", int(G.clone_group_size.astype(int).max()))
    return d.sort_values(["status", "gtdb_genus", "accession"], na_position="last")
def s2():
    rec = rd(f"{M}/gcf_records_c0.5.tsv"); st = rd(f"{base}/06_antismash/mibig_strict_per_region.tsv").drop_duplicates("gbk")
    d = rec.merge(st[["gbk"] + [c for c in st.columns if c not in rec.columns]], on="gbk", how="left"); d["strict"] = d["strict"].fillna("no hit")
    cat = d.strict.str.extract(r"^(known|short|weak|no hit)", expand=False)
    print("      MIBiG status of regions:", cat.value_counts().to_dict(), "of", len(d), "(Results 3: no hit 669, weak 106, short 77, known 116)")
    for p in ("T3PKS", "terpene", "RiPP-like", "ectoine"):
        s = d[d["product"] == p]; print("      %-10s n=%3d | no hit %3.0f%% | known %3.0f%%" % (p, len(s), 100 * (cat[s.index] == "no hit").mean(), 100 * (cat[s.index] == "known").mean()))
    mc = next((c for c in d.columns if "mibig" in c.lower() and "desc" not in c.lower() and "type" not in c.lower()), None)
    if mc: print("      most frequent MIBiG entries among known regions:", d.loc[cat == "known", mc].value_counts().head(8).to_dict())
    d = d.rename(columns={"family": "GCF_c0.5", "singleton": "singleton_c0.5"})
    for c in ("0.3", "0.7"):
        r = rd(f"{M}/gcf_records_c{c}.tsv")[["gbk", "family", "singleton"]].rename(columns={"family": f"GCF_c{c}", "singleton": f"singleton_c{c}"}); d = d.merge(r, on="gbk", how="left")
    return d
def s3():
    rec = rd(f"{M}/gcf_records_c0.5.tsv"); st = rd(f"{base}/06_antismash/mibig_strict_per_region.tsv").drop_duplicates("gbk")[["gbk", "strict"]]
    rec = rec.merge(st, on="gbk", how="left"); rec["strict"] = rec.strict.fillna("no hit"); ec = next((c for c in rec.columns if "edge" in c.lower()), None)
    rec["is_known"] = rec.strict.str.startswith("known"); rec["is_short"] = rec.strict.str.startswith("short"); rec["is_weak"] = rec.strict.str.startswith("weak")
    agg = dict(regions=("gbk", "size"), genomes=("genome", "nunique"), main_product=("product", lambda s: s.value_counts().index[0]), n_known=("is_known", "sum"), n_short=("is_short", "sum"), n_weak=("is_weak", "sum"))
    if ec: agg["pct_regions_on_contig_edge"] = (ec, lambda s: round(100 * s.mean(), 1))
    g = rec.groupby("family").agg(**agg)
    g["MIBiG_status"] = np.select([g.n_known >= 0.5 * g.regions, g.n_known > 0, g.n_short > 0, g.n_weak > 0], ["known", "partly known", "short match only", "weak only"], "no hit")
    print("      family status:", g.MIBiG_status.value_counts().to_dict(), "(Results 3: no hit 155, known 35, weak only 22, short match only 13, partly known 5)")
    return g.reset_index().sort_values("regions", ascending=False)
def s4():
    m = rd(f"{M}/gcf_matrix_c0.5_all.tsv", index_col=0); m.index.name = "accession"; return m.reset_index()
def s5():
    recs = {c: rd(f"{M}/gcf_records_c{c}.tsv") for c in ("0.3", "0.5", "0.7")}; rows = []
    def ari(a, b):
        ct = pd.crosstab(np.asarray(a), np.asarray(b)).values.astype(float); n = ct.sum(); c2 = lambda x: x * (x - 1) / 2
        s, sa, sb = c2(ct).sum(), c2(ct.sum(axis=1)).sum(), c2(ct.sum(axis=0)).sum(); e = sa * sb / c2(n); return (s - e) / ((sa + sb) / 2 - e)
    for c, r in recs.items():
        ec = next((x for x in r.columns if "edge" in x.lower()), None); tot = r.family.nunique() + int(r.family.isna().sum())
        rows += [("regions", c, len(r)), ("GCFs including single-region families", c, tot), ("GCFs with >= 2 regions", c, r.loc[~r.singleton, "family"].nunique()), ("single-region families", c, int(r.singleton.sum()))]
        if ec: rows += [("% of single-region families on a contig edge", c, round(100 * r.loc[r.singleton, ec].mean(), 1)), ("% of regions in multi-region families on a contig edge", c, round(100 * r.loc[~r.singleton, ec].mean(), 1))]
        print("      cutoff %s: GCFs %d (Results 3 says 300, 230, 161), single-region %d" % (c, tot, int(r.singleton.sum())))
    for x, y, cited in (("0.3", "0.5", 0.59), ("0.5", "0.7", 0.89), ("0.3", "0.7", 0.54)):
        m = recs[x][["gbk", "family", "singleton"]].merge(recs[y][["gbk", "family", "singleton"]], on="gbk"); lx = m.family_x.where(m.family_x.notna(), "s_" + m.gbk); ly = m.family_y.where(m.family_y.notna(), "s_" + m.gbk)
        v = ari(lx, ly); rows.append(("adjusted Rand index", "%s vs %s" % (x, y), round(v, 3))); print("      ARI %s vs %s = %.3f (text: %.2f)" % (x, y, v, cited))
    return pd.DataFrame(rows, columns=["metric", "cutoff_or_pair", "value"])
def s6s():
    for f in sorted(glob.glob(f"{M}/gcf_scatter_test*.tsv")):
        if f.endswith("_Aprime.tsv"): continue
        suf = os.path.basename(f)[len("gcf_scatter_test"):-4].strip("_"); add(("S6_GCF_tree_pattern" if not suf else "S6_tree_pattern_" + suf)[:31], lambda f=f: rd(f))
def s7():
    cg = cgt(); cg = cg[cg.in_A == "True"].copy(); cg["N50"] = pd.to_numeric(cg.Contig_N50, errors="coerce")
    r995 = set(rd(f"{M}/analysis_table_clonegroups_995.tsv", dtype=str).accession)
    t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick"); tips = {x.name: x for x in t.get_terminals()}; cache = {}; rng = np.random.default_rng(5)
    def pdist(a, b):
        k = (a, b) if a < b else (b, a)
        if k not in cache: cache[k] = t.distance(tips["usr_" + a], tips["usr_" + b])
        return cache[k]
    def stats(D, P, lab, nperm=4999):
        n = len(lab); iu = np.triu_indices(n, 1); a = pd.Series(D[iu]).rank().values; b = pd.Series(P[iu]).rank().values
        rho = np.corrcoef(a, b)[0, 1]; RP = np.zeros((n, n)); RP[iu] = b; RP = RP + RP.T
        hit = sum(np.corrcoef(a, RP[np.ix_(p, p)][iu])[0, 1] >= rho for p in (rng.permutation(n) for _ in range(nperm)))
        D2 = D ** 2; sst = D2[iu].sum() / n; ssw = 0.0
        for g in pd.unique(lab):
            idx = np.where(lab == g)[0]; k = len(idx)
            if k > 1: ssw += D2[np.ix_(idx, idx)][np.triu_indices(k, 1)].sum() / k
        R2 = (sst - ssw) / sst; ng = len(set(lab)); return rho, (hit + 1) / (nperm + 1), R2, 1 - (1 - R2) * (n - 1) / (n - ng)
    def one(label, cut, mpath, sub):
        sub = sub[sub.accession.map(lambda a: "usr_" + a in tips)]; rec = rd(f"{M}/gcf_records_c{cut}.tsv"); fp = rec.groupby("family")["product"].agg(lambda s: s.value_counts().index[0])
        X = rd(mpath, index_col=0).reindex(sub.accession.tolist()).dropna(how="all"); sub = sub.set_index("accession").loc[X.index].reset_index()
        keep = [f for f in X.columns if X[f].sum() >= 2 and fp.get(f, "") != "terpene-precursor"]
        B = X[keep].fillna(0).values.astype(int); inter = B @ B.T; sz = B.sum(axis=1); un = sz[:, None] + sz[None, :] - inter
        D = 1 - np.where(un > 0, inter / np.maximum(un, 1), 0); np.fill_diagonal(D, 0)
        r = sub.accession.tolist(); P = np.array([[0.0 if a == b else pdist(a, b) for b in r] for a in r]); rho, p, R2, R2a = stats(D, P, sub.genus.values)
        print("      %-50s genomes %3d | GCFs %3d | Mantel r %.3f (p %.3f) | genus R2 %.3f (adj %.3f)" % (label, len(r), len(keep), rho, p, R2, R2a))
        return {"setting": label, "bigscape_cutoff": cut, "genomes": len(r), "GCFs_analysed": len(keep), "Mantel_r": round(rho, 3), "Mantel_p_4999_permutations": round(p, 4), "PERMANOVA_R2_genus": round(R2, 3), "PERMANOVA_R2_adjusted": round(R2a, 3)}
    mats = sorted(glob.glob(f"{M}/gcf_matrix_c0.5_*.tsv")); print("      matrices at 0.5:", [os.path.basename(m) for m in mats]); main = f"{M}/gcf_matrix_c0.5_all.tsv"; comp = [m for m in mats if not m.endswith("_all.tsv")]
    plan = [("main: cutoff 0.5, all regions, 99.9% representatives", "0.5", main, cg), ("cutoff 0.3", "0.3", f"{M}/gcf_matrix_c0.3_all.tsv", cg), ("cutoff 0.7", "0.7", f"{M}/gcf_matrix_c0.7_all.tsv", cg), ("genomes with N50 >= 100 kb", "0.5", main, cg[cg.N50 >= 100000]), ("99.5% representatives", "0.5", main, cg[cg.accession.isin(r995)])]
    if comp: plan.insert(3, ("complete (non-edge) regions only", "0.5", comp[0], cg))
    rows = []
    for lab, cut, mp, sub in plan:
        try: rows.append(one(lab, cut, mp, sub))
        except Exception as e: print("      skipped %s: %s: %s" % (lab, type(e).__name__, str(e)[:120]))
    return pd.DataFrame(rows)
def s8a(): return pd.concat([rd(f) for f in sorted(glob.glob(f"{base}/13_figures/pangenome_summary_*.tsv"))], ignore_index=True)
def s8bc():
    res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]; gc = rd(f"{res}/Orthogroups/Orthogroups.GeneCount.tsv", index_col=0).drop(columns=["Total"])
    cg = cgt(); sp = dict(zip(cg.accession, cg.species)); cat_rows, sp_rows = [], []
    for name, f in (("99.9% (33 genomes)", "b_reps.tsv"), ("99.5% (17 genomes)", "b_reps_995.tsv")):
        reps = rd(f"{base}/10_pangenome/{f}", dtype=str); cols = [a for a in reps.accession if a in gc.columns]; A = gc[cols] > 0; A = A.loc[A.sum(axis=1) > 0]; n = len(cols); fq = A.sum(axis=1)
        c99, c95, c15 = np.ceil(.99 * n), np.ceil(.95 * n), np.ceil(.15 * n)
        for cat, m in (("core (>= 99%)", fq >= c99), ("soft core (95-99%)", (fq >= c95) & (fq < c99)), ("shell (15-95%)", (fq >= c15) & (fq < c95)), ("cloud (< 15%)", fq < c15)):
            cat_rows.append({"set": name, "category": cat, "orthogroups": int(m.sum()), "pct_of_orthogroups": round(100 * m.sum() / len(A), 1)})
        cat_rows.append({"set": name, "category": "all orthogroups", "orthogroups": len(A), "pct_of_orthogroups": 100.0})
        for s in sorted({sp.get(a, "") for a in cols} - {""}):
            mine = [a for a in cols if sp.get(a) == s]; other = [a for a in cols if sp.get(a) != s]; only = (A[mine].sum(axis=1) > 0) & (A[other].sum(axis=1) == 0)
            sp_rows.append({"set": name, "species": s, "genomes": len(mine), "species_only_orthogroups": int(only.sum()), "of_which_in_90pct_or_more_of_the_species_genomes": int((only & (A[mine].sum(axis=1) >= np.ceil(.9 * len(mine)))).sum())})
    print("      99.9% set categories:", [(r["category"], r["orthogroups"]) for r in cat_rows[:5]], "(summary: core 2,146, soft core 159, shell 3,845, cloud 1,846; 7,996 total)")
    sheets["S8c_species_restricted"] = pd.DataFrame(sp_rows); print("ok    %-27s %6d rows" % ("S8c_species_restricted", len(sp_rows))); return pd.DataFrame(cat_rows)
def s9a():
    ph = f"{base}/12_phages"; P = rd(f"{ph}/prophages.tsv")
    try: P = P.merge(rd(f"{ph}/prophage_idmap.tsv"), on="pid", how="left")
    except Exception as e: print("      no id map:", e)
    for f in sorted(glob.glob(f"{ph}/clusters_phabox_*.tsv")):
        tag = os.path.basename(f).replace("clusters_phabox_", "").replace(".tsv", ""); P = P.merge(pd.read_csv(f, sep="\t", header=None, names=["pid", "phabox_" + tag + "_cluster"]), on="pid", how="left"); print("      added", os.path.basename(f))
    try:
        I = rd(f"{ph}/integrase_by_prophage.tsv").rename(columns={"integrase_type": "integrase_type_annotation", "integrase_not_xer": "integrase_annotation_excluding_XerCD"})
        P = P.merge(I[["pid", "integrase_type_annotation", "integrase_annotation_excluding_XerCD"]], on="pid", how="left")
        for c in ("integrase_type_annotation", "integrase_annotation_excluding_XerCD"): P[c] = P[c].fillna(False).astype(bool)
        print("      prophages with an integrase-type annotation:", int(P.integrase_type_annotation.sum()), "(text: 29) | excluding XerC/D:", int(P.integrase_annotation_excluding_XerCD.sum()), "(text: 18)")
    except Exception as e: print("      no integrase file:", e)
    return P
def s9b():
    fs = glob.glob(f"{base}/**/mge_per_genome.tsv", recursive=True)
    if not fs: raise Exception("mge_per_genome.tsv not found")
    d = rd(fs[0]); print("      census file:", fs[0].replace(base, "~")); return d.drop(columns=[c for c in ("n_prophage_integrase",) if c in d.columns])
add("S1_genomes", s1); add("S2_BGC_regions", s2); add("S3_GCFs_c0.5", s3); add("S4_genome_x_GCF_c0.5", s4); add("S5_cluster_stability", s5); s6s()
add("S7_lineage_robustness", s7); add("S8a_pangenome_summary", s8a); add("S8b_orthogroup_categories", s8bc); add("S9a_prophages", s9a); add("S9b_prophage_census", s9b)
desc = {"S1_genomes": "All 150 GenBank assemblies: NCBI metadata, CheckM2, GTDB r232 classification, QC status, study-set flags, clone groups (ANI >= 99.9%), representative flags at 99.9% and 99.5%",
        "S2_BGC_regions": "968 antiSMASH regions: product class, contig-edge flag, MIBiG similarity and strict status, GCF at BiG-SCAPE cutoffs 0.3, 0.5 and 0.7",
        "S3_GCFs_c0.5": "GCFs at cutoff 0.5 (including single-region families): members, genomes, main product, MIBiG status",
        "S4_genome_x_GCF_c0.5": "Presence/absence of each GCF per genome (cutoff 0.5, all regions); other matrices are in tables_csv",
        "S5_cluster_stability": "Family counts and adjusted Rand index between BiG-SCAPE cutoffs",
        "S7_lineage_robustness": "Mantel (4,999 permutations) and PERMANOVA R2 by genus of GCF content under each sensitivity setting",
        "S8a_pangenome_summary": "Pangenome accumulation: Heaps' alpha under three ways of counting genome-specific genes",
        "S8b_orthogroup_categories": "Core, soft core, shell and cloud orthogroups for the 99.9% and 99.5% sets",
        "S8c_species_restricted": "Orthogroups restricted to one GTDB species",
        "S9a_prophages": "geNomad prophages with PhaBOX2 vOTU and AAI group assignments",
        "S9b_prophage_census": "Prophage and plasmid census per genome"}
desc.update({k: "Per-GCF test of whether carriers cluster on the tree (SES, FDR)" for k in sheets if k.startswith("S6_")})
with pd.ExcelWriter(f"{out}/Supplementary_Tables.xlsx", engine="openpyxl") as w:
    pd.DataFrame({"sheet": list(sheets), "content": [desc.get(k, "") for k in sheets], "rows": [len(v) for v in sheets.values()]}).to_excel(w, sheet_name="README", index=False)
    for k, df in sheets.items():
        df.to_excel(w, sheet_name=k[:31], index=False); w.sheets[k[:31]].freeze_panes = "A2"; df.to_csv(f"{out}/tables_csv/{k}.csv", index=False)
for f in glob.glob(f"{M}/gcf_matrix_*.tsv"): shutil.copy(f, f"{out}/tables_csv/")
print("wrote", f"{out}/Supplementary_Tables.xlsx", "| sheets:", len(sheets))
