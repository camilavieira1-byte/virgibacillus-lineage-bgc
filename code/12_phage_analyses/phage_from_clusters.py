import pandas as pd, numpy as np, os, re, itertools
base = "/PATH/TO/virgibacillus_v2"; ph = f"{base}/12_phages"
P = pd.read_csv(f"{ph}/prophages.tsv", sep="\t")
C = pd.read_csv(os.environ["CLUSTERS"], sep="\t", header=None, names=["pid", "c"], dtype=str).drop_duplicates("pid")
reps = set(pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str).accession)
G = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna("").set_index("accession")
H = pd.read_csv(f"{base}/05_ani/skani_edges.tsv", sep="\t")
acc = lambda s: re.search(r"(GCA_\d+\.\d)", s).group(1)
hostani = {frozenset((acc(r.Ref_file), acc(r.Query_file))): r.ANI for r in H.itertuples()}
gen, spp = G.genus.to_dict(), G.species.to_dict()
def hostrel(a, b):
    if a == b: return "same genome"
    ani = hostani.get(frozenset((a, b)))
    if ani is not None and ani >= 99.5: return "near-identical hosts (ANI>=99.5)"
    ga, gb, sa, sb = gen.get(a, ""), gen.get(b, ""), spp.get(a, ""), spp.get(b, "")
    if not (ga and gb): return "host not classified"
    if sa and sa == sb: return "same species, different strains"
    if ga == gb: return "same genus, different species"
    return "different genera"
def level(d):
    if len(d) == 1: return "singleton"
    accs = sorted(set(d.accession))
    if len(accs) == 1: return "several in one genome"
    labs = [a for a in accs if gen.get(a, "")]
    if len({gen[a] for a in labs}) > 1: return "shared across genera"
    if len({spp.get(a, "") for a in labs}) > 1: return "shared across species (same genus)"
    arr = np.array([hostani.get(frozenset(p), np.nan) for p in itertools.combinations(accs, 2)], float)
    return "near-identical hosts only" if (not np.isnan(arr).all() and np.nanmin(arr) >= 99.5) else "same species, different strains"
def show(Q, title):
    Qc = Q.merge(C, on="pid", how="left"); miss = int(Qc.c.isna().sum())
    Qc.loc[Qc.c.isna(), "c"] = ["unassigned_" + p for p in Qc[Qc.c.isna()].pid]
    sz = Qc.c.value_counts(); lv = pd.Series({c: level(d) for c, d in Qc.groupby("c")})
    print("\n==================== %s: %d prophages from %d genomes (not in the cluster table, treated as singletons: %d)" % (title, len(Qc), Qc.accession.nunique(), miss))
    print("clusters %d | singletons %d | with >= 2 prophages %d | largest %d" % (len(sz), int((sz == 1).sum()), int((sz >= 2).sum()), sz.max()))
    print(lv.value_counts().to_string())
    rel = []
    for c, d in Qc.groupby("c"):
        if len(d) > 1: rel += [hostrel(a, b) for a, b in itertools.combinations(d.accession.tolist(), 2)]
    print("\nprophage pairs within clusters: %d, by host relationship:" % len(rel)); print(pd.Series(rel, dtype=str).value_counts().to_string())
    for c in [c for c in lv.index if lv[c].startswith("shared across")][:10]:
        d = Qc[Qc.c == c]; print("  shared cluster: %d prophages in %d genomes, median length %.0f kb: %s" % (len(d), d.accession.nunique(), d.length.median() / 1e3, "; ".join(sorted({"%s | %s" % (gen.get(a, ""), spp.get(a, "")) for a in d.accession}))))
ne = P.in_seq_edge != 1; rp = P.accession.isin(reps)
show(P, "all prophages"); show(P[rp], "clone-group representatives only")
show(P[ne], "without prophages touching a contig end"); show(P[rp & ne], "representatives, without contig-end prophages")
