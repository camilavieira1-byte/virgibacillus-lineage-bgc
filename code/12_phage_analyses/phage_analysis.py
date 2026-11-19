import pandas as pd, numpy as np, os, re, glob, itertools
base = "/PATH/TO/virgibacillus_v2"; ph = f"{base}/12_phages"
P = pd.read_csv(f"{ph}/prophages.tsv", sep="\t")
E = pd.read_csv(f"{ph}/skani_edges.tsv", sep="\t")
key = lambda p: os.path.basename(p)[:-4]
E["a"], E["b"] = E.Ref_file.map(key), E.Query_file.map(key)
E["af"] = E[["Align_fraction_ref", "Align_fraction_query"]].max(axis=1)
G = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna("").set_index("accession")
H = pd.read_csv(f"{base}/05_ani/skani_edges.tsv", sep="\t")
acc = lambda s: re.search(r"(GCA_\d+\.\d)", s).group(1)
hostani = {frozenset((acc(r.Ref_file), acc(r.Query_file))): r.ANI for r in H.itertuples()}
gen, spp = G.genus.to_dict(), G.species.to_dict(); acc_of = dict(zip(P.pid, P.accession))

def cluster(ids, ani_thr, af_thr):
    ids = set(ids); par = {i: i for i in ids}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    for r in E[(E.ANI >= ani_thr) & (E.af >= af_thr) & E.a.isin(ids) & E.b.isin(ids)].itertuples(): par[f(r.a)] = f(r.b)
    return {i: f(i) for i in ids}
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

def analyze(Q, title):
    print("\n==================== %s: %d prophages from %d genomes" % (title, len(Q), Q.accession.nunique()))
    out = {}
    for name, (ani, af) in (("vOTU level (ANI>=95, AF>=85 of the shorter)", (95, 85)), ("related level (ANI>=80, AF>=50)", (80, 50))):
        cl = cluster(Q.pid, ani, af); Qc = Q.assign(c=Q.pid.map(cl)); sz = Qc.c.value_counts()
        lv = pd.Series({c: level(d) for c, d in Qc.groupby("c")})
        print("\n%s: %d clusters | singletons %d | clusters with >= 2 prophages %d | largest %d" % (name, len(sz), int((sz == 1).sum()), int((sz >= 2).sum()), sz.max()))
        print(lv.value_counts().to_string()); out[ani] = (Qc, lv)
    sel = E[(E.ANI >= 95) & (E.af >= 85) & E.a.isin(set(Q.pid)) & E.b.isin(set(Q.pid))]
    print("\nProphage pairs at vOTU level: %d, by host relationship (large clusters give many pairs):" % len(sel))
    print(pd.Series([hostrel(acc_of[a], acc_of[b]) for a, b in zip(sel.a, sel.b)]).value_counts().to_string())
    return out

res = analyze(P, "all prophages")
Qc, lv = res[95]
shared = [c for c in lv.index if lv[c] in ("shared across genera", "shared across species (same genus)")]
print("\nvOTUs shared across species/genera: %d" % len(shared))
for c in shared[:12]:
    d = Qc[Qc.c == c]; hs = sorted({"%s | %s" % (gen.get(a, ""), spp.get(a, "")) for a in d.accession})
    print("  cluster of %d prophages in %d genomes, length %.0f kb: %s" % (len(d), d.accession.nunique(), d.length.median() / 1e3, "; ".join(hs)))
analyze(P[P.in_seq_edge != 1], "excluding prophages that touch a contig end")

print("\n== positive control: prophages in near-identical host genome pairs (host ANI >= 99.5)")
byacc = P.groupby("accession").pid.apply(list); n = hit = 0
for pair, ani in hostani.items():
    if ani < 99.5: continue
    a, b = tuple(pair)
    if a not in byacc.index or b not in byacc.index: continue
    for pa in byacc[a]:
        n += 1; sub = E[((E.a == pa) & E.b.isin(byacc[b])) | ((E.b == pa) & E.a.isin(byacc[b]))]
        hit += int(((sub.ANI >= 95) & (sub.af >= 85)).any())
print("prophages with a counterpart in a near-identical host: %d | of which matched at vOTU level: %d (%.0f%%)" % (n, hit, 100 * hit / max(n, 1)))

rx = re.compile(r"integrase|tyrosine recombinase|Xer[CD]|site-specific.*recombinase", re.I); has = {}
for f in glob.glob(f"{base}/11_genomad/usr_*/usr_*_find_proviruses/*_provirus_genes.tsv"):
    a = os.path.basename(f)[4:].split("_provirus_genes")[0]
    try: d = pd.read_csv(f, sep="\t", usecols=["gene", "annotation_description"], dtype=str)
    except Exception: continue
    for g, ann in zip(d.gene, d.annotation_description):
        if isinstance(ann, str) and rx.search(ann): has[a + "__" + g.rsplit("_", 1)[0].replace("|", "_")] = True
P["integrase_like"] = P.pid.map(has).fillna(False)
print("\nprophages with >= 1 integrase/tyrosine-recombinase-type gene: %d of %d (%.0f%%)" % (int(P.integrase_like.sum()), len(P), 100 * P.integrase_like.mean()))
P["c95"] = P.pid.map(res[95][0].set_index("pid").c); P["c80"] = P.pid.map(res[80][0].set_index("pid").c)
P.to_csv(f"{ph}/prophages_clustered.tsv", sep="\t", index=False)
