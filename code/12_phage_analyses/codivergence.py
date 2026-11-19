import pandas as pd, numpy as np, itertools
from scipy.stats import spearmanr, rankdata, mannwhitneyu
from Bio import Phylo
rng = np.random.default_rng(11)
base = "/PATH/TO/virgibacillus_v2"; ph = f"{base}/12_phages"
P = pd.read_csv(f"{ph}/prophages.tsv", sep="\t").merge(pd.read_csv(f"{ph}/prophage_idmap.tsv", sep="\t"), on="pid")
C = pd.read_csv(f"{ph}/clusters_phabox_aai.tsv", sep="\t", header=None, names=["pid", "c"])
reps = set(pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str).accession)
A = pd.read_csv(f"{ph}/phabox_votu_aai/midfolder/self_results_aai.tsv", sep="\t"); A = A[A["query"] != A["target"]]
aai = {}
for r in A.itertuples(index=False):
    k = frozenset((r.query, r.target)); aai[k] = max(aai.get(k, 0), r.aai)
t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick"); tips = {x.name: x for x in t.get_terminals()}
G = pd.read_csv(f"{base}/09_matrix/analysis_table_genomes.tsv", sep="\t", dtype=str).fillna("").set_index("accession"); spp = G.species.to_dict()
cache = {}
def hd(a, b):
    if a == b: return 0.0
    k = frozenset((a, b))
    if k not in cache: cache[k] = t.distance(tips["usr_" + a], tips["usr_" + b])
    return cache[k]
def analyze(Q, title, nperm=9999):
    D = Q.merge(C, on="pid"); gl = []
    for c, d in D.groupby("c"):
        d = d.drop_duplicates("short")
        if len(d) < 3 or d.accession.nunique() < 3: continue
        hosts, ids = d.accession.tolist(), d.short.tolist(); n = len(ids)
        dm = np.array([[hd(a, b) for b in hosts] for a in hosts])
        pr = [(i, j, aai[frozenset((ids[i], ids[j]))]) for i, j in itertools.combinations(range(n), 2) if frozenset((ids[i], ids[j])) in aai and hosts[i] != hosts[j]]
        if len(pr) >= 3: gl.append((c, dm, np.array([p[0] for p in pr]), np.array([p[1] for p in pr]), np.array([p[2] for p in pr]), hosts))
    print("\n==================== %s" % title)
    if not gl: print("no usable groups"); return
    y = np.concatenate([g[4] for g in gl]); x = np.concatenate([g[1][g[2], g[3]] for g in gl]); ry = rankdata(y)
    obs = np.corrcoef(rankdata(x), ry)[0, 1]; null = []
    for _ in range(nperm):
        xp = np.concatenate([g[1][(pm := rng.permutation(len(g[5])))[g[2]], pm[g[3]]] for g in gl]); null.append(np.corrcoef(rankdata(xp), ry)[0, 1])
    null = np.array(null)
    print("groups used: %d | prophage pairs: %d | pairs per group: %s" % (len(gl), len(y), sorted([len(g[4]) for g in gl], reverse=True)[:8]))
    print("Spearman (prophage AAI vs host patristic distance): rho = %.3f | permutation p (rho <= obs) = %.4f | two-sided p = %.4f" % (
        obs, (np.sum(null <= obs) + 1) / (nperm + 1), (np.sum(np.abs(null) >= abs(obs)) + 1) / (nperm + 1)))
    sp = np.concatenate([[spp.get(g[5][i], "") == spp.get(g[5][j], "") and spp.get(g[5][i], "") != "" for i, j in zip(g[2], g[3])] for g in gl])
    print("pairs in the same host species: %d (median prophage AAI %.1f) | in different species: %d (median %.1f)" % (sp.sum(), np.median(y[sp]) if sp.any() else np.nan, (~sp).sum(), np.median(y[~sp]) if (~sp).any() else np.nan))
    if sp.any() and (~sp).any(): print("   Mann-Whitney same vs different species: p = %.3g" % mannwhitneyu(y[sp], y[~sp])[1])
    q = pd.qcut(x, 3, duplicates="drop"); print("median prophage AAI by host distance tertile:"); print(pd.Series(y).groupby(q, observed=True).agg(["size", "median"]).round(2).to_string())
    for g in gl:
        if len(g[4]) >= 6: print("   group %s: %d pairs, rho = %.2f" % (g[0], len(g[4]), spearmanr(g[1][g[2], g[3]], g[4])[0]))
rp = P.accession.isin(reps); ne = P.in_seq_edge != 1
analyze(P[rp], "clone-group representatives"); analyze(P[rp & ne], "representatives, without contig-end prophages")
