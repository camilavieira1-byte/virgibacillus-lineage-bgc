import pandas as pd, numpy as np, itertools
from scipy.stats import spearmanr, rankdata
from Bio import Phylo
rng = np.random.default_rng(21)
base = "/PATH/TO/virgibacillus_v2"; ph = f"{base}/12_phages"
P = pd.read_csv(f"{ph}/prophages.tsv", sep="\t").merge(pd.read_csv(f"{ph}/prophage_idmap.tsv", sep="\t"), on="pid")
C = pd.read_csv(f"{ph}/clusters_phabox_aai.tsv", sep="\t", header=None, names=["pid", "c"])
reps = set(pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str).accession)
A = pd.read_csv(f"{ph}/phabox_votu_aai/midfolder/self_results_aai.tsv", sep="\t"); A = A[A["query"] != A["target"]]
aai, sg = {}, {}
for r in A.itertuples(index=False):
    k = frozenset((r.query, r.target))
    if r.aai > aai.get(k, -1): aai[k] = r.aai; sg[k] = r.sgenes
t = Phylo.read(f"{base}/08_tree/virgi_iqtree.treefile", "newick"); tips = {x.name: x for x in t.get_terminals()}
cache = {}
def hd(a, b):
    if a == b: return 0.0
    k = frozenset((a, b))
    if k not in cache: cache[k] = t.distance(tips["usr_" + a], tips["usr_" + b])
    return cache[k]
def groups(Q, min_sg):
    D = Q.merge(C, on="pid"); gl = []
    for c, d in D.groupby("c"):
        d = d.drop_duplicates("short")
        if len(d) < 3 or d.accession.nunique() < 3: continue
        hosts, ids = d.accession.tolist(), d.short.tolist(); n = len(ids)
        dm = np.array([[hd(a, b) for b in hosts] for a in hosts])
        pr = [(i, j, aai[frozenset((ids[i], ids[j]))]) for i, j in itertools.combinations(range(n), 2)
              if frozenset((ids[i], ids[j])) in aai and sg[frozenset((ids[i], ids[j]))] >= min_sg and hosts[i] != hosts[j]]
        if len(pr) >= 3: gl.append((c, dm, np.array([p[0] for p in pr]), np.array([p[1] for p in pr]), np.array([p[2] for p in pr]), hosts, d.length.median() / 1e3))
    return gl
def pooled(gs, nperm):
    y = np.concatenate([g[4] for g in gs]); x = np.concatenate([g[1][g[2], g[3]] for g in gs]); ry = rankdata(y)
    obs = np.corrcoef(rankdata(x), ry)[0, 1]; null = []
    for _ in range(nperm):
        xp = np.concatenate([g[1][pm[g[2]], pm[g[3]]] for g in gs for pm in [rng.permutation(len(g[5]))]])
        null.append(np.corrcoef(rankdata(xp), ry)[0, 1])
    return obs, (np.sum(np.array(null) <= obs) + 1) / (nperm + 1), len(y)
def equalw(gs, nperm, minpairs=6):
    gs = [g for g in gs if len(g[4]) >= minpairs]
    if not gs: return np.nan, np.nan, 0
    ry = [rankdata(g[4]) for g in gs]
    def stat(perms):
        r = []
        for g, ryg, pm in zip(gs, ry, perms):
            xg = g[1][pm[g[2]], pm[g[3]]]
            r.append(np.corrcoef(rankdata(xg), ryg)[0, 1] if np.ptp(xg) > 0 and np.ptp(g[4]) > 0 else np.nan)
        return np.nanmean(r)
    obs = stat([np.arange(len(g[5])) for g in gs]); null = np.array([stat([rng.permutation(len(g[5])) for g in gs]) for _ in range(nperm)])
    return obs, (np.sum(null <= obs) + 1) / (nperm + 1), len(gs)
def run(Q, title, min_sg, nperm=4999):
    gl = groups(Q, min_sg); print("\n==================== %s | pairs need >= %d shared genes" % (title, min_sg))
    if len(gl) < 2: print("too few groups:", len(gl)); return
    o, p, n = pooled(gl, nperm); print("pooled, all groups (%d groups, %d pairs): rho = %.3f, p = %.4f" % (len(gl), n, o, p))
    big = max(gl, key=lambda g: len(g[4])); rest = [g for g in gl if g is not big]
    o, p, n = pooled(rest, nperm); print("pooled, without the largest group %s (%d pairs): rho = %.3f, p = %.4f" % (big[0], len(big[4]), o, p))
    o, p, k = equalw(gl, nperm); print("each group weighted equally (%d groups with >= 6 pairs): mean rho = %.3f, p = %.4f" % (k, o, p))
    print("per group:")
    for g in sorted(gl, key=lambda g: -len(g[4])):
        if len(g[4]) >= 6: print("   %-9s prophages %2d | pairs %3d | median length %3.0f kb | rho = %5.2f" % (g[0], len(g[5]), len(g[4]), g[6], spearmanr(g[1][g[2], g[3]], g[4])[0]))
rp = P.accession.isin(reps); ne = P.in_seq_edge != 1
for ms in (0, 5): run(P[rp], "clone-group representatives", ms)
run(P[rp & ne], "representatives, without contig-end prophages", 0)
