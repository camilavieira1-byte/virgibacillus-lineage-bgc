import pandas as pd, numpy as np, glob
base = "/PATH/TO/virgibacillus_v2"
res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]
gc = pd.read_csv(f"{res}/Orthogroups/Orthogroups.GeneCount.tsv", sep="\t", index_col=0).drop(columns=["Total"])
un = pd.read_csv(f"{res}/Orthogroups/Orthogroups_UnassignedGenes.tsv", sep="\t", index_col=0, dtype=str)
reps = pd.read_csv(f"{base}/10_pangenome/b_reps_995.tsv", sep="\t", dtype=str)
cols = [a for a in reps.accession if a in gc.columns]; N = len(cols)
print("Genomes:", N, "| unassigned-gene file columns found:", sum(c in un.columns for c in cols))

def lens(a):
    d, cur = {}, None
    for line in open(f"{base}/10_pangenome/proteins/{a}.faa"):
        if line.startswith(">"): cur = line[1:].split()[0]; d[cur] = 0
        else: d[cur] += len(line.strip())
    return d

A0 = (gc[cols] > 0); A0 = A0.loc[A0.sum(1) > 0].values.astype(bool)
single_len, single_genome = [], []
for gi, c in enumerate(cols):
    L = lens(c)
    for cell in un[c].dropna():
        for g in str(cell).split(", "):
            g = g.strip()
            if g: single_len.append(L.get(g, np.nan)); single_genome.append(gi)
single_len = np.array(single_len, float); single_genome = np.array(single_genome)
per = np.bincount(single_genome, minlength=N)
print("Unassigned (genome-specific) genes in these genomes: %d | per genome min %d, median %d, max %d" % (len(single_genome), per.min(), np.median(per), per.max()))
print("Of these, shorter than 100 aa: %.0f%% | shorter than 50 aa: %.0f%%" % (100 * np.nanmean(single_len < 100), 100 * np.nanmean(single_len < 50)))

def with_singletons(mask):
    g = single_genome[mask]; S = np.zeros((len(g), N), bool); S[np.arange(len(g)), g] = True
    return np.vstack([A0, S])
scen = {"orthogroups only": A0,
        "+ all singletons": with_singletons(np.ones(len(single_genome), bool)),
        "+ singletons >= 100 aa": with_singletons(np.nan_to_num(single_len, nan=0) >= 100)}

def heaps(A, NP=300, seed=3):
    rng = np.random.default_rng(seed); new = np.zeros((NP, N - 1)); Ns = np.arange(2, N + 1)
    for p in range(NP):
        seen = np.zeros(A.shape[0], bool)
        for k, g in enumerate(rng.permutation(N)):
            col = A[:, g]
            if k > 0: new[p, k - 1] = (col & ~seen).sum()
            seen |= col
    alpha = np.array([-np.polyfit(np.log(Ns), np.log(np.maximum(new[p], 1)), 1)[0] for p in range(NP)])
    return alpha, new

print("\n%-26s %9s %9s %10s %12s" % ("families counted", "families", "alpha", "95% range", "new in last"))
for k, A in scen.items():
    al, new = heaps(A)
    print("%-26s %9d %9.2f %10s %12.0f" % (k, A.shape[0], al.mean(), "%.2f-%.2f" % tuple(np.percentile(al, [2.5, 97.5])), new[:, -1].mean()))
print("\nalpha < 1: open pangenome; alpha > 1: closed. Core genome (families in all genomes) is the same in all three rows.")
