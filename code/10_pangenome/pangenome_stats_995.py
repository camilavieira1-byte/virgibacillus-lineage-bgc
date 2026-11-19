import pandas as pd, numpy as np, glob
base = "/PATH/TO/virgibacillus_v2"
res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]
gc = pd.read_csv(f"{res}/Orthogroups/Orthogroups.GeneCount.tsv", sep="\t", index_col=0).drop(columns=["Total"])
reps = pd.read_csv(f"{base}/10_pangenome/b_reps_995.tsv", sep="\t", dtype=str)
def col_for(a):
    if a in gc.columns: return a
    c = [x for x in gc.columns if x == a or a.startswith(x + ".") or x.startswith(a)]
    return c[0] if len(c) == 1 else None
reps["col"] = reps.accession.map(col_for)
print("Genome columns found for the representatives:", int(reps.col.notna().sum()), "of", len(reps))
reps = reps[reps.col.notna()]
cols = reps.col.tolist(); N = len(cols)
G = gc[cols]; X = (G > 0)
X = X.loc[X.sum(1) > 0]; freq = X.sum(1)
c99, c95, c15 = int(np.ceil(0.99 * N)), int(np.ceil(0.95 * N)), int(np.ceil(0.15 * N))
cat = pd.Series(pd.Categorical(np.select([freq >= c99, freq >= c95, freq >= c15], ["core", "soft core", "shell"], "cloud"), categories=["core", "soft core", "shell", "cloud"]), index=freq.index)
tab = cat.value_counts().reindex(["core", "soft core", "shell", "cloud"])
print("\nOrthogroups in the %d genomes: %d  (core >= %d genomes, soft core >= %d, shell >= %d)" % (N, len(X), c99, c95, c15))
print(pd.DataFrame({"orthogroups": tab, "percent": (100 * tab / len(X)).round(1)}).to_string())
print("core + soft core:", int((freq >= c95).sum()), "| single-copy in all genomes:", int((G == 1).all(axis=1).sum()),
      "| present in exactly one genome:", int((freq == 1).sum()))
pg = X.sum(0)
print("orthogroups per genome: min %d, median %d, max %d" % (pg.min(), pg.median(), pg.max()))
pd.DataFrame({"genomes_with_og": freq, "category": cat}).to_csv(f"{base}/10_pangenome/pangenome_categories_995.tsv", sep="\t")

rng = np.random.default_rng(3); A = X.values.astype(bool); NP = 300
new = np.zeros((NP, N - 1)); pan = np.zeros((NP, N)); core = np.zeros((NP, N))
for p in range(NP):
    seen = np.zeros(A.shape[0], bool); inall = np.ones(A.shape[0], bool)
    for k, g in enumerate(rng.permutation(N)):
        col = A[:, g]
        if k > 0: new[p, k - 1] = (col & ~seen).sum()
        seen |= col; inall &= col; pan[p, k] = seen.sum(); core[p, k] = inall.sum()
Ns = np.arange(2, N + 1)
alpha = np.array([-np.polyfit(np.log(Ns), np.log(np.maximum(new[p], 1)), 1)[0] for p in range(NP)])
print("\nHeaps' law (new orthogroups per added genome ~ k * N^-alpha; alpha < 1 means open): alpha = %.2f (95%% of orders: %.2f-%.2f)" % (
    alpha.mean(), np.percentile(alpha, 2.5), np.percentile(alpha, 97.5)))
print("Mean new orthogroups added by the last genome: %.0f" % new[:, -1].mean())
curve = pd.DataFrame({"N": np.arange(1, N + 1), "pan_mean": pan.mean(0), "core_mean": core.mean(0), "new_mean": np.r_[np.nan, new.mean(0)]})
curve.to_csv(f"{base}/10_pangenome/heaps_curves_995.tsv", sep="\t", index=False)
print(curve.iloc[sorted({0, 4, 9, N - 1})].round(1).to_string(index=False))

print("\n== orthogroups restricted to one GTDB species (among the %d genomes)" % N)
for sp, d in reps.groupby("species"):
    S = d.col.tolist(); O = [c for c in cols if c not in S]
    only = (X[S].sum(1) > 0) & (X[O].sum(1) == 0)
    diag = only & (X[S].mean(1) >= 0.9)
    print("%-30s genomes=%2d  species-only OGs=%4d  in >=90%% of its genomes=%4d" % (sp, len(S), int(only.sum()), int(diag.sum())))
