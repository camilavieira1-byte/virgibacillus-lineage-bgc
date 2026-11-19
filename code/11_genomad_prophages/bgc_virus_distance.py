import re, os, glob, argparse
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, mannwhitneyu, spearmanr

ap = argparse.ArgumentParser()
ap.add_argument("--min-contig", type=int, default=100000, help="keep BGCs on contigs >= this length (bp)")
ap.add_argument("--low", type=int, default=3, help="low prevalence = GCF with <= this many genomes")
ap.add_argument("--reps", default=None, help="file with representative genomes (GCA_xxx.x per line)")
ap.add_argument("--complete-only", action="store_true", help="drop BGCs on a contig edge")
args = ap.parse_args()

base = "/PATH/TO/virgibacillus_v2"
norm = lambda s: re.sub(r"\.\d+$", "", str(s))
WINDOWS_KB = [0, 10, 25, 50]

# ---- BGCs: ingroup, prevalence, optional dereplication ---------------------
a = pd.read_csv(f"{base}/06_antismash/bgc_mobility_prevalence.tsv", sep="\t")
tax = pd.read_csv(f"{base}/04_gtdbtk_b/gtdb_classification.tsv", sep="\t", usecols=["accession", "is_virgi"])
ingroup = set(tax.loc[tax["is_virgi"].astype(str) == "True", "accession"])
a = a[a["genome"].isin(ingroup)].copy()
a["edge"] = a["edge"].astype(str).map({"True": 1.0, "False": 0.0})

if args.reps:
    reps = {l.strip() for l in open(args.reps) if l.strip()}
    prev_rep = a[a["genome"].isin(reps)].groupby("fam")["genome"].nunique()
    a = a[a["genome"].isin(reps)].copy()
    a["prev"] = a["fam"].map(prev_rep).fillna(1)           # singletons = 1
    print("representatives used:", a["genome"].nunique())

# ---- contig lengths from the antiSMASH record GBKs (LOCUS line) ------------
lens = {}
for f in glob.glob(f"{base}/06_antismash/GCA_*/GCA_*.gbk"):
    if ".region" in os.path.basename(f): continue
    with open(f) as fh:
        for line in fh:
            if line.startswith("LOCUS"):
                q = line.split()
                if len(q) >= 3 and q[2].isdigit():
                    lens[norm(q[1])] = int(q[2])
a["clen"] = a["contig"].map(lens)
print("contig lengths found:", len(lens), "| BGCs without length:", a["clen"].isna().sum())
if a["clen"].notna().sum() == 0:
    raise SystemExit("No contig lengths found: check the 06_antismash layout (see the ls commands).")
a = a[a["clen"] >= args.min_contig].copy()
if args.complete_only: a = a[a["edge"] == 0].copy()
a["low"] = (a["prev"] <= args.low).astype(int)
print(f"BGCs kept (contig >= {args.min_contig} bp, complete_only={args.complete_only}): {len(a)} "
      f"| low-prevalence: {a['low'].sum()} | genomes: {a['genome'].nunique()}")

# ---- geNomad virus calls -----------------------------------------------------
vall = pd.read_csv(f"{base}/11_genomad/all_viruses.tsv", sep="\t")
vall["contig"] = vall["seq_name"].str.split("|").str[0].map(norm)
x = vall["seq_name"].str.extract(r"\|provirus_(\d+)_(\d+)")
vall["pstart"] = pd.to_numeric(x[0]).fillna(0)           # whole-contig virus -> 0..inf
vall["pend"] = pd.to_numeric(x[1]).fillna(np.inf)

def analyse(vir, label):
    print("\n" + "=" * 70 + f"\n{label}: {len(vir)} virus calls\n" + "=" * 70)
    vd = {c: g[["pstart", "pend"]].to_numpy() for c, g in vir.groupby("contig")}
    d = a.copy()
    def dist(r):
        v = vd.get(r["contig"])
        if v is None: return np.nan
        return np.maximum(0, np.maximum(v[:, 0] - r["end"], r["start"] - v[:, 1])).min()
    d["dist"] = d.apply(dist, axis=1)

    for strat, sub in [("all kept BGCs", d), ("non-edge only", d[d["edge"] == 0]),
                       ("edge only", d[d["edge"] == 1])]:
        print(f"\n-- {strat} (n={len(sub)})")
        for w in WINDOWS_KB:
            near = (sub["dist"] <= w * 1000)
            lo, hi = near[sub["low"] == 1], near[sub["low"] == 0]
            if len(lo) == 0 or len(hi) == 0:
                print(f"  <= {w:>2} kb: not enough BGCs in one group"); continue
            p = fisher_exact([[lo.sum(), (~lo).sum()], [hi.sum(), (~hi).sum()]])[1]
            print(f"  <= {w:>2} kb: low-prev {lo.sum()}/{len(lo)} ({100*lo.mean():.1f}%)  "
                  f"vs higher {hi.sum()}/{len(hi)} ({100*hi.mean():.1f}%)  Fisher p={p:.3g}")

    c = d[d["dist"].notna()].copy()
    c["ld"] = np.log10(c["dist"] + 1)
    print(f"\n-- continuous distance (BGCs on contigs with a virus call, n={len(c)})")
    if len(c) > 5 and c["low"].nunique() == 2:
        r, p = spearmanr(c["ld"], c["prev"])
        print(f"  Spearman log10(dist+1) vs genomes/GCF: rho={r:.3f} p={p:.3g}")
        lo, hi = c.loc[c["low"] == 1, "dist"], c.loc[c["low"] == 0, "dist"]
        print(f"  median distance: low-prev {lo.median():.0f} bp (n={len(lo)}) vs higher {hi.median():.0f} bp (n={len(hi)})"
              f"  Mann-Whitney p={mannwhitneyu(lo, hi).pvalue:.3g}")

    try:
        import statsmodels.formula.api as smf
        d["near25"] = (d["dist"] <= 25000).astype(int)
        d["lclen"] = np.log10(d["clen"])
        d["grp"] = pd.factorize(d["fam"].fillna(pd.Series("s" + d.index.astype(str), index=d.index)))[0]
        dd = d.dropna(subset=["edge", "lclen"])
        print(f"\n-- logit near(<=25 kb) ~ low_prev + edge + log10(contig length); events={dd['near25'].sum()}; cluster=family")
        m = smf.logit("near25 ~ low + edge + lclen", dd).fit(disp=0, cov_type="cluster",
                                                             cov_kwds={"groups": dd["grp"]})
        print(m.summary().tables[1])
    except Exception as e:
        print("\n(logit skipped:", e, ")")
    return d

res = analyse(vall, "ALL geNomad virus calls")
strict = vall[(vall["virus_score"] >= 0.8) & (vall["n_hallmarks"] >= 1)]
analyse(strict, "STRICT calls (virus_score >= 0.8 and >= 1 hallmark)")
res.to_csv(f"{base}/11_genomad/bgc_virus_distance.tsv", sep="\t", index=False)
