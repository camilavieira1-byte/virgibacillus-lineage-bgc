import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu, fisher_exact
base = "/PATH/TO/virgibacillus_v2"
cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str)
cg = cg[cg.in_A == "True"].copy()
for c in ["n_bgc", "n_edge", "Contig_N50"]: cg[c] = pd.to_numeric(cg[c])
cg["group"] = np.where(cg.in_B == "True", "Virgibacillus (B)", np.where(cg.genus.str.startswith("Virgibacillus"), "other Virgibacillus_X", "other genera"))
reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t")
reg = reg[reg.accession.isin(cg.accession)].copy()
nop = reg[reg["product"] != "terpene-precursor"]
cg["n_nopre"] = cg.accession.map(nop.groupby("accession").size()).fillna(0)
cg["n_nopre_complete"] = cg.accession.map(nop[~nop.contig_edge].groupby("accession").size()).fillna(0)

def summarize(d, label):
    print("\n=========== %s: %d representatives" % (label, len(d)))
    for g, x in d.groupby("group"):
        q = x.n_bgc.quantile([.25, .75]).values
        print("%-22s n=%3d | all regions: median %.0f (IQR %.0f-%.0f) | no terpene-precursor: median %.0f | complete, no precursor: median %.0f | median N50 %.0f kb" % (
            g, len(x), x.n_bgc.median(), q[0], q[1], x.n_nopre.median(), x.n_nopre_complete.median(), x.Contig_N50.median() / 1e3))
    b = d[d.group == "Virgibacillus (B)"]; o = d[d.group != "Virgibacillus (B)"]
    for col, name in (("n_bgc", "all regions"), ("n_nopre", "regions without terpene-precursor"), ("n_nopre_complete", "complete regions without terpene-precursor")):
        u, p = mannwhitneyu(b[col], o[col], alternative="two-sided")
        print("  B vs all others, %-44s Mann-Whitney p = %.3g | rank-biserial = %+.2f (positive: B higher)" % (name + ":", p, 2 * u / (len(b) * len(o)) - 1))
    tok = reg[reg.accession.isin(d.accession)].assign(tok=lambda x: x["product"].str.split("+")).explode("tok")
    pres = pd.crosstab(tok.accession, tok.tok).clip(upper=1).reindex(d.accession).fillna(0)
    isB = (d.group == "Virgibacillus (B)").values
    rows = []
    for k in [c for c in pres.columns if c != "terpene-precursor" and pres[c].sum() >= 8]:
        a = pres[k].values.astype(bool)
        t = [[a[isB].sum(), (~a[isB]).sum()], [a[~isB].sum(), (~a[~isB]).sum()]]
        rows.append((k, 100 * a[isB].mean(), 100 * a[~isB].mean(), fisher_exact(t)[1]))
    r = pd.DataFrame(rows, columns=["class", "pct_B", "pct_others", "p"]).sort_values("p")
    order = np.argsort(r.p.values); q = r.p.values[order] * len(r) / np.arange(1, len(r) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; r["fdr"] = np.nan; r.iloc[order, r.columns.get_loc("fdr")] = np.minimum(q, 1)
    print("\n  Share of representatives carrying each class (classes in >= 8 genomes):")
    print(r.round(3).to_string(index=False))

summarize(cg, "all representatives of the clade")
summarize(cg[cg.Contig_N50 >= 500e3], "N50 >= 500 kb only")
