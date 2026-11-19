import pandas as pd, numpy as np, glob, re
rng = np.random.default_rng(5)
base = "/PATH/TO/virgibacillus_v2"; gd = f"{base}/11_genomad"
mge = pd.read_csv(f"{gd}/mge_regions.tsv", sep="\t")
pro = mge[mge.kind == "prophage"]; vcs = set(zip(mge[mge.kind == "virus_contig"].accession, mge[mge.kind == "virus_contig"].contig))
pro_by = {k: d[["start", "end"]].values for k, d in pro.groupby(["accession", "contig"])}

def clen(fa):
    L, cur = {}, None
    for line in open(fa):
        if line.startswith(">"): cur = line[1:].split()[0]; L[cur] = 0
        elif cur: L[cur] += len(line.strip())
    return L

reg = pd.read_csv(f"{base}/06_antismash/bgc_regions.tsv", sep="\t")
def coords(acc, f):
    s = e = None; ctg = None
    for i, line in enumerate(open(f"{base}/06_antismash/{acc}/{f}")):
        if i == 0: ctg = line.split()[1]
        m = re.search(r"Orig\. start\s*::\s*(\d+)", line); s = int(m.group(1)) if m else s
        m = re.search(r"Orig\. end\s*::\s*(\d+)", line); e = int(m.group(1)) if m else e
        if (s is not None and e is not None) or i > 30: break
    return ctg, s, e
cc = [coords(a, f) for a, f in zip(reg.accession, reg.file)]
reg["ctg"] = [x[0] for x in cc]; reg["bs"] = [x[1] + 1 if x[1] is not None else np.nan for x in cc]; reg["be"] = [x[2] for x in cc]
print("BGC regions with coordinates: %d of %d" % (int(reg.bs.notna().sum()), len(reg)))

cg = pd.read_csv(f"{base}/09_matrix/analysis_table_clonegroups.tsv", sep="\t", dtype=str)
cg = cg[cg.in_A == "True"].copy(); cg["n50"] = pd.to_numeric(cg.Contig_N50); cg = cg[cg.n50 >= 500e3]
sub = reg[reg.accession.isin(cg.accession) & (~reg.contig_edge) & reg.bs.notna()].copy()
lens = {a: clen(f"{base}/02_isolates/fna/usr_{a}.fna") for a in cg.accession}
print("Contiguous representatives: %d | complete BGC regions analysed: %d" % (len(cg), len(sub)))

def hit(acc, ctg, s, e, pad):
    arr = pro_by.get((acc, ctg))
    if arr is None: return False
    return bool(((arr[:, 0] - pad <= e) & (arr[:, 1] + pad >= s)).any())
def pnull(acc, ctg, s, e, pad, n=1500):
    arr = pro_by.get((acc, ctg)); C = lens[acc].get(ctg)
    if arr is None or not C or C <= (e - s): return 0.0
    st = rng.integers(1, C - (e - s) + 1, n); en = st + (e - s)
    ov = np.zeros(n, bool)
    for a, b in arr: ov |= (a - pad <= en) & (b + pad >= st)
    return ov.mean()
for pad in (0, 10000):
    sub["obs"] = [hit(a, c, s, e, pad) for a, c, s, e in zip(sub.accession, sub.ctg, sub.bs, sub.be)]
    p_i = np.array([pnull(a, c, s, e, pad) for a, c, s, e in zip(sub.accession, sub.ctg, sub.bs, sub.be)])
    null = (rng.random((3000, len(p_i))) < p_i).sum(1); o = int(sub.obs.sum())
    print("\nBGC regions %s a prophage (pad = %d bp): observed %d | expected by chance %.1f | p(>=obs) = %.3f | p(<=obs) = %.3f" % (
        "overlapping" if pad == 0 else "within 10 kb of", pad, o, p_i.sum(), (np.sum(null >= o) + 1) / 3001, (np.sum(null <= o) + 1) / 3001))
sub["on_virus_contig"] = [(a, c) in vcs for a, c in zip(sub.accession, sub.ctg)]
print("BGC regions on a whole-contig virus call:", int(sub.on_virus_contig.sum()))
print(sub[sub.obs][["accession", "ctg", "product", "bs", "be"]].head(15).to_string(index=False))

print("\n================ accessory genome decomposition (33 Virgibacillus representatives)")
b = pd.read_csv(f"{base}/10_pangenome/b_reps.tsv", sep="\t", dtype=str)
res = sorted(glob.glob(f"{base}/10_pangenome/proteins/OrthoFinder/Results_*"))[-1]
og = pd.read_csv(f"{res}/Orthogroups/Orthogroups.tsv", sep="\t", index_col=0, dtype=str)
un = pd.read_csv(f"{res}/Orthogroups/Orthogroups_UnassignedGenes.tsv", sep="\t", index_col=0, dtype=str)
cat = pd.read_csv(f"{base}/10_pangenome/pangenome_categories.tsv", sep="\t", index_col=0).category
def gene_in_pro(acc):
    out = {}
    for line in open(f"{base}/10_pangenome/proteins/{acc}.faa"):
        if line.startswith(">"):
            p = line[1:].strip().split(" # "); gid = p[0].split()[0]; ctg = gid.rsplit("_", 1)[0]
            mid = (int(p[1]) + int(p[2])) / 2; arr = pro_by.get((acc, ctg))
            out[gid] = bool(arr is not None and ((arr[:, 0] <= mid) & (arr[:, 1] >= mid)).any())
    return out
rows = []
for acc in b.accession:
    inp = gene_in_pro(acc)
    for tab, label in ((og, None), (un, "genome-specific")):
        if acc not in tab.columns: continue
        for ogname, cell in tab[acc].dropna().items():
            c = label or cat.get(ogname)
            if c is None or (isinstance(c, float)): continue
            for gid in str(cell).split(", "):
                if gid.strip() in inp: rows.append((acc, ogname, c, inp[gid.strip()]))
d = pd.DataFrame(rows, columns=["accession", "og", "category", "in_prophage"])
t = d.groupby("category").agg(genes=("og", "size"), in_prophage=("in_prophage", "sum"))
t["pct_in_prophage"] = (100 * t.in_prophage / t.genes).round(2)
print(t.reindex(["core", "soft core", "shell", "cloud", "genome-specific"]).dropna().to_string())
print("\nTotal genes: %d | in prophages: %d (%.2f%%)" % (len(d), int(d.in_prophage.sum()), 100 * d.in_prophage.mean()))
acc_genes = d[d.category.isin(["shell", "cloud", "genome-specific"])]
print("Accessory genes (shell + cloud + genome-specific): %d | in prophages: %d (%.1f%%) | outside prophages: %d" % (
    len(acc_genes), int(acc_genes.in_prophage.sum()), 100 * acc_genes.in_prophage.mean(), int((~acc_genes.in_prophage).sum())))
